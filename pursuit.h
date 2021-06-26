#include <stack>
#include <queue>
#include <cassert>
#include <thread>
#include <mutex>
#include <condition_variable>

#include <fftw3.h>

#include "boundary.h"
#include "modal.h"

namespace EPP
{
    static std::default_random_engine generator;
    static std::recursive_mutex mutex;
    static std::condition_variable_any work_available;
    static std::condition_variable_any work_completed;
    static int work_outstanding;

    // these are essential constants that are read only
    // so safely shared by all threads
    struct worker_sample
    {
        const int measurements;
        const long events;
        const float *const data;
        const std::vector<bool> subset;
    };

    struct worker_output
    {
        enum worker_result
        {
            EPP_success,
            EPP_no_qualified,
            EPP_no_cluster,
            EPP_not_interesting,
            EPP_error
        } outcome;
        int X, Y;
        double edge_weight, balance_factor;
        volatile double best_score;
        std::vector<bool> in, out;
        long in_events, out_events;
        ClusterSeparatrix separatrix;
    } result;

    // abstract class representing a unit of work to be done
    // virtual functions let subclasses specialize tasks
    // handles work_completed and work_outstanding

    class Work
    {
    public:
        const struct worker_sample sample;

        // many threads can execute this in parallel
        virtual void parallel() noexcept
        {
            assert(("unimplemented", false));
        };
        // then only one thread at a time can run
        virtual void serial() noexcept
        {
            assert(("unimplemented", false));
        };

        ~Work()
        {
            std::unique_lock<std::recursive_mutex> lock(EPP::mutex);
            if (--work_outstanding == 0)
                work_completed.notify_all();
        };

    protected:
        explicit Work(const worker_sample &sample) noexcept : sample(sample)
        {
            std::unique_lock<std::recursive_mutex> lock(EPP::mutex);
            ++work_outstanding;
        };
    };

    volatile static bool kiss_of_death = false;
    static std::queue<Work *> work_list;

    // a generic worker thread. looks for work, does it, deletes it
    // virtual functions in the work object do all the real work
    class Worker
    {
    public:
        Worker() noexcept
        {
            std::unique_lock<std::recursive_mutex> lock(mutex);
            while (!kiss_of_death)
                if (work_list.empty())
                    work_available.wait(lock);
                else
                {
                    Work *work = work_list.front();
                    work_list.pop();
                    lock.unlock();
                    work->parallel();
                    lock.lock();
                    work->serial();
                    delete work;
                }
        };
    };

    // pursue a particular X, Y pair
    class PursueProjection : public Work
    {
        // fftw needs special alignment to take advantage of vector instructions
        class FFTData
        {
            float *data;

        public:
            FFTData()
            {
                data = nullptr;
            }

            ~FFTData();

            float *operator*() noexcept;

            inline float &operator[](const int i)
            {
                return data[i];
            }

            void zero() noexcept;
        };

        class Transform
        {
            void *DCT;
            void *IDCT;

        public:
            Transform() noexcept;

            ~Transform();

            void forward(FFTData &in, FFTData &out) noexcept;

            void reverse(FFTData &in, FFTData &out) noexcept;
        };

        static Transform transform;

        // this is filtering with a progressively wider gausian kernel
        void applyKernel(FFTData &cosine, FFTData &filtered, int pass) noexcept
        {
            double k[N + 1];
            double width = .001 * N * pass;
            for (int i = 0; i <= N; i++)
                k[i] = exp(-i * i * width * width);

            float *data = *cosine;
            float *smooth = *filtered;
            for (int i = 0; i <= N; i++)
            {
                for (int j = 0; j < i; j++)
                {
                    smooth[i + (N + 1) * j] = data[i + (N + 1) * j] * k[i] * k[j];
                    smooth[j + (N + 1) * i] = data[j + (N + 1) * i] * k[j] * k[i];
                }
                smooth[i + (N + 1) * i] = data[i + (N + 1) * i] * k[i] * k[i];
            }
        }

    public:
        int X, Y;
        worker_output::worker_result outcome;
        double best_score, best_balance_factor, best_edge_weight;
        ClusterSeparatrix separatrix;
        std::vector<bool> in;
        std::vector<bool> out;
        long graph_count, in_events, out_events;

        PursueProjection(
            const worker_sample sample,
            const int X,
            const int Y) noexcept
            : Work(sample), X(X), Y(Y){};

        ~PursueProjection() = default;
        ;

        virtual void parallel() noexcept;

        virtual void serial() noexcept;

        static worker_output *start(const int measurements, const long events, const float *const data, std::vector<bool> &subset) noexcept;
    };

    static std::vector<int> qualified_measurements;

    class QualifyMeasurement : public Work
    {
        class Scratch
        {
            float *data;
            long size;

        public:
            Scratch() noexcept
            {
                data = nullptr;
                size = 0;
            }

            ~Scratch()
            {
                delete[] data;
            }

            float *&reserve(long size) noexcept
            {
                if (this->size < size)
                {
                    delete[] data;
                    data = nullptr;
                }
                if (!data)
                {
                    data = new float[size];
                    this->size = size;
                }
                return data;
            }

            inline float &operator[](const int i) noexcept
            {
                return data[i];
            }
        };

    public:
        const int X;
        double KLDn = 0;
        double KLDe = 0;
        bool qualified = false;

        QualifyMeasurement(
            const worker_sample sample,
            const int X)
            : Work(sample), X(X){};

        virtual void parallel() noexcept;

        virtual void serial() noexcept;
    };

    // pursue a particular X, Y pair
    void PursueProjection::parallel() noexcept
    {
        if (Y < X)
            std::swap(X, Y);
        thread_local PursueProjection::FFTData weights;
        // compute the weights and sample statistics from the data for this subset
        long n = 0;
        weights.zero();
        double Sx = 0, Sy = 0, Sxx = 0, Sxy = 0, Syy = 0;
        for (long event = 0; event < sample.events; event++)
            if (sample.subset[event])
            {
                ++n;
                double x = sample.data[event * sample.measurements + X];
                double y = sample.data[event * sample.measurements + Y];

                int i = (int)(x * N);
                int j = (int)(y * N);
                double dx = x * N - i;
                double dy = y * N - j;
                weights[i + (N + 1) * j] += (float)(dx * dy);
                weights[i + 1 + (N + 1) * j] += (float)((1 - dx) * dy);
                weights[i + (N + 1) * j + (N + 1)] += (float)(dx * (1 - dy));
                weights[i + 1 + (N + 1) * j + (N + 1)] += (float)((1 - dx) * (1 - dy));

                Sx += x;
                Sy += y;
                Sxx += x * x;
                Sxy += x * y;
                Syy += y * y;
            }
        double Mx = Sx / (double)n; // means
        double My = Sy / (double)n;
        double Cxx = (Sxx - Sx * Mx) / (double)(n - 1); // covariance
        double Cxy = (Sxy - Sx * My) / (double)(n - 1);
        double Cyy = (Syy - Sy * My) / (double)(n - 1);

        // discrete cosine transform (FFT of real even function)
        thread_local FFTData cosine;
        transform.forward(weights, cosine);

        int clusters;
        int pass = 0;
        thread_local FFTData filtered;
        thread_local FFTData density;
        thread_local ModalClustering modal;
        do
        {
            // apply kernel to cosine transform
            applyKernel(cosine, filtered, ++pass);

            // inverse discrete cosine transform
            // gives a smoothed density estimator
            transform.reverse(filtered, density);

            // modal clustering
            clusters = modal.findClusters(*density);
        } while (clusters > 10);
        if (clusters < 2)
        {
            std::cout << "no cluster" << std::endl;
            outcome = worker_output::EPP_no_cluster;
            return;
        }

        // Kullback-Leibler Divergence
        double KLD = 0;
        double NQ = 0;
        double NP = 0;
        for (int i = 0; i <= N; i++)
            for (int j = 0; j <= N; j++)
            {
                double p = density[i + (N + 1) * j]; // density is *not* normalized
                NP += p;
                double x = i / (double)N - Mx;
                double y = j / (double)N - My;
                if (p <= 0)
                    continue;
                // Mahalanobis distance squared over 2 is unnormalized - ln Q
                double MD2 = (x * x / Cxx - 2 * x * y * Cxy / Cxx / Cyy + y * y / Cyy) / (1 - Cxy * Cxy / Cxx / Cyy) / 2;
                NQ += exp(-MD2);
                // unnormalized P ln(P/Q) = P * (ln P - ln Q) where P is density and Q is bivariant normal
                KLD += p * (log(p) + MD2);
            }

        // Normalize the density
        KLD /= NP;
        // subtract off normalization constants factored out of the sum above
        KLD -= log(NP / NQ);
        if (KLD < .16)
        {
            std::cout << "not interesting" << std::endl;
            outcome = worker_output::worker_result::EPP_not_interesting;
            return;
        }

        thread_local ClusterBoundary cluster_bounds;
        modal.getBoundary(*density, cluster_bounds);

        // compute the cluster weights
        auto cluster_map = cluster_bounds.getMap();
        long cluster_weight[clusters + 1];
        std::fill(cluster_weight, cluster_weight + clusters + 1, 0);
        for (long event = 0; event < sample.events; event++)
            if (sample.subset[event])
            {
                double x = sample.data[event * sample.measurements + X];
                double y = sample.data[event * sample.measurements + Y];
                short cluster = cluster_map->colorAt(x, y);
                ++cluster_weight[cluster];
            }

        // get the edges, which have their own weights
        auto edges = cluster_bounds.getEdges();
        assert(edges.size() <= max_booleans);

        // get the dual graph of the map
        auto graph = cluster_bounds.getDualGraph();

        // pile of graphs to consider
        std::stack<DualGraph> pile;
        pile.push(*graph);

        // find and score simple sub graphs
        best_score = std::numeric_limits<double>::infinity();
        booleans best_edges;
        booleans best_clusters;
        graph_count = 0;
        while (!pile.empty())
        {
            ++graph_count;
            DualGraph graph = pile.top();
            pile.pop();
            if (graph.isSimple())
            { // one edge, i.e., two populations
                booleans left_clusters = graph.left();
                long left_weight = 0;
                for (int i = 1; i <= clusters; i++)
                {
                    if (left_clusters & (1 << (i - 1)))
                        left_weight += cluster_weight[i];
                }
                booleans dual_edges = graph.edge();
                if (left_weight == 0 || left_weight == n) // empty cluster!
                {
                    // std::cout << "empty cluster" << std::endl;
                    continue;
                }
                double edge_weight = 0;
                for (int i = 0; i < edges.size(); i++)
                {
                    if (dual_edges & (1 << i))
                        edge_weight += edges[i].weight;
                }
                double P = (double)left_weight / (double)n;
                double balanced_weight = edge_weight / 4 / P / (1 - P);
                assert(balanced_weight > 0);

                // score this separatrix
                if (balanced_weight < best_score)
                {
                    best_score = balanced_weight;
                    best_edges = dual_edges;
                    best_clusters = left_clusters;
                    best_balance_factor = 4 * P * (1 - P);
                    best_edge_weight = edge_weight;
                }
            }
            else
            { // not simple so simplify it some, i.e., remove one dual edge at a time
                // and merge two adjacent subsets. that makes a bunch more graphs to look at
                std::vector<DualGraph> simplified = graph.simplify();
                for (const auto &graph : simplified)
                    pile.push(graph);
            }
        }
        if (best_score == std::numeric_limits<double>::infinity())
        {
            outcome = worker_output::EPP_no_cluster;
            return;
        }

        // thread_local ColoredBoundary<short, bool> subset_boundary;
        ColoredBoundary<short, bool> subset_boundary;
        subset_boundary.clear();
        for (int i = 0; i < edges.size(); i++)
        {
            if (best_edges & (1 << i))
            {
                ColoredEdge<short, short> edge = edges[i];
                bool lefty = best_clusters & (1 << (edge.widdershins - 1));
                subset_boundary.addEdge(edge.points, lefty, !lefty);
                // end points on the boundaries of data space are verticies
                ColoredPoint<short> point = edge.points[0];
                if (point.i == 0 || point.i == N || point.j == 0 || point.j == N)
                    subset_boundary.addVertex(point);
                point = edge.points[edge.points.size() - 1];
                if (point.i == 0 || point.i == N || point.j == 0 || point.j == N)
                    subset_boundary.addVertex(point);
            }
        }
        subset_boundary.setColorful(2);

        // create in/out subsets
        in.resize(sample.events);
        in.clear();
        in_events = 0;
        out.resize(sample.events);
        out.clear();
        out_events = 0;

        auto subset_map = subset_boundary.getMap();
        for (long event = 0; event < sample.events; event++)
            if (sample.subset[event])
            {
                double x = sample.data[event * sample.measurements + X];
                double y = sample.data[event * sample.measurements + Y];
                bool member = subset_map->colorAt(x, y);
                if (member)
                {
                    ++in_events;
                    in[event] = true;
                }
                else
                {
                    ++out_events;
                    out[event] = true;
                }
            }

        outcome = worker_output::worker_result::EPP_success;
        separatrix = subset_boundary.getEdges();

        // separatrix, in and out are the payload
    }

    void PursueProjection::serial() noexcept
    {
        std::cout << graph_count << " graphs considered" << std::endl;
        std::cout << "pursuit completed " << X << " vs " << Y << "  ";
        switch (outcome)
        {
        case worker_output::worker_result::EPP_error:
        {
            std::cout << "error";
            result.outcome = outcome;
            break;
        }
        case worker_output::worker_result::EPP_no_cluster:
        case worker_output::worker_result::EPP_not_interesting:
        {
            std::cout << "no split";
            result.outcome == outcome;
            break;
        }
        case worker_output::worker_result::EPP_success:
        {
            std::cout << best_score;
            if (result.best_score > best_score)
            {
                result.outcome = outcome;
                result.X = X;
                result.Y = Y;
                result.best_score = best_score;
                result.edge_weight = best_edge_weight;
                result.balance_factor = best_balance_factor;
                result.separatrix = separatrix;
                result.in = in;
                result.in_events = in_events;
                result.out = out;
                result.out_events = out_events;
            }
            break;
        }
        }
        std::cout << std::endl;
    }

    worker_output *PursueProjection::start(const int measurements, const long events, const float *const data, std::vector<bool> &subset) noexcept
    {
        worker_sample constants{measurements, events, data, subset};
        qualified_measurements.clear();
        result.outcome = worker_output::EPP_no_qualified;
        result.best_score = std::numeric_limits<double>::infinity();

        std::unique_lock<std::recursive_mutex> lock(mutex);
        for (int measurement = 0; measurement < constants.measurements; ++measurement)
            work_list.push(new QualifyMeasurement(constants, measurement));
        work_available.notify_all();

        return &result;
    }

    PursueProjection::Transform PursueProjection::transform;

    void QualifyMeasurement::parallel() noexcept
    {
        // get statistics for this measurement for this subset
        thread_local QualifyMeasurement::Scratch scratch;
        float *x = scratch.reserve(sample.events + 1);
        float *p = x;
        double Sx = 0, Sxx = 0;
        long n = 0;
        for (long event = 0; event < sample.events; event++)
            if (sample.subset[event])
            {
                float value = sample.data[event * sample.measurements + X];
                ++n;
                Sx += value;
                Sxx += value * value;
                *p++ = value;
            }
        const double Mx = Sx / n;
        const double sigma = sqrt((Sxx - Sx * Mx) / (double)(n - 1));

        // compute Kullback-Leibler Divergence
        std::sort(x, x + n);
        x[n] = 1;
        if (sigma > 0)
        {
            const double sqrt2 = sqrt(2);
            // normalization factors for truncated distributions
            double NQn = .5 * (erf((x[n] - Mx) / sigma / sqrt2) - erf((x[0] - Mx) / sigma / sqrt2));
            double NQe = exp(-x[0] / Mx) - exp(-x[n] / Mx);
            for (long i = 0, j; i < n; i = j)
            {
                j = i + 1;
                while ((x[j] - x[i]) < .001 && j < n)
                    j++;
                double P = (double)(j - i) / (double)n;
                double Qn = .5 * (erf((x[j] - Mx) / sigma / sqrt2) - erf((x[i] - Mx) / sigma / sqrt2)) / NQn;
                double Qe = (exp(-x[i] / Mx) - exp(-x[j] / Mx)) / NQe;
                KLDn += P * log(P / Qn);
                KLDe += P * log(P / Qe);
            }
        }
    }

    void QualifyMeasurement::serial() noexcept
    {
        qualified = KLDn > .16 && KLDe > .16;
        if (qualified)
        {
            // start pursuit on this measurement vs all the others found so far
            for (int Y : qualified_measurements)
                work_list.push(new PursueProjection(sample, X, Y));
            work_available.notify_all();

            qualified_measurements.push_back(X);
            std::cout << "dimension qualified " << X << std::endl;
        }
        else
            std::cout << "dimension disqualified " << X << std::endl;
    }

    PursueProjection::FFTData::~FFTData()
    {
        if (data)
            fftwf_free(data);
    }

    float *PursueProjection::FFTData::operator*() noexcept
    {
        if (!data)
            data = (float *)fftw_malloc(sizeof(float) * (N + 1) * (N + 1));
        return data;
    }

    void PursueProjection::FFTData::zero() noexcept
    {
        if (!data)
            data = (float *)fftw_malloc(sizeof(float) * (N + 1) * (N + 1));
        std::fill(data, data + (N + 1) * (N + 1), 0);
    }

    PursueProjection::Transform::Transform() noexcept
    {
        PursueProjection::FFTData in;
        PursueProjection::FFTData out;
        // FFTW planning is slow and not thread safe so we do it here
        DCT = (void *)fftwf_plan_r2r_2d((N + 1), (N + 1), *in, *out,
                                        FFTW_REDFT00, FFTW_REDFT00, 0);
        // actually they are the same in this case but leave it for now
        IDCT = (void *)fftwf_plan_r2r_2d((N + 1), (N + 1), *in, *out,
                                         FFTW_REDFT00, FFTW_REDFT00, 0);
        assert(DCT && IDCT);
    }

    PursueProjection::Transform::~Transform() noexcept
    {
        fftwf_destroy_plan((fftwf_plan)DCT);
        fftwf_destroy_plan((fftwf_plan)IDCT);
    }

    void PursueProjection::Transform::forward(FFTData &in, FFTData &out) noexcept
    {
        fftwf_execute_r2r((fftwf_plan)DCT, *in, *out);
    }

    void PursueProjection::Transform::reverse(FFTData &in, FFTData &out) noexcept
    {
        fftwf_execute_r2r((fftwf_plan)IDCT, *in, *out);
    }
}