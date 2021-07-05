#include <stack>
#include <queue>
#include <cassert>
#include <thread>
#include <mutex>
#include <condition_variable>

#include <fftw3.h>

#include "constants.h"
#include "client.h"
#include "boundary.h"
#include "modal.h"
// #include "pursuer.h"

namespace EPP
{
    static std::recursive_mutex mutex;
    static std::condition_variable_any work_available;
    static std::condition_variable_any work_completed;
    static int work_outstanding;

    std::shared_ptr<Result> _result;

    // abstract class representing a unit of work to be done
    // virtual functions let subclasses specialize tasks
    // handles work_completed and work_outstanding
    static std::vector<int> qualified_measurements;

    template <class ClientSample>
    class Work
    {
    public:
        const ClientSample sample;
        const Parameters parameters;

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
        explicit Work(
            const ClientSample &sample,
            const Parameters Parameters) noexcept
            : sample(sample), parameters(parameters)
        {
            std::unique_lock<std::recursive_mutex> lock(EPP::mutex);
            ++work_outstanding;
        };
    };

    volatile static bool kiss_of_death = false;

    // a generic worker thread. looks for work, does it, deletes it
    // virtual functions in the work object do all the real work
    template <class ClientSample>
    class Worker
    {
    public:
        static std::queue<Work<ClientSample> *> work_list;

        Worker() noexcept
        {
            std::unique_lock<std::recursive_mutex> lock(mutex);
            while (!kiss_of_death)
                if (work_list.empty())
                    work_available.wait(lock);
                else
                {
                    Work<ClientSample> *work = work_list.front();
                    work_list.pop();
                    lock.unlock();
                    work->parallel();
                    lock.lock();
                    work->serial();
                    delete work;
                }
        };
    };

    template <class ClientSample>
    std::queue<Work<ClientSample> *> Worker<ClientSample>::work_list;

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

    // pursue a particular X, Y pair
    template <class ClientSample>
    class PursueProjection : public Work<ClientSample>
    {
        // this is filtering with a progressively wider Gaussian kernel
        static void applyKernel(FFTData &cosine, FFTData &filtered, int pass) noexcept
        {
            const double pi = 3.14159265358979323846;
            double k[N + 1];
            double width = W * pass;
            for (int i = 0; i <= N; i++)
                k[i] = exp(-i * i * width * width * pi * pi * 2);

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
        Result::Status outcome;
        double best_score, best_balance_factor, best_edge_weight;
        ColoredEdge<short, bool> separatrix;
        std::vector<bool> in;
        std::vector<bool> out;
        long graph_count, in_events, out_events;
        int clusters;
        int pass;

        PursueProjection(
            const ClientSample sample,
            const Parameters parameters,
            const int X,
            const int Y) noexcept
            : Work<ClientSample>(sample, parameters), X(X), Y(Y){};

        ~PursueProjection() = default;

        virtual void parallel() noexcept;

        virtual void serial() noexcept;

        static std::shared_ptr<Result> start(const int measurements, const long events, const float *const data, std::vector<bool> &subset) noexcept;
    };

    template <class ClientSample>
    class QualifyMeasurement : public Work<ClientSample>
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
            const ClientSample sample,
            const Parameters parameters,
            const int X)
            : Work<ClientSample>(sample, parameters), X(X){};

        virtual void parallel() noexcept;

        virtual void serial() noexcept;
    };

    // pursue a particular X, Y pair
    template <class ClientSample>
    void PursueProjection<ClientSample>::parallel() noexcept
    {
        if (Y < X)
            std::swap(X, Y);
        pass = 0;
        clusters = 0;
        graph_count = 0;
        thread_local FFTData weights;
        // compute the weights and sample statistics from the data for this subset
        long n = 0;
        weights.zero();
        double Sx = 0, Sy = 0, Sxx = 0, Sxy = 0, Syy = 0;
        for (long event = 0; event < Work<ClientSample>::sample.events; event++)
            if (Work<ClientSample>::sample.subset[event])
            {
                ++n;
                double x = Work<ClientSample>::sample(event, X);
                double y = Work<ClientSample>::sample(event, Y);

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

        double KLD = 0;
        thread_local FFTData filtered;
        thread_local FFTData density;
        thread_local ModalClustering modal;
        thread_local ClusterBoundary cluster_bounds;
        std::vector<ClusterEdge> edges;
        do
        {
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
                outcome = Result::EPP_no_cluster;
                return;
            }

            // Kullback-Leibler Divergence
            if (KLD == 0)
            { // if it was complex enough to go around this is unlikely to be relavant the second time
                double NQ = 0;
                double NP = 0;
                for (int i = 0; i <= N; i++)
                    for (int j = 0; j <= N; j++)
                    {
                        double p = density[i + (N + 1) * j]; // density is *not* normalized
                        NP += p;
                        if (p <= 0)
                            continue;
                        double x = i / (double)N - Mx;
                        double y = j / (double)N - My;
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
                    outcome = Result::EPP_not_interesting;
                    return;
                }
            }

            modal.getBoundary(*density, cluster_bounds);
            // get the edges, which have their own weights
            edges = cluster_bounds.getEdges();
            // smooth some more if graph is too complex to process
        } while (edges.size() > max_booleans);

        // compute the cluster weights
        auto cluster_map = cluster_bounds.getMap();
        long cluster_weight[clusters + 1];
        std::fill(cluster_weight, cluster_weight + clusters + 1, 0);
        for (long event = 0; event < Work<ClientSample>::sample.events; event++)
            if (Work<ClientSample>::sample.subset[event])
            {
                double x = Work<ClientSample>::sample(event, X);
                double y = Work<ClientSample>::sample(event, Y);
                short cluster = cluster_map->colorAt(x, y);
                ++cluster_weight[cluster];
            }

        // get the dual graph of the map
        auto graph = cluster_bounds.getDualGraph();

        // pile of graphs to consider
        std::stack<DualGraph> pile;
        pile.push(*graph);

        // find and score simple sub graphs
        best_score = std::numeric_limits<double>::infinity();
        booleans best_edges;
        booleans best_clusters;
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
                if (left_weight == 0 || left_weight == n) // empty cluster!
                {
                    // std::cout << "empty cluster" << std::endl;
                    continue;
                }
                booleans dual_edges = graph.edge();
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
            outcome = Result::EPP_no_cluster;
            return;
        }

        thread_local ColoredBoundary<short, bool, booleans> subset_boundary;
        subset_boundary.clear();
        for (int i = 0; i < edges.size(); i++)
        {
            if (best_edges & (1 << i))
            {
                ColoredEdge<short, short> edge = edges[i];
                bool lefty = best_clusters & (1 << (edge.widdershins - 1));
                subset_boundary.addEdge(edge.points, lefty, !lefty);
                // end points on the boundaries of data space are vertices
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
        in.resize(Work<ClientSample>::sample.events);
        in.clear();
        in_events = 0;
        out.resize(Work<ClientSample>::sample.events);
        out.clear();
        out_events = 0;

        auto subset_map = subset_boundary.getMap();
        for (long event = 0; event < Work<ClientSample>::sample.events; event++)
            if (Work<ClientSample>::sample.subset[event])
            {
                double x = Work<ClientSample>::sample(event, X);
                double y = Work<ClientSample>::sample(event, Y);
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

        separatrix = subset_boundary.getEdges().at(0);
        outcome = Result::Status::EPP_success;
    }

    template <class ClientSample>
    void PursueProjection<ClientSample>::serial() noexcept
    {
        // std::cout << "pass " << pass << " clusters " << clusters << " graphs considered " << graph_count << std::endl;
        // std::cout << "pursuit completed " << X << " vs " << Y << "  ";
        switch (outcome)
        {
        case Result::EPP_error:
            _result->outcome = outcome;
            break;
        case Result::EPP_no_cluster:
        case Result::EPP_not_interesting:
            _result->outcome == outcome;
            break;
        case Result::EPP_success:
            if (_result->best_score > best_score)
            {
                _result->outcome = outcome;
                _result->X = X;
                _result->Y = Y;
                _result->separatrix.clear();
                ;
                _result->separatrix.reserve(separatrix.points.size());
                for (ColoredPoint<short> cp : separatrix.points)
                    _result->separatrix.push_back(Point(cp.i, cp.j));
                if (separatrix.widdershins)
                    std::reverse(_result->separatrix.begin(), _result->separatrix.end());
                _result->best_score = best_score;
                _result->edge_weight = best_edge_weight;
                _result->balance_factor = best_balance_factor;
                _result->in = in;
                _result->in_events = in_events;
                _result->out = out;
                _result->out_events = out_events;
                _result->pass = pass;
                _result->clusters = clusters;
                _result->graphs = graph_count;
                _result->qualified = qualified_measurements;
            }
            break;
        }
    }

    template <class ClientSample>
    void QualifyMeasurement<ClientSample>::parallel() noexcept
    {
        // get statistics for this measurement for this subset
        thread_local QualifyMeasurement::Scratch scratch;
        float *x = scratch.reserve(Work<ClientSample>::sample.events + 1);
        float *p = x;
        double Sx = 0, Sxx = 0;
        long n = 0;
        for (long event = 0; event < Work<ClientSample>::sample.events; event++)
            if (Work<ClientSample>::sample.subset[event])
            {
                float value = Work<ClientSample>::sample(event, X);
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

    template <class ClientSample>
    void QualifyMeasurement<ClientSample>::serial() noexcept
    {
        qualified = KLDn > .16 && KLDe > .16;
        if (qualified)
        {
            // start pursuit on this measurement vs all the others found so far
            for (int Y : qualified_measurements)
                Worker<ClientSample>::work_list.push(new PursueProjection<ClientSample>(Work<ClientSample>::sample, Work<ClientSample>::parameters, X, Y));
            work_available.notify_all();

            qualified_measurements.push_back(X);
            // std::cout << "dimension qualified " << X << std::endl;
        }
        // else
        // std::cout << "dimension disqualified " << X << std::endl;
    }

    FFTData::~FFTData()
    {
        if (data)
            fftwf_free(data);
    }

    float *FFTData::operator*() noexcept
    {
        if (!data)
            data = (float *)fftw_malloc(sizeof(float) * (N + 1) * (N + 1));
        return data;
    }

    void FFTData::zero() noexcept
    {
        if (!data)
            data = (float *)fftw_malloc(sizeof(float) * (N + 1) * (N + 1));
        std::fill(data, data + (N + 1) * (N + 1), 0);
    }

    Transform::Transform() noexcept
    {
        FFTData in;
        FFTData out;
        // FFTW planning is slow and not thread safe so we do it here
        DCT = (void *)fftwf_plan_r2r_2d((N + 1), (N + 1), *in, *out,
                                        FFTW_REDFT00, FFTW_REDFT00, 0);
        // actually they are the same in this case but leave it for now
        IDCT = (void *)fftwf_plan_r2r_2d((N + 1), (N + 1), *in, *out,
                                         FFTW_REDFT00, FFTW_REDFT00, 0);
        assert(DCT && IDCT);
    }

    Transform::~Transform() noexcept
    {
        fftwf_destroy_plan((fftwf_plan)DCT);
        fftwf_destroy_plan((fftwf_plan)IDCT);
    }

    void Transform::forward(FFTData &in, FFTData &out) noexcept
    {
        fftwf_execute_r2r((fftwf_plan)DCT, *in, *out);
    }

    void Transform::reverse(FFTData &in, FFTData &out) noexcept
    {
        fftwf_execute_r2r((fftwf_plan)IDCT, *in, *out);
    }
}