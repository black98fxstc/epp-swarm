#include <stack>
#include <queue>
#include <cassert>
#include <thread>
#include <random>
#include <mutex>
#include <condition_variable>
#include <iostream>
#include <fstream>

#include "constants.h"
#include "client.h"
#include "worker.h"
#include "boundary.h"
#include "transform.h"
#include "modal.h"

namespace EPP
{
    static Transform transform;

    // pursue a particular X, Y pair
    template <class ClientSample>
    class PursueProjection : public Work<ClientSample>
    {
        friend class MATLAB_Pursuer;
        friend class MATLAB_Local;
        friend class CloudPursuer;

    protected:
    public:
        static void start(
            Request<ClientSample> *request) noexcept;

        Candidate *const candidate;

        // this is filtering with a progressively narrower Gaussian
        // which gives a wider estimation kernel, i.e., more smoothing, lower resolution
        void applyKernel(FFTData &cosine, FFTData &filtered, int pass) noexcept
        {
            double k[N + 1];
            double width = this->parameters.W;
            for (int i = 1; i < pass; i++)
                width *= 1.5; // each pass increases width by 1/2
            for (int i = 0; i <= N; i++)
                k[i] = exp(-i * i * width * width * pi * pi * 2);

            float *data = *cosine;
            float *smooth = *filtered;
            for (int i = 0; i <= N; i++)
            {
                for (int j = 0; j < i; j++)
                {
                    smooth[i + (N + 1) * j] = (float)(data[i + (N + 1) * j] * k[i] * k[j]);
                    smooth[j + (N + 1) * i] = (float)(data[j + (N + 1) * i] * k[j] * k[i]);
                }
                smooth[i + (N + 1) * i] = (float)(data[i + (N + 1) * i] * k[i] * k[i]);
            }
        }

        PursueProjection(
            Request<ClientSample> *request,
            const int X,
            const int Y) noexcept
            : Work<ClientSample>(request), candidate(new Candidate(request->sample, X, Y)){};

        ~PursueProjection() = default;

        virtual void parallel() noexcept;

        virtual void serial() noexcept;
    };

    template <class ClientSample>
    class QualifyMeasurement : public Work<ClientSample>
    {
    protected:
    public:
        const Measurement X;
        double KLDn = 0;
        double KLDe = 0;
        bool qualified = false;

        QualifyMeasurement(
            Request<ClientSample> *request,
            const int X) noexcept
            : Work<ClientSample>(request), X(X){};

        virtual void parallel() noexcept;

        virtual void serial() noexcept;
    };

    // pursue a particular X, Y pair
    template <class ClientSample>
    void PursueProjection<ClientSample>::parallel() noexcept
    {
        thread_local FFTData weights;
        // compute the weights and sample statistics from the data for this subset
        Event n = 0;
        weights.zero();
        double Sx = 0, Sy = 0, Sxx = 0, Sxy = 0, Syy = 0;
        for (Event event = 0; event < this->sample.events; event++)
            if (this->subset->contains(event))
            {
                ++n;
                double x = this->sample(event, candidate->X);
                double y = this->sample(event, candidate->Y);

                int i = (int)(x * N);
                int j = (int)(y * N);
                double dx = x * N - i;
                double dy = y * N - j;
                weights[i + (N + 1) * j] += (float)((1 - dx) * (1 - dy));
                weights[i + 1 + (N + 1) * j] += (float)(dx * (1 - dy));
                weights[i + (N + 1) * j + (N + 1)] += (float)((1 - dx) * dy);
                weights[i + 1 + (N + 1) * j + (N + 1)] += (float)(dx * dy);

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
                applyKernel(cosine, filtered, ++candidate->pass);
                // inverse discrete cosine transform
                // gives a smoothed density estimator
                transform.reverse(filtered, density);
                // modal clustering
                candidate->clusters = modal.findClusters(*density, candidate->pass, this->parameters);
            } while (candidate->clusters > this->parameters.max_clusters);
            if (candidate->clusters < 2)
            {
                candidate->outcome = Status::EPP_no_cluster;
                return;
            }

            // Kullback-Leibler Divergence
            if (KLD == 0)
            { // if it was complex enough to go around this is unlikely to be relevant the second time
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
                if (KLD < this->parameters.kld.Normal2D)
                {
                    candidate->outcome = Status::EPP_not_interesting;
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
        unsigned long int *cluster_weight = nullptr;
        if (this->parameters.goal == Parameters::Goal::best_balance)
        {
            cluster_weight = new unsigned long int[candidate->clusters + 1];
            std::fill(cluster_weight, cluster_weight + candidate->clusters + 1, 0);
            for (Event event = 0; event < this->sample.events; event++)
                if (this->subset->contains(event))
                {
                    double x = this->sample(event, candidate->X);
                    double y = this->sample(event, candidate->Y);
                    short cluster = cluster_map->colorAt(x, y);
                    ++cluster_weight[cluster];
                }
        }

        // get the dual graph of the map
        auto graph = cluster_bounds.getDualGraph();

        // pile of graphs to consider
        std::stack<DualGraph> pile;
        pile.push(*graph);

        // find and score simple sub graphs
        double best_score = std::numeric_limits<double>::infinity();
        Booleans best_edges;
        Booleans best_clusters;
        double best_balance_factor;
        double best_edge_weight;
        unsigned long int best_in_weight, best_out_weight;
        while (!pile.empty())
        {
            ++candidate->graphs;
            DualGraph graph = pile.top();
            pile.pop();
            if (graph.isSimple()) // one edge, i.e., two populations
            {
                // because the mode is always in the first cluster
                // we choose the the node that includes it, i.e.,
                // the "in" set will always include the sample mode
                Booleans in_clusters = graph.left() & 1 ? graph.left() : graph.right();

                Booleans dual_edges = graph.edge();
                double edge_weight = 0;
                for (unsigned int i = 0; i < edges.size(); i++)
                {
                    if (dual_edges & (1 << i))
                        edge_weight += edges[i].weight;
                }
                edge_weight /= 8 * N * N; // approximates number of events within a border region of width W;
                double score = edge_weight;
                double balanced_factor = 0;
                unsigned long int in_weight = 0;
                if (this->parameters.goal == Parameters::Goal::best_balance)
                {
                    for (unsigned int i = 1; i <= candidate->clusters; i++)
                    {
                        if (in_clusters & (1 << (i - 1)))
                            in_weight += cluster_weight[i];
                    }
                    if (in_weight == 0 || in_weight == n) // empty cluster!
                    {
                        continue;
                    }
                    double P = (double)in_weight / (double)n;
                    balanced_factor = 4 * P * (1 - P);
                    score /= balanced_factor;
                }
                assert(score > 0);
                // score this separatrix
                if (score < best_score)
                {
                    best_score = score;
                    best_edges = dual_edges;
                    best_clusters = in_clusters;
                    best_in_weight = in_weight;
                    best_out_weight = (unsigned long)(n - in_weight);
                    best_balance_factor = balanced_factor;
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
        delete[] cluster_weight;
        if (best_score == std::numeric_limits<double>::infinity())
        {
            candidate->outcome = Status::EPP_no_cluster;
            return;
        }

        thread_local ColoredBoundary subset_boundary;
        subset_boundary.clear();
        for (unsigned int i = 0; i < edges.size(); i++)
        {
            if (best_edges & (1 << i))
            {
                ColoredEdge edge = edges[i];
                bool lefty = best_clusters & (1 << (edge.widdershins - 1));
                subset_boundary.addEdge(edge.points, !lefty, lefty);
                // end points on the boundaries of data space are vertices
                ColoredPoint point = edge.points[0];
                if (point.i == 0 || point.i == N || point.j == 0 || point.j == N)
                    subset_boundary.addVertex(point);
                point = edge.points[edge.points.size() - 1];
                if (point.i == 0 || point.i == N || point.j == 0 || point.j == N)
                    subset_boundary.addVertex(point);
            }
        }
        subset_boundary.setColorful(2);

        candidate->outcome = Status::EPP_success;

        ColoredEdge separatrix = subset_boundary.getEdges().at(0);
        candidate->separatrix.reserve(separatrix.points.size());
        for (ColoredPoint cp : separatrix.points)
            candidate->separatrix.push_back(cp);
        if (separatrix.widdershins)
            std::reverse(candidate->separatrix.begin(), candidate->separatrix.end());

        // create in/out subsets
        auto subset_map = subset_boundary.getMap();
        for (Event event = 0; event < this->sample.events; event++)
            if (this->subset->contains(event))
            {
                double x = this->sample(event, candidate->X);
                double y = this->sample(event, candidate->Y);
                bool member = subset_map->colorAt(x, y);
                if (member)
                {
                    ++candidate->in_events;
                    candidate->in.member(event, true);
                }
                else
                {
                    ++candidate->out_events;
                    candidate->out.member(event, true);
                }
            }

        candidate->score = best_score;
        candidate->edge_weight = best_edge_weight;
        candidate->balance_factor = best_balance_factor;
    }

    template <class ClientSample>
    void PursueProjection<ClientSample>::serial() noexcept
    {
        this->request->projections++;
        this->request->passes += candidate->pass;
        this->request->clusters += candidate->clusters;
        this->request->graphs += candidate->graphs;

        // keep the finalists in order, even failures get inserted so we return some error message
        int i = this->request->candidates.size();
        if (i < this->parameters.finalists)
            this->request->candidates.push_back(candidate);
        else if (*candidate < *this->request->candidates[--i])
            delete this->request->candidates[i];
        else
        {
            delete candidate;
            return;
        };
        for (; i > 0 && *candidate < *this->request->candidates[i - 1]; i--)
            this->request->candidates[i] = this->request->candidates[i - 1];
        this->request->candidates[i] = candidate;
    }

    template <class ClientSample>
    void QualifyMeasurement<ClientSample>::parallel() noexcept
    {
        thread_local struct Scratch
        {
            float *data;
            unsigned long int size;
        } scratch = {nullptr, 0};
        if (scratch.size < this->sample.events + 1)
        {
            delete[] scratch.data;
            scratch.data = new float[this->sample.events + 1];
        }

        // get statistics for this measurement for this subset
        float *x = scratch.data;
        float *p = x;
        double Sx = 0, Sxx = 0;
        long n = 0;
        for (Event event = 0; event < this->sample.events; event++)
            if (this->subset->contains(event))
            {
                float value = (float)this->sample(event, X);
                ++n;
                Sx += value;
                Sxx += value * value;
                *p++ = value;
            }
        const double Mx = Sx / (double)n;
        const double sigma = sqrt((Sxx - Sx * Mx) / (double)(n - 1));

        // compute Kullback-Leibler Divergence
        std::sort(x, x + n);
        x[n] = 1;
        if (sigma > 0)
        {
            // const double sqrt2 = sqrt(2);
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
        qualified = KLDn > this->parameters.kld.Normal1D && KLDe > this->parameters.kld.Exponential1D;
        if (qualified)
        {
            // start pursuit on this measurement vs all the others found so far
            for (int Y : this->request->qualified)
                Worker<ClientSample>::enqueue(
                    new PursueProjection<ClientSample>(this->request, X < Y ? X : Y, X < Y ? Y : X));

            this->request->qualified.push_back(X);
        }
    }

    template <class ClientSample>
    void PursueProjection<ClientSample>::start(
        Request<ClientSample> *request) noexcept
    {
        const SampleSubset<ClientSample> *subset = request->subset;
        const Parameters &parameters = request->parameters;
        for (Measurement measurement = 0; measurement < subset->sample.measurements; ++measurement)
            if (std::find(begin(parameters.censor), end(parameters.censor), measurement) == end(parameters.censor))
                Worker<ClientSample>::enqueue(
                    new QualifyMeasurement<ClientSample>(request, measurement));
    }
}