/*
 * Developer: Wayne Moore <wmoore@stanford.edu>
 * Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
 * License: BSD 3 clause
 */

#include <iostream>
#include <fstream>

#include "client.h"
#include "worker.h"
#include "boundary.h"
#include "modal.h"
#include "taxonomy.h"

namespace EPP
{

    // pursue a particular X, Y pair
    template <class ClientSample>
    class PursueProjection : public Work<ClientSample>
    {
        friend class MATLAB_Pursuer;
        friend class MATLAB_Local;
        friend class CloudPursuer;

    private:
        Transform &transform;
        const double *const *const kernel;

        static void applyKernel(float *cosine, float *filtered, const double *k) noexcept
        {
            for (int i = 0; i <= N; i++)
            {
                for (int j = 0; j < i; j++)
                {
                    filtered[i + (N + 1) * j] = (float)(cosine[i + (N + 1) * j] * k[i] * k[j]);
                    filtered[j + (N + 1) * i] = (float)(cosine[j + (N + 1) * i] * k[j] * k[i]);
                }
                filtered[i + (N + 1) * i] = (float)(cosine[i + (N + 1) * i] * k[i] * k[i]);
            }
        }

    public:
        Candidate *const candidate;

        static void start(
            Request<ClientSample> *request) noexcept;

        PursueProjection(
            Request<ClientSample> *request,
            Measurement X,
            Measurement Y) noexcept
            : Work<ClientSample>(request), transform(request->analysis->pursuer->transform),
              kernel(request->analysis->kernel), candidate(new Candidate(request->sample, X, Y)) {}

        virtual ~PursueProjection() = default;

        virtual void parallel() noexcept;

        virtual void serial() noexcept;
    };

    template <class ClientSample>
    class QualifyMeasurement : public Work<ClientSample>
    {
    public:
        const Measurement X;
        double KLDn = 0;
        double KLDe = 0;
        bool qualified = false;

        QualifyMeasurement(
            Request<ClientSample> *request,
            Measurement X) noexcept
            : Work<ClientSample>(request), X(X){};
        
        virtual ~QualifyMeasurement() = default;

        virtual void parallel() noexcept;

        virtual void serial() noexcept;
    };

    // pursue a particular X, Y pair
    template <class ClientSample>
    void PursueProjection<ClientSample>::parallel() noexcept
    {
        // compute the weights from the data for this subset
        thread_local Transform::Data weights(transform);
        weights.clear();
        Event n = 0;
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
                // in order to avoid the overhead of range checking, when x or y is 1 the
                // algorithm adds zeros out of bounds. For y = 1 this would be outside the
                // normal array bounds so (N+2)*(N+1)+1 memory must be allocated.
            }

        // discrete cosine transform (FFT of real even function)
        thread_local Transform::Data cosine(transform);
        transform.forward(weights, cosine);
        // transform.dump(weights, "weights.csv");

        std::vector<ColoredEdge> edges;
        thread_local Transform::Data filtered(transform);
        thread_local Transform::Data density(transform);
        thread_local Transform::Data variance(transform);
        thread_local ColoredBoundary cluster_bounds;
        thread_local ModalClustering modal;
        do
        {
            do
            {
                if (candidate->pass == max_passes)
                {
                    candidate->outcome = Status::EPP_error;
                    return;
                }
                // last density becomes this variance estimator
                std::swap(density.data, variance.data);
                // apply kernel to cosine transform
                applyKernel(cosine.data, filtered.data, this->kernel[++candidate->pass]);
                // inverse discrete cosine transform
                // gives a smoothed density estimator
                transform.reverse(filtered, density);
                // transform.dump(density, "density.csv");
                // modal clustering
                candidate->clusters = modal.findClusters(density.data, candidate->pass, this->parameters);
            } while (candidate->clusters > this->parameters.max_clusters);
            if (candidate->clusters < 2)
            {
                candidate->outcome = Status::EPP_no_cluster;
                return;
            }

            modal.getBoundary(density.data, cluster_bounds);
            // get the edges, which have their own weights
            edges = cluster_bounds.getEdges();
            // smooth some more if graph is too complex to process
        } while (edges.size() > max_booleans);

        // get the dual graph of the map
        ColoredGraph graph = cluster_bounds.getDualGraph();

        // Density Based Merging
        if (candidate->pass == 1)
        { // otherwise it was swapped in above
            applyKernel(cosine.data, filtered.data, this->kernel[candidate->pass - 1]);
            transform.reverse(filtered, variance);
            // transform.dump(variance, "variance.csv");
        }
        thread_local struct cluster_merge
        {
            float edge_max, edge_var;
            BitPosition i;

			bool operator<(
				const cluster_merge &other) const noexcept
			{ // larger taken first so sense inverted
				return edge_max > other.edge_max;
			};
        } merge[max_booleans];
        // find the highest point on each edge
        for (BitPosition i = 0; i < edges.size(); ++i)
        {
            const ColoredEdge &edge = edges[i];
            ColoredPoint point = edge.points[0];
            float edge_max = density[point.i + (N + 1) * point.j];
            for (Count j = 1; j < edge.points.size(); ++j)
            {
                ColoredPoint p = edge.points[j];
                float d = density[p.i + (N + 1) * p.j];
                if (d > edge_max)
                {
                    point = p;
                    edge_max = d;
                }
            }
            merge[i].edge_max = edge_max;
            merge[i].edge_var = variance[point.i + (N + 1) * point.j];
            merge[i].i = i;
        }
        std::sort(merge, merge + edges.size());
        // now test them in decreasong order of height
        for (BitPosition j = 0; j < edges.size(); ++j)
        {
            const ColoredEdge &edge = edges[merge[j].i];
            // can't apply it to half edges of border square
            if (edge.clockwise == 0 || edge.widdershins == 0)
                continue;

            double cluster_max, cluster_var;
            // the smaller of the maxima of the clusters the edge divides
            if (modal.maxima[edge.clockwise] < modal.maxima[edge.widdershins])
            {
                cluster_max = modal.maxima[edge.clockwise];
                ColoredPoint point = modal.center[edge.clockwise];
                cluster_var = variance[point.i + (N + 1) * point.j];
            }
            else
            {
                cluster_max = modal.maxima[edge.widdershins];
                ColoredPoint point = modal.center[edge.widdershins];
                cluster_var = variance[point.i + (N + 1) * point.j];
            }
            // formulas from DBM paper. 4N^2 normalizes the FFT
            // phi^2 is also gausian with sqrt(2) narrower kernel, but doesn't integrate to 1
            // compute with the normalized kernel and then 4 * radius^2 * pi corrects the integral
            double radius = this->parameters.N * this->parameters.kernelRadius(candidate->pass);
            double f_e = merge[j].edge_max / 4 / N / N / n;
            double sigma2_e = (merge[j].edge_var / 4 / N / N / 4 / radius / radius / pi / n - f_e * f_e) / (n - 1);
            double f_c = cluster_max / 4 / N / N / n;
            double sigma2_c = (cluster_var / 4 / N / N / 4 / radius / radius / pi / n - f_c * f_c) / (n - 1);
            // if the dip isn't significant, merge the two clusters and remove the edge
            if (sigma2_e > 0 && sigma2_c > 0 && f_c - std::sqrt(sigma2_c) < f_e + std::sqrt(sigma2_e))
            {
                // this may cause a lower one to become more significant
                if (modal.maxima[edge.clockwise] < modal.maxima[edge.widdershins])
                {
                    modal.maxima[edge.clockwise] = modal.maxima[edge.widdershins];
                    modal.center[edge.clockwise] = modal.center[edge.widdershins];
                }
                else
                {
                    modal.maxima[edge.widdershins] = modal.maxima[edge.clockwise];
                    modal.center[edge.widdershins] = modal.center[edge.clockwise];
                }

                graph = graph.simplify(merge[j].i);
                ++candidate->merges;
            }
        }
        // make sure there's anything left
        if (graph.isTrivial())
        {
            candidate->outcome = Status::EPP_not_significant;
            return;
        }

        // compute the cluster weights
        auto cluster_map = cluster_bounds.getMap();
        thread_local Event cluster_weight[max_booleans + 1];
        std::fill(cluster_weight, cluster_weight + candidate->clusters + 1, 0);
        for (Event event = 0; event < this->sample.events; event++)
            if (this->subset->contains(event))
            {
                double x = this->sample(event, candidate->X);
                double y = this->sample(event, candidate->Y);
                Color cluster = cluster_map->colorAt(x, y);
                ++cluster_weight[cluster];
            }

        // find and score simple sub graphs
        struct
        {
            double score = std::numeric_limits<double>::infinity();
            Booleans edges;
            Booleans clusters;
            double balance_factor;
            double edge_weight;
        } best;
        // pile of graphs to consider
        thread_local std::stack<ColoredGraph> pile;
        assert(pile.empty());
        pile.push(graph);
        while (!pile.empty())
        {
            ++candidate->graphs;
            ColoredGraph graph = pile.top();
            pile.pop();
            if (graph.isSimple()) // one edge, i.e., two populations
            {
                // because the mode is always in the first cluster
                // we choose the node that includes it, i.e.,
                // the "in" set will always include the sample mode
                Booleans in_clusters = graph.left() & 1 ? graph.left() : graph.right();

                Booleans dual_edges = graph.edge();
                double edge_weight = 0;
                for (size_t i = 0; i < edges.size(); i++)
                {
                    if (dual_edges & (1 << i))
                        edge_weight += edges[i].weight;
                }
                edge_weight /= 8 * N * N; // approximates number of events within a border region of width W;
                double score = edge_weight;
                double balance_factor = 0;
                unsigned long int in_weight = 0;
                for (BitPosition i = 1; i <= candidate->clusters; i++)
                {
                    if (in_clusters & (1 << (i - 1)))
                        in_weight += cluster_weight[i];
                }
                if (in_weight == 0 || in_weight == n) // empty cluster!
                {
                    continue;
                }
                if (this->parameters.goal == Parameters::Goal::best_balance)
                {
                    double P = (double)in_weight / (double)n;
                    balance_factor = 4 * P * (1 - P);
                    score /= balance_factor;
                }
                assert(score > 0);
                // score this separatrix
                if (score < best.score)
                {
                    best.score = score;
                    best.edges = dual_edges;
                    best.clusters = in_clusters;
                    best.balance_factor = balance_factor;
                    best.edge_weight = edge_weight;
                }
            }
            else
            { // not simple so simplify it some, i.e., remove one dual edge at a time
                // and merge two adjacent subsets. that makes a bunch more graphs to look at
                std::vector<ColoredGraph> simplified = graph.simplify();
                for (const auto &graph : simplified)
                    pile.push(graph);
            }
        }
        if (best.score == std::numeric_limits<double>::infinity())
        {
            candidate->outcome = Status::EPP_no_cluster;
            return;
        }

        // find the separatrix
        thread_local ColoredBoundary subset_boundary;
        subset_boundary.clear();
        thread_local std::vector<ColoredPoint> interior_vertex;
        assert(interior_vertex.empty());
        thread_local std::vector<BitPosition> half_edges;
        half_edges.clear();
        for (BitPosition i = 0; i < edges.size(); i++)
        {
            ColoredEdge &edge = edges[i];
            // skip half edges from vertex squares for now
            if (edge.clockwise == 0 || edge.widdershins == 0)
            {
                half_edges.push_back(i);
                continue;
            }
            if (best.edges & (1 << i))
            {
                bool lefty = best.clusters & (1 << (edge.widdershins - 1));
                subset_boundary.addEdge(edge.points, !lefty, lefty);
                // end points on the boundaries of data space are vertices
                ColoredPoint point = edge.points.front();
                if (point.i == 0 || point.i == N || point.j == 0 || point.j == N)
                    subset_boundary.addVertex(point);
                else
                {
                    auto position = std::find(interior_vertex.begin(), interior_vertex.end(), point);
                    if (position != interior_vertex.end())
                        interior_vertex.erase(position);
                    else
                        interior_vertex.push_back(point);
                }
                point = edge.points.back();
                if (point.i == 0 || point.i == N || point.j == 0 || point.j == N)
                    subset_boundary.addVertex(point);
                else
                {
                    auto position = std::find(interior_vertex.begin(), interior_vertex.end(), point);
                    if (position != interior_vertex.end())
                        interior_vertex.erase(position);
                    else
                        interior_vertex.push_back(point);
                }
            }
        }
        while (!interior_vertex.empty())
        {   // need to fill in using some of the half edges
            bool making_progress = false;
            ColoredPoint end_point = interior_vertex.back(); 
            // so it will be the last one pushed bellow if any
            // there are two half edges and we must always extend the
            // same one to completion or it can deadlock or switchback
            for (auto hep = half_edges.begin(); hep != half_edges.end(); ++hep)
            {
                ColoredEdge &edge = edges[*hep];
                if (edge.points.front() != end_point && edge.points.back() != end_point)
                    continue;
                // found the next connecting piece
                bool lefty;
                // since there are now only two, we can compute the missing color
                if (edge.widdershins == 0)
                    lefty = !(best.clusters & (1 << (edge.clockwise - 1)));
                else 
                    lefty = best.clusters & (1 << (edge.widdershins - 1));
                subset_boundary.addEdge(edge.points, !lefty, lefty);
                ColoredPoint point = edge.points.front();
                if (point.i == 0 || point.i == N || point.j == 0 || point.j == N)
                    subset_boundary.addVertex(point);
                else
                {
                    auto position = std::find(interior_vertex.begin(), interior_vertex.end(), point);
                    if (position != interior_vertex.end())
                        interior_vertex.erase(position);
                    else
                        interior_vertex.push_back(point);
                }
                point = edge.points.back();
                if (point.i == 0 || point.i == N || point.j == 0 || point.j == N)
                    subset_boundary.addVertex(point);
                else
                {
                    auto position = std::find(interior_vertex.begin(), interior_vertex.end(), point);
                    if (position != interior_vertex.end())
                        interior_vertex.erase(position);
                    else
                        interior_vertex.push_back(point);
                }
                // don't use this one again
                half_edges.erase(hep);
                making_progress = true;
                break;
            }
            assert(making_progress);
        }
        subset_boundary.setColorful(2);
        assert(!subset_boundary.empty());

        ColoredEdge separatrix = subset_boundary.getEdges().front();
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
        // and the check is in the mail
        candidate->outcome = Status::EPP_success;
        candidate->score = best.score;
        candidate->edge_weight = best.edge_weight;
        candidate->balance_factor = best.balance_factor;
    }

    template <class ClientSample>
    void PursueProjection<ClientSample>::serial() noexcept
    {
        this->request->projections++;
        this->request->passes += candidate->pass;
        this->request->clusters += candidate->clusters;
        this->request->graphs += candidate->graphs;
        this->request->merges += candidate->merges;

        // keep the finalists in order, even failures get inserted so we return some error message
        size_t i = this->request->candidates.size();
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

        if (this->request->success())
            this->request->status = EPP_success;
    }

    template <class ClientSample>
    void QualifyMeasurement<ClientSample>::parallel() noexcept
    {
        thread_local struct scratch
        {
            float *data;
            unsigned long int size;

            ~scratch () { delete[] data; }
        } scratch = {nullptr, 0};
        if (scratch.size < this->sample.events + 1)
        {
            delete[] scratch.data;
            scratch.size = this->sample.events + 1;
            scratch.data = new float[scratch.size];
        }

        // get statistics for this measurement for this subset
        float *x = scratch.data;
        float *p = x;
        double Sx = 0, Sxx = 0;
        Event n = 0, m = 0;
        for (Event event = 0; event < this->sample.events; event++)
            if (this->subset->contains(event))
            {
                float value = (float)this->sample(event, X);
                ++n;
                Sx += value;
                Sxx += value * value;
                *p++ = value;
            }
        std::sort(x, x + n);
        x[n] = 1;
        while (x[m] == 0)
            ++m; // for CyToF/exponential we censor true zeros
        if (m == n)
            return;
        const double mu = Sx / (double)n;
        const double sigma = sqrt((Sxx - Sx * mu) / (double)(n - 1));
        const double lambda = Sx / (double)(n - m);

        // compute Kullback-Leibler Divergence
        if (sigma > 0)
        {
            // normalization factors for truncated distribution
            double NQn = .5 * (erf((x[n] - mu) / sigma / sqrt2) - erf((x[0] - mu) / sigma / sqrt2));
            double NQe = exp(-x[0] / lambda) - exp(-x[n] / lambda);
            for (Event i = 0, j; i < n; i = j)
            {
                j = i + 1;
                while ((x[j] - x[i]) < .001 && j < n)
                    j++;

                double P = (double)(j - i) / (double)n;
                double Q = .5 * (erf((x[j] - mu) / sigma / sqrt2) - erf((x[i] - mu) / sigma / sqrt2)) / NQn;
                if (Q > 0)                  // catch underflow that causes infinite result
                    KLDn += P * log(P / Q); // I didn't think it was possible either

                if (i == 0 && m > 0)
                    P = (double)(j - m) / (double)(n - m);
                else
                    P = (double)(j - i) / (double)(n - m);
                Q = (exp(-x[i] / mu) - exp(-x[j] / mu)) / NQe;
                if (Q > 0)
                    KLDe += P * log(P / Q);
            }
        }
    }

    template <class ClientSample>
    void QualifyMeasurement<ClientSample>::serial() noexcept
    {
        qualified = this->KLDn > this->parameters.kld.Normal1D &&
                    this->KLDe > this->parameters.kld.Exponential1D;
        if (qualified)
        {
            // start pursuit on this measurement vs all the others found so far
            for (Measurement Y : this->request->qualified)
                Worker<ClientSample>::enqueue(
                    new PursueProjection<ClientSample>(this->request, X < Y ? X : Y, X < Y ? Y : X));

            this->request->qualified.push_back(X);
        }
        else
        {
            // remember the best two of the unqualified
            if (this->KLDn > this->request->fallback.X_KLD)
            {
                this->request->fallback.Y_KLD = this->request->fallback.X_KLD;
                this->request->fallback.Y = this->request->fallback.X;
                this->request->fallback.X_KLD = this->KLDn;
                this->request->fallback.X = this->X;
            }
            else if (this->KLDn > this->request->fallback.Y_KLD)
            {
                this->request->fallback.Y_KLD = this->KLDn;
                this->request->fallback.Y = this->X;
            }
        }

        if (--this->request->qualifying == 0)
        {
            // in case nothing else has been tried
            Measurement X, Y;
            switch (this->request->qualified.size())
            {
            case 0:
                X = this->request->fallback.X;
                Y = this->request->fallback.Y;
                break;
            case 1:
                X = this->request->qualified.front();
                Y = this->request->fallback.X;
                break;
            default:
                return;
            }
            Worker<ClientSample>::enqueue(
                new PursueProjection<ClientSample>(this->request, X < Y ? X : Y, X < Y ? Y : X));
        }
    }

    template <class ClientSample>
    void PursueProjection<ClientSample>::start(
        Request<ClientSample> *request) noexcept
    {
        SampleSubset<ClientSample> *subset = request->subset;
        for (Measurement measurement = 0; measurement < subset->sample.measurements; ++measurement)
            if (!request->analysis->censor(measurement))
                Worker<ClientSample>::enqueue(
                    new QualifyMeasurement<ClientSample>(request, measurement));
    }
}
