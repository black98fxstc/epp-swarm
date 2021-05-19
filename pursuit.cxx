#include <work.h>
#include <modal.h>
#include <fftw3.h>

#include <stack>

namespace EPP
{
    // pursue a particular X, Y pair
    void PursueProjection::parallel()
    {
        // compute the weights from the data for this subset
        long n = 0;
        weights.zero();
        const double divisor = 1.0 / N;
        for (long event = 0; event < sample.events; event++)
            if (sample.subset[event])
            {
                ++n;
                double x = sample.data[event * sample.measurments + X];
                double y = sample.data[event * sample.measurments + Y];
                int i, j;
                double dx = remquo(x, divisor, &i);
                double dy = remquo(y, divisor, &j);
                weights[i + (N + 1) * j] += (1 - dx) * (1 - dy);
                weights[i + 1 + (N + 1) * j] += dx * (1 - dy);
                weights[i + (N + 1) * j + (N + 1)] += (1 - dx) * dy;
                weights[i + 1 + (N + 1) * j + (N + 1)] += dx * dy;
            }

        // discrete cosine transform (FFT of real even function)
        transform.forward(weights, cosine);

        int clusters = 0;
        thread_local ModalClustering modal;
        do
        {
            // apply kernel to cosine transform
            // each application reduces the bandwidth further,
            // i.e., increases smoothing
            kernel.apply(cosine);

            // inverse discrete cosine transform
            // gives a smoothed density estimator
            transform.reverse(cosine, density);

            // modal clustering
            clusters = modal.findClusters(*density);
        } while (clusters > 10);

        thread_local ClusterBoundary cluster_bounds;
        modal.getBoundary(*density, cluster_bounds);

        // compute the cluster weights
        auto cluster_map = cluster_bounds.getMap();
        long cluster_weight[clusters + 1];
        for (long event = 0; event < sample.events; event++)
            if (sample.subset[event])
            {
                double x = sample.data[event * sample.measurments + X];
                double y = sample.data[event * sample.measurments + Y];
                short cluster = cluster_map->colorAt(x, y);
                ++cluster_weight[cluster];
            };

        // get the edges, which have their own weights
        auto edges = cluster_bounds.getEdges();

        // get the dual graph of the map
        auto graph = cluster_bounds.getDualGraph();

        // pile of graphs to consider
        std::stack<DualGraph> pile;
        pile.push(*graph);

        double best_score = std::numeric_limits<double>::infinity();
        booleans best_edges;
        booleans best_clusters;

        // find and score simple sub graphs
        while (!pile.empty())
        {
            DualGraph graph = pile.top();
            pile.pop();
            if (graph.isSimple())
            { // one edge, i.e., two populations
                booleans left_clusters = graph.left();
                double cluster_weight = 0;
                for (int i = 1; i <= clusters; i++)
                {
                    if (left_clusters & (1 << i))
                        cluster_weight += weights[i];
                }
                booleans dual_edges = graph.edge();
                double edge_weight = 0;
                for (int i = 0; i <= edges.size(); i++)
                {
                    if (dual_edges & (1 << i))
                        edge_weight += edges[i].weight;
                }
                double P = (double)cluster_weight / (double)n;
                double balanced_weight = 4 * P * (1 - P) * edge_weight;

                // score this separatrix
                if (balanced_weight < best_score)
                {
                    best_edges = dual_edges;
                    best_clusters = left_clusters;
                }
            }
            else
            { // not simple so simplify it some, i.e., remove one dual edge at a time
                // and merge two adjacent subsets. that makes a bunch more graphs to look at
                std::vector<DualGraph> simplified = graph.simplify();
                for (auto graph : simplified)
                    pile.push(graph);
            }
        }

        thread_local ColoredBoundary<short, bool> subset_boundary;
        subset_boundary.clear();
        for (int i = 0; i <= edges.size(); i++)
        {
            if (best_edges & (1 << i))
            {
                bool lefty = best_clusters & (1 << edges[i].clockwise);
                subset_boundary.addEdge(edges[i].points, lefty, !lefty, 0);
            }
        }
        subset_boundary.setColorful(2);

        // create in/out subsets
        auto subset_map = subset_boundary.getMap();
        for (long event = 0; event < sample.events; event++)
            if (sample.subset[event])
            {
                double x = sample.data[event * sample.measurments + X];
                double y = sample.data[event * sample.measurments + Y];
                short member = subset_map->colorAt(x, y);
                if (member)
                    in[event] = true;
                else
                    out[event] = true;
            };

        separatrix = subset_boundary.getEdges();

        // separatrix, in and out are the payload

        std::this_thread::sleep_for(std::chrono::milliseconds(binomial(generator)));
    }

    void PursueProjection::serial()
    {
        // see if this the best yet found
        // if no more to try produce result

        // make a boundry from the separatrix

        std::cout << "pursuit completed " << X << " vs " << Y << std::endl;
    };

    PursueProjection::Transform PursueProjection::transform;
    PursueProjection::Kernel PursueProjection::kernel;
    thread_local PursueProjection::FFTData PursueProjection::weights;
    thread_local PursueProjection::FFTData PursueProjection::cosine;
    thread_local PursueProjection::FFTData PursueProjection::density;

    void QualifyMeasurment::parallel()
    {
        // get statistics for this measurment for this subset
        double sum = 0, sum2 = 0;
        long n = 0;
        float *x = scratch.reserve(sample.events);
        float *p = x;
        for (long event = 0; event < sample.events; event++)
            if (sample.subset[event])
            {
                double value = sample.data[event * sample.measurments + X];
                ++n;
                sum += value;
                sum2 += value * value;
                *p++ = value;
            }
        const double mu = sum / n;
        const double sigma = sqrt((sum2 - sum * sum / n) / (n - 1));

        // compute Kulbach-Leibler Divergence
        std::sort(x, x + n);
        x[n] = 1;
        if (sigma > 0)
        {
            for (long i = 0, j; i < n; i = j)
            {
                j = i + 1;
                while ((x[j] - x[i]) < .001 && j < n)
                    j++;
                double p = (double)(j - i) / (double)n;
                double Q = x[j] - x[i];
                double Pn = erf((x[j] - mu) / sigma) - erf((x[i] - mu) / sigma);
                double Pe = exp(-x[i] / mu) - exp(-x[j] / mu);
                KLDn += p * log(Pn / Q);
                KLDe += p * log(Pe / Q);
            }
            // I need to look up the formulas again but it's something like this
        }

        std::this_thread::sleep_for(std::chrono::milliseconds(quality(generator)));
    };

    void QualifyMeasurment::serial()
    {
        if (coin_toss(generator))
            qualified = true;

        if (qualified)
        {
            // start pursuit on this measurement vs all the others found so far
            for (int Y : qualified_measurments)
                EPP::work_list.push(new EPP::PursueProjection(sample, X, Y));
            EPP::work_available.notify_all();

            qualified_measurments
                .push_back(X);
            std::cout << "dimension qualified " << X << std::endl;
        }
        else
            std::cout << "dimension disqualified " << X << std::endl;
    };

    thread_local QualifyMeasurment::Scratch QualifyMeasurment::scratch;

    PursueProjection::FFTData::~FFTData()
    {
        if (data)
            fftwf_free(data);
    }

    float *PursueProjection::FFTData::operator*()
    {
        if (!data)
            data = (float *)fftw_malloc(sizeof(float) * (N + 1) * (N + 1));
        return data;
    }

    PursueProjection::Transform::Transform()
    {
        // FFTW planning is slow and not thread safe so we do it here
        if (fftw_import_system_wisdom())
        {
            DCT = (void *)fftwf_plan_r2r_2d((N + 1), (N + 1), *weights, *cosine,
                                    FFTW_REDFT00, FFTW_REDFT00, 0);
            //  FFTW_WISDOM_ONLY);
            // actually they are the same in this case but leave it for now
            IDCT = (void *)fftwf_plan_r2r_2d((N + 1), (N + 1), *cosine, *density,
                                     FFTW_REDFT00, FFTW_REDFT00, 0);
            if (!DCT || !IDCT)
                throw std::runtime_error("can't initialize FFTW");
        }
        else
            throw std::runtime_error("can't initialize FFTW");
    };

    PursueProjection::Transform::~Transform()
    {
        fftwf_destroy_plan((fftwf_plan)DCT);
        fftwf_destroy_plan((fftwf_plan)IDCT);
    }

    void PursueProjection::Transform::forward(FFTData &in, FFTData &out)
    {
        fftwf_execute_r2r((fftwf_plan)DCT, *in, *out);
    }

    void PursueProjection::Transform::reverse(FFTData &in, FFTData &out)
    {
        fftwf_execute_r2r((fftwf_plan)IDCT, *in, *out);
    }
}