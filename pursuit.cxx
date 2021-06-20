#include <work.h>
#include <modal.h>
#include <fftw3.h>

#include <stack>

namespace EPP
{
    // pursue a particular X, Y pair
    void PursueProjection::parallel()
    {
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
		thread_local PursueProjection::FFTData cosine;
        transform.forward(weights, cosine);

        int clusters;
        int pass = 0;
		thread_local PursueProjection::FFTData filtered;
		thread_local PursueProjection::FFTData density;
        thread_local ModalClustering modal;
        do
        {
            // apply kernel to cosine transform
            kernel.apply(cosine, filtered, ++pass);

            // inverse discrete cosine transform
            // gives a smoothed density estimator
            transform.reverse(filtered, density);

            // modal clustering
            clusters = modal.findClusters(*density);
        } while (clusters > 10);

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

        // Normalize the density P, n for weights, (2N)^2 for discrete cosine transform
//        double NP = (double)(n * 4 * N * N);
        KLD /= NP;
        // subtract off normalization constants factored out of the sum above
        KLD -= log(NP / NQ);
        // OK now what do we do with it?

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

        // get the dual graph of the map
        auto graph = cluster_bounds.getDualGraph();

        // pile of graphs to consider
        std::stack<DualGraph> pile;
        pile.push(*graph);

        // find and score simple sub graphs
		double best_score = std::numeric_limits<double>::infinity();
		booleans best_edges;
		booleans best_clusters;
        long count = 0;
        while (!pile.empty())
        {
			++count;
            DualGraph graph = pile.top();
            pile.pop();
            if (graph.isSimple())
            { // one edge, i.e., two populations
                booleans left_clusters = graph.left();
                double left_weight = 0;
                for (int i = 1; i <= clusters; i++)
                {
                    if (left_clusters & (1 << i))
                        left_weight += cluster_weight[i];
                }
                booleans dual_edges = graph.edge();
                double edge_weight = 0;
                for (int i = 0; i <= edges.size(); i++)
                {
                    if (dual_edges & (1 << i))
                        edge_weight += edges[i].weight;
                }
                double P = (double)left_weight / (double)n;
                double balanced_weight = 4 * P * (1 - P) * edge_weight;

                // score this separatrix
                if (balanced_weight < best_score)
                {
                    best_edges = dual_edges;
                    best_clusters = left_clusters;
                }
            }
            else
            { 	// not simple so simplify it some, i.e., remove one dual edge at a time
                // and merge two adjacent subsets. that makes a bunch more graphs to look at
                std::vector<DualGraph> simplified = graph.simplify();
                for (const auto& graph : simplified)
                    pile.push(graph);
            }
        }
		std::cout << count << " graphs considered" << std::endl;

        thread_local ColoredBoundary<short, bool> subset_boundary;
        subset_boundary.clear();
        for (int i = 0; i < edges.size(); i++)
        {
            if (best_edges & (1 << i))
            {
                bool lefty = best_clusters & (1 << edges[i].widdershins);
                subset_boundary.addEdge(edges[i].points, lefty, !lefty);
            }
        }
        subset_boundary.setColorful(2);

        // create in/out subsets
        in.resize(n);
        in.clear();
        out.resize(n);
        out.clear();

        auto subset_map = subset_boundary.getMap();
		count = 0;
        for (long event = 0; event < sample.events; event++)
            if (sample.subset[event])
            {
            	++count;
                double x = sample.data[event * sample.measurements + X];
                double y = sample.data[event * sample.measurements + Y];
                bool member = subset_map->colorAt(x, y);
                if (member)
                    in[event] = true;
                else
                    out[event] = true;
            }

        separatrix = subset_boundary.getEdges();

        // separatrix, in and out are the payload
    }

    void PursueProjection::serial()
    {
        // see if this the best yet found
        // if no more to try produce result

        // make a boundary from the separatrix

        std::cout << "pursuit completed " << X << " vs " << Y << std::endl;
    }

    PursueProjection::Transform PursueProjection::transform;
    PursueProjection::Kernel PursueProjection::kernel;

    void QualifyMeasurement::parallel()
    {
        // get statistics for this measurement for this subset
		thread_local QualifyMeasurement::Scratch scratch;
        float *x = scratch.reserve(sample.events);
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

    void QualifyMeasurement::serial()
    {
        qualified = KLDn > .16 && KLDe > .16;
        if (qualified)
        {
            // start pursuit on this measurement vs all the others found so far
            for (int Y : qualified_measurements)
                EPP::work_list.push(new EPP::PursueProjection(sample, X, Y));
            EPP::work_available.notify_all();

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

    float *PursueProjection::FFTData::operator*()
    {
        if (!data)
            data = (float *)fftw_malloc(sizeof(float) * (N + 1) * (N + 1));
        return data;
    }

    void PursueProjection::FFTData::zero()
    {
        if (!data)
            data = (float *)fftw_malloc(sizeof(float) * (N + 1) * (N + 1));
        std::fill(data, data + (N + 1) * (N + 1), 0);
    }

    PursueProjection::Transform::Transform()
    {
        PursueProjection::FFTData in;
        PursueProjection::FFTData out;
        // FFTW planning is slow and not thread safe so we do it here
		DCT = (void *)fftwf_plan_r2r_2d((N + 1), (N + 1), *in, *out,
										FFTW_REDFT00, FFTW_REDFT00, 0);
		// actually they are the same in this case but leave it for now
		IDCT = (void *)fftwf_plan_r2r_2d((N + 1), (N + 1), *in, *out,
										 FFTW_REDFT00, FFTW_REDFT00, 0);
		if (!DCT || !IDCT)
			throw std::runtime_error("can't initialize FFTW");
    }

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