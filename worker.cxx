#include <cstring>
#include <cmath>
#include <string>
#include <iostream>
#include <exception>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <random>
#include <chrono>
#include <algorithm>
#include <stack>

#include <client.h>
#include <boundary.h>
#include <modal.h>

// testing stuff
std::default_random_engine generator;
std::binomial_distribution<int> binomial(2000, 0.5);
std::binomial_distribution<int> quality(500, 0.5);
std::binomial_distribution<int> coin_toss(1, 0.5);

namespace EPP
{
    std::recursive_mutex mutex;
    std::condition_variable_any work_available;
    std::condition_variable_any work_completed;
    int work_outstanding = 0;
    volatile bool kiss_of_death = false;

    // these are essential constants that are read only
    // so safely shared by all threads
    struct worker_sample
    {
        const int measurments;
        const long events;
        const float *const data;
        const std::vector<bool> subset;
    };

    // abstract class representing a unit of work to be done
    // virtual functions let subclasses specialize tasks
    // handles work_completed and work_ounstanding
    class Work
    {
    public:
        const struct worker_sample sample;

        // many threads can execute this in parallel
        virtual void parallel()
        {
            throw std::runtime_error("unimplemented");
        };
        // then only one thread at a time can run
        virtual void serial()
        {
            throw std::runtime_error("unimplemented");
        };

        ~Work()
        {
            std::unique_lock<std::recursive_mutex> lock(EPP::mutex);
            if (--work_outstanding == 0)
                work_completed.notify_all();
        };

    protected:
        Work(const worker_sample &sample) : sample(sample)
        {
            std::unique_lock<std::recursive_mutex> lock(EPP::mutex);
            ++work_outstanding;
        };
    };

    std::queue<Work *> work_list;

    // a generic worker thread. looks for work, does it, deletes it
    // virtual functions in the work object do all the real work
    class Worker
    {
    public:
        Worker()
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
                };
        };
    };

    // pursue a particular X, Y pair
    class PursueProjection : public Work
    {
        // fftw needs sepcial alignment to take advantage of vector instructions
        class FFTData
        {
            float *data;

        public:
            FFTData()
            {
                data = NULL;
            }

            ~FFTData()
            {
                if (data)
                    fftwf_free(data);
            }

            inline float *operator*()
            {
                if (!data)
                    data = (float *)fftw_malloc(sizeof(float) * (N + 1) * (N + 1));
                return data;
            }

            inline float &operator[](const int i)
            {
                return data[i];
            }

            inline void zero()
            {
                if (data)
                    std::fill(data, data + (N + 1) * (N + 1), 0);
            }
        };

        static thread_local FFTData weights;
        static thread_local FFTData cosine;
        static thread_local FFTData density;

        class Transform
        {
            fftwf_plan DCT;
            fftwf_plan IDCT;

        public:
            Transform()
            {
                // FFTW planning is slow and not thread safe so we do it here
                if (fftw_import_system_wisdom())
                {
                    DCT = fftwf_plan_r2r_2d((N + 1), (N + 1), *weights, *cosine,
                                            FFTW_REDFT00, FFTW_REDFT00, 0);
                    //  FFTW_WISDOM_ONLY);
                    // actually they are the same in this case but leave it for now
                    IDCT = fftwf_plan_r2r_2d((N + 1), (N + 1), *cosine, *density,
                                             FFTW_REDFT00, FFTW_REDFT00, 0);
                    if (!DCT || !IDCT)
                        throw std::runtime_error("can't initialize FFTW");
                }
                else
                    throw std::runtime_error("can't initialize FFTW");
            };

            ~Transform()
            {
                fftwf_destroy_plan(DCT);
                fftwf_destroy_plan(IDCT);
            }

            void forward(FFTData &in, FFTData &out)
            {
                fftwf_execute_r2r(DCT, *in, *out);
            }

            void reverse(FFTData &in, FFTData &out)
            {
                fftwf_execute_r2r(IDCT, *in, *out);
            }
        };

        static Transform transform;

        class Kernel
        {
            double k[N + 1];

        public:
            Kernel()
            {
                const double bw = N / 2;
                for (int i = 0; i <= N; i++)
                    k[i] = exp(-(i + .5) * (i - .5) / bw / bw);
            }

            void apply(FFTData &data)
            {
                for (int i = 0; i <= N; i++)
                {
                    for (int j = 0; j < i; j++)
                    {
                        data[i + (N + 1) * j] *= k[i] * k[j];
                        data[j + (N + 1) * i] *= k[j] * k[i];
                    }
                    data[i + (N + 1) * i] *= k[i] * k[i];
                }
            }
        };

        static Kernel kernel;

    public:
        const int X, Y;

        std::vector<ColoredEdge<short, bool>> separatrix;
        std::vector<bool> in;
        std::vector<bool> out;

        PursueProjection(
            const worker_sample sample,
            const int X,
            const int Y)
            : Work(sample), X(X), Y(Y){};

        ~PursueProjection(){};

        virtual void parallel()
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

                clusters = 5;
            } while (clusters > 10);

            thread_local ClusterBoundary cluster_bounds;
            modal.getBoundary(*density, cluster_bounds);

            // compute the cluster weights
            ClusterMap *cluster_map = cluster_bounds.getMap();
            short cluster_weight[clusters + 1];
            for (long event = 0; event < sample.events; event++)
                if (sample.subset[event])
                {
                    double x = sample.data[event * sample.measurments + X];
                    double y = sample.data[event * sample.measurments + Y];
                    short cluster = cluster_map->colorAt(x, y);
                    ++cluster_weight[cluster];
                };
            delete cluster_map;

            // get the edges, which have their own weights
            auto edges = cluster_bounds.getEdges();

            // get the dual graph of the map
            DualGraph *graph = cluster_bounds.getDualGraph();

            // pile of graphs to consider
            std::stack<DualGraph> pile;
            booleans best_edges;
            booleans left_clusters;
            pile.push(*graph);

            // find and score simple sub graphs
            while (!pile.empty())
            {
                DualGraph graph = pile.top();
                pile.pop();
                if (graph.isSimple())
                { // one edge, i.e., two populations
                    booleans clusters = graph.left();
                    double cluster_weight = 0;
                    for (int i = 1; i <= clusters; i++)
                    {
                        if (clusters & (1 << i))
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
                    best_edges = dual_edges;
                    left_clusters = clusters;
                }
                else
                { // not simple so simplify it some, i.e., remove one dual edge at a time
                    // and merge two adjacent subsets. that makes a bunch more graphs to look at
                    std::vector<DualGraph> simplified = graph.simplify();
                    for (auto graph : simplified)
                        pile.push(graph);
                }
            }
            delete graph;

            thread_local ColoredBoundary<short, bool> subset_boundary;
            subset_boundary.clear();
            for (int i = 0; i <= edges.size(); i++)
            {
                if (best_edges & (1 << i))
                {
                    bool lefty = left_clusters & (1 << edges[i].clockwise);
                    subset_boundary.addEdge(edges[i].points, lefty, !lefty, 0);
                }
            }

            // create in/out subsets
            ColoredMap<short, bool> *subset_map = subset_boundary.getMap();
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
            delete subset_map;

            separatrix = subset_boundary.getEdges();

            // separatrix, in and out are the payload

            std::this_thread::sleep_for(std::chrono::milliseconds(binomial(generator)));
        }

        virtual void
        serial()
        {
            // see if this the best yet found
            // if no more to try produce result

            // make a boundry from the separatrix

            std::cout << "pursuit completed " << X << " vs " << Y << std::endl;
        };
    };

    PursueProjection::Transform PursueProjection::transform;
    PursueProjection::Kernel PursueProjection::kernel;
    thread_local PursueProjection::FFTData PursueProjection::weights;
    thread_local PursueProjection::FFTData PursueProjection::cosine;
    thread_local PursueProjection::FFTData PursueProjection::density;

    std::vector<int> qualified_measurments;

    class QualifyMeasurment : public Work
    {
        class Scratch
        {
            float *data;
            int size;

        public:
            Scratch()
            {
                data = NULL;
                size = 0;
            }

            ~Scratch()
            {
                if (data)
                    delete[] data;
            }

            float *&reserve(int size)
            {
                if (this->size < size)
                {
                    delete[] data;
                    data = NULL;
                }
                if (!data)
                {
                    data = new float[size];
                    this->size = size;
                }
                return data;
            }

            inline float &operator[](const int i)
            {
                return data[i];
            }
        };

        static thread_local Scratch scratch;

    public:
        const int X;
        double KLDn = 0;
        double KLDe = 0;
        bool qualified = false;

        QualifyMeasurment(
            const worker_sample sample,
            const int X)
            : Work(sample), X(X){};

        virtual void parallel()
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

        virtual void serial()
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
    };

    thread_local QualifyMeasurment::Scratch QualifyMeasurment::scratch;
};

using json = nlohmann::json;

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " endpoint\n";
        return 1;
    }

    // AWS has to be initted *before* our constructors can run
    EPP::Init();
    try
    {
        EPP::Client client;

        // trivial ajax transaction

        json request;
        request["action"] = "Do something!";
        request["argument"] = "something to work on";

        json response = client.ajax(argv[1], request);

        std::cout << request.dump(4) << std::endl;
        std::cout << response.dump(4) << std::endl;

        // stage double data in memory to the server as floats

        double data[100][10];
        EPP::DefaultSample<double> one(10, 100, (double *)&data);
        client.stage(one);

        // remote workers fetch data as floats by the hash passed via JSON

        float small[100][10];
        EPP::DefaultSample<float> two(one.measurments, one.events, (float *)&small, one.get_key());
        client.fetch(two);

        // other data models

        double transpose[10][100];
        EPP::TransposeSample<double> three(10, 100, (double *)&transpose);
        client.stage(three);

        float **pointers;
        pointers = new float *[10];
        for (int i = 0; i < 10; i++)
            pointers[i] = new float[100];
        EPP::PointerSample<float> four(10, 100, pointers);
        client.stage(four);

        // subsets of samples are std::vector<bool>

        EPP::Subset first(one);
        first[1] = true;
        client.stage(first);

        EPP::Subset second(*first.sample, first.get_key());
        client.fetch(second);
        bool is_in_subset = second[1];

        // start some worker threads

        std::thread workers[10];
        for (int i = 0; i < 10; i++)
            workers[i] = std::thread([]() {
                EPP::Worker worker;
            });

        // start parallel projection pursuit
        {
            EPP::worker_sample constants{two.measurments, two.events, (const float *const)small, first};
            EPP::qualified_measurments.clear();
            std::unique_lock<std::recursive_mutex> lock(EPP::mutex);
            for (int measurment = 0; measurment < constants.measurments; ++measurment)
                EPP::work_list.push(new EPP::QualifyMeasurment(constants, measurment));
            EPP::work_available.notify_all();
        }
        // wait for everything to finish
        {
            std::unique_lock<std::recursive_mutex> lock(EPP::mutex);
            while (EPP::work_outstanding)
                EPP::work_completed.wait(lock);
        }

        // tell the workers to exit and wait for them to shut down
        EPP::kiss_of_death = true;
        {
            std::unique_lock<std::recursive_mutex> lock(EPP::mutex);
            EPP::work_available.notify_all();
        }
        for (int i = 0; i < 10; i++)
            workers[i].join();
    }
    catch (std::runtime_error e)
    {
        std::cout << e.what() << std::endl;
        return 1;
    }
    EPP::Finish();

    return 0;
};