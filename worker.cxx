#include <cstring>
#include <cmath>
#include <string>
#include <iostream>
#include <exception>
#include <client.h>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <random>
#include <chrono>
#include <algorithm>

// testing stuff
std::default_random_engine generator;
std::binomial_distribution<int> binomial(2000, 0.5);
std::binomial_distribution<int> quality(500, 0.5);
std::binomial_distribution<int> coin_toss(1, 0.5);

namespace EPP
{
    // resolution of the density estimator, best if it's a power of 2 for FFTW
    const int N = 1 << 9;

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

    // thread local storage, can be safely accessed only on
    // it's own thread, but saved to avoid storage allocation
    class worker_kit
    {
    private:
        int current_measurments;
        long current_events;
        float *scratch_area = NULL;
        float *weights_array = NULL;
        float *cosine_transform = NULL;
        float *densities_array = NULL;

        void check_size(const worker_sample &sample)
        {
            if (current_events < sample.events)
            {
                delete[] scratch_area;
                current_measurments = sample.measurments;
                current_events = sample.events;
            }
        }

    public:
        float *scratch(const worker_sample &sample)
        {
            check_size(sample);
            if (!scratch_area)
                scratch_area = new float[sample.events + 1];
            return scratch_area;
        };

        float *weights(const worker_sample &sample)
        {
            check_size(sample);
            if (!weights_array)
                weights_array = (float *)fftwf_malloc(sizeof(float) * N * N);
            return weights_array;
        };

        float *cosine(const worker_sample &sample)
        {
            check_size(sample);
            if (!cosine_transform)
                cosine_transform = (float *)fftwf_malloc(sizeof(float) * N * N);
            return cosine_transform;
        };

        float *density(const worker_sample &sample)
        {
            check_size(sample);
            if (!densities_array)
                densities_array = (float *)fftwf_malloc(sizeof(float) * N * N);
            return densities_array;
        };
    };

    // abstract class representing a unit of work to be done
    // virtual functions let subclasses specialize tasks
    // handles work_completed and work_ounstanding
    class Work
    {
    public:
        const struct worker_sample sample;

        // many threads can execute this in parallel
        virtual void parallel(worker_kit &kit)
        {
            throw std::runtime_error("unimplemented");
        };
        // then only one thread at a time can run
        virtual void serial(worker_kit &kit)
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
        // kit is thread local storage safe without syncronization
        worker_kit kit;

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
                    work->parallel(kit);
                    lock.lock();
                    work->serial(kit);
                    delete work;
                };
        };
    };

    // pursue a particular X, Y pair
    class PursueProjection : public Work
    {
        int _cluster;
        int &cluster(int &i, int &j)
        {
            return _cluster;
        };
        bool _contiguous;
        bool &contiguous(const int &i, const int &j)
        {
            return _contiguous;
        };
        struct vertex
        {
            float f;
            int i, j;
        } v[EPP::N * N], *pv = v;
        struct
        {
            bool operator()(vertex a, vertex b) const { return a.f > b.f; }
        } decreasing_density;

        void visit(
            int &result,
            int i,
            int j)
        {
            // this point is contiguous with a classified point
            contiguous(i, j) = true;

            // if this point has been assigned to a cluster
            if (cluster(i, j) > 0)
                // and our starting point has not
                if (result < 1)
                    // assign it to our cluster
                    result = cluster(i, j);
                else if (result != cluster(i, j))
                    // if we found something different it's really a boundary point
                    result = 0;
        };

    public:
        const int X, Y;

        virtual void parallel(worker_kit &kit)
        {
            // compute the weights from the data for this subset
            long n = 0;
            float *weights = kit.weights(sample);
            const double divisor = 1.0 / (N - 1.0);
            memset(weights, 0, N * N & sizeof(float));
            for (long event = 0; event < sample.events; event++)
                if (sample.subset[event])
                {
                    ++n;
                    double x = sample.data[event * sample.measurments + X];
                    double y = sample.data[event * sample.measurments + Y];
                    int i, j;
                    double dx = remquo(x, divisor, &i);
                    double dy = remquo(y, divisor, &j);
                    weights[i + N * j] += (1 - dx) * (1 - dy);
                    weights[i + 1 + N * j] += dx * (1 - dy);
                    weights[i + N * j + N] += (1 - dx) * dy;
                    weights[i + 1 + N * j + N] += dx * dy;
                }

            // discrete cosine transform (FFT of real even function)
            float *cosine = kit.cosine(sample);
            fftwf_execute_r2r(EPP::DCT, weights, cosine);

            int clusters = 0;
            do
            {
                // apply kernel to cosine transform
                // each application reduces the bandwidth further,
                // i.e., increases smoothing
                double bandwidth = 1;
                for (int i = 0; i < N; i++)
                {
                    double kernel = exp(-i * i * 2 / bandwidth);
                    cosine[i + N * i] *= kernel;
                    for (int j = 0; j < i; j++)
                    {
                        kernel = exp(-(i * i + j * j) / bandwidth);
                        cosine[i + N * j] *= kernel;
                        cosine[j + N * i] *= kernel;
                    } // missing some constant factors I don't remember
                }

                // inverse discrete cosine transform
                // gives a smoothed density estimator
                float *density = kit.density(sample);
                fftwf_execute_r2r(EPP::IDCT, cosine, density);

                // modal clustering

                // collect all the grid points
                pv = v;
                for (int i = 0; i < N; i++)
                    for (int j = 0; n < N; j++)
                    {
                        pv->f = density[i + N * j];
                        pv->i = i;
                        pv->j = j;
                        pv++;
                    }
                // get all comparisons out of the way early and efficiently
                std::sort(v, v + N * N, decreasing_density);
                // this should really be the significance threshold but this won't deadlock for now
                float threshold = v[N * N / 2].f;
                // for the points that are above threshold, i.e., cluster points
                for (pv = v; pv < v + N * N; pv++)
                {
                    if (pv->f < threshold)
                        break;
                    // visit the neighbors to see what clusters they belong to
                    int result = -1;
                    visit(result, pv->i - 1, pv->j);
                    visit(result, pv->i + 1, pv->j);
                    visit(result, pv->i, pv->j - 1);
                    visit(result, pv->i, pv->j + 1);
                    // if we didn't find one this is a new mode
                    if (result < 0)
                        cluster(pv->i, pv->j) = ++clusters;
                    else
                        cluster(pv->i, pv->j) = result;
                }
                // we don't trust these small densities so we take the rest
                // randomly so the border will grow approximately uniformly
                std::random_device rd;
                std::mt19937 g(rd());
                std::shuffle(pv, v + N * N, g);
                for (; pv < v + N * N; pv++)
                {
                    // find the next unassigned point that is contiguous with those already classified
                    vertex *pw;
                    for (pw = pv; pw < v + N * N; pw++)
                        if (contiguous(pw->i, pw->j))
                            break;
                    if (pw != pv)
                    {
                        struct vertex t;
                        t = *pv;
                        *pv = *pw;
                        *pw = t;
                    }
                    // visit the neighbors and then allocate it as a background point
                    int result = -1;
                    visit(result, pv->i - 1, pv->j);
                    visit(result, pv->i + 1, pv->j);
                    visit(result, pv->i, pv->j - 1);
                    visit(result, pv->i, pv->j + 1);
                    if (result < 0)
                        cluster(pv->i, pv->j) = 0;
                    else
                        cluster(pv->i, pv->j) = result;
                }

                clusters = 5;
            } while (clusters > 10);

            // find and score spearatrix

            std::this_thread::sleep_for(std::chrono::milliseconds(binomial(generator)));
        }

        virtual void
        serial(worker_kit &kit)
        {
            // see if this the best yet found
            // if no more to try produce result
            std::cout << "pursuit completed " << X << " vs " << Y << std::endl;
        };

        PursueProjection(
            const worker_sample sample,
            const int X,
            const int Y)
            : Work(sample), X(X), Y(Y){};
    };

    std::vector<int> qualified_measurments;

    class QualifyMeasurment : public Work
    {
    public:
        const int X;
        double KLDn = 0;
        double KLDe = 0;
        bool qualified = false;

        QualifyMeasurment(
            const worker_sample sample,
            const int X)
            : Work(sample), X(X){};

        virtual void parallel(worker_kit &kit)
        {
            // get statistics for this measurment for this subset
            double sum = 0, sum2 = 0;
            long n = 0;
            float *x = kit.scratch(sample);
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

        virtual void serial(worker_kit &kit)
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

        // FFTW planning is slow and not thread safe so we do it here
        if (fftw_import_system_wisdom())
        {
            float *in = (float *)fftw_malloc(sizeof(float) * EPP::N * EPP::N);
            float *out = (float *)fftw_malloc(sizeof(float) * EPP::N * EPP::N);
            EPP::DCT = fftwf_plan_r2r_2d(EPP::N, EPP::N, in, out,
                                         FFTW_REDFT10, FFTW_REDFT10, 0);
            //  FFTW_WISDOM_ONLY);
            EPP::IDCT = fftwf_plan_r2r_2d(EPP::N, EPP::N, in, out,
                                          FFTW_REDFT01, FFTW_REDFT01, 0);
            fftw_free(in);
            fftw_free(out);
            if (!EPP::DCT || !EPP::IDCT)
                throw std::runtime_error("can't initialize FFTW");
        }
        else
            throw std::runtime_error("can't initialize FFTW");

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
