#include <cstring>
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

namespace EPP
{
    std::recursive_mutex mutex;
    std::condition_variable_any work_available;
    std::condition_variable_any work_completed;
    int work_outstanding = 0;
    volatile bool kiss_of_death = false;

    // testing stuff
    std::default_random_engine generator;
    std::binomial_distribution<int> binomial(2000, 0.5);
    std::binomial_distribution<int> quality(500, 0.5);
    std::binomial_distribution<int> coin_toss(1, 0.5);

    // thread local storage
    struct work_kit
    {
        float *weights;
        float *densities;
    };

    // this data is read only so safely shared by threads
    struct worker_sample
    {
        const int measurments;
        const long events;
        const float *const data;
        const std::vector<bool> subset;
    };

    class Work
    {
    public:
        const struct worker_sample sample;

        virtual void parallel(work_kit &kit)
        {
            throw std::runtime_error("unimplemented");
        };

        virtual void serial(work_kit &kit)
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
                    work->parallel(kit);
                    lock.lock();
                    work->serial(kit);
                    delete work;
                };
        };

    private:
        struct work_kit kit;
    };

    class PursueProjection : public Work
    {
    public:
        const int X, Y;

        // many threads run this
        virtual void parallel(work_kit &kit)
        {
            // kit is thread local storage safe without syncronization
            // persue X vs Y
            std::this_thread::sleep_for(std::chrono::milliseconds(binomial(generator)));
        }
        // when thats done only one thread runs this
        virtual void serial(work_kit &kit)
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

    std::vector<int> qualified_dimensions;

    class QualifyDimension : public Work
    {
    public:
        const int X;
        double KLDn = 0;
        double KLDe = 0;
        bool qualified = false;

        QualifyDimension(
            const worker_sample sample,
            const int X)
            : Work(sample), X(X){};

        virtual void parallel(work_kit &kit)
        {
            double sum = 0, sum2 = 0;
            long n = 0;
            float x[sample.events + 1];
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

            // Kulbach-Leibler Divergence
            std::sort(x, x + sample.events);
            x[sample.events] = 1;
            if (sigma > 0)
            {
                long i = 0, j = 0;
                for (; i < sample.events; i = j)
                {
                    j = i + 1;
                    while (x[j] == x[i] && j < sample.events)
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

            std::this_thread::sleep_for(std::chrono::milliseconds(EPP::quality(EPP::generator)));
        };

        virtual void serial(work_kit &kit)
        {
            if (EPP::coin_toss(EPP::generator))
                qualified = true;

            if (qualified)
            {
                for (int Y : qualified_dimensions)
                {
                    EPP::PursueProjection *pursue = new EPP::PursueProjection(sample, X, Y);
                    EPP::work_list.push(pursue);
                };
                EPP::work_available.notify_all();

                qualified_dimensions.push_back(X);
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

        // start parallel projection pursuit
        {
            EPP::worker_sample constants{two.measurments, two.events, (const float *const)small, first};
            EPP::qualified_dimensions.clear();
            std::unique_lock<std::recursive_mutex> lock(EPP::mutex);
            for (int measurment = 0; measurment < constants.measurments; ++measurment)
            {
                EPP::QualifyDimension *qualify = new EPP::QualifyDimension(constants, measurment);
                EPP::work_list.push(qualify);
            };
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
