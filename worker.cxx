#include <work.h>
#include <modal.h>

#include <string>
#include <iostream>
#include <exception>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>

namespace EPP
{
    volatile static bool kiss_of_death = false;

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