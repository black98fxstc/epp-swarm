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

        // start some worker threads
        std::thread workers[10];
        for (int i = 0; i < 10; i++)
            workers[i] = std::thread([]() {
                EPP::Worker worker;
            });

        float *data = NULL;
        long data_size = 0;
        std::string last_data;
        while (!EPP::kiss_of_death)
        {
            json request;
            request["action"] = "Give me work";
            request["argument"] = "something to work on";

            json response = client.ajax(argv[1], request);

            std::cout << request.dump(4) << std::endl;
            std::cout << response.dump(4) << std::endl;

            // fake it for now
            response.clear();
            json smp;
            smp["measurments"] = 10;
            smp["events"] = 100;
            smp["key"] = "whatever";
            response["sample"] = smp;
            json sub;
            sub["key"] = "blah blah";
            response["subset"] = sub;
            //

            int measurments = response["sample"]["measurments"];
            long events = response["sample"]["events"];
            std::string sample_key = response["sample"]["key"];
            std::string subset_key = response["subset"]["key"];

            if (measurments * events < data_size)
            {
                delete[] data;
                data = new float[measurments * events];
                data_size = measurments * events;
            }

            EPP::DefaultSample<float> sample(measurments, events, data, sample_key);
            if (sample_key != last_data)
            {
                client.fetch(sample);
                last_data = sample_key;
            }

            EPP::Subset start(sample, subset_key);
            client.fetch(start);

            // start parallel projection pursuit
            {
                EPP::worker_sample constants{measurments, events, (const float *const)data, start};
                EPP::qualified_measurments.clear();
                std::unique_lock<std::recursive_mutex> lock(EPP::mutex);
                for (int measurment = 0; measurment < constants.measurments; ++measurment)
                    EPP::work_list.push(new EPP::QualifyMeasurment(constants, measurment));
                EPP::work_available.notify_all();
            };

            // wait for everything to finish
            {
                std::unique_lock<std::recursive_mutex> lock(EPP::mutex);
                while (EPP::work_outstanding)
                    EPP::work_completed.wait(lock);
            }

            // presumably report back to the dispatcher

            // only go around once for now
            EPP::kiss_of_death = true;
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