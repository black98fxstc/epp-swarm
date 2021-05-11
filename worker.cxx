#include <cstring>
#include <string>
#include <iostream>
#include <exception>
#include <client.h>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>

namespace EPP
{
    std::mutex avail, ready;
    std::condition_variable work_available;
    std::condition_variable work_completed;
    int work_outstanding = 0;
    volatile bool kiss_of_death = false;

    // thread local storage

    struct work_kit
    {
        float *weights;
        float *densities;
    };

    class Work
    {
    public:
        virtual void parallel(work_kit &kit)
        {
            throw std::runtime_error("unimplemented");
        };

        virtual void serial(work_kit &kit)
        {
            throw std::runtime_error("unimplemented");
        };

        Work()
        {
            std::unique_lock<std::mutex> lock(EPP::ready);
            ++work_outstanding;
        };
        ~Work()
        {
            std::unique_lock<std::mutex> lock(EPP::ready);
            if (--work_outstanding == 0)
                work_completed.notify_all();
        };
    };

    std::queue<Work> work_list;

    class Worker
    {
    public:
        Worker()
        {
            while (true)
            {
                std::unique_lock<std::mutex> lock(avail);
                while (work_list.empty())
                    work_available.wait(lock);
                if (kiss_of_death)
                    return;
                Work work = work_list.front();
                work_list.pop();
                lock.unlock();
                work.parallel(kit);
                lock.lock();
                work.serial(kit);
                std::cout << "work unit completed" << std::endl;
            };
        };

    private:
        struct work_kit kit;
    };

    class PursueProjection : public Work
    {
    public:
        // this data is read only so safely shared by threads
        const int X, Y;
        const float *data;
        const std::vector<bool> subset;

        // many threads run this
        virtual void parallel(work_kit &kit)
        {
            // kit is thread local storage safe without syncronization
            // persue X vs Y
        }
        // when thats done only one thread runs this
        virtual void serial(work_kit &kit){
            // see if this the best yet found
            // if no more to try produce result
        };

        PursueProjection(const int X, const int Y) : X(X), Y(Y){};
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

        // start parallel projection pursuits
        {
            std::unique_lock<std::mutex> lock(EPP::avail);
            for (int x = 0; x < 25; x++)
                for (int y = 0; y < x; y++)
                {
                    EPP::PursueProjection pursue(x, y);
                    EPP::work_list.push(pursue);
                };
            EPP::work_available.notify_all();
        }
        // wait for everything to finish
        {
            std::unique_lock<std::mutex> lock(EPP::ready);
            while (EPP::work_outstanding)
                EPP::work_completed.wait(lock);
        }

        // tell the workers to exit and wait for them to shut down

        EPP::kiss_of_death == true;
        {
            std::unique_lock<std::mutex> lock(EPP::avail);
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
}