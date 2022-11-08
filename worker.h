
/*
 * Developer: Wayne Moore <wmoore@stanford.edu> 
 * Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
 * License: BSD 3 clause
 */
namespace EPP
{
    // abstract class representing a unit of work to be done
    // virtual functions let subclasses specialize tasks

    template <class ClientSample>
    class Work
    {
    public:
        const ClientSample &sample;
        const SampleSubset<ClientSample> *subset;
        const Parameters &parameters;
        Request<ClientSample> *request;

        // many threads can execute this in parallel
        virtual void parallel() noexcept
        {
            assert(("unimplemented", false));
        };
        // then only one thread at a time can run
        virtual void serial() noexcept
        {
            assert(("unimplemented", false));
        };

        ~Work()
        {
            if (request->analysis->pursuer->decrement(request))
                request->analysis->pursuer->finish(request);
        };

    protected:
        explicit Work(
            Request<ClientSample> *request) noexcept
            : sample(request->sample), subset(request->subset), parameters(request->analysis->parameters), request(request)
        {
            request->analysis->pursuer->increment(request);
        };
    };

    // a generic worker thread. looks for work, does it, deletes it
    // virtual functions in the work object do all the real work
    template <class ClientSample>
    class Worker
    {
        friend class Work<ClientSample>;

    protected:
        static std::mutex serialize;
        static std::mutex mutex;
        static std::queue<Work<ClientSample> *> work_list;
        static std::condition_variable work_available;
        volatile static bool kiss_of_death;

        void work() noexcept
        {
            std::chrono::time_point<std::chrono::steady_clock> begin, end;

            Work<ClientSample> *work = dequeue();
            if (!work) // spurious wakeup
                return;

            begin = std::chrono::steady_clock::now();
            work->parallel();
            end = std::chrono::steady_clock::now();
            {
                std::lock_guard<std::mutex> lock(serialize);
                work->serial();
            }
            work->request->milliseconds += std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

            delete work;
        };

        Work<ClientSample> *dequeue() noexcept
        {
            std::lock_guard<std::mutex> lock(mutex);
            if (work_list.empty())
                return nullptr;
            Work<ClientSample> *work = work_list.front();
            work_list.pop();
            return work;
        }

        bool idle()
        {
            std::lock_guard<std::mutex> lock(mutex);
            return work_list.empty();
        }

        void wait()
        {
            std::unique_lock<std::mutex> lock(mutex);
            while (work_list.empty() && !kiss_of_death)
                work_available.wait(lock);
        }

    public:
        static void enqueue(
            Work<ClientSample> *work) noexcept
        {
            {
                std::lock_guard<std::mutex> lock(mutex);
                work_list.push(work);
            }
            work_available.notify_one();
        }

        static void kiss() noexcept
        {
            {
                std::lock_guard<std::mutex> lock(mutex);
                kiss_of_death = true;
            }
            work_available.notify_all();
        }

        static void revive() noexcept
        {
            std::lock_guard<std::mutex> lock(mutex);
            while (!work_list.empty())
            {
                delete work_list.front();
                work_list.pop();
            }
            kiss_of_death = false;
        }

        Worker(bool threaded = true) noexcept
        {
            if (threaded)
                while (!Worker<ClientSample>::kiss_of_death)
                    if (idle())
                        wait();
                    else
                        work();
            else
                while (!idle())
                    work();
        };
    };

    template <class ClientSample>
    std::mutex Worker<ClientSample>::mutex;

    template <class ClientSample>
    std::mutex Worker<ClientSample>::serialize;

    template <class ClientSample>
    volatile bool Worker<ClientSample>::kiss_of_death = false;

    template <class ClientSample>
    std::condition_variable Worker<ClientSample>::work_available;

    template <class ClientSample>
    std::queue<Work<ClientSample> *> Worker<ClientSample>::work_list;
}
