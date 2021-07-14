namespace EPP
{
    /**
     * Thread safe queues and so forth, to manage dividing up
     * the workload between threads. boilerplate basically
     * but the templates mean it has to be in the headers
     **/
    class WorkRequest : public Request
    {
    protected:
        static std::mutex mutex;
        static std::condition_variable completed;
        volatile unsigned int outstanding = 0;

    public:
        Result *working_result()
        {
            return _result.get();
        }

        void start()
        {
            std::unique_lock<std::mutex> lock(WorkRequest::mutex);
            ++outstanding;
        }

        void finish()
        {
            std::unique_lock<std::mutex> lock(WorkRequest::mutex);
            if (--outstanding == 0)
            {
                Request::finish();
                completed.notify_all();
            }
        }

        bool finished()
        {
            return outstanding == 0;
        };

        void wait()
        {
            std::unique_lock<std::mutex> lock(mutex);
            while (outstanding > 0)
                completed.wait(lock);
        };

        WorkRequest(
            Parameters parameters,
            Pursuer *pursuer) noexcept
            : Request(parameters, pursuer){};

        WorkRequest(
            const json &encoded,
            Pursuer *pursuer) noexcept
            : Request(pursuer){};
    };

    // abstract class representing a unit of work to be done
    // virtual functions let subclasses specialize tasks

    template <class ClientSample>
    class Work
    {
    public:
        const ClientSample sample;
        const Parameters parameters;
        WorkRequest *const request;

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
            request->finish();
        };

    protected:
        explicit Work(
            const ClientSample &sample,
            const Parameters parameters,
            WorkRequest *request) noexcept
            : sample(sample), parameters(parameters), request(request)
        {
            request->start();
        };
    };

    // a generic worker thread. looks for work, does it, deletes it
    // virtual functions in the work object do all the real work
    template <class ClientSample>
    class Worker
    {
    protected:
        static std::mutex serialize;
        static std::mutex mutex;
        static std::queue<Work<ClientSample> *> work_list;
        static std::condition_variable work_available;
        volatile static bool kiss_of_death;

        void work() noexcept
        {
            Work<ClientSample> *work = dequeue();
            work->parallel();
            {
                std::unique_lock<std::mutex> lock(serialize);
                work->serial();
            }
            delete work;
        };

        Work<ClientSample> *dequeue() noexcept
        {
            std::unique_lock<std::mutex> lock(mutex);
            Work<ClientSample> *work = work_list.front();
            work_list.pop();
            return work;
        }

        bool idle()
        {
            std::unique_lock<std::mutex> lock(mutex);
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
            std::unique_lock<std::mutex> lock(mutex);
            work_list.push(work);
            work_available.notify_one();
        }

        static void kiss() noexcept
        {
            std::unique_lock<std::mutex> lock(mutex);
            kiss_of_death = true;
            work_available.notify_all();
        }

        static void revive() noexcept
        {
            std::unique_lock<std::mutex> lock(mutex);
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

    static std::random_device random;

    std::mt19937_64 Request::generate(random());

    std::mutex WorkRequest::mutex;

    std::condition_variable WorkRequest::completed;

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