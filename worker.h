namespace EPP
{
    /**
     * Thread safe queues and so forth, to manage dividing up
     * the workload between threads. boilerplate basically
     * but the templates mean it has to be in the headers
     **/
    // class _WorkRequest : public _Request
    // {
    // protected:
    //     // static std::mutex mutex;
    //     // static std::condition_variable completed;
    //     // volatile unsigned int outstanding = 0;

    // public:
    //     // void start()
    //     // {
    //     //     std::unique_lock<std::mutex> lock(_WorkRequest::mutex);
    //     //     ++outstanding;
    //     // }

    //     // void finish()
    //     // {
    //     //     std::unique_lock<std::mutex> lock(_WorkRequest::mutex);
    //     //     if (--outstanding == 0)
    //     //     {
    //     //         completed.notify_all();
    //     //         _Request::finish();
    //     //     }
    //     // }

    //     // bool finished()
    //     // {
    //     //     return outstanding == 0 && Request::finished();
    //     // };

    //     // void wait()
    //     // {
    //     //     std::unique_lock<std::mutex> lock(mutex);
    //     //     while (outstanding > 0)
    //     //         completed.wait(lock);
    //     //     _Request::wait();
    //     // };
    // public:
    //     _WorkRequest(
    //         Parameters parameters,
    //         Pursuer *pursuer) noexcept
    //         : _Request(pursuer, parameters){};

    //     _WorkRequest(
    //         const json &encoded,
    //         Pursuer *pursuer) noexcept
    //         : _Request(pursuer){};
    // };

    // class WorkRequest : public Request
    // {
    // public:
    //     _Result *working_result()
    //     {
    //         auto i = (*this)->outstanding;
    //         return (*this)->_result.get();
    //     }

    //     WorkRequest(
    //         _Request *request) : Request(request)
    //     {
    //         request->outstanding = 0;
    //     };

    //     void start()
    //     {
    //         {
    //             std::unique_lock<std::mutex> lock((*this)->pursuer->mutex);
    //             ++(*this)->outstanding;
    //         }
    //     }

    //     void finish()
    //     {
    //         {
    //             std::unique_lock<std::mutex> lock((*this)->pursuer->mutex);
    //             --(*this)->outstanding;
    //         }
    //         if ((*this)->outstanding == 0)
    //             Request::finish();
    //     }

    //     void wait()
    //     {
    //         {
    //             std::unique_lock<std::mutex> lock((*this)->pursuer->mutex);
    //             while ((*this)->outstanding > 0)
    //                 (*this)->completed.wait(lock);
    //         }
    //         Request::wait();
    //     };
    // };

    // abstract class representing a unit of work to be done
    // virtual functions let subclasses specialize tasks

    template <class ClientSample>
    class Work
    {
    public:
        const ClientSample sample;
        const Parameters parameters;
        Request request;

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
            --request;
        };

    protected:
        explicit Work(
            const ClientSample &sample,
            const Parameters parameters,
            Request request) noexcept
            : sample(sample), parameters(parameters), request(request)
        {
            ++request;
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

    // std::mutex _WorkRequest::mutex;

    // std::condition_variable _WorkRequest::completed;

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