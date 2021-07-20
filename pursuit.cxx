#include <memory>

#include "constants.h"
#include "boundary.h"
#include "modal.h"
// #include "pursuer.h"
#include "pursuit.h"
#include "client.h"

namespace EPP
{
    template <class ClientSample>
    void Pursuer<ClientSample>::start(
        const ClientSample sample,
        const Parameters parameters) noexcept
    {
        qualified_measurements.clear();

        _result = std::shared_ptr<Result>(new Result);
        _result->outcome = Result::EPP_no_qualified;
        _result->best_score = std::numeric_limits<double>::infinity();
        _result->begin = std::chrono::steady_clock::now();

        std::unique_lock<std::recursive_mutex> lock(mutex);
        for (int measurement = 0; measurement < sample.measurements; ++measurement)
            work_list.push(new QualifyMeasurement(sample, parameters, result));
        work_available.notify_all();

        if (threads == 0)
            while (!work_list.empty())
            {
                Work *work = work_list.front();
                work_list.pop();
                lock.unlock();
                work->parallel();
                lock.lock();
                work->serial();
                delete work;
            }
    };

    template <class ClientSample>
    bool Pursuer<ClientSample>::finished() noexcept
    {
        std::unique_lock<std::recursive_mutex> lock(EPP::mutex);
        return !EPP::work_outstanding;
    };

    template <class ClientSample>
    void Pursuer<ClientSample>::wait() noexcept
    {
        std::unique_lock<std::recursive_mutex> lock(EPP::mutex);
        while (EPP::work_outstanding)
            EPP::work_completed.wait(lock);
    };

    template <class ClientSample>
    std::shared_ptr<Result> Pursuer<ClientSample>::result() noexcept
    {
        wait();
        _result->end = std::chrono::steady_clock::now();
        _result->milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(_result->end - _result->begin);
        return _result;
    };

    template <class ClientSample>
    std::shared_ptr<Result> pursue(
        const int measurements,
        const long events,
        const float *const data,
        std::vector<bool> &subset) noexcept;

    template <class ClientSample>
    Pursuer<ClientSample>::Pursuer() noexcept
        : Pursuer(std::thread::hardware_concurrency()){};

    template <class ClientSample>
    Pursuer<ClientSample>::Pursuer(int threads) noexcept
        : threads(threads < 0 ? std::thread::hardware_concurrency() : threads)
    {
        // start some worker threads
        workers = new std::thread[threads];
        for (int i = 0; i < threads; i++)
            workers[i] = std::thread(
                []()
                { EPP::Worker worker; });
    };

    template <class ClientSample>
    Pursuer<ClientSample>::~Pursuer()
    {
        // tell the workers to exit and wait for them to shut down
        kiss_of_death = true;
        {
            std::unique_lock<std::recursive_mutex> lock(EPP::mutex);
            work_available.notify_all();
        }
        for (int i = 0; i < threads; i++)
            workers[i].join();
        delete[] workers;
        _result.reset();
    };
};