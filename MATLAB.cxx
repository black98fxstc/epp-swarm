#include "constants.h"
#include "client.h"
#include "boundary.h"
#include "modal.h"
#include "pursuit.h"

namespace EPP
{
    void MATLAB_Pursuer::start(
        MATLAB_Sample sample,
        Parameters parameters) noexcept
    {
        qualified_measurements.clear();

        _result = std::shared_ptr<Result>(new Result);
        _result->outcome = Result::EPP_no_qualified;
        _result->best_score = std::numeric_limits<double>::infinity();
        _result->begin = std::chrono::steady_clock::now();

        std::unique_lock<std::recursive_mutex> lock(mutex);
        for (int measurement = 0; measurement < sample.measurements; ++measurement)
            Worker<MATLAB_Sample>::work_list.push(new QualifyMeasurement<MATLAB_Sample>(sample, parameters, measurement));
        work_available.notify_all();

        if (threads == 0)
            while (!Worker<MATLAB_Sample>::work_list.empty())
            {
                Work<MATLAB_Sample> *work = Worker<MATLAB_Sample>::work_list.front();
                Worker<MATLAB_Sample>::work_list.pop();
                lock.unlock();
                work->parallel();
                lock.lock();
                work->serial();
                delete work;
            }
    };

    void MATLAB_Pursuer::start(
        const int measurements,
        const long events,
        float *data,
        std::vector<bool> &subset) noexcept
    {
        MATLAB_Sample sample{measurements, events, data, subset};
        start(sample, Default);
    };

    void MATLAB_Pursuer::start(
        const int measurements,
        const long events,
        float *data) noexcept
    {
        MATLAB_Sample sample{measurements, events, data};
        start(sample, Default);
    };

    bool MATLAB_Pursuer::finished() noexcept
    {
        std::unique_lock<std::recursive_mutex> lock(EPP::mutex);
        return !EPP::work_outstanding;
    };

    void MATLAB_Pursuer::wait() noexcept
    {
        std::unique_lock<std::recursive_mutex> lock(EPP::mutex);
        while (EPP::work_outstanding)
            EPP::work_completed.wait(lock);
    };

    std::shared_ptr<Result> MATLAB_Pursuer::result() noexcept
    {
        wait();
        _result->end = std::chrono::steady_clock::now();
        _result->milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(_result->end - _result->begin);
        return _result;
    };

    std::shared_ptr<Result> MATLAB_Pursuer::pursue(
        const int measurements,
        const long events,
        float *data,
        std::vector<bool> &subset) noexcept
    {
        start(measurements, events, data, subset);
        return result();
    };

    MATLAB_Pursuer::MATLAB_Pursuer(int threads) noexcept
        : threads(threads < 0 ? std::thread::hardware_concurrency() : threads)
    {
        // start some worker threads
        kiss_of_death = false;
        workers = new std::thread[threads];
        for (int i = 0; i < threads; i++)
            workers[i] = std::thread(
                []()
                { EPP::Worker<MATLAB_Sample> worker; });
    };

    MATLAB_Pursuer::MATLAB_Pursuer() noexcept
        : MATLAB_Pursuer(std::thread::hardware_concurrency()){};

    MATLAB_Pursuer::~MATLAB_Pursuer()
    {
        // tell the workers to exit and wait for them to shut down
        EPP::kiss_of_death = true;
        {
            std::unique_lock<std::recursive_mutex> lock(EPP::mutex);
            EPP::work_available.notify_all();
        }
        for (int i = 0; i < threads; i++)
            workers[i].join();
        delete[] workers;
        _result.reset();
    };
}
