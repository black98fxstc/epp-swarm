#include "constants.h"
#include "client.h"
#include "boundary.h"
#include "modal.h"
#include "pursuit.h"

namespace EPP
{
    std::unique_ptr<Request> MATLAB_Pursuer::start(
        MATLAB_Sample sample,
        Parameters parameters) noexcept
    {
        std::unique_ptr<WorkRequest> request = std::unique_ptr<WorkRequest>(new WorkRequest(parameters));
        for (int measurement = 0; measurement < sample.measurements; ++measurement)
            if (parameters.censor.empty() || !parameters.censor.at(measurement))
                Worker<MATLAB_Sample>::enqueue(
                    new QualifyMeasurement<MATLAB_Sample>(sample, parameters, request.get(), measurement));

        if (threads == 0)
        {
            Worker<MATLAB_Sample> worker(false);
            // std::unique_lock<std::recursive_mutex> lock(Worker<MATLAB_Sample>::mutex);
            // while (!Worker<MATLAB_Sample>::work_list.empty())
            // {
            //     Work<MATLAB_Sample> *work = Worker<MATLAB_Sample>::work_list.front();
            //     Worker<MATLAB_Sample>::work_list.pop();
            //     lock.unlock();
            //     work->parallel();
            //     lock.lock();
            //     work->serial();
            //     delete work;
            // }
        }

        return request;
    }

    std::unique_ptr<Request> MATLAB_Pursuer::start(
        const unsigned short int measurements,
        const unsigned long int events,
        const float *const data,
        std::vector<bool> &subset) noexcept
    {
        MATLAB_Sample sample(measurements, events, data, subset);
        return start(sample, Default);
    }

    std::unique_ptr<Request> MATLAB_Pursuer::start(
        const unsigned short int measurements,
        const unsigned long int events,
        const float *const data) noexcept
    {
        MATLAB_Sample sample(measurements, events, data);
        return start(sample, Default);
    }

    // bool MATLAB_Pursuer::finished() noexcept
    // {
    //     std::unique_lock<std::recursive_mutex> lock(EPP::mutex);
    //     return !EPP::work_outstanding;
    // }

    // void MATLAB_Pursuer::wait() noexcept
    // {
    //     std::unique_lock<std::recursive_mutex> lock(EPP::mutex);
    //     while (EPP::work_outstanding)
    //         EPP::work_completed.wait(lock);
    // }

    // std::shared_ptr<Result> MATLAB_Pursuer::result() noexcept
    // {
    //     wait();
    //     _result->end = std::chrono::steady_clock::now();
    //     _result->milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(_result->end - _result->begin);
    //     return _result;
    // }

    std::shared_ptr<Result> MATLAB_Pursuer::pursue(
        const unsigned short int measurements,
        const unsigned long int events,
        const float *const data,
        std::vector<bool> &subset) noexcept
    {
        return start(measurements, events, data, subset)->result();
    }

    std::shared_ptr<Result> MATLAB_Pursuer::pursue(
        const unsigned short int measurements,
        const unsigned long int events,
        const float *const data) noexcept
    {
        return start(measurements, events, data)->result();
    };

    MATLAB_Pursuer::MATLAB_Pursuer(int threads) noexcept
        : threads(threads < 0 ? std::thread::hardware_concurrency() : threads)
    {
        Worker<MATLAB_Sample>::revive();
        // start some worker threads
        workers = new std::thread[threads];
        for (int i = 0; i < threads; i++)
            workers[i] = std::thread(
                []()
                { EPP::Worker<MATLAB_Sample> worker; });
    }

    MATLAB_Pursuer::MATLAB_Pursuer() noexcept
        : MATLAB_Pursuer(std::thread::hardware_concurrency())
        {
            Worker<MATLAB_Sample>::revive();
        };

    MATLAB_Pursuer::~MATLAB_Pursuer()
    {
        // tell the workers to exit and wait for them to shut down
        Worker<MATLAB_Sample>::kiss();
        for (int i = 0; i < threads; i++)
            workers[i].join();
        delete[] workers;
    }
}
