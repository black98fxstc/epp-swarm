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
        PursueProjection<MATLAB_Sample>::start(sample, parameters, request);

        if (threads == 0)
            Worker<MATLAB_Sample> worker(false);

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
        : MATLAB_Pursuer(std::thread::hardware_concurrency()){};

    MATLAB_Pursuer::~MATLAB_Pursuer()
    {
        // tell the workers to exit and wait for them to shut down
        Worker<MATLAB_Sample>::kiss();
        for (int i = 0; i < threads; i++)
            workers[i].join();
        delete[] workers;
    }
}
