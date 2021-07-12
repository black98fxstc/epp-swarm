#include "pursuit.h"

namespace EPP
{
    std::unique_ptr<Request> MATLAB_Pursuer::start(
        const unsigned short int measurements,
        const unsigned long int events,
        const float *const data,
        std::vector<bool> &subset) noexcept
    {
        MATLAB_Sample sample(measurements, events, data, subset);
        return SamplePursuer<MATLAB_Sample>::start(sample);
    }

    std::unique_ptr<Request> MATLAB_Pursuer::start(
        const unsigned short int measurements,
        const unsigned long int events,
        const float *const data) noexcept
    {
        MATLAB_Sample sample(measurements, events, data);
        return SamplePursuer<MATLAB_Sample>::start(sample);
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
    }

    std::unique_ptr<Request> MATLAB_Local::start(
        const MATLAB_Sample sample,
        const Parameters parameters) noexcept
    {
        std::unique_ptr<WorkRequest> request = std::unique_ptr<WorkRequest>(new WorkRequest(parameters, this));
        PursueProjection<MATLAB_Sample>::start(sample, parameters, request);

        if (workers.size() == 0)
            Worker<MATLAB_Sample> worker(false);

        return request;
    }

    MATLAB_Local::MATLAB_Local(int threads) noexcept
        : MATLAB_Pursuer(threads < 0 ? std::thread::hardware_concurrency() : threads)
    {
        Worker<MATLAB_Sample>::revive();
        for (int i = 0; i < workers.size(); i++)
            workers[i] = std::thread(
                []()
                { EPP::Worker<MATLAB_Sample> worker; });
    }

    MATLAB_Local::MATLAB_Local() noexcept
        : MATLAB_Local(std::thread::hardware_concurrency()){}

    MATLAB_Local::~MATLAB_Local()
    {
        Worker<MATLAB_Sample>::kiss();
        for (int i = 0; i < workers.size(); i++)
            workers[i].join();
    }

    std::unique_ptr<Request> MATLAB_Remote::start(
        const MATLAB_Sample sample,
        const Parameters parameters) noexcept
    {
        std::unique_ptr<WorkRequest> request = std::unique_ptr<WorkRequest>(new WorkRequest(parameters, this));
        json encoded = (json)*(request.get());
        // send it in it's way
        return request;
    }

    void MATLAB_Remote::finish(
        std::string incoming) noexcept
    {
        json payload;
        // from somewhere
        Key request_key; // from JSON
        Request *request = requests.find(request_key)->second;
        Result *result = request->result().get();
        *result = payload;
        request->finish();
    }

    MATLAB_Remote::MATLAB_Remote() noexcept : MATLAB_Pursuer(0) {}

    MATLAB_Remote::~MATLAB_Remote() = default;
}
