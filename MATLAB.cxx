#include "pursuit.h"

namespace EPP
{
    std::unique_ptr<Request> MATLAB_Pursuer::start(
        const MATLAB_Sample sample,
        const Parameters parameters) noexcept
    {
        std::unique_ptr<WorkRequest> request = std::unique_ptr<WorkRequest>(new WorkRequest(parameters, this));
        PursueProjection<MATLAB_Sample>::start(sample, parameters, request);

        if (workers.size() == 0)
            Worker<MATLAB_Sample> worker(false);

        return request;
    };

    /**
     * MATLAB convenience routines
     * always use the pursuer parameters
     **/
    std::unique_ptr<Request> MATLAB_Pursuer::start(
        const unsigned short int measurements,
        const unsigned long int events,
        const float *const data,
        Subset &subset) noexcept
    {
        MATLAB_Sample sample(measurements, events, data, subset);
        return start(sample, parameters);
    }

    std::unique_ptr<Request> MATLAB_Pursuer::start(
        const unsigned short int measurements,
        const unsigned long int events,
        const float *const data) noexcept
    {
        MATLAB_Sample sample(measurements, events, data);
        return start(sample, parameters);
    }

    std::shared_ptr<Result> MATLAB_Pursuer::pursue(
        const unsigned short int measurements,
        const unsigned long int events,
        const float *const data,
        Subset &subset) noexcept
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

    /**
     * in process client for MATLAB
     **/

    MATLAB_Local::MATLAB_Local(
        Parameters parameters,
        int threads) noexcept
        : MATLAB_Pursuer(parameters, threads < 0 ? std::thread::hardware_concurrency() : threads)
    {
        Worker<MATLAB_Sample>::revive();
        for (unsigned int i = 0; i < workers.size(); i++)
            workers[i] = std::thread(
                []()
                { EPP::Worker<MATLAB_Sample> worker; });
    }

    MATLAB_Local::~MATLAB_Local()
    {
        Worker<MATLAB_Sample>::kiss();
        for (unsigned int i = 0; i < workers.size(); i++)
            workers[i].join();
    }

    /**
     * MATLAB client accessing remote instance
     **/

    std::unique_ptr<Request> MATLAB_Remote::start(
        MATLAB_Sample sample,
        const Parameters parameters) noexcept
    {
        std::unique_ptr<WorkRequest> request = std::unique_ptr<WorkRequest>(new WorkRequest(parameters, this));
        json encoded = (json) * (request.get());
        Remote::out(encoded);
        return request;
    }

    void MATLAB_Remote::finish(
        const json &encoded)
    {
        // from somewhere
        Key request_key; // from JSON
        Request *request = requests.find(request_key)->second;
        Result *result = request->result().get();
        *result = encoded;
        request->finish();
    }

    MATLAB_Remote::MATLAB_Remote(
        Parameters parameters) noexcept
        : MATLAB_Pursuer(parameters, 0){};

    MATLAB_Remote::~MATLAB_Remote() = default;
}
