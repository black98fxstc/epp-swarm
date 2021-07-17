#include "pursuit.h"

namespace EPP
{
    Request MATLAB_Pursuer::start(
        const MATLAB_Sample sample,
        const Parameters parameters) noexcept
    {
        Request request = Request(this, parameters);
        PursueProjection<MATLAB_Sample>::start(sample, parameters, request);

        if (workers.size() == 0)
            Worker<MATLAB_Sample> worker(false);

        return request;
    };

    /**
     * MATLAB convenience routines
     * always use the pursuer parameters
     **/
    Request MATLAB_Pursuer::start(
        const unsigned short int measurements,
        const unsigned long int events,
        const float *const data,
        Subset &subset) noexcept
    {
        MATLAB_Sample sample(measurements, events, data, subset);
        return start(sample, parameters);
    }

    Request MATLAB_Pursuer::start(
        const unsigned short int measurements,
        const unsigned long int events,
        const float *const data) noexcept
    {
        MATLAB_Sample sample(measurements, events, data);
        return start(sample, parameters);
    }

    Result MATLAB_Pursuer::pursue(
        const unsigned short int measurements,
        const unsigned long int events,
        const float *const data,
        Subset &subset) noexcept
    {
        return start(measurements, events, data, subset).result();
    }

    Result MATLAB_Pursuer::pursue(
        const unsigned short int measurements,
        const unsigned long int events,
        const float *const data) noexcept
    {
        return start(measurements, events, data).result();
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

    Request MATLAB_Remote::start(
        MATLAB_Sample sample,
        const Parameters parameters) noexcept
    {
        Request request = Request(this, parameters);
        json encoded = (json) * (request.get());

        Remote::out(encoded);

        return request;
    }

    MATLAB_Remote::MATLAB_Remote(
        Parameters parameters) noexcept
        : MATLAB_Pursuer(parameters, 0){};

    MATLAB_Remote::~MATLAB_Remote() = default;
}
