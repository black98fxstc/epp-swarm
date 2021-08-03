#include "pursuit.h"

namespace EPP
{
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

    MATLAB_Remote::MATLAB_Remote(
        const Parameters &parameters) noexcept
        : MATLAB_Pursuer(parameters, 0){};

    MATLAB_Remote::~MATLAB_Remote() = default;
}
