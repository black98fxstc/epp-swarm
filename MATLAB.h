#include "client.h"

namespace EPP
{
    typedef TransposeSample<float> MATLAB_Sample;

    class MATLAB_Pursuer : public Pursuer<MATLAB_Sample>
    {
    public:
    protected:
        MATLAB_Pursuer() = delete;

        MATLAB_Pursuer(
            const Parameters &parameters,
            int threads) noexcept
            : Pursuer<MATLAB_Sample>(parameters, threads){};

        ~MATLAB_Pursuer() = default;
    };

    class MATLAB_Local : public MATLAB_Pursuer
    {
    public:
        int getThreads() const noexcept
        {
            return workers.size();
        };

        MATLAB_Local(
            const Parameters &parameters,
            int threads) noexcept
            : MATLAB_Pursuer(parameters, threads < 0 ? std::thread::hardware_concurrency() : threads)
        {
            Worker<MATLAB_Sample>::revive();
            for (unsigned int i = 0; i < workers.size(); i++)
                workers[i] = std::thread(
                    []()
                    { EPP::Worker<MATLAB_Sample> worker; });
        }

        ~MATLAB_Local()
        {
            Worker<MATLAB_Sample>::kiss();
            for (unsigned int i = 0; i < workers.size(); i++)
                workers[i].join();
        }
    };

    class MATLAB_Remote : Remote, public MATLAB_Pursuer
    {
    public:
        void finish(
            const json &encoded);

        MATLAB_Remote(
            const Parameters &parameters) noexcept;

        ~MATLAB_Remote();
    };
    /**
     * in process client for MATLAB
     **/

    /**
     * MATLAB client accessing remote instance
     **/

    MATLAB_Remote::MATLAB_Remote(
        const Parameters &parameters) noexcept
        : MATLAB_Pursuer(parameters, 0){};

    MATLAB_Remote::~MATLAB_Remote() = default;
}
