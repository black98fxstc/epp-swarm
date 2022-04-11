#include "client.h"

namespace EPP
{
    typedef TransposeSample<float> MATLAB_Sample;

    /**
     * in process client for MATLAB
     **/

    class MATLAB_Local : public Pursuer<MATLAB_Sample>
    {
    public:
        int getThreads() const noexcept
        {
            return (int)workers.size();
        };

        MATLAB_Local(
            const Parameters &parameters,
            int threads = -1)
            : Pursuer<MATLAB_Sample>(parameters, threads){};
    };

    /**
     * MATLAB client accessing remote instance
     **/

    class MATLAB_Remote : Remote, public Pursuer<MATLAB_Sample>
    {
    public:
        MATLAB_Remote(
            const Parameters &parameters) noexcept
            : Pursuer<MATLAB_Sample>(parameters, 0){};
    };
}
