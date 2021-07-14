#include <iostream>
#include <sstream>
#include <exception>

#include "pursuit.h"

namespace EPP
{
    typedef DefaultSample<float> CloudSample;

    class CloudPursuer : public SamplePursuer<CloudSample>, Server
    {
    public:
        std::unique_ptr<Request> start(
            const CloudSample sample,
            const Parameters parameters) noexcept
        {
            std::unique_ptr<WorkRequest> request = std::unique_ptr<WorkRequest>(new WorkRequest(parameters, this));
            PursueProjection<CloudSample>::start(sample, parameters, request);

            return request;
        }

        void start(const json &encoded)
        {
            CloudSample sample(encoded);
            Parameters parameters(encoded);
            sample.wait();
            sample.subset.wait();
            start(sample, parameters);
        }

        void finish(Request *request) noexcept
        {
            Result *result = request->_result.get();
            json *encoded; //= *result;
            // send it out the wire
            Pursuer::finish(request);
        };

        CloudPursuer(Parameters parameters) noexcept
            : SamplePursuer<CloudSample>(parameters, std::thread::hardware_concurrency())
        {
            Worker<CloudSample>::revive();
            for (unsigned int i = 0; i < workers.size(); i++)
                workers[i] = std::thread(
                    []()
                    { Worker<CloudSample> worker; });
        }

        ~CloudPursuer()
        {
            Worker<CloudSample>::kiss();
            for (unsigned int i = 0; i < workers.size(); i++)
                workers[i].join();
        };
    };
}

using namespace EPP;

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " endpoint\n";
        return 1;
    }
    
    // set up the network

    CloudPursuer pursuer(Default);

    return 0;
};