#include <iostream>
#include <sstream>
#include <exception>

#include "pursuit.h"

namespace EPP
{
    /**
     * remote worker instance
     **/
    typedef DefaultSample<float> CloudSample;

    class CloudPursuer : Remote, public SamplePursuer<CloudSample>
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
            unsigned short measurements; // from json
            unsigned long int events;
            Key key1, key2;

            Subset subset(events, key1);
            CloudSample sample(measurements, events, subset, key2);
            Parameters parameters(encoded);
            sample.wait();  // blob fault and wait if necessary
            subset.wait();  // may not need to reload recent data
            start(sample, parameters);
        }

        void finish(Request *request) noexcept
        {
            Result *result = request->_result.get();
            json encoded; //= *result;
            // send it out the wire
            out(encoded);
            Pursuer::finish(request);
        };

        void finish(const json &encoded)
        {
            Pursuer::finish(encoded);
        };

        json in()
        {
            return Remote::in();
        }

        CloudPursuer(Parameters parameters) noexcept
            : SamplePursuer<CloudSample>(parameters)
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
    Parameters parameters;
    
    // set up the network

    CloudPursuer pursuer(parameters);
    while (true)
    {
        json encoded = pursuer.in();
        Remote::Service service = Remote::Service::request;
        switch (service)
        {
            case Remote::Service::request:
            pursuer.start(encoded);
            break;

            case Remote::Service::result:
            pursuer.finish(encoded);
            break;
        }
    }

    return 0;
};