#include <iostream>
#include <sstream>
#include <exception>

#include "pursuit.h"

namespace EPP
{
    /**
     * remote worker instance
     **/

    void CloudPursuer::start(const json &encoded)
    {
        // unsigned short measurements; // from json
        // unsigned long int events;
        // Key key1, key2;

        // SampleSubset<CloudSample> subset(events, key1);
        // CloudSample sample(measurements, events, subset, key2);
        // Parameters parameters(encoded);
        // sample.wait(); // blob fault and wait if necessary
        // subset.wait(); // may not need to reload recent data
        // start(sample, parameters);
    }

    void CloudPursuer::finish(ClientRequest<CloudSample> *request) noexcept
    {
        // _Result *result = request;
        json encoded; //= *result;
        // send it out the wire
        out(encoded);
        // Pursuer::finish(request);
    };

    void CloudPursuer::finish(const json &encoded)
    {
        Pursuer::finish(encoded);
    };

    json CloudPursuer::remote()
    {
        return Remote::in();
    }

    CloudPursuer::CloudPursuer(const Parameters &parameters) noexcept
        : SamplePursuer<CloudSample>(parameters)
    {
        Worker<CloudSample>::revive();
        for (unsigned int i = 0; i < workers.size(); i++)
            workers[i] = std::thread(
                []()
                { Worker<CloudSample> worker; });
    }

    CloudPursuer::~CloudPursuer()
    {
        Worker<CloudSample>::kiss();
        for (unsigned int i = 0; i < workers.size(); i++)
            workers[i].join();
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
        json encoded = pursuer.remote();
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