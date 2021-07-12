#include <iostream>
#include <exception>

#include "pursuit.h"

namespace EPP
{
    typedef DefaultSample<float> CloudSample;

    class CloudPursuer : public SamplePursuer<CloudSample>
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

        void finish(json *encoded) noexcept
        {
            // send it out the wire
        };

        void finish(Request *request) noexcept
        {
            Result *result = request->_result.get();
            json *encoded ;//= *result;
            finish(encoded);
        };

        // void finish(Request *request) noexcept
        // {
        //     std::unique_lock<std::mutex> lock(mutex);
        //     requests.erase(request->key());
        //     completed.notify_all();
        // }
        // };

        CloudPursuer() noexcept
            : SamplePursuer<CloudSample>(std::thread::hardware_concurrency())
        {
            Worker<CloudSample>::revive();
            for (int i = 0; i < workers.size(); i++)
                workers[i] = std::thread(
                    []()
                    { EPP::Worker<CloudSample> worker; });
        }

        ~CloudPursuer()
        {
            Worker<CloudSample>::kiss();
            for (int i = 0; i < workers.size(); i++)
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

    CloudPursuer pursuer;
    while (true)
    {
        json payload;
        CloudSample sample(payload);
        Parameters parameters(payload);
        pursuer.start(sample, parameters);
    }

    return 0;
};