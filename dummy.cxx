/**
 * dummy EPP client for testing
 * */
#include <client.h>

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " endpoint\n";
        return 1;
    }

    // AWS has to be initted *before* our constructors can run
    EPP::Init();
    try
    {
        EPP::Client client;

        // get some data from somewhere? CSV?        
        double data[100][10];

        // stage it
        EPP::DefaultSample<double> sample(10, 100, (double *)&data);
        client.stage(sample);
        // make a subset of everything and stage that
        EPP::Subset everything(sample);
        for (long i = 0; i < sample.events; i++)
            everything[i] = true;
        client.stage(everything);

        // make up a request object
        json request;
        request["action"] = "EPP!";
        json smp;
        smp["measurments"] = sample.measurments;
        smp["events"] = sample.events;
        smp["key"] = sample.get_key();
        request["sample"] = smp;
        json sub;
        sub["key"] = everything.get_key();
        request["subset"] = sub;

        json response = client.ajax(argv[1], request);

        std::cout << request.dump(4) << std::endl;
        std::cout << response.dump(4) << std::endl;
    }
    catch (std::runtime_error e)
    {
        std::cout << e.what() << std::endl;
        return 1;
    }
    EPP::Finish();

    return 0;
};