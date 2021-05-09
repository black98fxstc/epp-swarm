#include <cstring>
#include <string>
#include <iostream>
#include <client.h>

using json = nlohmann::json;

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " endpoint\n";
        return 1;
    }
    EPP::Client client;

    // trivial ajax transaction

    json request;
    request["action"] = "Do something!";
    request["argument"] = "something to work on";

    json response = client.ajax(argv[1], request);

    std::cout << request.dump(4) << std::endl;
    std::cout << response.dump(4) << std::endl;

    // stage double data in memory to the server as floats

    double data[100][10];
    EPP::DefaultSample<double> one(10, 100, (double*)&data);
    client.stage(one);

    // remote workers fetch data as floats by the hash passed via JSON

    EPP::DefaultSample<float> two(one.measurments, one.events, one.get_key());
    client.fetch(two);

    double transpose[10][100];
    EPP::TransposeSample<double> three(10, 100, (double*)&transpose);
    client.stage(three);

    float **pointers;
    EPP::PointerSample<float> four(10, 100, pointers);
    client.stage(four);

    EPP::Subset first(one);
    first[1] = true;
    client.stage(first);
    client.fetch(first);

    return 0;
}