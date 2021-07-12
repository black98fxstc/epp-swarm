#include <iostream>
#include <fstream>
#include <thread>

#include "client.h"

int main(int argc, char *argv[])
{
    if (argc < 4 || argc > 5)
    {
        std::cout << "Usage: " << argv[0] << "<measurements> <events> <csv-file> [<threads>]\n";
        return 1;
    }

    try
    {
        // program arguments
        int measurements = std::stoi(argv[1]);
        unsigned long events = std::stol(argv[2]);
        int threads = std::thread::hardware_concurrency();
        if (argc == 5)
            threads = std::stoi(argv[4]);
        if (threads < 0)
            threads = std::thread::hardware_concurrency();

        // get the data file
        float *data = new float[measurements * events];
        std::ifstream datafile(argv[3], std::ios::in);
        std::string line;
        std::string value;
        std::getline(datafile, line);
        for (unsigned long int i = 0; i < events; i++)
        {
            std::getline(datafile, line);
            std::stringstream sstr(line, std::ios::in);
            for (int j = 0; j < measurements; j++)
            {
                std::getline(sstr, value, ',');
                data[i + events * j] = std::stof(value);
            }
        }
        datafile.close();

        EPP::MATLAB_Local pursuer(threads);                  // reusable, you can do many start/result calls
        EPP::MATLAB_Sample sample(measurements, events, data); // default constructor does range check
        EPP::Parameters parameters = EPP::Default;             // this is the default
        parameters.finalists = 6;
        parameters.W = .006;             // works well on Eliver and Cytek
        parameters.sigma = 3.0;          // 3 to 5 maybe 6 are reasonable
                                         // less than three probably very noisy
        // parameters.censor.resize(measurements); // if empty censoring is disabled
        // parameters.censor[5] = true; // censor measurment 5
        std::unique_ptr<EPP::Request> request = pursuer.start(sample, parameters);
        if (!request->finished()) // optional, used when you want to do something else while it runs
            request->wait();
        std::shared_ptr<EPP::Result> result = request->result();
        // result = pursuer.pursue(measurements, events, data); // convenience routine takes default parameters

        if (result->outcome() != EPP::Status::EPP_success)
            std::cout << "oops" << std::endl;
        else
        {
            // suitable for FlowJo and GatingML
            std::vector<EPP::Point> in_polygon = result->winner().in_polygon(parameters.W);

            std::cout << "projections " << result->projections << " avg passes " << (double)result->passes / (double)result->projections << " clusters " << (double)result->clusters / (double)result->projections << " graphs " << (double)result->graphs / (double)result->projections << " ms " << result->milliseconds.count() << std::endl;
            std::cout << "best score " << result->winner().X << " " << result->winner().Y << "  " << result->winner().score << std::endl;
        }

        delete[] data;
    }
    catch (std::runtime_error e)
    {
        std::cout << e.what() << std::endl;
        return 1;
    }

    return 0;
}