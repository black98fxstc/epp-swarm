#include <iostream>
#include <fstream>
#include <thread>

#include "client.h"
#include "pursuit.h"
#include "MATLAB.h"

int main(int argc, char *argv[])
{
    if (argc < 4 || argc > 6)
    {
        std::cout << "Usage: " << argv[0] << " <measurements> <events> <csv-file> [<parameter-file>|default [<threads>]]\n";
        return 1;
    }

    try
    {
        // program arguments
        EPP::Measurement measurements = std::stoi(argv[1]);
        EPP::Event events = std::stol(argv[2]);
        int threads = std::thread::hardware_concurrency();
        if (argc > 5)
            threads = std::stoi(argv[5]);
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

        // get the parameters
        EPP::Parameters parameters = EPP::Default;
        if (argc > 4 && std::strcmp(argv[4], "default"))
        {
            std::ifstream paramfile(argv[4], std::ios::in);
            parameters = json::parse(paramfile);
            paramfile.close();
        }

        EPP::MATLAB_Local pursuer(parameters, threads);
        const EPP::MATLAB_Sample sample(measurements, events, data);
        EPP::SampleSubset<EPP::MATLAB_Sample> *subset = new EPP::SampleSubset<EPP::MATLAB_Sample>(sample);

        EPP::Analysis<EPP::MATLAB_Sample> *analysis = pursuer.analyze(sample, subset, parameters);
        while (!analysis->complete())
            analysis->wait();

        std::cout << subset->tree().dump() << std::endl;

        delete analysis;
        delete subset;
        delete[] data;
    }
    catch (std::runtime_error e)
    {
        std::cout << e.what() << std::endl;
        return 1;
    }

    return 0;
}