#include <iostream>
#include <fstream>
#include <thread>

#include "client.h"
#include "pursuit.h"
#include "MATLAB.h"

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
        EPP::Measurment measurements = std::stoi(argv[1]);
        EPP::Event events = std::stol(argv[2]);
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

        EPP::Parameters parameters = EPP::Default;
        parameters.finalists = 6;
        parameters.recursive = true;
        parameters.W = .01; // .006 works well on Eliver and Cytek
        // parameters.sigma = 3.0;          // 3 to 5 maybe 6 are reasonable
        // less than three probably very noisy
        // parameters.censor.resize(measurements); // if empty censoring is disabled
        // parameters.censor[5] = true; // censor measurment 5
        EPP::MATLAB_Local pursuer(parameters, threads);
        const EPP::MATLAB_Sample sample(measurements, events, data);
        EPP::SampleSubset<EPP::MATLAB_Sample> *subset = new EPP::SampleSubset<EPP::MATLAB_Sample>(sample);

        EPP::Analysis<EPP::MATLAB_Sample> *analysis = pursuer.analyze(sample, subset, parameters);
        int i = 0;
        while (true)
            if (i < analysis->size())
            {
                const EPP::Lysis *lysis = (*analysis)(i++);

                std::cout << "projections " << lysis->projections << " avg passes " << (double)lysis->passes / (double)lysis->projections << " clusters " << (double)lysis->clusters / (double)lysis->projections << " graphs " << (double)lysis->graphs / (double)lysis->projections << " ms " << lysis->milliseconds.count() << std::endl;
                if (lysis->success())
                    std::cout << "best score " << lysis->winner().X << " " << lysis->winner().Y << "  " << lysis->winner().score << std::endl;
                else
                    std::cout << "no split" << std::endl;
            }
            else if (analysis->complete())
                break;
            else
                analysis->wait();
        std::cout << "compute " << analysis->compute_time.count() << "   clock " << analysis->milliseconds.count() << std::endl;

        delete[] data;
    }
    catch (std::runtime_error e)
    {
        std::cout << e.what() << std::endl;
        return 1;
    }

    return 0;
}