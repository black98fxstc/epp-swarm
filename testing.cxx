#include <iostream>
#include <fstream>
#include <thread>

#include "client.h"
#include "pursuit.h"
#include "MATLAB.h"

int main(int argc, char *argv[])
{
    if (argc < 4 || argc > 7)
    {
        std::cout << "Usage: " << argv[0] << " <measurements> <events> <csv-file> [<parameter-file>|default [<gating-file>|- [<threads>]]]\n";
        return 1;
    }

    try
    {
        // program arguments
        EPP::Measurement measurements = std::stoi(argv[1]);
        EPP::Event events = std::stol(argv[2]);
        int threads = std::thread::hardware_concurrency();
        if (argc > 6)
            threads = std::stoi(argv[6]);
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
        };

        EPP::MATLAB_Local pursuer(parameters, threads);
        const EPP::MATLAB_Sample sample(measurements, events, data);
        EPP::SampleSubset<EPP::MATLAB_Sample> *subset = new EPP::SampleSubset<EPP::MATLAB_Sample>(sample);

        EPP::Analysis<EPP::MATLAB_Sample> *analysis = pursuer.analyze(sample, subset, parameters);
        unsigned int i = 0;
        // report results as they come in (optional)
        while (!analysis->complete())
            if (i < analysis->size())
            {
                const EPP::Lysis *lysis = (*analysis)(i++);

                std::cerr << "projections " << lysis->projections << " avg passes " << (double)lysis->passes / (double)lysis->projections << " clusters " << (double)lysis->clusters / (double)lysis->projections << " graphs " << (double)lysis->graphs / (double)lysis->projections << " ms " << lysis->milliseconds.count() << std::endl;
                if (lysis->success())
                    std::cerr << "best score " << lysis->winner().X << " " << lysis->winner().Y << "  " << lysis->winner().score << std::endl;
                else
                    std::cerr << "no split" << std::endl;
            }
            else
                analysis->wait();

        std::cerr << "total projections " << analysis->projections << " passes " << analysis->passes << " clusters " << analysis->clusters << " graphs " << analysis->graphs << std::endl;
        std::cerr << "avg passes " << (double)analysis->passes / (double)analysis->projections << " clusters " << (double)analysis->clusters / (double)analysis->projections << " graphs " << (double)analysis->graphs / (double)analysis->projections << std::endl;
        std::cerr << analysis->types.size() << " types in " << analysis->size() << " subsets found    compute " << analysis->compute_time.count() << " clock " << analysis->milliseconds.count() << " ms" << std::endl;

        if (argc > 5 && std::strcmp(argv[5], "-"))
        {
            std::ofstream treefile(argv[5], std::ios::out);
            treefile << subset->tree().dump();
            treefile.close();
        }
        else
            std::cout << subset->tree().dump(2) << std::endl;

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