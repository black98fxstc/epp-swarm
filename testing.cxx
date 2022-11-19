
/*
 * Developer: Wayne Moore <wmoore@stanford.edu>
 * Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
 * License: BSD 3 clause
 */
#include <iostream>
#include <fstream>

#include "pursuit.h"
#include "MATLAB.h"

int main(int argc, char *argv[])
{
    if (argc < 4 || argc > 9)
    {
        std::cout << "Usage: " << argv[0] << " <measurements> <events> <data-csv> [<parameters-json>|default [<result-json|none [ <class-csv>|none [<threads>]]]]\n";
        return 1;
    }

    try
    {
        // program arguments
        EPP::Measurement measurements = (EPP::Measurement)std::stoi(argv[1]);
        std::vector<std::string> markers(measurements);
        EPP::Event events = std::stol(argv[2]);
        int threads = std::thread::hardware_concurrency();
        if (argc > 7)
            threads = std::stoi(argv[7]);
        if (threads < 0)
            threads = std::thread::hardware_concurrency();

        // get the data file
        float *data = new float[measurements * events];
        std::ifstream datafile(argv[3], std::ios::in);
        std::string line;
        std::string value;
        std::getline(datafile, line);
        std::stringstream sstr(line, std::ios::in);
        for (int j = 0; j < measurements; j++)
        {
            std::getline(sstr, value, ',');
            while (std::isspace(value.front()))
                value.erase(value.begin());
            markers[j] = value;
        }
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
        parameters.min_relative = .005;

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

                std::cerr << "projections " << lysis->projections << " avg passes " << (double)lysis->passes / (double)lysis->projections << " clusters " << (double)lysis->clusters / (double)lysis->projections << " graphs " << (double)lysis->graphs / (double)lysis->projections << " merges " << (double)lysis->merges / (double)lysis->projections << " ms " << lysis->milliseconds.count() << std::endl;
                if (lysis->success())
                    std::cerr << "best score " << lysis->winner().X << " " << lysis->winner().Y << "  " << lysis->winner().score << std::endl;
                else
                    std::cerr << "no split" << std::endl;
            }
            else
                analysis->wait();

        std::cerr << "total projections " << analysis->projections << " passes " << analysis->passes << " clusters " << analysis->clusters << " graphs " << analysis->graphs << " merges " << analysis->merges << std::endl;
        std::cerr << "avg passes " << (double)analysis->passes / (double)analysis->projections << " clusters " << (double)analysis->clusters / (double)analysis->projections << " graphs " << (double)analysis->graphs / (double)analysis->projections << " merges " << (double)analysis->merges / (double)analysis->projections << std::endl;
        std::cerr << analysis->types() << " types in " << analysis->size() << " subsets found    compute " << analysis->compute_time.count() << " clock " << analysis->milliseconds.count() << " ms" << std::endl;

        json result;
        result["version"] = "0.1";
        result["gating"] = analysis->gating();
        result["taxonomy"] = (json)*analysis->classify();
        std::vector<EPP::Taxon *> phenogram = analysis->phenogram();

        if (argc > 5 && std::strcmp(argv[5], "none"))
        {
            std::ofstream resultfile(argv[5], std::ios::out);
            resultfile << result.dump();
            resultfile.close();
        }

        if (argc > 6 && std::strcmp(argv[6], "none"))
        {
            std::ofstream classfile(argv[6], std::ios::out);
            for (EPP::Unique *u = analysis->classification; u != analysis->classification + events; ++u)
                classfile << *u << std::endl;
            classfile.close();
        }

        std::cout << std::endl;
        std::cout << EPP::Taxonomy::ascii(phenogram, markers);
        std::cout << std::endl;

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
