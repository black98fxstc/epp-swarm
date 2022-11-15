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
        }

        // set up the machinery
        EPP::MATLAB_Local pursuer(parameters, threads);
        const EPP::MATLAB_Sample sample(measurements, events, data);
        EPP::SampleSubset<EPP::MATLAB_Sample> *subset = new EPP::SampleSubset<EPP::MATLAB_Sample>(sample);

        // run the analysis until completes while reporting progress
        EPP::Analysis<EPP::MATLAB_Sample> *analysis = pursuer.analyze(sample, subset, parameters);
        EPP::Count i = 0;
        while (!analysis->complete())
            if (i < analysis->size())
                if ((*analysis)(i++)->success())
                    std::cerr << "+" << std::flush;
                else
                    std::cerr << "-" << std::flush;
            else
                analysis->wait();
        std::cerr << std::endl;

        // save the gating tree
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
