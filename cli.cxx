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
        std::cout << "Usage: " << argv[0] << " <measurements> <events> <data-csv> [<parameters-json>|default [<gating-json>|- [<taxonomy-json>|none [ <classify-csv>|none [<threads>]]]]]\n";
        return 1;
    }

    try
    {
        // program arguments
        EPP::Measurement measurements = (EPP::Measurement)std::stoi(argv[1]);
        std::vector<std::string> markers(measurements);
        EPP::Event events = std::stol(argv[2]);
        int threads = std::thread::hardware_concurrency();
        if (argc > 8)
            threads = std::stoi(argv[8]);
        if (threads < 0)
            threads = std::thread::hardware_concurrency();

        // get the parameters
        EPP::Parameters parameters = EPP::Default;
        if (argc > 4 && std::strcmp(argv[4], "default"))
        {
            std::ifstream paramfile(argv[4], std::ios::in);
            parameters = json::parse(paramfile);
            paramfile.close();
        };

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

        EPP::MATLAB_Local pursuer(parameters, threads);
        const EPP::MATLAB_Sample sample(measurements, events, data);
        EPP::SampleSubset<EPP::MATLAB_Sample> subset(sample);

        auto analysis = pursuer.analyze(sample, subset, parameters);
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
            std::ofstream gatefile(argv[5], std::ios::out);
            gatefile << subset.gating().dump();
            gatefile.close();
        }
        else
            std::cout << subset.gating().dump(2) << std::endl;

        // save the taxonomy
        if (argc > 6 && std::strcmp(argv[6], "none"))
        {
            std::ofstream taxonfile(argv[6], std::ios::out);
            json taxonomy = (json)*analysis->taxonomy();
            taxonfile << taxonomy.dump();
            taxonfile.close();
        }

        // save the classification and mahalanoabis vectors as one file
        if (argc > 7 && std::strcmp(argv[7], "none"))
        {
            std::ofstream classfile(argv[7], std::ios::out);
            for (EPP::Event event = 0; event < sample.events; ++event)
                classfile << analysis->classification[event] << "," << analysis->mahalanobis[event] << std::endl;
            classfile.close();

            // TODO evolve the next 3 lines of code so that the locations
            // of phenogram.html and template.html can be provided 
            // as additional command line arguments.  This also requires
            // refactoring taxonomy.h Phenogram::toHtml to accept
            // location of template.html.  Without these new arguments
            // these file locations are the invoker's current
            // working directory.   If template.html is not found this process
            // still exits normally and phenogram.html contains only data.
            // This internal dependency was super easy for FlowJoBridge to
            // handle, but other invoker implementations will likely find 
            // its inelegance annoying.
            std::ofstream phenofile("phenogram.html", std::ios::out);
            EPP::Phenogram::toHtml(analysis->taxonomy(), markers, phenofile);
            phenofile.close();

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
