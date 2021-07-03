#include <iostream>
#include <fstream>
#include <sstream>

#include "constants.h"
#include "boundary.h"
#include "modal.h"
#include "pursuit.h"
#include "pursuer.h"

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
        long events = std::stol(argv[2]);
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
        for (long i = 0; i < events; i++)
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

        // initial subset is everything in range
        // there is no bounds checking in the algorithm
        std::vector<bool> start(events);
        for (long i = 0; i < events; i++)
        {
            bool in_range = true;
            for (int j = 0; j < measurements; j++)
            {
                float value = data[i + events * j];
                if (value < 0)
                    in_range = false;
                if (value > 1)
                    in_range = false;
            }
            start[i] = in_range;
        }

        EPP::Pursuer pursuer(threads);  // reusable, you can do many start/result calls
        pursuer.start(measurements, events, data, start);
        if (!pursuer.finished()) // optional, used when you want to do something else while it runs
            pursuer.wait();
        std::shared_ptr<EPP::Result> result = pursuer.result();

        if (result->outcome != EPP::Result::EPP_success)
            std::cout << "oops" << std::endl;
        else
        {
            std::cout << "pass " << result->pass << " clusters " << result->clusters 
                << " graphs " << result->graphs << " ms " << result->milliseconds.count() << std::endl;
            std::cout << "best score " << result->X << " " << result->Y << "  " << result->best_score << std::endl;
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