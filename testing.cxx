#include <iostream>
#include <fstream>
#include <fftw3.h>

#include <work.h>

#include "pursuit.cxx"

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " endpoint\n";
        return 1;
    }

    try
    {
        int measurments = 10;
        long events = 55000;
        int threads = 1;

        // get some data from somewhere? CSV?
        float data[measurments * events];
        std::ifstream datafile("/home/wmoore/Downloads/sample55k.csv", std::ios::in);
        std::string line;
        std::string value;
        std::getline(datafile, line);
        for (long i = 0; i < events; i++)
        {
            std::getline(datafile, line);
            std::stringstream sstr(line, std::ios::in);
            for (int j = 0; j < measurments; j++)
            {
                std::getline(sstr, value, ',');
                data[i * measurments + j] = std::stof(value);
            }
        }
        datafile.close();

        // start some worker threads
        std::thread workers[threads];
        for (int i = 0; i < threads; i++)
            workers[i] = std::thread([]()
                                     { EPP::Worker worker; });

        while (!EPP::kiss_of_death)
        {
            EPP::DefaultSample<float> sample(measurments, events, data);
            EPP::Subset start(sample);
            for (long i = 0; i < events; i++)
            {
                bool in_range = true;
                for (int j = 0; j < measurments; j++)
                {
                    float value = data[i * measurments + j];
                    if (value < 0)
                        in_range = false;
                    if (value > 1)
                        in_range = false;
                }
                start[i] = in_range;
            }

            // start parallel projection pursuit
            {
//                EPP::worker_sample constants{measurments, events, (const float *const)data, start};
//                EPP::qualified_measurements.clear();
				EPP::PursueProjection::start(sample, data, start);
//                std::unique_lock<std::recursive_mutex> lock(EPP::mutex);
//                for (int measurment = 0; measurment < constants.measurements; ++measurment)
//                    EPP::work_list.push(new EPP::QualifyMeasurement(constants, measurment));
//                EPP::work_available.notify_all();
            }

            // wait for everything to finish
            {
                std::unique_lock<std::recursive_mutex> lock(EPP::mutex);
                while (EPP::work_outstanding)
                    EPP::work_completed.wait(lock);
            }

            // presumably report back to the dispatcher

            // only go around once for now
            EPP::kiss_of_death = true;
        }

        // tell the workers to exit and wait for them to shut down
        EPP::kiss_of_death = true;
        {
            std::unique_lock<std::recursive_mutex> lock(EPP::mutex);
            EPP::work_available.notify_all();
        }
        for (int i = 0; i < threads; i++)
            workers[i].join();
    }
    catch (std::runtime_error e)
    {
        std::cout << e.what() << std::endl;
        return 1;
    }

    return 0;
}