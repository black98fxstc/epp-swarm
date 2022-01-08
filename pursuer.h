#ifndef _EPP_PURSUER_H
#define _EPP_PURSUER_H 1

#include <vector>
#include <chrono>
#include <thread>

namespace EPP
{
    // // in process client
    // struct Point
    // {
    //     short i, j;

    //     Point (short i, short j) noexcept : i(i), j(j) {};
    // };

    struct Result
    {
        std::vector<int> qualified;
        std::vector<bool> in, out;
        std::vector<Point> separatrix;
        std::chrono::time_point<std::chrono::steady_clock> begin, end;
        std::chrono::milliseconds milliseconds;
        double edge_weight, balance_factor, best_score;
        long in_events, out_events;
        int X, Y;
        int pass, clusters, graphs;
        enum Status
        {
            EPP_success,
            EPP_no_qualified,
            EPP_no_cluster,
            EPP_not_interesting,
            EPP_error
        } outcome;
    };

    class Pursuer
    {
        int threads;
        std::thread *workers;

    public:
        void start(
            const int measurements, 
            const long events, 
            const float *const data, 
            std::vector<bool> &subset) noexcept;
        bool finished() noexcept;
        void wait() noexcept;
        std::shared_ptr<Result> result() noexcept;
        std::shared_ptr<Result> pursue(
            const int measurements, 
            const long events, 
            const float *const data, 
            std::vector<bool> &subset) noexcept;
        Pursuer() noexcept;
        Pursuer(int threads) noexcept;
        ~Pursuer();
    };
}
#endif /* _EPP_PURSUER_H */