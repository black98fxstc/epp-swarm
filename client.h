#ifndef _EPP_CLIENT_H
#define _EPP_CLIENT_H 1

#include <ios>
#include <sstream>
#include <vector>
#include <queue>
#include <chrono>
#include <thread>

namespace EPP
{
    typedef uint32_t epp_word;
    typedef unsigned char epp_hash[32];

    /**
     * These classes define the client interface to EPP and the necessary data structures
     * 
     */

    class Sample
    {
    public:
        std::vector<bool> subset;
        const long events;
        const int measurements;

        Sample(int measurements,
               long events,
               std::vector<bool> subset) noexcept
            : measurements(measurements), events(events), subset(subset){};

        Sample(int measurements,
               long events) noexcept
            : measurements(measurements), events(events), subset(events)
        {
            std::fill(subset.begin(), subset.end(), true);
        };

    private:
        // these are virtual because our friend stream won't know which variant it will be
        virtual epp_word get_word(int measurement, long event) const noexcept
        {
            return (epp_word)0;
        }
        virtual void put_word(int measurement, long event, epp_word data) noexcept {};
        friend class SampleStream;
    };

    template <typename _float>
    class DefaultSample : public Sample
    {
    public:
        inline double operator()(long event, int measurement) const noexcept
        {
            return (double)data[measurements * event + measurement];
        };

        DefaultSample(const int measurements,
                      const long events,
                      _float *data) noexcept
            : Sample(measurements, events), data(data)
        {
            for (long event = 0; event < events; event++)
                for (int measurement = 0; measurement < measurements; measurement++)
                    if (data[measurements * event + measurement] < 0 || data[measurements * event + measurement] > 1)
                        subset[event] = false;
        };

        DefaultSample(const int measurements,
                      const long events,
                      _float *data,
                      std::vector<bool> subset) noexcept
            : Sample(measurements, events, subset), data(data){};

    protected:
        epp_word get_word(int measurement, long event) const noexcept
        {
            float f = data[measurements * event + measurement];
            return *(epp_word *)&f;
        };

        void put_word(int measurement, long event, epp_word value) noexcept
        {
            float f = *(float *)&value;
            data[measurements * event + measurement] = (_float)f;
        };

    private:
        const _float *data;
    };

    template <typename _float>
    class TransposeSample : public Sample
    {
    public:
        inline double operator()(long event, int measurement) const noexcept
        {
            return (double)data[events * measurement + event];
        };

        TransposeSample(const int measurements,
                        const long events,
                        _float *data) noexcept
            : Sample(measurements, events), data(data)
        {
            for (long event = 0; event < events; event++)
                for (int measurement = 0; measurement < measurements; measurement++)
                    if (data[events * measurement + event] < 0 || data[events * measurement + event] > 1)
                        subset[event] = false;
        };

        TransposeSample(const int measurements,
                        const long events,
                        _float *data,
                        std::vector<bool> subset) noexcept
            : Sample(measurements, events, subset), data(data){};

    protected:
        epp_word get_word(int measurement, long event) const noexcept
        {
            float f = data[events * measurement + event];
            return *(epp_word *)&f;
        };

        void put_word(int measurement, long event, epp_word value) noexcept
        {
            float f = *(float *)&value;
            data[events * measurement + event] = (_float)f;
        };

    private:
        _float *data;
    };

    template <typename _float>
    class PointerSample : public Sample
    {
    public:
        inline double operator()(long event, int measurement) const noexcept
        {
            return (double)data[measurement][event];
        };

        PointerSample(const int measurements,
                      const long events,
                      const _float **const data) noexcept
            : Sample(measurements, events), data(data)
        {
            for (long event = 0; event < events; event++)
                for (int measurement = 0; measurement < measurements; measurement++)
                    if (data[measurement][event] < 0 || data[measurement][event] > 1)
                        subset[event] = false;
        };

        PointerSample(const int measurements,
                      const long events,
                      _float *data,
                      std::vector<bool> subset) noexcept
            : Sample(measurements, events, subset), data(data){};

    protected:
        epp_word get_word(int measurement, long event) const noexcept
        {
            float f = data[measurement][event];
            return *(epp_word *)&f;
        };

        void put_word(int measurement, long event, epp_word value) noexcept
        {
            float f = *(float *)&value;
            data[measurement][event] = (_float)f;
        };

    private:
        const _float **data;
    };

    class Subset : public std::vector<bool>
    {
    public:
        explicit Subset(Sample &sample);
        Sample *sample;

    private:
        friend class SubsetStream;
    };

    struct Parameters
    {
        static const int N = 1 << 8; // resolution of points and boundaries
                                     // optimized when there are lots of small factors
        double W = .005;             // standard deviation of kernel, this combination gives
                                     // points and features a precision of two significant figures
        enum Goal
        {
            best_separation, // which objective function
            best_balance
        } goal = best_balance;

        int max_clusters = 12;

        int finalists = 1; // remember this many of the best candidates

        struct KLD
        {
            double Normal2D = .16; // KLD threshold for informative cases
            double Normal1D = .16;
            double Exponential1D = .16;
        } kld;

        std::vector<bool> censor; // omit measurments from consideration

                                   // these control the density threshold for a new cluster
        double sigma = 5;          // significance of difference from zero, probably to high
        double A = pi * 4 * W * W; // area of threshold spot, equivalent to */-2W

        Parameters(
            Goal goal = best_balance,
            KLD kld = {.16, .16, .16},
            double W = .005,
            double sigma = 5,
            double A = pi * 4 * .005 * .005)
            : goal(goal), kld(kld), W(W), sigma(sigma), A(pi * 4 * W * W), 
            censor(0), finalists(1), max_clusters(12){};
    };

    const Parameters Default;

    struct Point
    {
        short i, j;

        inline double x() const noexcept { return (double)i / (double)Parameters::N; };
        inline double y() const noexcept { return (double)j / (double)Parameters::N; };

        Point(short i, short j) noexcept : i(i), j(j){};
    };

    enum Status
    {
        EPP_success,
        EPP_no_qualified,
        EPP_no_cluster,
        EPP_not_interesting,
        EPP_error
    };

    struct Candidate
    {
        std::vector<Point> separatrix;
        std::vector<bool> in, out;
        double score, edge_weight, balance_factor;
        long in_events, out_events;
        int X, Y, pass, clusters, graphs;
        enum Status outcome;

        bool operator<(const Candidate &other) const noexcept
        {
            return score < other.score;
        }

        Candidate(
            const int X,
            const int Y)
            : X(X < Y ? X : Y), Y(X < Y ? Y : X),
              outcome(Status::EPP_error),
              score(std::numeric_limits<double>::infinity()),
              pass(0), clusters(0), graphs(0){};
    };

    struct Result
    {
        std::vector<Candidate> candidates;
        std::vector<short> qualified;
        std::chrono::milliseconds milliseconds;
        int projections, passes, clusters, graphs;

        Candidate winner() const noexcept
        {
            return candidates[0];
        }

        enum Status outcome ()
        {
            return winner().outcome;
        };

    protected:
        std::chrono::time_point<std::chrono::steady_clock> begin, end;
        friend class MATLAB_Pursuer;
    };

    template <class ClientSample>
    class Pursuer
    {
    public:
        void start(
            const ClientSample sample,
            const Parameters parameters) noexcept;
        void start(
            const ClientSample sample) noexcept;
        bool finished() noexcept;
        void wait() noexcept;
        std::shared_ptr<Result> result() noexcept;
        std::shared_ptr<Result> pursue(
            const ClientSample sample,
            const Parameters parameters) noexcept;
        std::shared_ptr<Result> pursue(
            const ClientSample sample) noexcept;

    protected:
        Pursuer() noexcept = default;
        ~Pursuer() = default;
    };

    typedef TransposeSample<float> MATLAB_Sample;

    class MATLAB_Pursuer : public Pursuer<MATLAB_Sample>
    {
        int threads;
        std::thread *workers;

    public:
        void start(
            const MATLAB_Sample sample,
            const Parameters parameters) noexcept;
        void start(
            const MATLAB_Sample sample) noexcept;
        void start(
            const int measurements,
            const long events,
            float *data,
            std::vector<bool> &subset) noexcept;
        void start(
            const int measurements,
            const long events,
            float *data) noexcept;
        bool finished() noexcept;
        void wait() noexcept;
        std::shared_ptr<Result> result() noexcept;
        std::shared_ptr<Result> pursue(
            const MATLAB_Sample sample,
            const Parameters parameters) noexcept;
        std::shared_ptr<Result> pursue(
            const MATLAB_Sample sample) noexcept;
        std::shared_ptr<Result> pursue(
            const int measurements,
            const long events,
            float *data,
            std::vector<bool> &subset) noexcept;
        std::shared_ptr<Result> pursue(
            const int measurements,
            const long events,
            float *data) noexcept;
        MATLAB_Pursuer() noexcept;
        MATLAB_Pursuer(int threads) noexcept;
        ~MATLAB_Pursuer();
    };

    /**
     * utility classes for searializing sample and subset as streams
     */

    class SampleStream : public std::iostream
    {
    protected:
        class sample_buffer : public std::streambuf
        {

        public:
            explicit sample_buffer(Sample &sample);
            virtual ~sample_buffer();

        protected:
            virtual std::streambuf::int_type underflow();
            virtual std::streambuf::int_type overflow(std::streambuf::int_type value);
            virtual std::streambuf::int_type sync();

        private:
            Sample *sample;
            epp_word *buffer;
            long next_event;
        };

    public:
        explicit SampleStream(Sample &sample);
    };

    class SubsetStream : public std::iostream
    {
    protected:
        class subset_buffer : public std::streambuf
        {
        public:
            explicit subset_buffer(Subset &subset);
            virtual ~subset_buffer();
            virtual std::streambuf::int_type underflow();
            virtual std::streambuf::int_type overflow(std::streambuf::int_type value);
            virtual std::streambuf::int_type sync();

        private:
            Subset *subset;
            uint8_t *buffer;
            long next_event;
            friend class SubsetStream;
        };

    public:
        explicit SubsetStream(Subset &subset);
    };
}
#endif /* _EPP_CLIENT_H */