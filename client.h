#ifndef _EPP_CLIENT_H
#define _EPP_CLIENT_H 1

#include <ios>
#include <sstream>
#include <vector>
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
        inline double operator()(long event, int measurment) const noexcept
        {
            return (double)data[measurements * event + measurment];
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
        epp_word get_word(int measurement, long event)
        {
            float f = data[measurements * event + measurement];
            return *(epp_word *)&f;
        };

        void put_word(int measurement, long event, epp_word value)
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
        epp_word get_word(int measurement, long event)
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
        epp_word get_word(int measurement, long event)
        {
            float f = data[measurement][event];
            return *(epp_word *)&f;
        };

        void put_word(int measurement, long event, epp_word value)
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
        static const int N = 1 << 8;    // resolution of points and boundaries
                                        // optimized when there are lots of small factors
        enum Goal
        {
            best_separation,            // which objective function
            best_balance
        } goal;
        
        struct KLD
        {
            double Normal2D = .16,      // KLD tests for significance 
            Normal1D = .16, 
            Exponential1D = .16;
        } kld;

        std::vector<bool> censor;

        double W = .01;                 // standard deviation of kernel
        double sigma = 5;               // significance of difference from zero, probably to high
        double A = pi * W * W;          // area of threshold spot, equivalent to */-W probably to low

        Parameters(
            Goal goal = best_balance,
            KLD kld = {.16, .16, .16},
            double W = .01,
            double sigma = 5,
            double A = pi * .01 * .01)
            : goal(goal), kld(kld), W(W), sigma(sigma), A(pi * W * W), censor(0) {};
    };

    const Parameters Default;

    struct Point
    {
        short i, j;

        inline double x() { return (double) i / (double) Parameters::N; };
        inline double y() { return (double) j / (double) Parameters::N; };

        Point(short i, short j) noexcept : i(i), j(j){};
    };

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
        int projections, passes, clusters, graphs;
        enum Status
        {
            EPP_success,
            EPP_no_qualified,
            EPP_no_cluster,
            EPP_not_interesting,
            EPP_error
        } outcome;
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