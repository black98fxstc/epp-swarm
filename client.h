#ifndef _EPP_CLIENT_H
#define _EPP_CLIENT_H 1

#include <ios>
#include <chrono>
#include <sstream>
#include <algorithm>
#include <random>
#include <vector>
#include <queue>
#include <chrono>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <unordered_map>
#include <cassert>
#include <cstring>
#include <cmath>

#include "constants.h"
#include <nlohmann/json.hpp>
// using json = nlohmann::json;
typedef void *json;

namespace EPP
{
    typedef std::uint32_t epp_word;
    static std::random_device random;
    static std::mt19937_64 generate(random());

    class Key; // forward references needed
    struct Meta;
    // class Subset;
    class Pursuer;
    // class _Result;
    // typedef std::shared_ptr<const _Result> Result;

    /**
     * Provides a content based associative memory service
     * of blobs shared between clients and servers
     **/

    struct KeyHash
    {
        std::size_t operator()(Key const &key) const noexcept;
    };

    struct Key
    {
        friend class Blob;
        friend class Sample;
        // friend class Subset;
        friend class _Request;
        friend class Remote;

        template <class ClientSample>
        friend class SamplePursuer;

    private:
        union
        { // 256 bit key
            std::uint8_t bytes[32];
            std::uint_fast64_t random[4];
        };

        static std::unordered_map<Key, std::weak_ptr<Meta>, KeyHash> metadata;

        static void vacuum() noexcept;

        std::shared_ptr<Meta> meta() const noexcept;

        // const Key &operator=(const Key &that) const noexcept
        // {
        //     std::move(that.bytes, that.bytes + 32, bytes);
        //     return *this;
        // };

    public:
        explicit operator json() const noexcept;

        bool operator==(
            const Key &other) const noexcept
        {
            return std::memcmp(bytes, other.bytes, 32) == 0;
        }

        Key &operator=(Key &&that) noexcept
        {
            if (!(*this == that))
                std::move(that.bytes, that.bytes + 32, bytes);
            return *this;
        }

        Key &operator=(const json &encoded);

        explicit Key(const json &encoded)
        {
            *this = encoded;
        };

        Key(const Key &key)
        {
            std::move(key.bytes, key.bytes + 32, bytes);
        };

        Key(std::istream &stream);

        Key()
        {
            std::move(no_key, no_key + 32, bytes);
        };
    };

    const Key NoKey;

    struct Meta
    {
        std::iostream *stream = nullptr;
        std::uint8_t *_buffer = nullptr;
        unsigned long int size = 0;
        bool valid = false;
        bool fault = false;

        std::uint8_t *buffer(unsigned long int size)
        {
            if (_buffer == nullptr)
            {
                _buffer = new std::uint8_t[size];
                valid = false;
            }
            return _buffer;
        }

        ~Meta()
        {
            delete[] _buffer;
        }
    };

    /**
     * a blob is identified by it's content. it's key is a hash of it's stream.
     * the contant may be valid or invalid. if it's invalid a blob fault may
     * be generated. the client can then test the validity later or wait for 
     * it to become valid
     **/
    class Blob
    {
        friend class Remote;

    private:
        static std::mutex mutex;
        static std::condition_variable wakeup;

    protected:
        Key _key;
        std::shared_ptr<Meta> meta;

        std::istream istream();

        std::ostream ostream();

        class Handler
        {
        public:
            void startFault(
                Key key){};
        };

        static Handler *handler;

    public:
        bool valid();

        bool fault();

        void wait();

    protected:
        static void content(
            Meta *meta);

        Blob(
            const Key &key);

        Blob();
    };

    /**
     * These classes define the client interface to EPP and the necessary data structures
     * First those needed to specify an analysis request
     */

    class Sample : public Blob
    {
    public:
        const unsigned long int events;
        const unsigned short int measurements;

        const Key &key() noexcept;

        // virtual operator Subset() const noexcept;

        unsigned long int size()
        {
            return sizeof(epp_word) * measurements * events;
        }

        // explicit operator json() noexcept;

        Sample &operator=(const json &encoded);

        // explicit Sample(const json &encoded) : Blob()
        // {
        //     *this = encoded;
        // };

    protected:
        Sample(unsigned short int measurements,
               unsigned long int events,
               Key key = NoKey) noexcept
            : measurements(measurements), events(events), Blob(key){};

        Sample() = default;

    private:
        // these are virtual because our friend stream won't know which variant it will be
        virtual epp_word get_word(unsigned short int measurement, unsigned long event) const noexcept
        {
            return (epp_word)0;
        }
        virtual void put_word(unsigned short measurement, long event, epp_word data) const noexcept {};
        friend class SampleStream;
    };

    class Subset : public std::vector<bool>, public Blob
    {
        friend class SubsetStream;

    private:
    public:
        const Sample *sample;

        Key key();

        operator json() noexcept;

        // Subset& operator=(const Subset &other) = delete;

        // Subset &operator=(const json &encoded);

        // explicit Subset(const json &encoded) : std::vector<bool>(0), Blob()
        // {
        //     *this = encoded;
        // };

        // Subset(
        //     Sample &sample,
        //     Key key = NoKey) : sample(&sample), std::vector<bool>(sample.events, true), Blob(key){};

        // Subset(
        //     Sample *sample,
        //     Key key = NoKey) : sample(sample), std::vector<bool>(sample->events, true), Blob(key){};

        // Subset(
        //     Sample &sample,
        //     std::vector<bool> included,
        //     Key key = NoKey) : sample(&sample), std::vector<bool>(included), Blob(key){};

        // Subset(
        //     const Sample *sample) : sample(sample){};
    protected:
        Subset(std::vector<bool> &subset) : std::vector<bool>(subset){};

        Subset() = default;
    };

    template <class ClientSample>
    class SampleSubset : public Subset
    {
    public:
        ClientSample *sample;

        SampleSubset(
            const ClientSample *sample)
            : sample(sample), Subset(){};

        SampleSubset(
            ClientSample *sample,
            std::vector<bool> &subset)
            : sample(sample), Subset(subset){};
    };

    // template <class ClientSample>
    // class AbstractSample : public Sample
    // {
    // protected:
    //     // inline const double operator()(unsigned long int event, unsigned short int measurement) const noexcept;

    //     AbstractSample(unsigned short int measurements,
    //            unsigned long int events,
    //            Key key = NoKey) noexcept
    //         : Sample(measurements, events, key){};

    // public:
    //     operator SampleSubset<ClientSample>() const noexcept
    //     {
    //         std::vector<bool> in_range(events);
    //         for (long int event = 0; event < events; event++)
    //             for (unsigned short int measurement = 0; measurement < measurements; measurement++)
    //                 // if ((*this)(event, measurement) < 0 || (*this(event, measurement) > 1)
    //                     in_range[event] = false;
    //         SampleSubset<ClientSample> subset(this, in_range);
    //         return this;
    //     };
    // };

    /**
     * EPP is an algorithim not an implementation, and it can be applied to the clients in memory data
     * therefore it's defined as a template based on the users preferred data model
     **/

    template <typename _float>
    class DefaultSample : public Sample
    {
    public:
        inline const double operator()(unsigned long int event, unsigned short int measurement) const noexcept
        {
            return (double)data[Sample::measurements * event + measurement];
        };

        operator SampleSubset<DefaultSample>() const noexcept
        {
            std::vector<bool> in_range(events);
            for (long int event = 0; event < events; event++)
                for (unsigned short int measurement = 0; measurement < measurements; measurement++)
                    if (data[measurements * event + measurement] < 0 || data[measurements * event + measurement] > 1)
                        in_range[event] = false;
            SampleSubset<DefaultSample> subset(this, in_range);
            return this;
        };

        DefaultSample(unsigned short int measurements,
                      unsigned long int events,
                      SampleSubset<DefaultSample> subset,
                      Key key) noexcept
            : Sample(measurements, events, subset, key), data(nullptr){};

        DefaultSample(const unsigned short int measurements,
                      const unsigned long int events,
                      const _float *const data,
                      SampleSubset<DefaultSample> subset) noexcept
            : Sample(measurements, events, subset, NoKey), data(data){};

        explicit operator json() const noexcept
        {
            return nullptr;
        };

        DefaultSample &operator=(const json &encoded)
        {
            return *this;
        }

        DefaultSample(const json &encoded) : data(nullptr)
        {
            *this = encoded;
        }

    protected:
        epp_word get_word(unsigned short int measurement, unsigned long int event) const noexcept
        {
            float f = data[Sample::measurements * event + measurement];
            return *(epp_word *)&f;
        };

        void put_word(unsigned short int measurement, unsigned long int event, epp_word value) noexcept
        {
            float f = *(float *)&value;
            data[Sample::measurements * event + measurement] = (_float)f;
        };

    private:
        const _float *const data;
    };

    template <typename _float>
    class TransposeSample : public Sample
    {
    public:
        inline double operator()(unsigned long int event, unsigned short int measurement) const noexcept
        {
            return (double)data[Sample::events * measurement + event];
        };

        operator SampleSubset<TransposeSample>()
        {
            std::vector<bool> in_range(events, true);
            for (unsigned long int event = 0; event < events; event++)
                for (unsigned short int measurement = 0; measurement < measurements; measurement++)
                    if (data[events * measurement + event] < 0 || data[events * measurement + event] > 1)
                        in_range[event] = false;
            SampleSubset<TransposeSample> *subset = new SampleSubset<TransposeSample>(this, in_range);
            return *subset;
        };

        TransposeSample(unsigned short int measurements,
                        unsigned long int events,
                        SampleSubset<TransposeSample> subset,
                        Key key) noexcept
            : Sample(measurements, events, subset, key), data(nullptr){};

        TransposeSample(const unsigned short int measurements,
                        const unsigned long int events,
                        const _float *const data) noexcept
            : Sample(measurements, events, NoKey), data(data){};

    protected:
        epp_word get_word(unsigned short int measurement, unsigned long int event) const noexcept
        {
            float f = data[Sample::events * measurement + event];
            return *(epp_word *)&f;
        };

        void put_word(unsigned short int measurement, unsigned long int event, epp_word value) noexcept
        {
            float f = *(float *)&value;
            data[Sample::events * measurement + event] = (_float)f;
        };

    private:
        const _float *const data;
    };

    template <typename _float>
    class PointerSample : public Sample
    {
    public:
        inline double operator()(unsigned long int event, unsigned short int measurement) const noexcept
        {
            return (double)data[measurement][event];
        };

        operator SampleSubset<PointerSample>()
        {
            std::vector<bool> in_range(events);
            for (unsigned long int event = 0; event < events; event++)
                for (unsigned short int measurement = 0; measurement < measurements; measurement++)
                    if (data[measurement][event] < 0 || data[measurement][event] > 1)
                        in_range[event] = false;
            SampleSubset<PointerSample> subset(this, in_range);
            return this;
        };

        PointerSample(unsigned short int measurements,
                      unsigned long int events,
                      SampleSubset<PointerSample> subset,
                      Key key) noexcept
            : Sample(measurements, events, subset, key), data(nullptr){};

        PointerSample(const unsigned short int measurements,
                      const unsigned long int events,
                      const _float *const *const data,
                      SampleSubset<PointerSample> subset) noexcept
            : Sample(measurements, events, subset, NoKey), data(data){};

    protected:
        epp_word get_word(unsigned short int measurement, unsigned long int event) const noexcept
        {
            float f = data[measurement][event];
            return *(epp_word *)&f;
        };

        void put_word(unsigned short int measurement, unsigned long int event, epp_word value) noexcept
        {
            float f = *(float *)&value;
            data[measurement][event] = (_float)f;
        };

    private:
        const _float *const *const data;
    };

    struct Parameters
    {
        // N = 256 gives points and features a precision of roughly two significant figures

        static const unsigned short N; // resolution of points and boundaries
                                       // optimized when there are lots of small factors

        double W = 1 / (double)N; // standard deviation of kernel,
                                  // this is the highest achievable resolution, in practice a higher
                                  // value might be used for application reasons or just performance

        double sigma = 5; // controls the density threshold for starting a new cluster

        enum Goal
        {                    // which objective function
            best_separation, // lowest edge weight
            best_balance     // edge weight biased towards more even splits
        } goal = best_balance;

        int finalists = 1; // remember this many of the best candidates

        struct KLD // KLD threshold for informative cases
        {
            double Normal2D = .16;      // is this population worth splitting?
            double Normal1D = .16;      // is the measurement just normal
            double Exponential1D = .16; // is this an exponential tail (CyToF)

            KLD(
                double Normal2D = .16,
                double Normal1D = .16,
                double Exponential1D = .16)
            noexcept
                : Normal2D(Normal2D), Normal1D(Normal1D), Exponential1D(Exponential1D){};
        };

        // const static KLD KLD_Default;
        KLD kld{.16, .16, .16};

        std::vector<bool> censor; // omit measurements from consideration

        // algorithm tweaks

        unsigned int max_clusters = 12; // most clusters the graph logic should handle

        bool suppress_in_out = false; // don't bother with in and out sets

        explicit operator json() const noexcept;

        Parameters &operator=(const json &encoded);

        Parameters(const json &encoded)
        {
            *this = encoded;
        };

        Parameters(
            Goal goal = best_balance,
            KLD kld = {.16, .16, .16},
            double sigma = 4,
            double W = 1 / (double)N)
            : goal(goal), kld(kld), W(W), sigma(sigma),
              censor(0), finalists(1), max_clusters(12),
              suppress_in_out(false){};
    };

    const Parameters Default;

    /**
     * a request is handled asynbronously and possibly remotely
     * a client can test whether it is finished and wait for it to finish
     * actually returned as a shared_ptr
     **/
    class _Request;
    class Result;
    class Request : protected std::shared_ptr<_Request>
    {
        friend class Pursuer;
        friend class MATLAB_Pursuer;
        friend class MATLAB_Local;
        friend class MATLAB_Remote;
        friend class CloudPursuer;

        template <class ClientSample>
        friend class Work;

    private:
    protected:
        const Key key() const noexcept;

        void finish();

        Request &operator++();

        Request &operator--();

        Request(
            Pursuer *pursuer,
            Parameters parameters);

        // Request(_Request *request) : std::shared_ptr<_Request>(request)
        // {
        //     request->_finished = false;
        // };

    public:
        bool finished() const noexcept;

        void wait() const noexcept;

        Result result() const noexcept;

        Request(
            _Request *request) : std::shared_ptr<_Request>(request){};

        // _Result *working()
        // {
        //     return nullptr;
        //     // return (*this)->working_result;
        // }
    };

    class _Result;
    class Result : public std::shared_ptr<const _Result>
    {
    public:
        Result(
            const _Result *result) : std::shared_ptr<const _Result>(result){};

        Result() = default;
    };

    /**
     * structures defining the result of an ananlysis
     **/
    enum Status
    {
        EPP_success,
        EPP_no_qualified,
        EPP_no_cluster,
        EPP_not_interesting,
        EPP_error
    };

    struct Point
    {
        short i, j;

        inline double x() const noexcept { return (double)i / (double)Parameters::N; };
        inline double y() const noexcept { return (double)j / (double)Parameters::N; };

        inline bool operator==(const Point &other) const noexcept
        {
            return i == other.i && j == other.j;
        }

        inline bool operator!=(const Point &other) const noexcept
        {
            return !(*this == other);
        }

        Point(short i, short j) noexcept : i(i), j(j){};

        Point() : i(0), j(0){};
    };

    typedef std::vector<Point> Polygon;

    struct Candidate
    {
        Polygon separatrix;
        // Subset in, out;
        double score, edge_weight, balance_factor;
        unsigned long int in_events, out_events;
        unsigned int pass, clusters, graphs;
        unsigned short int X, Y;
        enum Status outcome;

    private:
        static void close_clockwise(
            Polygon &polygon) noexcept;

        void simplify(
            const double tolerance,
            Polygon &simplified,
            const unsigned short int lo,
            const unsigned short int hi) const noexcept;

    public:
        bool operator<(const Candidate &other) const noexcept
        {
            if (score < other.score)
                return true;
            if (score > other.score)
                return false;
            return outcome < other.outcome;
        }

        Polygon simplify(
            const double tolerance) const noexcept;

        Polygon in_polygon() const noexcept;

        Polygon in_polygon(
            double tolerance) const noexcept;

        Polygon out_polygon() const noexcept;

        Polygon out_polygon(
            double tolerance) const noexcept;

        operator json();

        Candidate &operator=(const json &encoded);

        // Candidate(const json &encoded)
        // {
        //     *this = encoded;
        // }

        Candidate(
            unsigned short int X,
            unsigned short int Y)
            : X(X < Y ? X : Y), Y(X < Y ? Y : X),
              outcome(Status::EPP_error),
              score(std::numeric_limits<double>::infinity()),
              pass(0), clusters(0), graphs(0){};

        Candidate() : Candidate(0, 1){};
    };

    struct _Result
    {
        const Key &key;
        const Parameters &parameters;
        std::vector<Candidate> candidates;
        std::vector<unsigned short int> qualified;
        std::chrono::milliseconds milliseconds;
        unsigned int projections, passes, clusters, graphs;

        const Candidate &winner() const noexcept
        {
            return candidates[0];
        }

        enum Status outcome() const noexcept
        {
            return winner().outcome;
        };

        explicit operator json();

        // _Result &operator=(const json &encoded);

        // explicit _Result(const json &encoded)
        // {
        //     *this = encoded;
        // };

        explicit _Result(
            const Key &key,
            const Parameters &parameters)
            : key(key), parameters(parameters), projections(0),
              passes(0), clusters(0), graphs(0)
        {
            candidates.reserve(parameters.finalists);
        };
    };

    /**
     * Classes defining the actual clients
     **/
    class Pursuer
    {
        friend class Request;
        friend class Remote;

        template <class ClientSample>
        friend class Work;
    
        template <class ClientSample>
        friend class SamplePursuer;

    protected:
        std::unordered_map<const Key, _Request *, KeyHash> requests;
        std::vector<std::thread> workers;
        std::mutex mutex;

        void start(
            _Request *request) noexcept;

        void increment(
            _Request *request) noexcept;

        void decrement(
            _Request *request) noexcept;

        void finish(
            _Request *request) noexcept;

        void wait(
            _Request *request) noexcept;

        void start(
            const json &encoded){};

        void finish(
            const json &encoded);

        Pursuer(
            const Parameters &parameters) noexcept
            : parameters(parameters){};

        Pursuer(
            const Parameters &parameters,
            int threads) noexcept
            : parameters(parameters), workers(threads){};

    public:
        const Parameters &parameters;

        bool finished()
        {
            std::unique_lock<std::mutex> lock(mutex);
            for (auto it = requests.begin(); it != requests.end(); it++)
            {
                Request request = it->second;
                if (request.finished())
                    return true;
            }
            return false;
        };

        // void wait()
        // {
        //     std::unique_lock<std::mutex> lock(mutex);
        //     if (!finished())
        //         completed.wait(lock);
        // }

        // void waitAll()
        // {
        //     std::unique_lock<std::mutex> lock(mutex);
        //     if (!requests.empty())
        //         completed.wait(lock);
        // }
    };

    template <class ClientSample>
    class ClientRequest;

    template <class ClientSample>
    class SamplePursuer : public Pursuer
    {
    public:
        // this is the virtual function that is required for every client
        // it's the only one called by Pursuer but than means it must be virtual
        // because we can't link templates we have to instantiate them individually
        // the others are convenience methods for specific clients
        virtual ClientRequest<ClientSample> *start(
            const SampleSubset<ClientSample> *subset,
            const Parameters &parameters) noexcept
        {
            Key key;
            for (auto &random_bits : key.random)
                random_bits = generate();
            _Result *result = new _Result(key, parameters);
            ClientRequest<ClientSample> *request = new ClientRequest<ClientSample>(this, key, subset, parameters, result);

            SamplePursuer<ClientSample>::start(request);

            return request;
        };

        // ClientRequest<ClientSample> *start(
        //     const ClientSample &sample,
        //     const Parameters &parameters) noexcept
        //     {
        //         SampleSubset<ClientSample> subset(&sample);
        //         ClientRequest<ClientSample> *request = SamplePursuer<ClientSample>::start(subset, parameters);
        //         return request;
        //     };

        // default is to use the pursuer's parameters
        ClientRequest<ClientSample> *start(
            const ClientSample &sample) noexcept
        {
            return start(sample, parameters);
        };

        void start(
            ClientRequest<ClientSample> *request) noexcept;

        // Result pursue(
        //     const ClientSample &sample) noexcept
        // {
        //     return start(sample, parameters).result();
        // };

        // Request start(
        //     _Request *request) noexcept;

        // void finish(
        //     _Request<ClientSample> *request) noexcept;
        // ;

    protected:
        SamplePursuer() noexcept = default;

        SamplePursuer(
            const Parameters &parameters,
            int threads)
            : Pursuer(parameters, threads){};

        SamplePursuer(
            const Parameters &parameters = Default)
            : Pursuer(parameters){};

        ~SamplePursuer() = default;
        ;
    };

    class _Request
    {
        friend class Request;
        friend class Pursuer;

        template <class ClientSample>
        friend class Work;

         template <class ClientSample>
        class ClientRequest;

        template <class ClientSample>
        friend class QualifyMeasurement;

        template <class ClientSample>
        friend class PursueProjection;

    private:
        std::chrono::time_point<std::chrono::steady_clock> begin, end;

    protected:
        const Key _key;
        Pursuer *const pursuer;
        _Result *const _result;
        const Result result;

        std::condition_variable completed;
        volatile unsigned int outstanding = 0;
        volatile bool _finished;

    public:
        const Parameters &parameters;

        _Request(
            Pursuer *const pursuer,
            const Key &key,
            const Parameters &parameters,
            _Result *result)
            : pursuer(pursuer), _key(key), parameters(parameters), _result(result), result(result){};

        _Request() = default;
    };

    template <class ClientSample>
    class ClientRequest : public _Request
    {
        friend class Request;
        friend class SamplePursuer<ClientSample>;

    protected:
    public:
        const SampleSubset<ClientSample> *subset;

        ClientRequest(
            SamplePursuer<ClientSample> *const pursuer,
            const Key &key,
            const SampleSubset<ClientSample> *subset,
            const Parameters &parameters,
            _Result *const result) noexcept
            : _Request(pursuer, key, parameters, result), subset(subset){};

        explicit operator json() const noexcept;

        _Request &operator=(const json &encoded);

        // _Request(const json &encoded)
        // {
        //     *this = encoded;
        // };
    };

    typedef TransposeSample<float> MATLAB_Sample;

    class MATLAB_Pursuer : public SamplePursuer<MATLAB_Sample>
    {
    public:
        // Request start( // this one does the heavy lifting
        //     MATLAB_Sample &sample,
        //     const Parameters &parameters) noexcept
        //     {
        //         SampleSubset<MATLAB_Sample> subset = sample;
        //         ClientRequest<MATLAB_Sample> *request = SamplePursuer<MATLAB_Sample>::start(subset, parameters);
        //         Request shared(request);
        //         return shared;
        //     };

        // _Request *start( // thiese are convenience routines
        //     const unsigned short int measurements,
        //     const unsigned long int events,
        //     const float *const data,
        //     SampleSubset &subset) noexcept;
        // Request start(
        //     const unsigned short int measurements,
        //     const unsigned long int events,
        //     const float *const data) noexcept;
        // Result pursue(
        //     const unsigned short int measurements,
        //     const unsigned long int events,
        //     const float *const data,
        //     SampleSubset &subset) noexcept;
        // Result pursue(
        //     const unsigned short int measurements,
        //     const unsigned long int events,
        //     const float *const data) noexcept;

    protected:
        MATLAB_Pursuer() = delete;

        MATLAB_Pursuer(
            Parameters parameters,
            int threads) noexcept
            : SamplePursuer<MATLAB_Sample>(parameters, threads){};
        ~MATLAB_Pursuer() = default;
    };

    class MATLAB_Local : public MATLAB_Pursuer
    {
    public:
        int getThreads() const noexcept
        {
            return workers.size();
        };

        MATLAB_Local(
            Parameters parameters = Default,
            int threads = std::thread::hardware_concurrency()) noexcept;
        ~MATLAB_Local();
    };

    /**
     * communications services assuming control and data streams
     * are already set up. blob faults are handled here.
     * requests and results are safely queued for the pursuer
     **/

    class Remote : Blob::Handler
    {
        friend class Blob;

    public:
        enum Service
        {
            request,
            result,
            fault,
            content
        };

        // thread safe send/receive one json message

        void out(const json &encoded); // does not block

        json in(); // blocks calling thread

    protected:
        void copy(
            std::istream *in,
            std::ostream *out,
            unsigned long int count);

        void startFault(
            Key key);

        void transmit();

        void receive();

        Remote();

        ~Remote();

    private:
        std::mutex mutex;
        std::queue<json> incoming;
        std::queue<json> outgoing;
        std::condition_variable wake_in;
        std::condition_variable wake_out;
        std::iostream *remote_control;
        std::iostream *remote_data;
        std::thread receiver;
        std::thread transmitter;
        volatile bool on_the_air = true;
    };

    class MATLAB_Remote : Remote, public MATLAB_Pursuer
    {
    public:
        // Request start(
        //     MATLAB_Sample *sample,
        //     const Parameters &parameters) noexcept;

        void finish(
            const json &encoded);

        MATLAB_Remote(
            const Parameters &parameters) noexcept;

        ~MATLAB_Remote();
    };
    /**
     * remote worker instance
     **/
    typedef DefaultSample<float> CloudSample;

    class CloudPursuer : Remote, public SamplePursuer<CloudSample>
    {
    public:
        ClientRequest<CloudSample> *start(
            CloudSample *sample,
            const Parameters &parameters) noexcept;

        void start(const json &encoded);

        void finish(ClientRequest<CloudSample> *request) noexcept;

        void finish(const json &encoded);

        json remote();

        CloudPursuer(
            const Parameters &parameters) noexcept;

        ~CloudPursuer();
    };

    /**
     * utility classes for serializing sample and subset as streams
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