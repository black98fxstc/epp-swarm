#ifndef _EPP_CLIENT_H
#define _EPP_CLIENT_H 1

#include <ios>
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
// #include <nlohmann/json.hpp>

namespace EPP
{
    // using json = nlohmann::json;
    typedef void *json;
    typedef std::uint32_t epp_word;
    static std::random_device random;

    class Key; // forward references needed
    struct Meta;
    class Pursuer;
    class _Result;
    typedef std::shared_ptr<_Result> Result;

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
        friend class _Request;
        friend class Remote;

    private:
        union
        { // 256 bit key
            std::uint8_t bytes[32];
            std::uint_fast64_t random[4];
        };

        static std::unordered_map<Key, std::weak_ptr<Meta>, KeyHash> metadata;

        static void vacuum() noexcept;

        std::shared_ptr<Meta> meta() const noexcept;

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

        Key &operator=(const Key &that)
        {
            std::move(that.bytes, that.bytes + 32, bytes);
            return *this;
        };

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

    class Subset : public std::vector<bool>, public Blob
    {
        friend class SubsetStream;

    public:
        Key key();

        Subset(
            unsigned long int events,
            Key key = NoKey) : std::vector<bool>(events, true), Blob(key){};

        Subset() = default;
    };

    class Sample : public Blob
    {
    public:
        Subset subset;
        unsigned long int events;
        unsigned short int measurements;

        Key key();

        unsigned long int size()
        {
            return sizeof(epp_word) * measurements * events;
        }

        explicit operator json() const noexcept;

        Sample &operator=(const json &encoded);

        explicit Sample(const json &encoded) : subset(0), Blob()
        {
            *this = encoded;
        };

        Sample(unsigned short int measurements,
               unsigned long int events,
               Subset subset,
               Key key = NoKey) noexcept
            : measurements(measurements), events(events), subset(subset), Blob(key){};

        Sample(unsigned short int measurements,
               unsigned long int events,
               Key key = NoKey) noexcept
            : measurements(measurements), events(events), subset(events), Blob(key){};

    protected:
        Sample() = default;

    private:
        // these are virtual because our friend stream won't know which variant it will be
        virtual epp_word get_word(unsigned short int measurement, unsigned long event) const noexcept
        {
            return (epp_word)0;
        }
        virtual void put_word(unsigned short measurement, long event, epp_word data) noexcept {};
        friend class SampleStream;
    };

    /**
     * EPP is an algorithim not an implementation, and it can be applied to the clients in memory data
     * therefore it's defined as a template based on the users preferred data model
     **/

    template <typename _float>
    class DefaultSample : public Sample
    {
    public:
        inline double operator()(unsigned long int event, unsigned short int measurement) const noexcept
        {
            return (double)data[measurements * event + measurement];
        };

        DefaultSample(unsigned short int measurements,
                      unsigned long int events,
                      Subset subset,
                      Key key) noexcept
            : Sample(measurements, events, subset, key), data(nullptr){};

        DefaultSample(const unsigned short int measurements,
                      const unsigned long int events,
                      const _float *const data) noexcept
            : Sample(measurements, events), data(data)
        {
            for (long int event = 0; event < events; event++)
                for (unsigned short int measurement = 0; measurement < measurements; measurement++)
                    if (data[measurements * event + measurement] < 0 || data[measurements * event + measurement] > 1)
                        subset[event] = false;
        };

        DefaultSample(const unsigned short int measurements,
                      const unsigned long int events,
                      const _float *const data,
                      Subset subset) noexcept
            : Sample(measurements, events, subset), data(data){};

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
            float f = data[measurements * event + measurement];
            return *(epp_word *)&f;
        };

        void put_word(unsigned short int measurement, unsigned long int event, epp_word value) noexcept
        {
            float f = *(float *)&value;
            data[measurements * event + measurement] = (_float)f;
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
            return (double)data[events * measurement + event];
        };

        TransposeSample(unsigned short int measurements,
                        unsigned long int events,
                        Subset subset,
                        Key key) noexcept
            : Sample(measurements, events, subset, key), data(nullptr){};

        TransposeSample(const unsigned short int measurements,
                        const unsigned long int events,
                        const _float *const data) noexcept
            : Sample(measurements, events), data(data)
        {
            for (unsigned long int event = 0; event < events; event++)
                for (unsigned short int measurement = 0; measurement < measurements; measurement++)
                    if (data[events * measurement + event] < 0 || data[events * measurement + event] > 1)
                        subset[event] = false;
        };

        TransposeSample(const unsigned short int measurements,
                        const unsigned long int events,
                        const _float *const data,
                        Subset subset) noexcept
            : Sample(measurements, events, subset), data(data){};

    protected:
        epp_word get_word(unsigned short int measurement, unsigned long int event) const noexcept
        {
            float f = data[events * measurement + event];
            return *(epp_word *)&f;
        };

        void put_word(unsigned short int measurement, unsigned long int event, epp_word value) noexcept
        {
            float f = *(float *)&value;
            data[events * measurement + event] = (_float)f;
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

        PointerSample(unsigned short int measurements,
                      unsigned long int events,
                      Subset subset,
                      Key key) noexcept
            : Sample(measurements, events, subset, key), data(nullptr){};

        PointerSample(const unsigned short int measurements,
                      const unsigned long int events,
                      const _float *const *const data) noexcept
            : Sample(measurements, events), data(data)
        {
            for (unsigned long int event = 0; event < events; event++)
                for (unsigned short int measurement = 0; measurement < measurements; measurement++)
                    if (data[measurement][event] < 0 || data[measurement][event] > 1)
                        subset[event] = false;
        };

        PointerSample(const unsigned short int measurements,
                      const unsigned long int events,
                      const _float *const *const data,
                      Subset subset) noexcept
            : Sample(measurements, events, subset), data(data){};

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

        static const unsigned short N = 1 << 8; // resolution of points and boundaries
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
     * actually returned as a
     **/

    class _Request
    {
        friend class Request;
        friend class Pursuer;

    private:
        static std::mt19937_64 generate;
        std::chrono::time_point<std::chrono::steady_clock> begin, end;

    protected:
        std::condition_variable completed;
        Pursuer *const pursuer;
        volatile unsigned int outstanding = 0;
        volatile bool _finished;
        Result final_result;

        _Request(
            Pursuer *pursuer,
            Parameters parameters) noexcept;

    public:
        _Result *working_result;

        explicit operator json() const noexcept;

        _Request &operator=(const json &encoded);
    };

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

        Request(
            Pursuer *pursuer,
            Parameters parameters);

        Request(_Request *request) : std::shared_ptr<_Request>(request)
        {
            request->_finished = false;
        };

    public:
        bool finished() const noexcept;

        void wait() const noexcept;

        Result result() const noexcept;

        _Result *working()
        {
            return (*this)->working_result;
        }

        Request &operator++();

        Request &operator--();
    };

    // typedef std::shared_ptr<_Request> Request;

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
    };

    typedef std::vector<Point> Polygon;

    struct Candidate
    {
        Polygon separatrix;
        Subset in, out;
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

        explicit operator json();

        Candidate &operator=(const json &encoded);

        Candidate(const json &encoded)
        {
            *this = encoded;
        }

        Candidate(
            unsigned short int X,
            unsigned short int Y)
            : X(X < Y ? X : Y), Y(X < Y ? Y : X),
              outcome(Status::EPP_error),
              score(std::numeric_limits<double>::infinity()),
              pass(0), clusters(0), graphs(0){};
    };

    struct _Result
    {
        Key key;
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

        explicit operator json() const noexcept;

        _Result &operator=(const json &encoded);

        explicit _Result(const json &encoded)
        {
            *this = encoded;
        };

        explicit _Result(
            Parameters parameters)
            : projections(0),
              passes(0), clusters(0), graphs(0)
        {
            candidates.reserve(parameters.finalists);
        };
    };

    typedef std::shared_ptr<_Result> Result;

    /**
     * Classes defining the actual clients
     **/
    class Pursuer
    {
        friend class _Request;
        friend class Request;
        friend class Remote;

        template <class ClientSample>
        friend class SamplePursuer;

    protected:
        std::unordered_map<const Key, Request, KeyHash> requests;
        std::vector<std::thread> workers;
        std::mutex mutex;

        Request start(
            _Request *request) noexcept;

        void finish(
            _Request *request) noexcept;;

        void wait(
            _Request *request) noexcept;;

        void start(
            const json &encoded){};

        void finish(
            const json &encoded);

        Pursuer(
            Parameters parameters) noexcept
            : parameters(parameters){};

        Pursuer(
            Parameters parameters,
            int threads) noexcept
            : parameters(parameters), workers(threads){};

    public:
        Parameters parameters;

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
    class SamplePursuer : public Pursuer
    {
    public:
        // this is the virtual function that is required for every client
        // because we can't link templates we have to instantiate them individually
        // the others are convenience methods for specific clients
        virtual Request start(
            const ClientSample sample,
            const Parameters parameters) noexcept = 0;

        // default is to use the pursuer's parameters
        Request start(
            const ClientSample sample) noexcept
        {
            return start(sample, parameters);
        };

        Result pursue(
            const ClientSample sample) noexcept
        {
            return start(sample, parameters).result();
        };

    protected:
        SamplePursuer() noexcept = default;

        SamplePursuer(
            Parameters parameters,
            int threads)
            : Pursuer(parameters, threads){};

        SamplePursuer(
            Parameters parameters = Default)
            : Pursuer(parameters){};

        ~SamplePursuer() = default;
        ;
    };

    typedef TransposeSample<float> MATLAB_Sample;

    class MATLAB_Pursuer : public SamplePursuer<MATLAB_Sample>
    {
    public:
        Request start( // this one does the heavy lifting
            const MATLAB_Sample sample,
            const Parameters parameters) noexcept;

        Request start( // thiese are convenience routines
            const unsigned short int measurements,
            const unsigned long int events,
            const float *const data,
            Subset &subset) noexcept;
        Request start(
            const unsigned short int measurements,
            const unsigned long int events,
            const float *const data) noexcept;
        Result pursue(
            const unsigned short int measurements,
            const unsigned long int events,
            const float *const data,
            Subset &subset) noexcept;
        Result pursue(
            const unsigned short int measurements,
            const unsigned long int events,
            const float *const data) noexcept;

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

        void out(const json &encoded) // does not block
        {
            std::unique_lock<std::mutex> lock(mutex);
            outgoing.push(encoded);
            wake_out.notify_one();
        }

        json in() // blocks calling thread
        {
            std::unique_lock<std::mutex> lock(mutex);
            while (incoming.empty())
                wake_in.wait(lock);
            json encoded = incoming.front();
            incoming.pop();
            return encoded;
        }

    protected:
        void copy(
            std::istream *in,
            std::ostream *out,
            unsigned long int count)
        {
            char buffer[8192];
            while (count > 0)
            {
                unsigned long int chunk = count;
                if (chunk > 8192)
                    chunk = 8192;
                in->read(buffer, chunk);
                long int n = in->gcount();
                if (n > 0)
                    out->write(buffer, n);
                count -= n;
            }
        };

        void startFault(
            Key key)
        {
            json encoded;
            out(encoded);
        };

        void transmit()
        {
            // don't need to hold the lock while we do I/O
            json encoded;
            {
                std::unique_lock<std::mutex> lock(mutex);
                if (outgoing.empty())
                {
                    wake_out.wait(lock);
                    return; // make sure we're still on the air
                }
                else
                {
                    encoded = outgoing.front();
                    outgoing.pop();
                }
            }
            Service service = request; // from json
            switch (service)
            {
            case request:
            case result:
            case fault:
            {
                // serialize encoded and send on the control channel
                std::string serialized("<our json message>");
                remote_control->write(serialized.data(), serialized.size());
                break;
            }
            case content:
            {
                // send the blob on the data channel using stream in meta
                Key blob_key; // from json
                Meta *meta = blob_key.meta().get();
                std::iostream *blob = meta->stream;
                blob->clear();
                blob->seekp(0);

                // serialize encoded and send on the control channel
                std::string serialized("Whatever");
                remote_control->write(serialized.data(), serialized.size());
                // send the blob content on the data channel
                copy(blob, remote_data, meta->size);
                remote_data->flush();
                break;
            }
            }
            remote_control->flush();
        };

        void receive()
        {
            json encoded;              // deserialize from control channel
            Service service = request; // from json
            switch (service)
            {
            case request:
            case result:
            { // kick upstairs
                std::unique_lock<std::mutex> lock(mutex);
                incoming.push(encoded);
                wake_in.notify_one();
            }
            case content:
            {
                Key blob_key; // from json
                Meta *meta = blob_key.meta().get();
                std::ostream *content = meta->stream;
                content->clear();
                content->seekp(0);
                copy(remote_data, content, meta->size);
                content->flush();
                Blob::content(meta); // wakeup anyone waiting for the data
                break;
            }
            case fault:
            { // change json from fault to content service and
                // put it on the output queue
                out(encoded);
                break;
            }
            }
        };

        Remote()
        {
            Blob::handler = this;
            transmitter = std::thread(
                [this]()
                {
                    while (on_the_air)
                        transmit();
                });
            receiver = std::thread(
                [this]()
                {
                    while (on_the_air)
                        receive();
                });
        };

        ~Remote()
        {
            {
                std::unique_lock<std::mutex> lock(mutex);
                on_the_air = false;
                wake_in.notify_all();
                wake_out.notify_all();
            }
            transmitter.join();
            receiver.join();
        };

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
        Request start(
            const MATLAB_Sample sample,
            const Parameters parameters) noexcept;

        void finish(
            const json &encoded);

        MATLAB_Remote(
            Parameters parameters) noexcept;

        ~MATLAB_Remote();
    };
    /**
     * remote worker instance
     **/
    typedef DefaultSample<float> CloudSample;

    class CloudPursuer : Remote, public SamplePursuer<CloudSample>
    {
    public:
        Request start(
            const CloudSample sample,
            const Parameters parameters) noexcept;

        void start(const json &encoded);

        void finish(_Request *request) noexcept;

        void finish(const json &encoded);

        json remote();

        CloudPursuer(Parameters parameters) noexcept;

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