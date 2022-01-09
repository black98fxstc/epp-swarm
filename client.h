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
using json = nlohmann::json;
// typedef void *json;

namespace EPP
{
    typedef std::uint32_t epp_word;
    static std::random_device random;
    typedef uint32_t Booleans;

    /**
     * Exhaustive Projection Pursuit Client
     * 
     * Input parameters and output status
     */
    struct Parameters
    {
        // N = 256 gives points and features a precision of roughly two significant figures

        static const unsigned short N; // resolution of points and boundaries
                                       // optimized when there are lots of small factors

        double W = sqrt2 / (double)N; // standard deviation of kernel,
                                      // this is the highest achievable resolution, i.e., the resolution
                                      // along the diagonal. it works well but in practice a higher
                                      // value might be used for application reasons or just performance

        double sigma = 4; // controls the density threshold for starting a new cluster

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

        KLD kld{.16, .16, .16};

        std::vector<bool> censor; // omit measurements from consideration

        // algorithm tweaks

        bool recursive = false;

        unsigned int min_events = 0; // minimum events to try to split, max sigma squared

        unsigned int max_clusters = 12; // most clusters the graph logic should handle

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
            double W = sqrt2 / (double)N)
            : goal(goal), kld(kld), W(W), sigma(sigma),
              censor(0), finalists(1), max_clusters(12){};
    };

    const Parameters Default;

    enum Status
    {
        EPP_success,
        EPP_no_qualified,
        EPP_no_cluster,
        EPP_not_interesting,
        EPP_error
    };

    /**
     * Samples and Subsets
     */
    typedef uint64_t Event;
    typedef uint16_t Measurment;

    class Sample
    {
        friend class SampleStream;

    public:
        const Event events;
        const Measurment measurements;

    protected:
        Sample(Measurment measurements,
               Event events) noexcept
            : measurements(measurements), events(events){};

        Sample() = default;

    private:
        // these are virtual because our friend stream won't know which variant it will be
        virtual epp_word get_word(Measurment measurement, unsigned long event) const noexcept
        {
            return (epp_word)0;
        }
        virtual void put_word(unsigned short measurement, long event, epp_word data) const noexcept {};
    };

    template <class ClientSample>
    class SampleSubset;

    class Subset
    {
        friend class SubsetStream;

    public:
        const Sample &sample;

        bool contains(Event event) const noexcept
        {
            return data[event / 8] & 1 << event % 8;
        };

        void member(Event event, bool membership = false)
        {
            int q = event / 8;
            uint8_t b = 1 << event % 8;
            if (membership)
                data[q] |= b;
            else
                data[q] &= ~b;
        };

        Subset(const Sample &sample)
            : sample(sample)
        {
            data = new uint8_t[(sample.events + 7) / 8];
        };

        Subset(const Sample &sample, bool membership)
            : Subset(sample)
        {
            size_t q = (sample.events + 7) / 8;
            if (membership)
                std::memset(data, -1, q);
            else
                std::memset(data, 0, q);
        };

        Subset(const Sample &sample,
               const Subset &other) : Subset(sample)
        {
            std::memcpy(data, other.data, (sample.events + 7) / 8);
        };

        ~Subset()
        {
            delete[] data;
        };

        // Subset() = default;

    protected:
        uint8_t *data;
    };

    /**
     * structures defining the result of an ananlysis
     **/
    typedef int16_t Coordinate;

    struct Point
    {
        Coordinate i, j;

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

    class Candidate
    {
    public:
        Polygon separatrix;
        Subset in, out;
        Event in_events, out_events;

        double score, edge_weight, balance_factor;
        unsigned int pass, clusters, graphs;
        Measurment X, Y;
        enum Status outcome;

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
            const Sample &sample,
            Measurment X,
            Measurment Y)
            : X(X), Y(Y), in(sample, false), out(sample, false),
              outcome(Status::EPP_error),
              score(std::numeric_limits<double>::infinity()),
              pass(0), clusters(0), graphs(0){};

    private:
        static void close_clockwise(
            Polygon &polygon) noexcept;

        void simplify(
            const double tolerance,
            Polygon &simplified,
            const unsigned short int lo,
            const unsigned short int hi) const noexcept;
    };

    class Lysis
    {
    public:
        std::vector<Measurment> qualified;
        std::vector<Candidate *> candidates;
        std::chrono::milliseconds milliseconds;
        unsigned int projections, passes, clusters, graphs;

        const Candidate &winner() const noexcept
        {
            return *candidates[0];
        };

        enum Status outcome() const noexcept
        {
            return winner().outcome;
        };

        bool success() const noexcept
        {
            return winner().outcome == EPP_success;
        }

        Polygon separatrix() const noexcept
        {
            return winner().separatrix;
        };

        Polygon in_polygon() const noexcept
        {
            return winner().in_polygon();
        };

        Polygon in_polygon(
            double tolerance) const noexcept
        {
            return winner().in_polygon(tolerance);
        };

        Polygon out_polygon() const noexcept
        {
            return winner().out_polygon();
        };

        Polygon out_polygon(
            double tolerance) const noexcept
        {
            return winner().out_polygon(tolerance);
        };

        Measurment X()
        {
            return winner().X;
        };

        Measurment Y()
        {
            return winner().Y;
        };

        explicit operator json();

        // _Result &operator=(const json &encoded);

        // explicit _Result(const json &encoded)
        // {
        //     *this = encoded;
        // };

        Lysis(
            const Parameters &parameters)
            : projections(0), passes(0), clusters(0), graphs(0)
        {
            candidates.reserve(parameters.finalists);
        };

        ~Lysis()
        {
            for (Candidate *&candidate : candidates)
                delete candidate;
        };
    };


    /**
     * Provides a content based associative memory service
     * of blobs shared between clients and servers
     **/

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

    struct Key
    {
        friend class Blob;
        friend class Sample;
        friend class Remote;

        template <class ClientSample>
        friend class SamplePursuer;

        template <class ClientSample>
        friend class Request;

    protected:
        union
        { // 256 bit key
            std::uint8_t bytes[32];
            std::uint_fast64_t random[4];
        };

        static std::unordered_map<Key, std::weak_ptr<Meta>, Key> metadata;

        static void vacuum() noexcept;

        std::shared_ptr<Meta> meta() const noexcept;

        const Key &operator=(const Key &that) noexcept
        {
            std::memcpy(bytes, that.bytes, 32);
            return *this;
        };

    public:
        static std::mt19937_64 generate;

        std::size_t operator()(Key const &key) const noexcept
        {                                  // relies on the fact that the
            return *(std::size_t *)(&key); // key is already a good hash
        }

        bool operator==(
            const Key &other) const noexcept
        {
            return std::memcmp(bytes, other.bytes, 32) == 0;
        }

        Key &operator=(Key &&that) noexcept
        {
            if (!(*this == that))
                std::memcpy(bytes, that.bytes, 32);
            return *this;
        }

        Key &operator=(Key &that) noexcept
        {
            if (!(*this == that))
                std::memcpy(bytes, that.bytes, 32);
            return *this;
        }

        Key(const Key &key)
        {
            std::move(key.bytes, key.bytes + 32, bytes);
        };

        Key(std::mt19937_64 &generate)
        {
            for (auto &random_bits : random)
                random_bits = generate();
        };

        Key()
        {
            std::move(no_key, no_key + 32, bytes);
        };

        explicit operator json() const noexcept;

        Key &operator=(const json &encoded);

        explicit Key(const json &encoded)
        {
            *this = encoded;
        };

        Key(std::istream &stream);
    };

    const Key NoKey;

    class Blob
    {
        friend class Remote;

        // private:
        //     static std::mutex mutex;
        //     static std::condition_variable wakeup;

    protected:
        Key key;
        // size_t size();
        // std::shared_ptr<Meta> meta;

        // std::istream istream();

        // std::ostream ostream();

        // class Handler
        // {
        // public:
        //     void startFault(
        //         Key key){};
        // };

        // static Handler *handler;

    public:
        // bool valid();

        // bool fault();

        // void wait();

    protected:
        // static void content(
        //     Meta *meta);

        Blob(
            const Key &key);

        Blob();
    };

    /**
     * Templates depending on the clients data model
     */
    static int subset_count = 0;
    template <class ClientSample>
    class SampleSubset : public Subset, protected Blob
    {
        friend class SubsetStream;

    public:
        Event events;
        Measurment X, Y;
        Polygon polygon, simplified;
        SampleSubset *const parent;
        std::vector<SampleSubset *> children;

        SampleSubset(
            const ClientSample &sample)
            : Subset(sample, true), parent(nullptr)
        {
            for (Event event = 0; event < sample.events; event++)
                for (Measurment measurement = 0; measurement < sample.measurements; measurement++)
                    if (sample(event, measurement) < 0 || sample(event, measurement) > 1)
                        member(event, false);
        };

        operator json()
        {
            json subset;
            subset["ID"] = ++subset_count;
            subset["X"] = X;
            subset["Y"] = Y;
            subset["events"] = events;
            json polygon;
            int i = 0;
            for (auto &point : this->polygon)
            {
                json vertex;
                vertex[0] = point.x();
                vertex[1] = point.y();
                polygon[i++] = vertex;
            };
            subset["polygon"] = polygon;

            if (this->children.size() > 0)
            {
                json children;
                i = 0;
                for (SampleSubset *child : this->children)
                    children[i++] = (json)*child;
                subset["children"] = children;
            }
            return subset;
        };

        SampleSubset(
            const ClientSample &sample,
            SampleSubset *parent,
            const Subset &subset) : Subset(sample, subset), parent(parent){};

        ~SampleSubset()
        {
            for(auto &child : children)
                delete child;
        };
    };

    template <class ClientSample>
    class Analysis;

    template <class ClientSample>
    class Pursuer;

    template <class ClientSample>
    class Work;

    template <class ClientSample>
    class Request : public Lysis
    {
        friend class SampleStream;
        friend class Analysis<ClientSample>;
        friend class Pursuer<ClientSample>;
        friend class Work<ClientSample>;

    public:
        Analysis<ClientSample> *const analysis;
        const ClientSample &sample;
        SampleSubset<ClientSample> *subset;
        const Parameters &parameters;

    private:
        volatile bool finished;
        std::chrono::time_point<std::chrono::steady_clock> begin, end;
        Key key;
        volatile unsigned int outstanding = 0;

        Request(
            Analysis<ClientSample> *const analysis,
            const ClientSample &sample,
            SampleSubset<ClientSample> *subset,
            const Parameters &parameters) noexcept
            : analysis(analysis), sample(sample), subset(subset), parameters(parameters), Lysis(parameters){};

    protected:
        explicit operator json() const noexcept;
    };

    template <class ClientSample>
    class PursueProjection;

    template <class ClientSample>
    class Pursuer
    {
        friend class Work<ClientSample>;
        friend class Analysis<ClientSample>;

    public:
        const Parameters parameters;

        Analysis<ClientSample> *analyze(
            const ClientSample &sample,
            SampleSubset<ClientSample> *subset,
            const Parameters &parameters) noexcept
        {
            Analysis<ClientSample> *analysis = new Analysis<ClientSample>(this, sample, parameters);
            analysis->lyse(subset);
            return analysis;
        };

    protected:
        std::unordered_map<const Key, Request<ClientSample> *, Key> requests;
        std::vector<std::thread> workers;
        std::mutex mutex;

        void start(
            Request<ClientSample> *request) noexcept
        {
            Key key(Key::generate);
            request->key = key;
            request->finished = false;
            request->begin = std::chrono::steady_clock::now();

            PursueProjection<ClientSample>::start(request);

            std::unique_lock<std::mutex> lock(mutex);
            bool inserted = requests.insert(std::pair<const Key, Request<ClientSample> *>(request->key, request)).second;
            assert(inserted);
        }

        void increment(
            Request<ClientSample> *request) noexcept
        {
            std::unique_lock<std::mutex> lock(mutex);
            ++request->outstanding;
        }

        void decrement(
            Request<ClientSample> *request) noexcept
        {
            std::unique_lock<std::mutex> lock(mutex);
            --request->outstanding;
        }

        void finish(
            Request<ClientSample> *request) noexcept
        {
            request->end = std::chrono::steady_clock::now();
            request->milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(request->end - request->begin);
            request->analysis->finish(request);

            std::unique_lock<std::mutex> lock(mutex);

            auto it = requests.find(request->key);
            assert(it != requests.end());
            requests.erase(it);

            request->finished = true;
        };

        void wait(
            Request<ClientSample> *request) noexcept
        {
            std::unique_lock<std::mutex> lock(mutex);
            while (!request->finished)
                request->completed.wait(lock);
        };

        void start(
            const json &encoded){};

        void finish(
            const json &encoded)
        {
            Key request_key; // from JSON
            Request<ClientSample> *request = requests.find(request_key)->second;
            // _Result *result = request.working();
            // *result = encoded;
            // request.finish();
        }

        Pursuer(
            const Parameters &parameters) noexcept
            : parameters(parameters){};

        Pursuer(
            const Parameters &parameters,
            int threads) noexcept
            : parameters(parameters), workers(threads){};

    public:
        bool finished()
        {
            std::unique_lock<std::mutex> lock(mutex);
            for (auto it = requests.begin(); it != requests.end(); it++)
            {
                Request<ClientSample> request = it->second;
                if (request.finished())
                    return true;
            }
            return false;
        };
    };

    template <class ClientSample>
    class Analysis : public std::vector<Lysis *>
    {
        friend class Pursuer<ClientSample>;

    public:
        Pursuer<ClientSample> *const pursuer;
        const ClientSample &sample;
        const Parameters parameters;
        std::chrono::milliseconds milliseconds;
        std::chrono::milliseconds compute_time;
        unsigned int projections = 0, passes = 0, clusters = 0, graphs = 0;

        const Lysis *operator()(int i) const noexcept
        {
            return lysis[i];
        }

        size_type size() noexcept
        {
            return lysis.size();
        }

        bool complete()
        {
            std::unique_lock<std::mutex> lock(mutex);
            return lysis.size() == requests;
        }

        void wait() noexcept
        {
            std::unique_lock<std::mutex> lock(mutex);
            if (lysis.size() < requests)
                progress.wait(lock);
        };

    protected:
        std::mutex mutex;
        std::condition_variable progress;
        volatile int requests = 0;
        std::chrono::time_point<std::chrono::steady_clock> begin, end;
        std::vector<Request<ClientSample> *> lysis;

        void lyse(SampleSubset<ClientSample> *subset)
        {
            Request<ClientSample> *request = new Request<ClientSample>(this, this->sample, subset, parameters);
            pursuer->start(request);

            std::unique_lock<std::mutex> lock(mutex);
            ++requests;
        }

        void finish(
            Request<ClientSample> *request)
        {
            this->compute_time += request->milliseconds;
            this->projections += request->projections++;
            this->passes += request->passes;
            this->clusters += request->clusters;
            this->graphs += request->graphs;
            if (request->success() && request->analysis->parameters.recursive)
            {
                int threshold = std::max(
                    (unsigned int)(request->analysis->parameters.sigma * request->analysis->parameters.sigma),
                    request->analysis->parameters.max_clusters);
                if (request->winner().in_events > threshold)
                {
                    SampleSubset<ClientSample> *child = new SampleSubset<ClientSample>(this->sample, request->subset, request->winner().in);
                    lyse(child);
                    child->X = request->X();
                    child->Y = request->Y();
                    child->events = request->winner().in_events;
                    child->polygon = request->in_polygon(parameters.W);
                    request->subset->children.push_back(child);
                }
                if (request->winner().out_events > threshold)
                {
                    SampleSubset<ClientSample> *child = new SampleSubset<ClientSample>(this->sample, request->subset, request->winner().out);
                    lyse(child);
                    child->X = request->X();
                    child->Y = request->Y();
                    child->events = request->winner().out_events;
                    child->polygon = request->out_polygon(parameters.W);
                    request->subset->children.push_back(child);
                }
            }

            std::unique_lock<std::mutex> lock(mutex);
            lysis.push_back(request);
            progress.notify_all();
            end = std::chrono::steady_clock::now();
            milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
        };

        Analysis(
            Pursuer<ClientSample> *pursuer,
            const ClientSample &sample,
            const Parameters &parameters) : pursuer(pursuer), sample(sample), parameters(parameters)
        {
            begin = std::chrono::steady_clock::now();
        };

        Analysis(
            Pursuer<ClientSample> *pursuer,
            const ClientSample &sample)
            : pursuer(pursuer), sample(sample), parameters(Default)
        {
            begin = std::chrono::steady_clock::now();
        };
    };

    /**
     * EPP is an algorithim not an implementation, and it can be applied to the clients in memory data
     * therefore it's defined as a template based on the users preferred data model
     **/

    // distributed single precision measurment value
    class Value
    {
        union
        {
            float floating;
            uint32_t binary;
        };

        inline Value get(uint8_t *&ptr) const noexcept
        {
            Value v;
            v.binary |= *ptr++ << 24; // big endian
            v.binary |= *ptr++ << 16;
            v.binary |= *ptr++ << 8;
            v.binary |= *ptr++;
            return v;
        };

        inline void put(uint8_t *&ptr, Value v) const noexcept
        {
            *ptr++ = v.binary >> 24;
            *ptr++ = v.binary >> 16;
            *ptr++ = v.binary >> 8;
            *ptr++ = v.binary;
        }

        Value(float value) : floating(value){};
        Value(double value) : floating((float)value){};
        Value() : binary(0){};
    };

    template <typename _float>
    class DefaultSample : public Sample
    {
    public:
        inline const double operator()(unsigned long int event, unsigned short int measurement) const noexcept
        {
            return (double)data[Sample::measurements * event + measurement];
        };

        DefaultSample(Measurment measurements,
                      Event events,
                      SampleSubset<DefaultSample> subset,
                      Key key) noexcept
            : Sample(measurements, events, subset, key), data(nullptr){};

        DefaultSample(const Measurment measurements,
                      const Event events,
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
        epp_word get_word(Measurment measurement, Event event) const noexcept
        {
            float f = data[Sample::measurements * event + measurement];
            return *(epp_word *)&f;
        };

        void put_word(Measurment measurement, Event event, epp_word value) noexcept
        {
            float f = *(float *)&value;
            data[Sample::measurements * event + measurement] = (_float)f;
        };

    private:
        const _float *const data;
    };

    template <typename _float>
    class TransposeSample : public Sample, protected Blob
    {
    public:
        inline double operator()(unsigned long int event, Measurment measurement) const noexcept
        {
            return (double)data[Sample::events * measurement + event];
        };

        TransposeSample(Measurment measurements,
                        Event events,
                        SampleSubset<TransposeSample> subset,
                        Key key) noexcept
            : Sample(measurements, events, subset, key), data(nullptr){};

        TransposeSample(const Measurment measurements,
                        const Event events,
                        const _float *const data) noexcept
            : Sample(measurements, events), data(data){};

    protected:
        epp_word get_word(Measurment measurement, Event event) const noexcept
        {
            float f = data[Sample::events * measurement + event];
            return *(epp_word *)&f;
        };

        void put_word(Measurment measurement, Event event, epp_word value) noexcept
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
        inline double operator()(unsigned long int event, Measurment measurement) const noexcept
        {
            return (double)data[measurement][event];
        };

        PointerSample(Measurment measurements,
                      Event events,
                      SampleSubset<PointerSample> subset,
                      Key key) noexcept
            : Sample(measurements, events, subset, key), data(nullptr){};

        PointerSample(const unsigned short int measurements,
                      const Event events,
                      const _float *const *const data,
                      SampleSubset<PointerSample> subset) noexcept
            : Sample(measurements, events, subset, NoKey), data(data){};

    protected:
        epp_word get_word(Measurment measurement, Event event) const noexcept
        {
            float f = data[measurement][event];
            return *(epp_word *)&f;
        };

        void put_word(unsigned short int measurement, Event event, epp_word value) noexcept
        {
            float f = *(float *)&value;
            data[measurement][event] = (_float)f;
        };

    private:
        const _float *const *const data;
    };

    /**
     * communications services assuming control and data streams
     * are already set up. blob faults are handled here.
     * requests and results are safely queued for the pursuer
     **/

    class Remote //: Blob::Handler
    {
        // friend class Blob;

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
            Event count);

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
    /**
     * remote worker instance
     **/
    typedef DefaultSample<float> CloudSample;

    class CloudPursuer : Remote, public Pursuer<CloudSample>
    {
    public:
        void start(const json &encoded);

        void finish(Request<CloudSample> *request) noexcept;

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