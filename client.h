/*
 * Developer: Wayne Moore <wmoore@stanford.edu>
 * Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
 * License: BSD 3 clause
 */
#ifndef _EPP_CLIENT_H
#define _EPP_CLIENT_H 1

#include <chrono>
#include <random>
#include <thread>
#include <algorithm>
#include <vector>
#include <stack>
#include <queue>
#include <unordered_map>
#include <mutex>
#include <condition_variable>
#include <iomanip>
#include <cassert>

#include <nlohmann/json.hpp>
using json = nlohmann::json;
#include "constants.h"
#include "metadata.h"
#include "transform.h"

namespace EPP
{
    typedef uint32_t Event;
    typedef uint16_t Measurement;
    typedef uint32_t Count;
    typedef uint32_t Unique;
    typedef int16_t Coordinate;

    typedef std::uint32_t epp_word;

    /**
     * Exhaustive Projection Pursuit Client
     *
     * Input parameters and output status
     */
    class Parameters
    {
    public:
        // N = 256 gives points and features a precision of roughly two significant figures

        static const unsigned short N; // resolution of points and boundaries
                                       // optimized when there are lots of small factors

        double W = (double)2 / (double)N; // standard deviation of kernel,
                                          // this is the highest achievable resolution for DBM

        enum Goal
        {                    // which objective function
            best_separation, // lowest edge weight
            best_balance     // edge weight biased towards more even splits
        } goal = best_balance;

        struct KLD // KLD threshold for informative cases
        {
            double Normal2D = .24;      // is this population worth splitting?
            double Normal1D = .04;      // is the measurement just normal
            double Exponential1D = .40; // is this an exponential tail (CyToF)

            KLD(
                double Normal2D = .24,
                double Normal1D = .04,
                double Exponential1D = .40) noexcept
                : Normal2D(Normal2D), Normal1D(Normal1D), Exponential1D(Exponential1D){};
        };

        KLD kld{.24, .04, .40};

        std::vector<Measurement> censor; // omit measurements from consideration

        // algorithm tweaks

        bool recursive = true; // restart process on the two subsets

        Count finalists = 1; // remember this many of the best candidates

        Event min_events = 0;    // minimum events to try to split, max sigma squared
        double min_relative = 0; // minimum fraction of total events to try to split

        // implementation details, not intended for general use

        double sigma = 3;         // threshold for starting a new cluste
        Count max_clusters = 12;  // most clusters the graph logic should handle
        double tolerance = .01;   // default tolerance for polygon simplification
        double balance_power = 1; // exponent of balance factor

        double kernelRadius(unsigned int pass) const noexcept
        {
            return this->W * pow(sqrt2, pass) / sqrt2;
        }

        bool isCensored(Measurement measurement) const noexcept
        {
            return std::find(this->censor.begin(), this->censor.end(), measurement) != this->censor.end();
        }

        explicit operator json() const noexcept;

        Parameters &operator=(const json &encoded);

        Parameters(const json &encoded) { *this = encoded; };

        Parameters(
            Goal goal = best_balance,
            KLD kld = {.24, .04, .40},
            double W = sqrt2 / (double)N)
            : W(W), goal(goal), kld(kld), censor(0), finalists(1),
              sigma(3), max_clusters(12), tolerance(.01), balance_power(1){};
    };

    const Parameters Default;

    const static std::vector<std::string> Goal_strings{
        "best_separation", "best_balance"};

    enum Status
    {
        EPP_success,
        EPP_characterized,
        EPP_no_qualified,
        EPP_no_cluster,
        EPP_not_interesting,
        EPP_threshold,
        EPP_error
    };

    const static std::vector<std::string> Status_strings{
        "success", "characterized", "no_qualified", "no_cluster", "not_interesting", "error"};

    /**
     * Samples and Subsets
     */
    class Sample
    {
        friend class SampleStream;

    public:
        const Event events;
        const Measurement measurements;

        virtual ~Sample() = default;

    protected:
        Sample(Measurement measurements,
               Event events) noexcept
            : events(events), measurements(measurements){};

        // these are virtual because our friend stream won't know which variant it will be
        virtual epp_word get_word(Measurement measurement, Event event) const noexcept = 0;

        virtual void put_word(Measurement measurement, Event event, epp_word data) noexcept = 0;
    };

    class Subset
    {
        friend class SubsetStream;

    public:
        const Sample &sample;

        bool contains(Event event) const noexcept
        {
            return this->data[event / 8] & 1 << event % 8;
        };

        void member(Event event, bool membership = false)
        {
            if (membership)
                this->data[event / 8] |= 1 << event % 8;
            else
                this->data[event / 8] &= ~(1 << event % 8);
        };

        Subset(const Sample &sample)
            : sample(sample), data(new uint8_t[((size_t)sample.events + 7) / 8]){};

        Subset(const Sample &sample, bool membership)
            : Subset(sample)
        {
            if (membership)
                std::memset(this->data, -1, (size_t)(sample.events + 7) / 8);
            else
                std::memset(this->data, 0, (size_t)(sample.events + 7) / 8);
        };

        Subset(const Sample &sample,
               const Subset &other) : Subset(sample)
        {
            std::memcpy(this->data, other.data, (size_t)(sample.events + 7) / 8);
        };

        virtual ~Subset()
        {
            delete[] data;
        };

    protected:
        uint8_t *const data;
    };

    /**
     * structures defining the result of an ananlysis
     **/
    struct Point
    {
        Coordinate i, j;

        inline double x() const noexcept { return (double)i / (double)Parameters::N; };
        inline double y() const noexcept { return (double)j / (double)Parameters::N; };

        inline bool operator==(const Point &other) const noexcept
        {
            return this->i == other.i && this->j == other.j;
        }

        inline bool operator!=(const Point &other) const noexcept
        {
            return !(*this == other);
        }

        Point(short i, short j) noexcept : i(i), j(j){};

        Point() noexcept : i(0), j(0){};

        virtual ~Point() = default;
    };

    class Polygon : public std::vector<Point>
    {
    public:
        operator json() const noexcept;
    };

    class Candidate
    {
    public:
        Polygon separatrix;
        Subset in, out;
        Event in_events, out_events;
        double score, edge_weight, balance_factor;
        unsigned int pass, clusters, graphs, merges;
        Measurement X, Y;
        enum Status outcome;

        bool operator<(const Candidate &other) const noexcept
        {
            if (this->score < other.score)
                return true;
            if (this->score > other.score)
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
            Measurement X,
            Measurement Y)
            : in(sample, false), out(sample, false), in_events(0), out_events(0),
              score(std::numeric_limits<double>::infinity()),
              pass(0), clusters(0), graphs(0), merges(0),
              X(X), Y(Y), outcome(Status::EPP_error){};

    private:
        void close_clockwise(
            Polygon &polygon) const noexcept;

        void simplify(
            const double tolerance,
            Polygon &simplified,
            const size_t lo,
            const size_t hi) const noexcept;
    };

    class Lysis
    {
    public:
        std::vector<Measurement> qualified;
        std::vector<Candidate *> candidates;
        std::vector<double> markers;
        std::vector<double> covariance;
        std::vector<double> invcovariance;
        std::vector<Lysis *> children;
        std::chrono::milliseconds milliseconds = std::chrono::milliseconds::zero();
        Lysis *parent;
        Event events;
        double divergence;
        unsigned int projections, passes, clusters, graphs, merges;
        Unique ID;
        Unique taxon;
        bool in_set;

        const Candidate &winner() const noexcept
        {
            return *this->candidates[0];
        };

        enum Status outcome() const noexcept
        {
            if (candidates.size() > 0)
                return this->winner().outcome;
            else
                return EPP_no_cluster;
        };

        bool success() const noexcept
        {
            return this->outcome() == EPP_success;
        }

        Polygon separatrix() const noexcept
        {
            return this->winner().separatrix;
        };

        const Subset &in() const noexcept
        {
            return this->winner().in;
        };

        Event in_events() const noexcept
        {
            return this->winner().in_events;
        };

        Polygon in_polygon() const noexcept
        {
            return this->winner().in_polygon();
        };

        Polygon in_polygon(
            double tolerance) const noexcept
        {
            return this->winner().in_polygon(tolerance);
        };

        const Subset &out() const noexcept
        {
            return this->winner().out;
        };

        Event out_events() const noexcept
        {
            return this->winner().out_events;
        };

        Polygon out_polygon() const noexcept
        {
            return this->winner().out_polygon();
        };

        Polygon out_polygon(
            double tolerance) const noexcept
        {
            return this->winner().out_polygon(tolerance);
        };

        Measurement X() const noexcept
        {
            return this->winner().X;
        }

        Measurement Y() const noexcept
        {
            return this->winner().Y;
        }

        json gating(double tolerance = Default.tolerance) const noexcept;

        explicit operator json() const noexcept;

        virtual ~Lysis()
        {
            for (Candidate *candidate : this->candidates)
                delete candidate;
        };

    protected:
        Lysis(
            const Parameters &parameters,
            Event events,
            Measurement measurements)
            : markers(measurements, 0), events(events), divergence(0), parent(nullptr),
              projections(0), passes(0), clusters(0),
              graphs(0), merges(0)
        {
            this->candidates.reserve(parameters.finalists);
        }
    };

    class Taxon
    {
        friend class Taxonomy;

    public:
        std::vector<Taxon *> subtaxa;
        std::vector<double> markers;
        std::vector<double> covariance;
        std::vector<bool> sibling;
        Taxon *supertaxon;
        Lysis *subset;
        double divergence, dissimilarity, depth;
        Event population;
        int rank, height;
        Unique ID;

        bool isSpecific() const noexcept { return this->subtaxa.empty(); }
        bool isGeneric() const noexcept { return !this->subtaxa.empty(); }

        explicit operator json() const noexcept;

        Taxon(Lysis *subset);
        Taxon(Taxon *red, Taxon *blue, bool merge = false);

        bool operator<(const Taxon &that) { return this->population < that.population; }

    private:
        double walk(std::vector<Taxon *> &phenogram);
    };

    class Similarity
    {
    public:
        double dissimilarity;
        Taxon *red, *blue;

        bool operator<(const Similarity &that) const noexcept { return this->dissimilarity > that.dissimilarity; }

        Similarity(
            Taxon *red,
            Taxon *blue);
    };

    class Phenogram : public std::vector<Taxon *>
    {
    public:
        explicit operator json() const noexcept;

        void toHtml(
            std::vector<std::string> markers,
            std::ofstream &html) noexcept;
    };

    class Taxonomy : public std::vector<Taxon *>
    {
        friend class Taxon;

    public:
        static double cityBlockDistance(
            std::vector<double> &red,
            std::vector<double> &blue) noexcept;

        static Taxon *classify(std::vector<Taxon *> &taxonomy) noexcept;

        static Phenogram phenogram(std::vector<Taxon *> &taxonomy);

        static std::string ascii(std::vector<Taxon *> &phenogram,
                                 std::vector<std::string> markers);

        static std::string ascii(std::vector<Taxon *> &phenogram);
    };

    /**
     * Templates depending on the clients data model
     */
    template <class ClientSample>
    class SampleSubset : public Subset, protected Blob
    {
        friend class SubsetStream;

    public:
        Event events;
        Measurement X, Y;
        Polygon polygon;
        const SampleSubset *const parent;
        std::vector<const SampleSubset *> children;

        SampleSubset(
            const ClientSample &sample)
            : Subset(sample, true), events(sample.events), parent(nullptr)
        {
            for (Event event = 0; event < this->sample.events; event++)
                for (Measurement measurement = 0; measurement < this->sample.measurements; measurement++)
                    if (sample(event, measurement) < 0 || sample(event, measurement) > 1)
                    {
                        this->member(event, false);
                        --events;
                    };
        };

        json tree() const noexcept
        {
            static int subset_count = 0;
            json subset;
            subset["ID"] = ++subset_count;
            if (parent)
            {
                subset["X"] = X;
                subset["Y"] = Y;

                json polygon;
                for (auto &point : this->polygon)
                {
                    json vertex;
                    vertex[0] = point.x();
                    vertex[1] = point.y();
                    polygon += vertex;
                };
                subset["polygon"] = polygon;
            }
            subset["events"] = events;

            if (this->children.size() > 0)
            {
                json children;
                for (const SampleSubset *child : this->children)
                    children += child->tree();
                subset["children"] = children;
            }
            return subset;
        };

        SampleSubset(
            const ClientSample &sample,
            const SampleSubset *parent,
            const Subset &subset) noexcept : Subset(sample, subset), parent(parent){};

        virtual ~SampleSubset()
        {
            for (auto &child : this->children)
                delete child;
        };
    };

    template <class ClientSample>
    class Work;

    template <class ClientSample>
    class Worker;

    template <class ClientSample>
    class Pursuer;

    template <class ClientSample>
    class Analysis;

    template <class ClientSample>
    class PursueProjection;

    template <class ClientSample>
    class QualifyMeasurement;

    template <class ClientSample>
    class CharacterizeSubset;

    template <class ClientSample>
    class Request : public Lysis
    {
        friend class SampleStream;
        friend class Analysis<ClientSample>;
        friend class Pursuer<ClientSample>;
        friend class Work<ClientSample>;
        friend class PursueProjection<ClientSample>;
        friend class QualifyMeasurement<ClientSample>;

    public:
        Analysis<ClientSample> *const analysis;
        const ClientSample &sample;
        SampleSubset<ClientSample> *subset;

    private:
        volatile bool finished;
        Status status;
        std::chrono::time_point<std::chrono::steady_clock> begin, end;
        Key key;
        volatile unsigned int outstanding = 0;
        struct
        {
            Measurement X = 0, Y = 1;
            double X_KLD = 0;
            double Y_KLD = 0;
        } fallback;
        Measurement qualifying = 0;

        Request(
            Analysis<ClientSample> *const analysis,
            const ClientSample &sample,
            SampleSubset<ClientSample> *subset,
            const Parameters &parameters) noexcept
            : Lysis(parameters, subset->events, sample.measurements),
              analysis(analysis), sample(sample), subset(subset)
        {
            this->qualifying = analysis->qualifying;
            ID = analysis->unique();
        };

        virtual ~Request() = default;

    protected:
        explicit operator json() const noexcept;
    };

    /**
     * Pursuer orchestrates the asynchronous worker threads and returns results to an Analysis
     **/
    template <class ClientSample>
    class Pursuer
    {
        friend class Work<ClientSample>;
        friend class Analysis<ClientSample>;
        friend class PursueProjection<ClientSample>;

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

        Analysis<ClientSample> *analyze(
            const ClientSample &sample,
            SampleSubset<ClientSample> *subset) noexcept
        {
            return analyse(sample, subset, this->parameters);
        }

    protected:
        std::unordered_map<const Key, Request<ClientSample> *, Key> requests;
        std::vector<std::thread> workers;
        std::recursive_mutex mutex;
        std::mt19937_64 generate; // not thread safe
        Transform transform;      // because all the threads must use the same FFT size

        void start(
            Request<ClientSample> *request) noexcept
        {
            std::lock_guard<std::recursive_mutex> lock(mutex);
            request->finished = false;
            Key key(generate);
            request->key = key;
            bool inserted = requests.insert(std::pair<const Key, Request<ClientSample> *>(request->key, request)).second;
            assert(inserted);
        }

        void increment(
            Request<ClientSample> *request) noexcept
        {
            std::lock_guard<std::recursive_mutex> lock(mutex);
            ++request->outstanding;
        }

        bool decrement(
            Request<ClientSample> *request) noexcept
        {
            std::lock_guard<std::recursive_mutex> lock(mutex);
            --request->outstanding;
            return request->outstanding == 0;
        }

        void finish(
            Request<ClientSample> *request) noexcept
        {
            {
                std::lock_guard<std::recursive_mutex> lock(mutex);

                auto it = requests.find(request->key);
                assert(it != requests.end());
                requests.erase(it);

                request->finished = true;
            }

            request->analysis->finish(request);
        };

        // void start(
        //     const json &encoded){};

        // void finish(
        //     const json &encoded)
        // {
        //     Key request_key; // from JSON
        //     Request<ClientSample> *request = requests.find(request_key)->second;
        // }

        Pursuer(
            const Parameters &parameters,
            int threads) noexcept
            : parameters(parameters),
              workers(threads < 0 ? std::thread::hardware_concurrency() : threads),
              transform(parameters.N, this->mutex)
        {
            std::random_device random;
            generate.seed(random());
            Worker<ClientSample>::revive();
            for (size_t i = 0; i < workers.size(); i++)
                workers[i] = std::thread(
                    []()
                    { Worker<ClientSample> worker; });
        }

        ~Pursuer()
        {
            Worker<ClientSample>::kiss();
            for (size_t i = 0; i < workers.size(); i++)
                workers[i].join();
        }
    };

    /**
     * An Analysis tries to recursively split a subset using a Pursuer instance
     * then collects and marshals the results
     **/
    template <class ClientSample>
    class Analysis
    {
        friend class Pursuer<ClientSample>;
        friend class Request<ClientSample>;
        friend class Taxonomy;

    public:
        Pursuer<ClientSample> *const pursuer;
        const ClientSample &sample;
        const Parameters parameters;
        std::chrono::milliseconds milliseconds;
        std::chrono::milliseconds compute_time = std::chrono::milliseconds::zero();
        const double *const *const kernel; // constant for each analysis
        Unique *classification;
        float *mahalanobis;
        Count projections = 0, passes = 0, clusters = 0, graphs = 0, merges = 0;

        Unique unique()
        {
            std::lock_guard<std::mutex> lock(mutex);
            return ++this->uniques;
        }

        const Lysis *operator()(int i) const noexcept
        {
            return this->lysis[i];
        }

        json gating()
        {
            return this->lysis.front()->gating(this->parameters.tolerance);
        }

        Taxon *classify()
        {
            if (taxonomy.back()->subtaxa.size() > 0)
                return taxonomy.back();

            std::lock_guard<std::mutex> lock(mutex);
            Taxonomy::classify(this->taxonomy);
            for (Taxon *tax : this->taxonomy)
                if (tax->isGeneric())
                    tax->ID = ++this->uniques;
            for (Taxon *tax : this->taxonomy)
                if (tax->subset)
                    tax->subset->taxon = tax->ID;

            return this->taxonomy.back();
        }

        Phenogram phenogram()
        {
            return Taxonomy::phenogram(this->taxonomy);
        }

        Count types() noexcept
        {
            return (Count)this->_type.size();
        }

        const Lysis *type(int i) const noexcept
        {
            return this->_type[i];
        }

        Count size() noexcept
        {
            return (Count)this->lysis.size();
        }

        bool complete() noexcept
        {
            std::lock_guard<std::mutex> lock(mutex);
            return lysis.size() == this->requests;
        }

        bool censor(Measurement measurement) const noexcept
        {
            if (measurement < sample.measurements)
                return this->censored[measurement];
            else
                return true;
        }

        void wait() noexcept
        {
            std::unique_lock<std::mutex> lock(mutex);
            if (lysis.size() < this->requests)
                this->progress.wait(lock);
        }

        ~Analysis()
        {
            delete[] this->censored;
            delete[] this->classification;
            delete[] this->mahalanobis;
            for (Count i = 0; i <= max_passes; ++i)
                delete[] this->kernel[i];
            delete[] this->kernel;
            for (auto &ly : this->lysis)
                delete ly;
            for (auto &tax : this->taxonomy)
                delete tax;
        }

    protected:
        std::mutex mutex;
        std::condition_variable progress;
        std::chrono::time_point<std::chrono::steady_clock> begin, end;
        std::vector<Request<ClientSample> *> lysis;
        std::vector<Request<ClientSample> *> _type;
        std::vector<Taxon *> taxonomy;
        bool *censored;
        Event threshold;
        Measurement qualifying = 0;
        volatile Count requests = 0;
        Unique uniques = 0;

        Lysis *lyse(SampleSubset<ClientSample> *subset)
        {
            Request<ClientSample> *request = new Request<ClientSample>(this, this->sample, subset, this->parameters);
            {
                std::lock_guard<std::mutex> lock(mutex);
                ++requests;
            }
            this->pursuer->start(request);

            if (subset->events > this->threshold)
            {
                request->status = EPP_no_cluster;
                PursueProjection<ClientSample>::start(request);
            }
            else
            {
                request->status = EPP_threshold;
                this->pursuer->finish(request);
            }

            return request;
        }

        void finish(
            Request<ClientSample> *request)
        {
            Lysis *ly;
            SampleSubset<ClientSample> *child;
            switch (request->status)
            {
            case EPP_success:
            {
                child = new SampleSubset<ClientSample>(this->sample, request->subset, request->in());
                child->X = request->X();
                child->Y = request->Y();
                child->events = request->in_events();
                child->polygon = request->in_polygon(parameters.tolerance);
                request->subset->children.push_back(child);
                if (request->analysis->parameters.recursive)
                {
                    ly = lyse(child);
                    ly->parent = request;
                    ly->events = child->events;
                    ly->in_set = true;
                    request->children.push_back(ly);
                }

                child = new SampleSubset<ClientSample>(this->sample, request->subset, request->out());
                child->X = request->X();
                child->Y = request->Y();
                child->events = request->out_events();
                child->polygon = request->out_polygon(parameters.tolerance);
                request->subset->children.push_back(child);
                if (request->analysis->parameters.recursive)
                {
                    ly = lyse(child);
                    ly->parent = request;
                    ly->events = child->events;
                    ly->in_set = false;
                    request->children.push_back(ly);
                }

                std::lock_guard<std::mutex> lock(mutex);
                lysis.push_back(request);
                break;
            }

            case EPP_no_cluster:
            case EPP_threshold:
            {
                pursuer->start(request);

                Worker<ClientSample>::enqueue(
                    new CharacterizeSubset<ClientSample>(request));

                break;
            }

            case EPP_characterized:
            {
                std::lock_guard<std::mutex> lock(mutex);
                lysis.push_back(request);
                _type.push_back(request);

                Taxon *tax = new Taxon(request);
                tax->ID = ++this->uniques;
                taxonomy.push_back(tax);
                break;
            }

            default:
                assert(("oops", false));
            }

            std::lock_guard<std::mutex> lock(mutex);
            switch (request->status)
            {
            case EPP_success:
            case EPP_characterized:
            {
                this->compute_time += request->milliseconds;
                this->projections += request->projections;
                this->passes += request->passes;
                this->clusters += request->clusters;
                this->graphs += request->graphs;
                this->merges += request->merges;
                this->end = std::chrono::steady_clock::now();
                this->milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
                break;
            }

            case EPP_no_cluster:
            case EPP_threshold:
            {
                request->status = EPP_characterized;
                break;
            }
            }
            this->progress.notify_all();
        };

        static const double *const *initKernel(const Parameters &parameters) noexcept
        {
            // precompute all the kernel coefficients once
            double **kernel = new double *[max_passes + 1];
            for (Count pass = 0; pass <= max_passes; ++pass)
            {
                double *k = kernel[pass] = new double[N + 1];
                double radius = parameters.kernelRadius(pass);
                for (int i = 0; i <= N; i++)
                    k[i] = exp(-i * i * radius * radius * pi * pi * 2);
            }
            return kernel;
        }

        Analysis(
            Pursuer<ClientSample> *pursuer,
            const ClientSample &sample,
            const Parameters &parameters) : pursuer(pursuer), sample(sample), parameters(parameters), kernel(initKernel(parameters))
        {
            this->begin = std::chrono::steady_clock::now();
            this->classification = new Unique[sample.events];
            this->mahalanobis = new float[sample.events];
            std::fill(this->classification, this->classification + sample.events, 0);
            this->censored = new bool[sample.measurements];
            this->threshold = std::max(
                std::max(
                    (Event)(parameters.min_relative * sample.events), // relative to current sample
                    parameters.min_events),                           // absolute event count
                (Event)(parameters.sigma * parameters.sigma));        // algorithim limit
            for (Measurement measurement = 0; measurement < sample.measurements; ++measurement)
                if (parameters.isCensored(measurement))
                    this->censored[measurement] = true;
                else
                {
                    ++this->qualifying;
                    this->censored[measurement] = false;
                }
        };
    };

    const unsigned short Parameters::N = 1 << 8; // resolution of points and boundaries

    static size_t find_string(
        std::string string,
        const std::vector<std::string> strings)
    {
        for (size_t i = 0; i < strings.size(); i++)
            if (string == strings[i])
                return i;
        throw std::runtime_error("enumeration string did not match");
    };

    Parameters::operator json() const noexcept
    {
        json parameters;
        parameters["N"] = this->N;
        parameters["W"] = this->W;
        parameters["goal"] = Goal_strings[this->goal];
        parameters["finalists"] = this->finalists;
        json kld;
        kld["Normal2D"] = this->kld.Normal2D;
        kld["Normal1D"] = this->kld.Normal1D;
        kld["Exponential1D"] = this->kld.Exponential1D;
        parameters["KLD"] = kld;
        json censor = json::array();
        for (size_t i = 0; i < this->censor.size(); i++)
            censor[i] = this->censor.at(i);
        parameters["censor"] = censor;
        parameters["min_events"] = this->min_events;
        parameters["min_relative"] = this->min_relative;
        parameters["balance_power"] = this->balance_power;
        parameters["sigma"] = this->sigma;
        parameters["max_clusters"] = this->max_clusters;
        parameters["tolerance"] = this->tolerance;
        return parameters;
    }

    Parameters &Parameters::operator=(const json &encoded)
    {
        this->W = encoded.value("W", Default.W);
        if (encoded.contains("goal"))
        {
            std::string goal = encoded["goal"];
            this->goal = (Goal)find_string(goal, Goal_strings);
        }
        else
            this->goal = Default.goal;
        this->finalists = encoded.value("finalists", Default.finalists);
        if (encoded.contains("KLD"))
        {
            json kld = encoded["KLD"];
            this->kld.Normal2D = kld.value("Normal2D", Default.kld.Normal2D);
            this->kld.Normal1D = kld.value("Normal1D", Default.kld.Normal1D);
            this->kld.Exponential1D = kld.value("Exponential1D", Default.kld.Exponential1D);
        }
        this->recursive = encoded.value("recursive", Default.recursive);
        this->sigma = encoded.value("sigma", Default.sigma);
        this->min_events = encoded.value("min_events", Default.min_events);
        this->min_relative = encoded.value("min_relative", Default.min_relative);
        this->balance_power = encoded.value("balance_power", Default.balance_power);
        this->max_clusters = encoded.value("max_clusters", Default.max_clusters);
        this->tolerance = encoded.value("tolerance", Default.tolerance);
        if (encoded.contains("censor"))
        {
            json censor = encoded["censor"];
            for (size_t i = 0; i < censor.size(); i++)
                this->censor.push_back(censor.at(i));
        }
        return *this;
    }

    Polygon::operator json() const noexcept
    {
        json polygon;
        for (const Point &pt : *this)
        {
            json point;
            point[0] = pt.x();
            point[1] = pt.y();
            polygon += point;
        }
        return polygon;
    }

    Candidate::operator json()
    {
        json candidate;
        candidate["separatrix"] = (json)this->separatrix;
        candidate["score"] = this->score;
        candidate["edge_weight"] = this->edge_weight;
        candidate["balance_factor"] = this->balance_factor;
        candidate["in_events"] = this->in_events;
        candidate["out_events"] = this->out_events;
        candidate["pass"] = this->pass;
        candidate["clusters"] = this->clusters;
        candidate["graphs"] = this->graphs;
        candidate["X"] = this->X;
        candidate["Y"] = this->Y;
        candidate["outcome"] = Status_strings[this->outcome];
        return candidate;
    }

    Candidate &Candidate::operator=(const json &encoded)
    {
        json separatrix = encoded["separatrix"];
        // Polygon polygon(separatrix.size());
        // for (size_t i = 0; i < separatrix.size(); i++)
        // {
        //     json point = separatrix[i];
        //     polygon[i] = Point(point[0], point[1]);
        // }
        // this->separatrix = polygon;
        this->score = encoded["score"];
        this->edge_weight = encoded["edge_weight"];
        this->balance_factor = encoded["balance_factor"];
        this->in_events = encoded["in_events"];
        this->out_events = encoded["out_events"];
        this->pass = encoded["pass"];
        this->clusters = encoded["clusters"];
        this->graphs = encoded["graphs"];
        this->X = encoded["X"];
        this->Y = encoded["Y"];
        this->outcome = (Status)find_string(encoded["outcome"], Status_strings);
        return *this;
    }

    Lysis::operator json() const noexcept
    {
        json lysis;
        lysis["events"] = this->events;
        lysis["ID"] = this->ID;
        if (this->taxon)
            lysis["taxon"] = this->taxon;
        if (this->parent)
        {
            lysis["X"] = this->X();
            lysis["Y"] = this->Y();
        }
        if (this->children.size() > 0)
        {
            json children;
            for (const Lysis *child : this->children)
                children += (json)*child;
            lysis["children"] = children;
        }
        else
            lysis["divergence"] = this->divergence;

        return lysis;
    }

    json Lysis::gating(double tolerance) const noexcept
    {
        json gates;
        gates["ID"] = this->ID;
        gates["events"] = this->events;
        if (this->parent)
        {
            gates["X"] = parent->X();
            gates["Y"] = parent->Y();
            if (this->in_set)
                gates["polygon"] = parent->in_polygon(tolerance);
            else
                gates["polygon"] = parent->out_polygon(tolerance);
        }
        if (this->children.size())
        {
            json children;
            for (Lysis *ly : this->children)
                children += ly->gating(tolerance);
            gates["children"] = children;
        }
        else
            gates["taxon"] = this->taxon;
        return gates;
    }
}

#include "sample.h"
#include "stream.h"
#include "polygon.h"

#endif /* _EPP_CLIENT_H */
