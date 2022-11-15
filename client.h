/*
 * Developer: Wayne Moore <wmoore@stanford.edu>
 * Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
 * License: BSD 3 clause
 */
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

#include <nlohmann/json.hpp>
using json = nlohmann::json;
// typedef void *json;
#include "constants.h"
#include "metadata.h"

namespace EPP
{
    typedef uint32_t Event;
    typedef uint16_t Measurement;
    typedef uint32_t Count;
    typedef uint32_t Unique;
    typedef int16_t Coordinate;

    typedef std::uint32_t epp_word;
    static std::random_device random;

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
                double Exponential1D = .40)
            noexcept
                : Normal2D(Normal2D), Normal1D(Normal1D), Exponential1D(Exponential1D){};
        };

        KLD kld{.24, .04, .40};

        std::vector<Measurement> censor; // omit measurements from consideration

        // algorithm tweaks

        bool recursive = true; // restart process on the two subsets

        int finalists = 1; // remember this many of the best candidates

        unsigned int min_events = 0; // minimum events to try to split, max sigma squared
        double min_relative = 0;     // minimum fraction of total events to try to split
        double balance_power = 1;    // exponent of balance factor

        // implementation details, not intended for general use

        double sigma = 3;               // threshold for starting a new cluste
        unsigned int max_clusters = 12; // most clusters the graph logic should handle
        double tolerance = .01;         // default tolerance for polygon simplification

        double kernelWidth(unsigned int pass) const noexcept
        {
            return this->W * pow(sqrt2, pass); // kernel increases by sqrt(2) each pass so
        }                                      // FFT calculation can be reused

        explicit operator json() const noexcept;

        Parameters &operator=(const json &encoded);

        Parameters(const json &encoded) { *this = encoded; };

        Parameters(
            Goal goal = best_balance,
            KLD kld = {.12, .04, .40},
            double W = sqrt2 / (double)N)
            : goal(goal), kld(kld), W(W), sigma(3), tolerance(.01),
              censor(0), finalists(1), max_clusters(12), balance_power(1){};
    };

    const Parameters Default;

    enum Status
    {
        EPP_success,
        EPP_characterized,
        EPP_no_qualified,
        EPP_no_cluster,
        EPP_not_interesting,
        EPP_error
    };

    /**
     * Samples and Subsets
     */
    class Sample
    {
        friend class SampleStream;

    public:
        const Event events;
        const Measurement measurements;

    protected:
        Sample(Measurement measurements,
               Event events) noexcept
            : measurements(measurements), events(events){};

    private:
        // these are virtual because our friend stream won't know which variant it will be
        virtual epp_word get_word(Measurement measurement, Event event) const noexcept
        {
            return (epp_word)0;
        }
        virtual void put_word(Measurement measurement, Event event, epp_word data) const noexcept {};
    };

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
            if (membership)
                data[event / 8] |= 1 << event % 8;
            else
                data[event / 8] &= ~(1 << event % 8);
        };

        Subset(const Sample &sample)
            : sample(sample), data(new uint8_t[((size_t)sample.events + 7) / 8]){};

        Subset(const Sample &sample, bool membership)
            : Subset(sample)
        {
            if (membership)
                std::memset(data, -1, (size_t)(sample.events + 7) / 8);
            else
                std::memset(data, 0, (size_t)(sample.events + 7) / 8);
        };

        Subset(const Sample &sample,
               const Subset &other) : Subset(sample)
        {
            std::memcpy(data, other.data, (size_t)(sample.events + 7) / 8);
        };

        ~Subset()
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
            return i == other.i && j == other.j;
        }

        inline Point &operator=(const Point &other)
        {
            this->i = other.i;
            this->j = other.j;
            return *this;
        }

        inline bool operator!=(const Point &other) const noexcept
        {
            return !(*this == other);
        }

        Point(short i, short j) noexcept : i(i), j(j){};

        Point() : i(0), j(0){};
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
            Measurement X,
            Measurement Y)
            : X(X), Y(Y), in(sample, false), out(sample, false),
              in_events(0), out_events(0), outcome(Status::EPP_error),
              score(std::numeric_limits<double>::infinity()),
              pass(0), clusters(0), graphs(0), merges(0){};

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
        Event events;
        Unique ID;
        Lysis *parent;
        std::vector<Lysis *> children;
        Unique taxon = 0;
        std::chrono::milliseconds milliseconds = std::chrono::milliseconds::zero();
        unsigned int projections, passes, clusters, graphs, merges;
        bool in_set = true;

        const Candidate &winner() const noexcept
        {
            return *candidates[0];
        };

        enum Status outcome() const noexcept
        {
            if (candidates.size() > 0)
                return winner().outcome;
            else
                return EPP_no_cluster;
        };

        bool success() const noexcept
        {
            return outcome() == EPP_success;
        }

        Polygon separatrix() const noexcept
        {
            return winner().separatrix;
        };

        const Subset &in() const noexcept
        {
            return winner().in;
        };

        Event in_events() const noexcept
        {
            return winner().in_events;
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

        const Subset &out() const noexcept
        {
            return winner().out;
        };

        Event out_events() const noexcept
        {
            return winner().out_events;
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

        Measurement X() const noexcept
        {
            return winner().X;
        }

        Measurement Y() const noexcept
        {
            return winner().Y;
        }

        json gating(double tolerance = Default.tolerance) const noexcept;

        explicit operator json() const noexcept;

        ~Lysis()
        {
            for (Candidate *candidate : candidates)
                delete candidate;
        };

    protected:
        Lysis(
            const Parameters &parameters,
            Event events,
            Measurement measurements)
            : events(events), markers(measurements, 0), parent(nullptr),
              projections(0), passes(0), clusters(0),
              graphs(0), merges(0)
        {
            candidates.reserve(parameters.finalists);
        }
    };
    
    class Taxon
    {
        friend class Taxonomy;

    public:
        std::vector<Taxon *> subtaxa;
        std::vector<double> markers;
        std::vector<bool> connect;
        Taxon *supertaxon;
        Lysis *subset;
        double dissimilarity, depth;
        Event population;
        int rank, height;
        Unique ID;

        bool isSpecific() const noexcept { return subtaxa.empty(); }
        bool isGeneric() const noexcept { return !subtaxa.empty(); }

        explicit operator json() const noexcept;

        Taxon(Lysis *subset);
        Taxon(Event population, std::vector<double> &markers) : population(population), markers(markers) {}
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

        bool operator<(const Similarity &that) { return this->dissimilarity > that.dissimilarity; }

        Similarity(
            Taxon *red,
            Taxon *blue);
    };

    class Taxonomy : public std::vector<Taxon *>
    {
        friend class Taxon;

    public:
        static double cityBlockDistance(
            std::vector<double> &red,
            std::vector<double> &blue) noexcept;

        static Taxon *classify(std::vector<Taxon *> &taxonomy) noexcept;

        static std::vector<Taxon *> phenogram(std::vector<Taxon *> &taxonomy);

        static std::string Taxonomy::ascii(std::vector<Taxon *> &phenogram,
                                           std::vector<std::string> markers);

        static std::string Taxonomy::ascii(std::vector<Taxon *> &phenogram);
    };

    class Phenogram : public std::vector<Taxon *>
    {
        explicit operator json() const noexcept;
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
            : Subset(sample, true), parent(nullptr), events(sample.events)
        {
            for (Event event = 0; event < sample.events; event++)
                for (Measurement measurement = 0; measurement < sample.measurements; measurement++)
                    if (sample(event, measurement) < 0 || sample(event, measurement) > 1)
                    {
                        member(event, false);
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
            const Subset &subset) : Subset(sample, subset), parent(parent){};

        ~SampleSubset()
        {
            for (auto &child : children)
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
            Measurement X, Y;
            double X_KLD = std::numeric_limits<double>::infinity();
            double Y_KLD = std::numeric_limits<double>::infinity();
        } fallback;
        Measurement qualifying = 0;

        Request(
            Analysis<ClientSample> *const analysis,
            const ClientSample &sample,
            SampleSubset<ClientSample> *subset,
            const Parameters &parameters) noexcept
            : analysis(analysis), sample(sample), subset(subset), Lysis(parameters, subset->events, sample.measurements)
        {
            qualifying = analysis->qualifying;
            ID = ++analysis->uniques;
        };

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
        std::recursive_mutex mutex;
        static std::mt19937_64 generate; // not thread safe

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

        void start(
            const json &encoded){};

        void finish(
            const json &encoded)
        {
            Key request_key; // from JSON
            Request<ClientSample> *request = requests.find(request_key)->second;
        }

        Pursuer(
            const Parameters &parameters,
            int threads) noexcept
            : parameters(parameters), workers(threads < 0 ? std::thread::hardware_concurrency() : threads)
        {
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

    template <class ClientSample>
    std::mt19937_64 Pursuer<ClientSample>::generate(EPP::random());

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
        Count projections = 0, passes = 0, clusters = 0, graphs = 0, merges = 0;

        const Lysis *operator()(int i) const noexcept
        {
            return lysis[i];
        }

        json gating()
        {
            if (!lysis_unique)
                for (Lysis *ly : lysis)
                    ly->ID = ++uniques;
            lysis_unique = true;
            return lysis.front()->gating(parameters.tolerance);
        }

        Taxon *classify()
        {
            Taxonomy::classify(taxonomy);
            if (!lysis_unique)
                for (Lysis *ly : lysis)
                    ly->ID = ++uniques;
            lysis_unique = true;
            if (!taxon_unique)
            {
                for (Taxon *tax : taxonomy)
                    tax->ID = ++uniques;
                for (Taxon *tax : taxonomy)
                    if (tax->subset)
                        tax->subset->taxon = tax->ID;
            }
            taxon_unique = true;
            return taxonomy.back();
        }

        std::vector<Taxon *> phenogram()
        {
            return Taxonomy::phenogram(taxonomy);
        }

        Count types() noexcept { return _types; }

        Count size() noexcept
        {
            return (Count)lysis.size();
        }

        bool complete() noexcept
        {
            std::lock_guard<std::mutex> lock(mutex);
            return lysis.size() == requests;
        }

        bool censor(Measurement measurement) const noexcept
        {
            if (measurement < sample.measurements)
                return censored[measurement];
            else
                return true;
        }

        void wait() noexcept
        {
            std::unique_lock<std::mutex> lock(mutex);
            if (lysis.size() < requests)
                progress.wait(lock);
        }

        ~Analysis()
        {
            for (auto &ly : lysis)
                delete ly;
            for (auto &tax : taxonomy)
                delete tax;
        }

    protected:
        std::mutex mutex;
        std::condition_variable progress;
        std::chrono::time_point<std::chrono::steady_clock> begin, end;
        std::vector<Request<ClientSample> *> lysis;
        std::vector<Taxon *> taxonomy;
        bool *censored;
        Measurement qualifying = 0;
        volatile Count requests = 0, _types = 0;
        Unique uniques = 0;
        bool lysis_unique = false, taxon_unique = false;

        Lysis *lyse(SampleSubset<ClientSample> *subset)
        {
            Request<ClientSample> *request = new Request<ClientSample>(this, this->sample, subset, parameters);
            request->status = EPP_no_cluster;
            pursuer->start(request);

            PursueProjection<ClientSample>::start(request);

            std::lock_guard<std::mutex> lock(mutex);
            ++requests;

            return request;
        }

        Lysis *characterize(Request<ClientSample> *request)
        {
            request->status = EPP_characterized;
            pursuer->start(request);

            Worker<ClientSample>::enqueue(
                new CharacterizeSubset<ClientSample>(request));

            return request;
        }

        Lysis *characterize(SampleSubset<ClientSample> *subset)
        {
            Request<ClientSample> *request = new Request<ClientSample>(this, this->sample, subset, parameters);

            std::lock_guard<std::mutex> lock(mutex);
            ++requests;

            return characterize(request);
        }

        void finish(
            Request<ClientSample> *request)
        {
            Lysis *ly;
            switch (request->status)
            {
            case EPP_success:
            {
                unsigned int threshold = std::max(
                    std::max(
                        (unsigned int)(request->analysis->parameters.min_relative * this->sample.events),       // relative to current sample
                        request->analysis->parameters.min_events),                                              // absolute event count
                    (unsigned int)(request->analysis->parameters.sigma * request->analysis->parameters.sigma)); // algorithim limit

                SampleSubset<ClientSample> *child = new SampleSubset<ClientSample>(this->sample, request->subset, request->in());
                child->X = request->X();
                child->Y = request->Y();
                child->events = request->in_events();
                child->polygon = request->in_polygon(parameters.tolerance);
                request->subset->children.push_back(child);
                if (request->analysis->parameters.recursive && request->in_events() > threshold)
                {
                    ly = lyse(child);
                    ly->parent = request;
                    ly->events = request->in_events();
                }
                else
                {
                    ly = characterize(child);
                    ly->parent = request;
                }
                ly->in_set = true;
                request->children.push_back(ly);

                child = new SampleSubset<ClientSample>(this->sample, request->subset, request->out());
                child->X = request->X();
                child->Y = request->Y();
                child->events = request->out_events();
                child->polygon = request->out_polygon(parameters.tolerance);
                request->subset->children.push_back(child);
                if (request->analysis->parameters.recursive && request->out_events() > threshold)
                {
                    ly = lyse(child);
                    ly->events = request->out_events();
                    ly->parent = request;
                }
                else
                {
                    ly = characterize(child);
                    ly->parent = request;
                }
                ly->in_set = false;
                request->children.push_back(ly);

                std::lock_guard<std::mutex> lock(mutex);
                lysis.push_back(request);
                break;
            }

            case EPP_no_cluster:
            {
                characterize(request);
                break;
            }

            case EPP_characterized:
            {
                std::lock_guard<std::mutex> lock(mutex);
                lysis.push_back(request);
                taxonomy.push_back(new Taxon(request));
                ++_types;
                break;
            }

            default:
                assert(("oops", false));
            }

            std::lock_guard<std::mutex> lock(mutex);
            compute_time += request->milliseconds;
            projections += request->projections;
            passes += request->passes;
            clusters += request->clusters;
            graphs += request->graphs;
            merges += request->merges;
            end = std::chrono::steady_clock::now();
            milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

            progress.notify_all();
        };

        Analysis(
            Pursuer<ClientSample> *pursuer,
            const ClientSample &sample,
            const Parameters &parameters) : pursuer(pursuer), sample(sample), parameters(parameters)
        {
            begin = std::chrono::steady_clock::now();
            censored = new bool[sample.measurements];
            const std::vector<Measurement> *c = &parameters.censor;
            for (Measurement measurement = 0; measurement < sample.measurements; ++measurement)
                if (!(std::find(c->begin(), c->end(), measurement) != c->end()))
                {
                    ++qualifying;
                    censored[measurement] = false;
                }
                else
                    censored[measurement] = true;
        };
    };
}

#include "sample.h"

#endif /* _EPP_CLIENT_H */
