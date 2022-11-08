/*
 * Developer: Wayne Moore <wmoore@stanford.edu>
 * Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
 * License: BSD 3 clause
 */
#ifndef _EPP_TAXONOMY_H
#define _EPP_TAXONOMY_H 1

namespace EPP
{
    class Similarity
    {
    public:
        double dissimilarity;
        Taxon *red, *blue;

        bool operator<(const Similarity &other) { return this->dissimilarity > other.dissimilarity; }

        Similarity(
            Taxon *red,
            Taxon *blue);
    };

    class Taxonomy
    {
    public:
        static double cityBlockDistance(
            std::vector<double> &red,
            std::vector<double> &blue)
        {
            double d = 0;
            for (auto rp = red.begin(), bp = blue.begin(); rp != red.end(); ++rp, ++bp)
                d += abs(*rp - *bp);
            return d;
        }

        static Taxon *classify(std::vector<Taxon> &taxonomy)
        {
            if (taxonomy.back().subtaxa.size() > 0)
                return &taxonomy.back();

            std::deque<Similarity> similarities;
            std::forward_list<Taxon *> unclassified;
            // prime the pump with the EPP found subsets
            for (Taxon &taxp : taxonomy)
            {
                for (Taxon *taxq : unclassified)
                    similarities.emplace(similarities.begin(), &taxp, taxq);
                unclassified.push_front(&taxp);
            }
            while (!similarities.empty())
            {
                // find the most similar unclassified taxa
                std::make_heap(similarities.begin(), similarities.end());
                std::pop_heap(similarities.begin(), similarities.end());
                Similarity similar = similarities.back();
                similarities.pop_back();
                // create the new naxon
                taxonomy.emplace_back(similar.red, similar.blue);
                Taxon *taxp = &taxonomy.back();
                // remove any similarites mooted by this
                for (auto sp = similarities.begin(); sp != similarities.end();)
                    if (sp->red == similar.red || sp->red == similar.blue || sp->blue == similar.red || sp->blue == similar.blue)
                        similarities.erase(sp);
                    else
                        ++sp;
                // remove the now classified taxa
                unclassified.remove_if([&](Taxon *taxq)
                                       { return taxq == similar.red || taxq == similar.blue; });
                // compute the new similarities and list the new taxon as unclassified
                for (Taxon *taxq : unclassified)
                    similarities.emplace(similarities.begin(), taxp, taxq);
                if (!similarities.empty())
                    unclassified.push_front(taxp);
            }
            assert(unclassified.empty());
            return &taxonomy.back();
        }
    };

    Taxon::Taxon(Lysis *subset) : supertaxon(nullptr), subtaxa(0), subset(subset)
    {
        population = subset->events;
        markers = subset->markers;
    }

    Taxon::Taxon(Taxon *red, Taxon *blue) : supertaxon(nullptr), subset(nullptr)
    {
        population = red->population + blue->population;
        double p = ((double)red->population) / ((double)(red->population + blue->population));
        for (auto rp = red->markers.begin(), bp = blue->markers.begin(); rp < red->markers.end(); ++rp, ++bp)
            markers.push_back(p * *rp + (1 - p) * *bp);
        subtaxa.push_back(red);
        subtaxa.push_back(blue);
        red->supertaxon = this;
        blue->supertaxon = this;
        red->dissimilarity = Taxonomy::cityBlockDistance(this->markers, red->markers);
        blue->dissimilarity = Taxonomy::cityBlockDistance(this->markers, blue->markers);
    }

    Similarity::Similarity(
        Taxon *red,
        Taxon *blue) : red(red), blue(blue),
                       dissimilarity(Taxonomy::cityBlockDistance(red->markers, blue->markers)) {}

    template <class ClientSample>
    class CharacterizeSubset : public Work<ClientSample>
    {
        Event &events;
        std::vector<double> &markers;

    public:
        CharacterizeSubset(
            Request<ClientSample> *request) noexcept
            : Work<ClientSample>(request), events(request->subset->events), markers(request->markers) {}

        virtual void parallel() noexcept;

        virtual void serial() noexcept {}
    };

    template <class ClientSample>
    void CharacterizeSubset<ClientSample>::parallel() noexcept
    {
        double *data = markers.data();
        for (Event event = 0; event < sample.events; event++)
            if (subset->contains(event))
            {
                for (Measurement measurement = 0; measurement < sample.measurements; ++measurement)
                    data[measurement] += this->sample(event, measurement);
                ++events;
            }
        for (Measurement measurement = 0; measurement < sample.measurements; ++measurement)
            if (this->request->analysis->censor(measurement))
                data[measurement] = 0;
            else
                data[measurement] /= events;
    }
}

#endif /* _EPP_TAXONOMY_H */
