/*
 * Developer: Wayne Moore <wmoore@stanford.edu>
 * Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
 * License: BSD 3 clause
 */
#ifndef _EPP_TAXONOMY_H
#define _EPP_TAXONOMY_H 1

namespace EPP
{
    class Taxon
    {
    public:
        Event population;
        std::vector<double> markers;
        Taxon *supertaxon;
        double dissimilarity;
        std::vector<Taxon *> subtaxa;
        Lysis *subset;
        Unique ID;

        bool isSpecific() const noexcept { return subtaxa.empty(); }
        bool isGeneric() const noexcept { return !subtaxa.empty(); }

        explicit operator json() const noexcept;

        Taxon(Lysis *subset);
        Taxon(Event population, std::vector<double> &markers) : population(population), markers(markers) {}
        Taxon(Taxon *red, Taxon *blue, bool merge = false);
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
    public:
        static double cityBlockDistance(
            std::vector<double> &red,
            std::vector<double> &blue)
        {
            double d = 0;
            auto rp = red.begin();
            auto bp = blue.begin();
            while (rp != red.end())
                d += abs(*rp++ - *bp++);
            return d;
        }

        static Taxon *classify(std::vector<Taxon *> &taxonomy, Event true_pop)
        {
            if (taxonomy.back()->subtaxa.size() > 0)
                return taxonomy.back();

            std::deque<Similarity> similarities;
            std::forward_list<Taxon *> unclassified;
            // prime the pump with the EPP found subsets
            for (Taxon *taxp : taxonomy)
            {
                for (Taxon *taxq : unclassified)
                    similarities.emplace_back(taxp, taxq);
                unclassified.push_front(taxp);
            }
            double least_dissimilar = 0;
            while (!similarities.empty())
            {
                // find the most similar unclassified taxa
                std::make_heap(similarities.begin(), similarities.end());
                std::pop_heap(similarities.begin(), similarities.end());
                Similarity sim = similarities.back();
                double dissimilarity = sim.dissimilarity;
                Taxon *red = sim.red;
                Taxon *blue = sim.blue;
                similarities.pop_back();
                // create the new naxon
                Taxon *taxp;
                if (dissimilarity < least_dissimilar)
                    taxp = new Taxon(red, blue, true);
                else
                {
                    taxp = new Taxon(red, blue, false);
                    least_dissimilar = dissimilarity;
                }
                taxonomy.push_back(taxp);
                // remove any similarites mooted by this
                for (auto sp = similarities.begin(); sp != similarities.end();)
                    if (sp->red == red || sp->red == blue || sp->blue == red || sp->blue == blue)
                        sp = similarities.erase(sp);
                    else
                        ++sp;
                // remove the now classified taxa
                unclassified.remove_if([red, blue](Taxon *taxq)
                                       { return taxq == red || taxq == blue; });
                // compute the new similarities and list the new taxon as unclassified
                for (Taxon *taxq : unclassified)
                    similarities.emplace_back(taxp, taxq);
                unclassified.push_front(taxp);
            }
            return taxonomy.back();
        }
    };

    Taxon::Taxon(Lysis *subset)
        : population(subset->events), markers(subset->markers),
          supertaxon(nullptr), dissimilarity(std::numeric_limits<double>::quiet_NaN()),
          subtaxa(0), subset(subset) {}

    Taxon::Taxon(
        Taxon *red,
        Taxon *blue,
        bool merge)
        : population(0), markers(red->markers.size(), 0),
          supertaxon(nullptr), dissimilarity(std::numeric_limits<double>::quiet_NaN()),
          subtaxa(0), subset(nullptr)
    {
        if (merge)
        {
            if (red->isSpecific())
                subtaxa.push_back(red);
            else
                for (Taxon *tax : red->subtaxa)
                    subtaxa.push_back(tax);
            if (blue->isSpecific())
                subtaxa.push_back(blue);
            else
                for (Taxon *tax : blue->subtaxa)
                    subtaxa.push_back(tax);
        }
        else
        {
            subtaxa.push_back(red);
            subtaxa.push_back(blue);
        }

        for (Taxon *tax : subtaxa)
            population += tax->population;
        for (Taxon *tax : subtaxa)
        {
            double p = ((double)tax->population) / ((double)population);
            auto my_marker = markers.begin();
            auto tax_marker = tax->markers.begin();
            while (my_marker != markers.end())
                *my_marker++ += p * *tax_marker++;
        }
        for (Taxon *tax : subtaxa)
        {
            tax->supertaxon = this;
            tax->dissimilarity = Taxonomy::cityBlockDistance(markers, tax->markers);
        }
    }

    Taxon::operator json() const noexcept
    {
        json taxon;
        taxon["population"] = this->population;
        json markers;
        for (size_t i = 0; i < this->markers.size(); ++i)
            markers[i] = this->markers[i];
        taxon["markers"] = markers;
        taxon["ID"] = ID;
        taxon["dissimilarity"] = this->dissimilarity;
        if (this->supertaxon)
            taxon["supertaxon"] = this->supertaxon->ID;
        json subtaxa;
        for (size_t i = 0; i < this->subtaxa.size(); ++i)
        {
            Taxon *tax = this->subtaxa.at(i);
            subtaxa[i] = (json)*tax;
        }
        if (subtaxa.size() > 0)
            taxon["subtaxa"] = subtaxa;
        if (this->subset)
            taxon["gating"] = this->subset->ID;
        return taxon;
    }

    Similarity::Similarity(
        Taxon *red,
        Taxon *blue)
        : red(red), blue(blue),
          dissimilarity(Taxonomy::cityBlockDistance(red->markers, blue->markers)) {}

    template <class ClientSample>
    class CharacterizeSubset : public Work<ClientSample>
    {
        Event &events;
        std::vector<double> &markers;

    public:
        CharacterizeSubset(
            Request<ClientSample> *request) noexcept
            : Work<ClientSample>(request), events(request->events), markers(request->markers) {}

        virtual void parallel() noexcept;

        virtual void serial() noexcept {}
    };

    template <class ClientSample>
    void CharacterizeSubset<ClientSample>::parallel() noexcept
    {
        double *data = this->markers.data();
        for (Measurement measurement = 0; measurement < this->sample.measurements; ++measurement)
            if (!this->request->analysis->censor(measurement))
                for (Event event = 0; event < this->sample.events; ++event)
                    if (subset->contains(event))
                        data[measurement] += this->sample(event, measurement);
        for (Measurement measurement = 0; measurement < this->sample.measurements; ++measurement)
            data[measurement] /= events;
    }
}

#endif /* _EPP_TAXONOMY_H */
