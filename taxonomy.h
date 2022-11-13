/*
 * Developer: Wayne Moore <wmoore@stanford.edu>
 * Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
 * License: BSD 3 clause
 */
#ifndef _EPP_TAXONOMY_H
#define _EPP_TAXONOMY_H 1

#include <string>

namespace EPP
{
    class Taxon
    {
        friend class Taxonomy;

    public:
        std::vector<Taxon *> subtaxa;
        std::vector<double> markers;
        std::vector<bool> connect;
        Event population;
        Taxon *supertaxon;
        Lysis *subset;
        double dissimilarity, depth;
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
            std::vector<double> &blue)
        {
            double d = 0;
            auto rp = red.begin();
            auto bp = blue.begin();
            while (rp != red.end())
                d += abs(*rp++ - *bp++);
            return d;
        }

        static Taxon *classify(std::vector<Taxon *> &taxonomy)
        {
            if (taxonomy.back()->subtaxa.size() > 0)
                return taxonomy.back();

            std::deque<Similarity> similarities;
            std::forward_list<Taxon *> unclassified;
            // prime the pump with the subset found subsets
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
                similarities.pop_back();
                // create the new naxon
                Taxon *taxp;
                double dissimilarity = sim.dissimilarity;
                if (dissimilarity < least_dissimilar)
                    taxp = new Taxon(sim.red, sim.blue, true);
                else
                {
                    taxp = new Taxon(sim.red, sim.blue, false);
                    least_dissimilar = dissimilarity;
                }
                taxonomy.push_back(taxp);
                // remove any similarites mooted by this
                for (auto sp = similarities.begin(); sp != similarities.end();)
                    if (sp->red == sim.red || sp->red == sim.blue || sp->blue == sim.red || sp->blue == sim.blue)
                        sp = similarities.erase(sp);
                    else
                        ++sp;
                // remove the now classified taxa
                unclassified.remove_if([sim](Taxon *taxq)
                                       { return taxq == sim.red || taxq == sim.blue; });
                // compute the new similarities and list the new taxon as unclassified
                for (Taxon *taxq : unclassified)
                    similarities.emplace_back(taxp, taxq);
                unclassified.push_front(taxp);
            }
            return taxonomy.back();
        }

        static std::vector<Taxon *> phenogram(std::vector<Taxon *> &taxonomy)
        {
            std::vector<Taxon *> phenogram;
            phenogram.reserve(taxonomy.size());
            
            double depth = taxonomy.back()->walk(phenogram);
            // normalize the 
            for (Taxon *tax : taxonomy)
                tax->depth /= depth;
            auto one = phenogram.begin();
            auto two = ++phenogram.begin();
            std::vector<bool> connect(taxonomy.back()->height + 1, false);
            (*one)->connect = connect;
            for (int links = 0; two != phenogram.end(); ++one, ++two)
            {
                if ((*two)->rank > (*one)->rank)
                {
                    if ((*one)->rank > 0)
                        connect[(*one)->rank - 1] = true;
                    for (int r = (*one)->rank; r < connect.size(); ++r)
                        connect[r] = false;
                }
                else if ((*two)->rank < (*one)->rank)
                    if ((*two)->rank > 0)
                        connect[(*two)->rank - 1] = true;
                (*two)->connect = connect;
            }
            std::reverse(phenogram.begin(), phenogram.end());
            return phenogram;
        }

        static std::string Taxonomy::ascii(std::vector<Taxon *> &phenogram);
    };

    std::string Taxonomy::ascii(std::vector<Taxon *> &phenogram)
    {
        std::vector<char> mark(phenogram.front()->markers.size());
        std::vector<int> pos(100);
        EPP::Event population = phenogram.front()->population;
        std::string line;
        for (EPP::Taxon *tax : phenogram)
        {
            for (int i = 0; i < mark.size(); ++i)
                mark[i] = (char)('0' + (int)(10 * tax->markers[i]));

            pos[tax->rank] = ((int)(64 * tax->depth));
            int p = 0;
            for (int r = 0; r < tax->rank - 1; ++r)
            {
                while (p++ < pos[r])
                    line.push_back(' ');
                if (tax->connect[r])
                    line.push_back('|');
                else
                    line.push_back(' ');
            }
            if (tax->rank > 0)
            {
                while (p++ < pos[tax->rank - 1])
                    line.push_back(' ');
                line.push_back('+');
            }
            while (p++ < pos[tax->rank])
                line.push_back('=');
            line.push_back('x');
            while (p++ < 67)
                switch (tax->rank % 3)
                {
                case 0:
                    line.push_back('-');
                    break;
                case 1:
                    line.push_back('+');
                    break;
                case 2:
                    line.push_back('~');
                    break;
                }
            line.push_back('-');
            if (tax->isSpecific())
                line.push_back('*');
            else
                line.push_back(' ');
            for (char &c : mark)
                line.push_back(c);
            line.push_back(' ');
            int pop = (int)(40.0 * tax->population / population);
            char c = 'X';
            if (pop < 1)
            {
                c = 'x';
                pop = (int)(400.0 * tax->population / population);
            }
            while (pop-- > 0)
                line.push_back(c);
            line.push_back('\n');
        }
        return line;
    }

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
        std::sort(this->subtaxa.begin(), this->subtaxa.end(),
                  [](Taxon *a, Taxon *b)
                  { return *a < *b; });
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

    double Taxon::walk(std::vector<Taxon *> &phenogram)
    {
        if (this->supertaxon)
        {
            this->depth = this->supertaxon->depth + this->dissimilarity;
            this->rank = this->supertaxon->rank + 1;
        }
        else
        {
            this->depth = 0;
            this->rank = 0;
        }
        double maxd = this->depth;
        this->height = 0;
        if (this->isGeneric())
            for (Taxon *tax : this->subtaxa)
            {
                maxd = std::max(maxd, tax->walk(phenogram));
                this->height = std::max(this->height, tax->height + 1);
            }
        phenogram.push_back(this);
        return maxd;
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
