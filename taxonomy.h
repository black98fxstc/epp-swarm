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
    double Taxonomy::cityBlockDistance(
        std::vector<double> &red,
        std::vector<double> &blue) noexcept
    {
        double d = 0;
        auto rp = red.begin();
        auto bp = blue.begin();
        while (rp != red.end())
            d += std::abs(*rp++ - *bp++);
        return d;
    }

    Taxon* Taxonomy::classify(std::vector<Taxon *> &taxonomy) noexcept
    {
        std::deque<Similarity> similarities;
        std::forward_list<Taxon *> unclassified;
        // prime the pump with the characterized subset
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
            // if the disimilarity goes down the two were closer to 
            // each other than to anthing else so merge them
            if (sim.dissimilarity < least_dissimilar)
                taxp = new Taxon(sim.red, sim.blue, true);  // merge
            else
            {
                taxp = new Taxon(sim.red, sim.blue, false); // new taxon
                least_dissimilar = sim.dissimilarity;
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

    Phenogram Taxonomy::phenogram(std::vector<Taxon *> &taxonomy)
    {
        Phenogram phenogram;
        phenogram.reserve(taxonomy.size());

        // post order walk of the taxonomy tree to find the rank, hight and depth
        // of each taxon and to fill out the vector in (reverse) order
        double depth = taxonomy.back()->walk(phenogram);
        // normalize the depth
        for (Taxon *tax : taxonomy)
            tax->depth /= depth;
        // locate the connectors
        auto one = phenogram.begin();
        auto two = ++phenogram.begin();
        std::vector<bool> connect(taxonomy.back()->height + 1, false);
        (*one)->connect = connect;
        for (; two != phenogram.end(); ++one, ++two)
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
        // we want the root first in the returned vector
        std::reverse(phenogram.begin(), phenogram.end());
        return phenogram;
    }

    std::string Taxonomy::ascii(
        std::vector<Taxon *> &phenogram)
    {
        std::vector<char> mark(phenogram.front()->markers.size());
        std::vector<int> pos(phenogram.front()->height + 1);
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
            line.push_back('X');
            char c;
            if (tax->isSpecific())
                c = '+';
            else
                switch (tax->rank % 2)
                {
                case 0:
                    c = '-';
                    break;
                case 1:
                    c = '~';
                    break;
                }

            while (p++ < 67)
                line.push_back(c);
            if (tax->isSpecific())
                line.push_back('*');
            else
                line.push_back(' ');
            for (char &c : mark)
                line.push_back(c);
            line.push_back(' ');
            int pop = (int)(40.0 * tax->population / population);
            c = 'X';
            if (pop < 1)
            {
                c = 'o';
                pop = (int)(400.0 * tax->population / population);
            }
            while (pop-- > 0)
                line.push_back(c);
            line.push_back('\n');
        }
        return line;
    }
    
    std::string Taxonomy::ascii(
        std::vector<Taxon *> &phenogram,
        std::vector<std::string> markers)
    {
        std::string line;
        int p = 0, m = (int)markers.size();
        for (int i = 0; i < markers.size(); ++i)
        {
            p = 0;
            while (++p < 68)
                line.push_back(' ');
            for (int j = m - i; j > 0; --j)
                line.push_back(' ');
            line.push_back('/');
            for (int j = i; j > 0; --j)
                line.push_back('/');
            line.push_back(' ');
            line += markers[i];
            line.push_back('\n');
        }
        return line + ascii(phenogram);
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
        taxon["ID"] = this->ID;
        taxon["dissimilarity"] = this->dissimilarity;
        taxon["depth"] = this->depth;
        taxon["rank"] = this->rank;
        taxon["height"] = this->height;
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

    Phenogram::operator json() const noexcept
    {
        json phenogram;
        for (Taxon *tax : *this)
        {
            json taxon;
            taxon["population"] = tax->population;
            json markers;
            for (double m : tax->markers)
                markers += m;
            taxon["markers"] = markers;
            taxon["ID"] = tax->ID;
            if (tax->supertaxon)
            {
                taxon["supertaxon"] = tax->supertaxon->ID;
                taxon["dissimilarity"] = tax->dissimilarity;
            }
            taxon["depth"] = tax->depth;
            taxon["rank"] = tax->rank;
            taxon["height"] = tax->height;
            json connect;
            for (int i = 0; i < tax->rank; ++i)
                if (tax->connect.at(i))
                    connect += i;
            taxon["connect"] = connect;
            json subtaxa;
            for (Count i = 0; i < tax->subtaxa.size(); ++i)
                subtaxa[i] = tax->subtaxa.at(i)->ID;
            if (subtaxa.size() > 0)
                taxon["subtaxa"] = subtaxa;
            if (tax->subset)
                taxon["gating"] = tax->subset->ID;
            phenogram += taxon;
        }
        return phenogram;
    }

    Similarity::Similarity(
        Taxon *red,
        Taxon *blue)
        : red(red), blue(blue),
          dissimilarity(Taxonomy::cityBlockDistance(red->markers, blue->markers)) {}

    template <class ClientSample>
    class CharacterizeSubset : public Work<ClientSample>
    {
        std::vector<double> &markers;
        Unique *classification;
        Event &events;

    public:
        CharacterizeSubset(
            Request<ClientSample> *request) noexcept
            : Work<ClientSample>(request), events(request->events), markers(request->markers),
            classification(request->analysis->classification) {}

        virtual void parallel() noexcept;

        virtual void serial() noexcept {}
    };

    template <class ClientSample>
    void CharacterizeSubset<ClientSample>::parallel() noexcept
    {
        double *data = this->markers.data();
        for (Event event = 0; event < this->sample.events; ++event)
            if (this->subset->contains(event))
            {
                for (Measurement measurement = 0; measurement < this->sample.measurements; ++measurement)
                    data[measurement] += this->sample(event, measurement);
                classification[event] = this->request->ID;
            }
        for (Measurement measurement = 0; measurement < this->sample.measurements; ++measurement)
            if (this->request->analysis->censor(measurement))
                data[measurement] = 0;
            else
                data[measurement] /= events;
    }
}

#endif /* _EPP_TAXONOMY_H */
