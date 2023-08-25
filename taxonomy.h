/*
 * Developer: Wayne Moore <wmoore@stanford.edu>
 * Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
 * License: BSD 3 clause
 */
#ifndef _EPP_TAXONOMY_H
#define _EPP_TAXONOMY_H 1

#include <string>

#include "client.h"
#include "cholesky.h"

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

    double Taxon::walk(Taxon *supertaxon, std::vector<Taxon *> &phenogram)
    {
        if (supertaxon)
        {
            this->depth = supertaxon->depth + this->dissimilarity;
            this->rank = supertaxon->rank + 1;
        }
        else
        {
            this->depth = 0;
            this->rank = 0;
        }
        double maxd = this->depth;
        for (Taxon *tax : this->subtaxa)
        {
            maxd = std::max(maxd, tax->walk(this, phenogram));
        }
        phenogram.push_back(this);
        return maxd;
    }

    Taxon *Taxonomy::classify(std::vector<const Lysis *> &types)
    {
        std::deque<Similarity> similarities;
        std::forward_list<Taxon *> unclassified;
        Taxon *taxp = nullptr;
        // prime the pump with the characterized subsets
        for (const Lysis *subset : types)
        {
            taxp = new Taxon(subset);
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
            // if the disimilarity goes down the two were closer to
            // each other than to anthing else so merge them
            if (sim.dissimilarity < least_dissimilar)
                taxp = new Taxon(sim.red, sim.blue, true); // merge
            else
            {
                taxp = new Taxon(sim.red, sim.blue, false); // new taxon
                least_dissimilar = sim.dissimilarity;
            }
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
        // post order walk of the taxonomy tree to find the rank and depth
        // of each taxon and to fill out the vector in (reverse) order
        Taxon *root = taxp;
        std::vector<Taxon *> phenogram;
        assert(root->types == types.size());
        phenogram.reserve(root->taxa);
        double depth = root->walk(nullptr, phenogram);
        Count idx = 0;
        Count uniques = 0;
        for (Taxon *taxp : phenogram)
        {
            taxp->ID = ++uniques;
            // normalize the depth
            taxp->depth *= root->height / depth;
            // find the mid point
            if (taxp->isSpecific())
                taxp->index = idx++;
            else
            {
                double sum = 0;
                for (Taxon *taxq : taxp->subtaxa)
                {
                    taxq->supertaxon = taxp->ID;
                    sum += taxq->index;
                }
                taxp->index = sum / taxp->subtaxa.size();
            }
        }
        assert(idx == root->types);
        // locate the siblings
        auto one = phenogram.rbegin();
        auto two = ++phenogram.rbegin();
        std::vector<bool> sibling(taxp->height + 1, false);
        (*one)->sibling = sibling;
        for (; two != phenogram.rend(); ++one, ++two)
        {
            if ((*two)->rank > (*one)->rank)
            {
                if ((*one)->rank > 0)
                    sibling[(*one)->rank - 1] = true;
                for (Count r = (*one)->rank; r <= (*two)->rank; ++r)
                    sibling[r] = false;
            }
            else if ((*two)->rank == (*one)->rank)
                sibling[(*two)->rank - 1] = true;
            (*two)->sibling = sibling;
        }
        return root;
    }

    Taxon::Taxon(const Lysis *subset)
        : subtaxa(0), markers(subset->markers), covariance(subset->covariance), subset(subset->ID),
          divergence(subset->divergence), dissimilarity(std::numeric_limits<double>::quiet_NaN()),
          population(subset->events), height(0), taxa(1), types(1) {}

    Taxon::Taxon(
        Taxon *red,
        Taxon *blue,
        bool merge)
        : subtaxa(0), markers(red->markers.size(), 0), covariance(0), subset(0),
          divergence(std::numeric_limits<double>::quiet_NaN()), dissimilarity(std::numeric_limits<double>::quiet_NaN()),
          population(0), height(0), taxa(1), types(0)
    {
        if (merge)
        {
            if (red->isSpecific())
                subtaxa.push_back(red);
            else
            {
                for (Taxon *tax : red->subtaxa)
                    subtaxa.push_back(tax);
                red->subtaxa.clear();
                delete red;
            }
            if (blue->isSpecific())
                subtaxa.push_back(blue);
            else
            {
                for (Taxon *tax : blue->subtaxa)
                    subtaxa.push_back(tax);
                blue->subtaxa.clear();
                delete blue;
            }
        }
        else
        {
            subtaxa.push_back(red);
            subtaxa.push_back(blue);
        }

        for (Taxon *tax : subtaxa)
            population += tax->population;

        double *weight = new double[subtaxa.size()];
        double sum = 0;
        for (size_t i = 0; i < subtaxa.size(); ++i)
            sum += weight[i] = std::pow((double)subtaxa[i]->population, 1.0 / (double)markers.size());
        for (size_t i = 0; i < subtaxa.size(); ++i)
            for (size_t j = 0; j < markers.size(); ++j)
                markers[j] += weight[i] / sum * subtaxa[i]->markers[j];
        delete[] weight;
        for (Taxon *tax : subtaxa)
        {
            tax->supertaxon = this->ID;
            tax->dissimilarity = Taxonomy::cityBlockDistance(markers, tax->markers);
            height = std::max(height, tax->height + 1);
            taxa += tax->taxa;
            types += tax->types;
        }
        std::sort(this->subtaxa.begin(), this->subtaxa.end(),
                  [](Taxon *a, Taxon *b)
                  { return *b < *a; });
    }

    Taxon::~Taxon()
    {
        for (Taxon *taxp : subtaxa)
            delete taxp;
    }

    Taxon::operator json() const noexcept
    {
        json taxon;
        taxon["population"] = this->population;
        json markers = json::array();
        for (double m : this->markers)
            markers += m;
        taxon["markers"] = markers;
        taxon["ID"] = this->ID;
        if (this->supertaxon)
        {
            taxon["supertaxon"] = this->supertaxon;
            taxon["dissimilarity"] = this->dissimilarity;
        }
        taxon["depth"] = this->depth;
        taxon["rank"] = this->rank;
        taxon["height"] = this->height;
        taxon["index"] = this->index;
        taxon["taxa"] = this->taxa;
        taxon["types"] = this->types;
        json sibling = json::array();
        for (Count i = 0; i < this->rank; ++i)
            if (this->sibling.at(i))
                sibling += i;
        taxon["sibling"] = sibling;
        if (this->isSpecific())
        {
            taxon["gating"] = this->subset;
            taxon["divergence"] = this->divergence;
            json covariance = json::array();
            for (double c : this->covariance)
                covariance += c;
            taxon["covariance"] = covariance;
        }
        else
        {
            json subtaxa = json::array();
            for (Taxon *subtax : this->subtaxa)
                subtaxa += (json)*subtax;
            taxon["subtaxa"] = subtaxa;
        }
        return taxon;
    }

    Similarity::Similarity(
        Taxon *red,
        Taxon *blue)
        : dissimilarity(Taxonomy::cityBlockDistance(red->markers, blue->markers)),
          red(red), blue(blue) {}

    template <class ClientSample>
    class CharacterizeSubset : public Work<ClientSample>
    {
        Unique *classification;
        float *mahalanobis;
        const Event &events;
        const Measurement &measurements;

    public:
        CharacterizeSubset(
            Request<ClientSample> *request) noexcept
            : Work<ClientSample>(request), classification(request->analysis->classification.data()),
              mahalanobis(request->analysis->mahalanobis.data()), events(request->sample.events),
              measurements(request->sample.measurements) {}

        virtual ~CharacterizeSubset() = default;

        virtual void parallel() noexcept;

        virtual void serial() noexcept {}
    };

    template <class ClientSample>
    void CharacterizeSubset<ClientSample>::parallel() noexcept
    {
        // compute the mean vector for the phenotype
        double *mean = this->request->markers.data();
        for (Event event = 0; event < this->events; ++event)
            if (this->subset->contains(event))
                for (Measurement i = 0; i < this->measurements; ++i)
                    mean[i] += this->sample(event, i);
        for (Measurement i = 0; i < this->measurements; ++i)
            mean[i] /= this->request->events;

        // compute the covariance matrix in lower trangular form
        this->request->covariance.resize((this->measurements * (this->measurements + 1)) / 2, 0);
        double *cov = this->request->covariance.data();
        for (Event event = 0; event < this->events; ++event)
            if (this->subset->contains(event))
                for (Measurement i = 0, k = 0; i < this->measurements; ++i)
                    for (Measurement j = 0; j <= i; ++j, ++k)
                        cov[k] += (this->sample(event, i) - mean[i]) * (this->sample(event, j) - mean[j]);
        for (Measurement i = 0, k = 0; i < this->measurements; ++i)
            for (Measurement j = 0; j <= i; ++j, ++k)
                cov[k] /= this->request->events - 1;

        // remove censored dimensions
        for (Measurement i = 0; i < this->measurements; ++i)
            if (this->request->analysis->censor(i))
                mean[i] = 0;
        for (Measurement i = 0, k = 0; i < this->measurements; ++i)
            if (this->request->analysis->censor(i))
            {
                for (Measurement j = 0; j < i; ++j, ++k)
                    cov[k] = 0;
                cov[k++] = 1;
            }
            else
                for (Measurement j = 0; j <= i; ++j, ++k)
                    if (this->request->analysis->censor(j))
                        cov[k] = 0;

        // compute the inverse of the covariance also in lower trangular form
        this->request->invcovariance.resize((this->measurements * (this->measurements + 1)) / 2, 0);
        double *inv = this->request->invcovariance.data();
        double *work = new double[(this->measurements * (this->measurements + 1)) / 2];
        int ifault, nullty;
        syminv(cov, measurements, inv, work, &nullty, &ifault);
        assert(ifault == 0);
        delete[] work;
        if (nullty != 0)    // can't invert
            return;

        // compute the Mahalanobis distance and Kullback-Leibler divergence
        double KLD = 0, NQ = 0;
        // exclude censored
        for (Measurement i = 0; i < this->measurements; ++i)
            if (this->request->analysis->censor(i))
                inv[((i + 1) * (i + 2)) / 2 - 1] = 0;
        // d2 is chi squared
        double maha_mean = (double)(this->measurements - this->parameters.censor.size());
        double maha_sd = sqrt(maha_mean);
        for (Event event = 0; event < this->events; ++event)
            if (this->subset->contains(event))
            {
                double d2 = 0;
                for (Measurement i = 0, k = 0; i < this->measurements; ++i)
                    if (!this->request->analysis->censor(i))
                    {
                        for (Measurement j = 0; j < i; ++j, ++k)
                            if (!this->request->analysis->censor(j))
                                d2 += 2 * inv[k] * (this->sample(event, i) - mean[i]) * (this->sample(event, j) - mean[j]);
                        d2 += inv[k++] * (this->sample(event, i) - mean[i]) * (this->sample(event, i) - mean[i]);
                    }
                    else
                        k += i + 1;
                this->classification[event] = this->request->ID;
                this->mahalanobis[event] = (float)((1 + erf((d2 - maha_mean) / maha_sd / sqrt2)) / 2);
                KLD += d2 / 2;
                NQ += exp(-d2 / 2);
            }
        KLD /= this->request->events;
        KLD -= log(this->request->events / NQ);
        assert(KLD > 0);
        this->request->divergence = KLD;
    }

    void Phenogram::toHtml(
        Taxon *taxonomy,
        std::vector<std::string> &markers,
        std::ofstream &html)
    {
        std::ifstream crutch("template.html", std::ios::in);
        std::string line;
        while (std::getline(crutch, line))
        {
            if (line == "        <!--markers-->")
                break;
            html << line << std::endl;
        }
        json marks = json::array();
        for (const std::string marker : markers)
            marks += marker;
        html << marks << std::endl;
        while (std::getline(crutch, line))
        {
            if (line == "        <!--taxonomy-->")
                break;
            html << line << std::endl;
        }
        html << (json)*taxonomy << std::endl;
        while (std::getline(crutch, line))
            html << line << std::endl;
        crutch.close();
    }
}
#endif /* _EPP_TAXONOMY_H */
