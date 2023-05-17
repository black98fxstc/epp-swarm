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
            sum += weight[i] = std::pow(subtaxa[i]->population, 1.0 / (double)markers.size());
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
            : Work<ClientSample>(request), classification(request->analysis->classification),
              mahalanobis(request->analysis->mahalanobis), events(request->sample.events),
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
            if (this->request->analysis->censor(i))
                mean[i] = 0;
            else
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
                    else
                        cov[k] /= this->request->events - 1;

        // compute the inverse of the covariance also in lower trangular form
        this->request->invcovariance.resize((this->measurements * (this->measurements + 1)) / 2, 0);
        double *inv = this->request->invcovariance.data();
        double *work = new double[(this->measurements * (this->measurements + 1)) / 2];
        int ifault, nullty;
        syminv(cov, measurements, inv, work, &nullty, &ifault);
        assert(ifault == 0);
        assert(nullty == 0);
        for (Measurement i = 0; i < this->measurements; ++i)
            if (this->request->analysis->censor(i))
                inv[((i + 1) * (i + 2)) / 2 - 1] = 0;
        delete[] work;

        // compute the Mahalanobis distance and Kullback-Leibler divergence
        double KLD = 0, NQ = 0;
        // d2 is chi squared
        double maha_mean = (double)(this->measurements - this->parameters.censor.size());
        double maha_sd = sqrt(maha_mean);
        for (Event event = 0; event < this->events; ++event)
            if (this->subset->contains(event))
            {
                double d2 = 0;
                for (Measurement i = 0, k = 0; i < this->measurements; ++i)
                {
                    for (Measurement j = 0; j < i; ++j, ++k)
                        d2 += 2 * inv[k] * (this->sample(event, i) - mean[i]) * (this->sample(event, j) - mean[j]);
                    d2 += inv[k++] * (this->sample(event, i) - mean[i]) * (this->sample(event, i) - mean[i]);
                }
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

    std::string Taxonomy::ascii(
        std::vector<Taxon *> &phenogram)
    {
        std::vector<char> mark(phenogram.front()->markers.size());
        std::vector<int> pos(phenogram.front()->height + 1);
        EPP::Event population = phenogram.front()->population;
        std::string line;
        for (EPP::Taxon *tax : phenogram)
        {
            for (size_t i = 0; i < mark.size(); ++i)
                mark[i] = (char)('0' + (int)(10 * tax->markers[i]));

            pos[tax->rank] = ((int)(64 * tax->depth));
            int p = 0;
            for (Count r = 0; r + 1 < tax->rank; ++r)
            {
                while (p++ < pos[r])
                    line.push_back(' ');
                if (tax->sibling[r])
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
        for (size_t i = 0; i < markers.size(); ++i)
        {
            p = 0;
            while (++p < 68)
                line.push_back(' ');
            for (size_t j = m - i; j > 0; --j)
                line.push_back(' ');
            line.push_back('/');
            for (size_t j = i; j > 0; --j)
                line.push_back('/');
            line.push_back(' ');
            line += markers[i];
            line.push_back('\n');
        }
        return line + ascii(phenogram);
    }

    void Phenogram::toHtml2(
        Taxon *taxonomy,
        std::vector<std::string> &markers,
        std::ofstream &html)
    {
        std::ifstream crutch("test.html", std::ios::in);
        std::string line;
        while (std::getline(crutch, line))
        {
            if (line == "<!--break-->")
                break;
            html << line << std::endl;
        }
        crutch.close();

        json marks = json::array();
        for (const std::string marker : markers)
            marks += marker;

        html << "var markers = " << marks << ";" << std::endl;
        html << "var taxonomy = " << (json)*taxonomy << ";" << std::endl;

        html << "render(taxonomy,markers)";
        html << "</script>";
        html << "</body>";
        html << "</html>" << std::endl;
    }

    void Phenogram::toHtml(
        Taxon *taxonomy,
        std::vector<std::string> markers,
        std::ofstream &html)
    {
        html << "<html>";
        html << "<head>";
        html << "<style>";
        html << ".tooltip { position: relative; display: inline-block; } ";
        html << ".tooltip .tooltiptext { visibility: hidden; width: 120px; background: white; position: absolute; z-index: 1; } ";
        html << ".tooltip:hover .tooltiptext { visibility: visible; } ";
        html << "</style>";
        html << "<script src=\"https://code.jquery.com/jquery-3.6.4.min.js\" integrity=\"sha256-oP6HI9z1XaZNBrJURtCoUT5SUnxFr8s3BzRl+cbzUq8=\" crossorigin=\"anonymous\"></script>";
        html << "</head>";
        html << "<body style=\"font-family: monospace;\">" << std::endl;

        size_t marker_length = 0;
        int p = 0, m = (int)markers.size();
        while (++p < 66)
            html << "&nbsp;";
        for (size_t j = m; j > 0; --j)
            html << "&nbsp;";
        html << "&nbsp;&nbsp;<input type=\"checkbox\" checked id=\"all_or_none\">"
             << "All&nbsp;/&nbsp;None"
             << "</input><br/>" << std::endl;
        for (size_t i = 0; i < markers.size(); ++i)
        {
            p = 0;
            while (++p < 66)
                html << "&nbsp;";
            for (size_t j = m - i; j > 0; --j)
                html << "&nbsp;";
            html << "&#x2571;";
            for (size_t j = i; j > 0; --j)
                html << "&#x2571;";
            html << "&nbsp;<input type=\"checkbox\" class=\"allreagents\" checked id=\"reagent" << i << "label\">" << markers[i] << "</input><br/>" << std::endl;

            if (markers[i].length() > marker_length)
                marker_length = markers[i].length();
        }

        std::vector<double> mark(markers.size());
        double m_min = 1, m_max = 0;
        std::stack<Taxon *> stack;
        stack.push(taxonomy);
        while (!stack.empty())
        {
            Taxon *tax = stack.top();
            stack.pop();
            std::for_each(tax->subtaxa.begin(), tax->subtaxa.end(), [&](Taxon *tax)
                          { stack.push(tax); });
            for (size_t i = 0; i < mark.size(); ++i)
            {
                if (tax->markers[i] < m_min)
                    m_min = tax->markers[i];
                if (tax->markers[i] > m_max)
                    m_max = tax->markers[i];
            }
        }

        std::vector<int> pos(taxonomy->height + 1);
        EPP::Event population = taxonomy->population;
        stack.push(taxonomy);
        while (!stack.empty())
        {
            Taxon *tax = stack.top();
            stack.pop();
            std::for_each(tax->subtaxa.begin(), tax->subtaxa.end(), [&](Taxon *tax)
                          { stack.push(tax); });

            for (size_t i = 0; i < mark.size(); ++i)
                mark[i] = (tax->markers[i] - m_min) / (m_max - m_min);

            pos[tax->rank] = ((int)(64 * tax->depth));
            int p = 0;
            for (Count r = 0; r + 1 < tax->rank; ++r)
            {
                while (p++ < pos[r])
                    html << "&nbsp;";
                if (tax->sibling[r])
                    html << "&#x2503;";
                else
                    html << "&nbsp;";
            }
            if (tax->rank > 0)
            {
                while (p++ < pos[tax->rank - 1])
                    html << "&nbsp;";
                if (tax->sibling[tax->rank - 1])
                    html << "&#x2523;";
                else
                    html << "&#x2517;";
            }
            while (p++ < pos[tax->rank])
                html << "&#x2501;";
            if (tax->isSpecific())
                html << "&#x257E;";
            else
            {
                html << "&#x2513;<br/>" << std::endl;
                continue;
            }

            while (p++ < 65)
                html << "&#x2504;";
            if (tax->isSpecific())
                html << "&nbsp;";
            else
                html << "&nbsp;";

            html << "<span class=\"tooltip\">";
            for (size_t i = 0; i < markers.size(); ++i)
            {
                html << "<span class=\"reagent" << i << "value\" style=\"background-color: hsl(";
                html << (int)(240 * (1 - mark[i]));
                html << " 100% 50%); white-space: nowrap;\">&nbsp;</span>";
            }
            html << "<span class=\"tooltiptext\">";
            for (size_t i = 0; i < markers.size(); ++i)
            {
                html << "<span style=\"background-color: hsl(";
                html << (int)(240 * (1 - mark[i]));
                html << " 100% 50%);\">&nbsp;" << markers[i];
                for (size_t j = markers[i].length(); j < marker_length; ++j)
                    html << "&nbsp;";
                html << "&nbsp;</span><br/>";
            }
            html << "</span>";
            html << "</span>";
            html << "&nbsp;";
            html << "<span class=\"tooltip\">";
            html << "<span class=\"tooltiptext\" style=\"width: 40em;\">";
            char buffer[128];
            std::snprintf(buffer, 128, "%5.2g%% %8ld / %ld",
                          100.0 * tax->population / population, tax->population, population);
            html << buffer << "</span>";
            int pop = (int)(40.0 * tax->population / population);
            std::string c = "&#x2501;";
            if (pop < 1)
            {
                c = "&#x254C;";
                pop = (int)(400.0 * tax->population / population);
            }
            for (int i = 0; i < pop; ++i)
                html << c;
            for (int i = pop; i < 40; ++i)
                html << "&nbsp;";
            html << "</span>";
            html << "<br/>" << std::endl;
        }
        html << "Shorter horizontal distance is more similar. ";
        html << "Warmer colors are higher expression. ";
        html << "Population size, dashed lines are 1/10 of solid lines";
        html << "<script>" << std::endl;
        html << "function reagentChanged(event){if ($(\"#reagent\" + Number(event.data)  + \"label\").prop(\"checked\")) $(\".reagent\" + Number(event.data) + \"value\").css(\"visibility\", \"\");";
        html << "else $(\".reagent\" + Number(event.data) + \"value\").css(\"visibility\", \"hidden\"); };";
        html << "function allOrNone(event){if ($(\"#all_or_none\").prop(\"checked\")) $(\"input.allreagents\").prop(\"checked\",true).trigger(\"change\");";
        html << "else $(\"input.allreagents\").prop(\"checked\",false).trigger(\"change\"); };";
        html << "$(\"#all_or_none\").on(\"change\",allOrNone);";
        for (size_t i = 0; i < markers.size(); ++i)
            html << "$(\"#reagent" << i << "label\").on(\"change\", " << i << ",reagentChanged);";
        html << "</script>";
        html << "</body>";
        html << "</html>" << std::endl;
    }
}

#endif /* _EPP_TAXONOMY_H */
