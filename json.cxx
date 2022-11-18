
/*
 * Developer: Wayne Moore <wmoore@stanford.edu>
 * Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
 * License: BSD 3 clause
 */
#include <iomanip>
#include <exception>

#include "client.h"

namespace EPP
{
    const static std::vector<std::string> Status_strings{
        "success", "characterized", "no_qualified", "no_cluster", "not_interesting", "error"};
    const static std::vector<std::string> Goal_strings{
        "best_separation", "best_balance"};

    size_t find_string(
        std::string string,
        const std::vector<std::string> strings)
    {
        for (size_t i = 0; i < strings.size(); i++)
            if (string == strings[i])
                return i;
        throw std::runtime_error("enumeration string did not match");
    };

    Key::operator json() const noexcept
    {
        std::ostringstream hex;
        hex << std::hex << std::setw(2) << std::setfill('0');
        for (const std::uint8_t &value : this->bytes)
            hex << value;
        return hex.str();
    }

    Key &Key::operator=(const json &encoded)
    {
        // std::istringstream hex((std::string)encoded); problematic nlohmann cast
        std::string hex_encoded = encoded;
        std::istringstream hex(hex_encoded);
        hex >> std::hex >> std::setw(2);
        for (std::uint8_t &value : this->bytes)
            hex >> value;
        return *this;
    }

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

    Parameters &Parameters::operator=(const json &encoded) noexcept
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

    Polygon::operator json () const noexcept
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

    Lysis:: operator json() const noexcept
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
