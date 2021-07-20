#include <iomanip>
#include <exception>

#include "client.h"

namespace EPP
{
    const static std::vector<std::string> Status_strings 
    {
        "success", "no_qualified", "no_cluster", "not_interesting", "error"
    };
    const static std::vector<std::string> Goal_strings = 
    { 
        "best_separation", "best_balance" 
    };

    int find_string(
        std::string string,
        const std::vector<std::string> strings)
    {
        for (int i = 0; i < strings.size(); i++)
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
        std::istringstream hex((std::string)encoded);
        hex >> std::hex >> std::setw(2);
        for (std::uint8_t &value : this->bytes)
            hex >> value;
        return *this;
    }

	Sample::operator json() noexcept
	{
        json encoded;
        encoded["measurements"] = this->measurements;
        encoded["events"] = this->events;
        encoded["key"] = (json)key();
		return encoded;
	}

	Sample &Sample::operator=(const json &encoded)
	{
        // this->measurements = encoded["measurements"];
        // this->events = encoded["events"];
        // this->_key = encoded["key"];
		return *this;
	}

	Subset::operator json() noexcept
	{
        json encoded;
        encoded["key"] = (json)key();
		return encoded;
	}

	// Subset &Subset::operator=(const json &encoded)
	// {
    //     this->_key = encoded["key"];
	// 	return *this;
	// }

	Parameters::operator json() const noexcept
	{
        json parameters;
        parameters["N"] = this->N;
        parameters["W"] = this->W;
        parameters["sigma"] = this->sigma;
        parameters["goal"] = Goal_strings[this->goal];
        parameters["finalists"] = this->finalists;
        json kld;
        kld["Normal2D"] = this->kld.Normal2D;
        kld["Normal1D"] = this->kld.Normal1D;
        kld["Exponential1D"] = this->kld.Exponential1D;
        parameters["KLD"] = kld;
        json censor;
        for (int i = 0; i < this->censor.size(); i++)
            censor[i] = this->censor[i];
        parameters["censor"] = censor;
        parameters["suppress_in_out"] = this->suppress_in_out;
		return parameters;
	}

	Parameters &Parameters::operator=(const json &encoded)
	{
        this->W = encoded["W"];
        this->sigma = encoded["sigma"];
        this->goal = (Goal)find_string(encoded["goal"], Goal_strings);
        this->finalists = encoded["finalists"];
        this->kld.Normal2D = encoded["KLD"]["Normal2D]"];
        this->kld.Normal1D = encoded["KLD"]["Normal1D]"];
        this->kld.Exponential1D = encoded["KLD"]["Exponential1D]"];
        int n = encoded["censor"].size();
        if (n > 0)
        {
            this->censor.reserve(n);
            for (int i = 0; i < n; i++)
                this->censor.push_back(encoded["censor"][i]);
        }
        this->suppress_in_out = encoded["suppress_in_out"];
		return *this;
	}

    // _Request::operator json() const noexcept
    // {
    //     json request;
    //     return nullptr;
    // }

    // _Request &_Request::operator=(const json &encoded)
    // {
    //     return *this;
    // }

    Candidate::operator json()
    {
        json candidate;
        json separatrix;
        for (int i = 0; i < this->separatrix.size(); i++)
        {
            json point;
            point[0] = this->separatrix[i].i;
            point[1] = this->separatrix[i].j;
            separatrix[i] = point;
        }
        candidate["separatrix"] = separatrix;
        // candidate["in"] = this->in;
        // candidate["out"] = this->out;
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
        Polygon polygon(separatrix.size());
        for (int i = 0; i < separatrix.size(); i++)
        {
            json point = separatrix[i];
            polygon[i] = Point(point[0], point[1]);
        }
        this->separatrix = polygon;
        // this->in = encoded["in"];
        // this->out = encoded["out"];
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

    _Result::operator json()
    {
        json result;
        result["key"] = (json)this->key;
        json candidates;
        for (int i = 0; i < this->candidates.size(); i++)
            candidates[i] = this->candidates[i];
        result["candidates"]= candidates;
        json qualified;
        for (int i = 0; i < this->qualified.size(); i++)
            qualified[i] = this->qualified[i];
        result["qualified"] = qualified;
        result["milliseconds"] = this->milliseconds.count();
        result["projections"] = this->projections;
        result["passes"] = this->passes;
        result["clusters"] = this->clusters;
        result["graphs"] = this->graphs;
        return result;
    }

    _Result &_Result::operator=(const json &encoded)
    {
        this->key = encoded["key"];
        json candidates = encoded["candidates"];
        std::vector<Candidate> vector(candidates.size());
        for (int i = 0; i < candidates.size(); i++)
            vector[i] = candidates[i];
        this->candidates = vector;
        json qualified = encoded["qualified"];
        this->qualified.clear();
        for (json &measurement : qualified)
            this->qualified.push_back(measurement);
        this->milliseconds = std::chrono::duration<int, std::milli>(encoded["milliseconds"]);
        this->projections = encoded["projections"];
        this->passes = encoded["passes"];
        this->clusters = encoded["clusters"];
        this->graphs = encoded["graphs"];
        return *this;
    }
}
