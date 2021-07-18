#include "client.h"

namespace EPP
{
    Key::operator json() const noexcept
    {
        return nullptr;
    }

    Key &Key::operator=(const json &encoded)
    {
        return *this;
    }

	Sample::operator json()  noexcept
	{
		return nullptr;
	}

	Sample &Sample::operator=(const json &encoded)
	{
		return *this;
	}

	Parameters::operator json() const noexcept
	{
		return nullptr;
	}

	Parameters &Parameters::operator=(const json &encoded)
	{
		return *this;
	}

    Candidate::operator json()
    {
        return nullptr;
    }

    Candidate &Candidate::operator=(const json &encoded)
    {
        return *this;
    }

    _Request::operator json() const noexcept
    {
        return nullptr;
    }

    _Request &_Request::operator=(const json &encoded)
    {
        return *this;
    }

    _Result::operator json()
    {
        return nullptr;
    }

    _Result &_Result::operator=(const json &encoded)
    {
        return *this;
    }
}
