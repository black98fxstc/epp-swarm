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

	Sample::operator json() const noexcept
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

    Candidate::operator json() const noexcept
    {
        return nullptr;
    }

    Candidate &Candidate::operator=(const json &encoded)
    {
        return *this;
    }

    Result::operator json() const noexcept
    {
        return nullptr;
    }

    Result &Result::operator=(const json &encoded)
    {
        return *this;
    }
}
