#include <type_traits>

#include "client.h"

const int QUANTUM = 1000;

namespace EPP
{
    SampleStream::sample_buffer::sample_buffer(Sample &sample)
        : std::streambuf(), sample(&sample)
    {
        buffer = new epp_word[sample.measurements * QUANTUM];
        next_event = 0;
        setg(0, 0, 0);
        setp((char *)buffer, (char *)(buffer + sample.measurements * QUANTUM));
    }

    SampleStream::sample_buffer::~sample_buffer()
    {
        delete[] buffer;
    }

    std::streambuf::int_type SampleStream::sample_buffer::underflow()
    {
        long count = sample->events - next_event;
        if (count > QUANTUM)
            count = QUANTUM;

        if (count == 0)
            return traits_type::eof();

        epp_word *ptr = buffer;
        for (int i = 0; i < count; i++, next_event++)
            for (int measurment = 0; measurment < sample->measurements; measurment++)
                *ptr++ = sample->get_word(measurment, next_event);
        setg((char *)buffer, (char *)buffer, (char *)ptr);
        return traits_type::to_int_type(*gptr());
    }

    std::streambuf::int_type SampleStream::sample_buffer::overflow(std::streambuf::int_type value)
    {
        size_t write = pptr() - pbase();
        if (write)
        {
            long count = write / (sizeof(epp_word) * sample->measurements);
            epp_word *ptr = buffer;
            for (int i = 0; i < count; i++, next_event++)
                for (int measurment = 0; measurment < sample->measurements; measurment++)
                    sample->put_word(measurment, next_event, *ptr++);
        }
        setp((char *)buffer, (char *)(buffer + sample->measurements * QUANTUM));

        if (next_event == sample->events)
            return traits_type::eof();
        else
            return traits_type::not_eof(value);
    };

    int SampleStream::sample_buffer::sync()
    {
        std::streambuf::int_type result = this->overflow(traits_type::eof());
        return traits_type::eq_int_type(result, traits_type::eof()) ? -1 : 0;
    }

    SampleStream::SampleStream(Sample &sample) : std::iostream(new sample_buffer(sample)){};

    SubsetStream::subset_buffer::subset_buffer(Subset &subset)
        : subset(&subset)
    {
        buffer = new uint8_t[QUANTUM];
        next_event = 0;
        setg(0, 0, 0);
        setp((char *)buffer, (char *)(buffer + QUANTUM));
    }

    SubsetStream::subset_buffer::~subset_buffer()
    {
        delete[] buffer;
    }

    std::streambuf::int_type SubsetStream::subset_buffer::underflow()
    {
        long count = subset->size() - next_event;
        if (count > QUANTUM * 8)
            count = QUANTUM * 8;

        if (count == 0)
            return traits_type::eof();

        uint8_t *ptr = buffer;
        while (count > 8)
        {
            uint8_t data = 0;
            for (int bit = 1; bit < 1 << 8; bit <<= 1)
                if (subset->at(next_event++))
                    data |= bit;
            *ptr++ = data;
            count -= 8;
        }
        if (count > 0)
        {
            uint8_t data = 0;
            for (int bit = 1; bit < 1 << count; bit <<= 1)
                if (subset->at(next_event++))
                    data |= bit;
            *ptr++ = data;
            count = 0;
        }

        setg((char *)buffer, (char *)buffer, (char *)ptr);
        return traits_type::to_int_type(*gptr());
    }

    std::streambuf::int_type SubsetStream::subset_buffer::overflow(std::streambuf::int_type value)
    {
        long int write = pptr() - pbase();
        if (write)
        {
            long int count = subset->size() - next_event;
            if (count > write * 8)
                count = write * 8;

            uint8_t *ptr = buffer;
            while (count > 8)
            {
                uint8_t data = *ptr++;
                for (int bit = 1; bit < 1 << 8; bit <<= 1)
                    subset->at(next_event++) = data & bit;
                count -= 8;
            }
            if (count > 0)
            {
                uint8_t data = *ptr++;
                for (int bit = 1; bit < 1 << count; bit <<= 1)
                    subset->at(next_event++) = data & bit;
                count = 0;
            }
        }
        setp((char *)buffer, (char *)(buffer + QUANTUM));

        if (next_event == subset->size())
            return traits_type::eof();
        else
            return traits_type::not_eof(value);
    };

    int SubsetStream::subset_buffer::sync()
    {
        std::streambuf::int_type result = this->overflow(traits_type::eof());
        return traits_type::eq_int_type(result, traits_type::eof()) ? -1 : 0;
    }

    SubsetStream::SubsetStream(Subset &subset) : std::iostream(new subset_buffer(subset)){};

    // Subset::Subset(unsigned long int events)
    //     : std::vector<bool>(sample.events){};
}
