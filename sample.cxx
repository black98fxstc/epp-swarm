#include <client.h>
#include <type_traits>

const int QUANTUM = 1000;

namespace EPP
{
    SampleStream::sample_buffer::sample_buffer(Sample &sample)
        : sample(&sample)
    {
        buffer = new epp_word[sample.measurments * QUANTUM];
        next_event = 0;
        setg(0, 0, 0);
        setp((char *)buffer, (char *)(buffer + sample.measurments * QUANTUM));
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
            for (int measurment = 0; measurment < sample->measurments; measurment++)
                *ptr++ = sample->get_word(measurment, next_event);
        setg((char *)buffer, (char *)buffer, (char *)ptr);
        return traits_type::to_int_type(*gptr());
    }

    std::streambuf::int_type SampleStream::sample_buffer::overflow(std::streambuf::int_type value)
    {
        size_t write = pptr() - pbase();
        if (write)
        {
            long count = write / (sizeof(epp_word) * sample->measurments);
            epp_word *ptr = buffer;
            for (int i = 0; i < count; i++, next_event++)
                for (int measurment = 0; measurment < sample->measurments; measurment++)
                    sample->put_word(measurment, next_event, *ptr++);
        }
        setp((char *)buffer, (char *)(buffer + sample->measurments * QUANTUM));

        if (next_event == sample->events)
            return traits_type::eof();
        else
            return traits_type::not_eof(value);
    };

    std::string Sample::get_key()
    {
        if (!hash)
        {
            SHA256_CTX sha256;
            SHA256_Init(&sha256);

            for (long event = 0; event < events; event++)
                for (int measurment = 0; measurment < measurments; measurment++)
                {
                    epp_word word = get_word(measurment, event);
                    SHA256_Update(&sha256, &word, sizeof(word));
                }

            SHA256_Final(hash, &sha256);
        }

        if (key.length() == 0)
        {
            key = std::string();
            key.resize(2 * sizeof(hash_t));
            for (int i = 0, j = 0; i < sizeof(hash_t); i++)
            {
                int byte = hash[i];
                int nibble = byte & 0x0f;
                if (nibble < 10)
                    key[j++] = '0' + nibble;
                else
                    key[j++] = 'A' + nibble - 10;
                nibble = (byte >> 4) && 0x0f;
                if (nibble < 10)
                    key[j++] = '0' + nibble;
                else
                    key[j++] = 'A' + nibble - 10;
            }
        }
        return key;
    };

    Sample::Sample(const int measurments,
                   const long events)
        : measurments(measurments), events(events){};

    Sample::Sample(const int measurments,
                   const long events,
                   std::string key)
        : measurments(measurments), events(events), key(key){};
}
