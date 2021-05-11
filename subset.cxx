#include <client.h>
#include <type_traits>

const int QUANTUM = 1000;

namespace EPP
{
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
        long count = subset->sample->events - next_event;
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
        size_t write = pptr() - pbase();
        if (write)
        {
            long count = subset->sample->events - next_event;
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

        if (next_event == subset->sample->events)
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

    std::string Subset::get_key()
    {
        hash_t hash;
        SHA256_CTX sha256;
        SHA256_Init(&sha256);

        long count = sample->events;
        long next_event = 0;
        while (count > 8)
        {
            uint8_t data = 0;
            for (int bit = 1; bit < 1 << 8; bit <<= 1)
                if (this->at(next_event++))
                    data |= bit;
            count -= 8;
            SHA256_Update(&sha256, &data, sizeof(data));
        }
        if (count > 0)
        {
            uint8_t data = 0;
            for (int bit = 1; bit < 1 << count; bit <<= 1)
                if (this->at(next_event++))
                    data |= bit;
            count = 0;
            SHA256_Update(&sha256, &data, sizeof(data));
        }
        SHA256_Final(hash, &sha256);

        std::string key = sample->get_key();
        key.resize(4 * sizeof(hash_t));
        for (int i = 0, j = 2 * sizeof(hash_t); i < sizeof(hash_t); i++)
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
        return key;
    };

    Subset::Subset(Sample &sample)
        : sample(&sample), std::vector<bool>(sample.events){};

    Subset::Subset(Sample &sample, std::string key)
        : sample(&sample), key(key), std::vector<bool>(sample.events){};
}