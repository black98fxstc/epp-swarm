#ifndef _EPP_CLIENT_H
#define _EPP_CLIENT_H  1

#include <ios>
#include <sstream>
#include <memory>

#include <aws/core/Aws.h>
#include <aws/core/utils/logging/LogLevel.h>
#include <aws/s3/S3Client.h>

#include <curl/curl.h>
#include <nlohmann/json.hpp>
#include <openssl/sha.h>

using json = nlohmann::json;

namespace EPP
{
    typedef uint32_t epp_word;
    typedef unsigned char hash_t[SHA256_DIGEST_LENGTH];

    // resolution of the density estimator
    // FFT is fastest when N has lots of small prime factors
    const int N = 1 << 6;

    class Sample
    {
    public:
        const int measurements;
        const long events;
        std::string get_key();

    protected:
        Sample(int measurements,
               long events);
        Sample(int measurements,
               long events,
               std::string key);
        std::string key;
        hash_t hash;

    private:
        virtual epp_word get_word(int measurement, long event) = 0;
        virtual void put_word(int measurement, long event, epp_word data) = 0;
        friend class SampleStream;
    };

    template <typename _float>
    class DefaultSample : public Sample
    {
    public:
        DefaultSample(const int measurements,
                      const long events,
                      _float *data) : Sample(measurements, events), data(data){};
        DefaultSample(const int measurements,
                      const long events,
                      _float *data,
                      std::string key) : Sample(measurements, events, key), data(data){};

    protected:
        epp_word get_word(int measurement, long event)
        {
            float f = data[measurements * event + measurement];
            return *(epp_word *)&f;
        };

        void put_word(int measurement, long event, epp_word value)
        {
            float f = *(float *)&value;
            data[measurements * event + measurement] = (_float)f;
        };

    private:
        _float *data;
    };

    template <typename _float>
    class TransposeSample : public Sample
    {
    public:
        TransposeSample(const int measurements,
                        const long events,
                        _float *data) : Sample(measurements, events), data(data){};
        TransposeSample(const int measurements,
                        const long events,
                        _float *data,
                        std::string key) : Sample(measurements, events, key), data(data){};

    protected:
        epp_word get_word(int measurement, long event)
        {
            float f = data[events * measurement + event];
            return *(epp_word *)&f;
        };

        void put_word(int measurement, long event, epp_word value)
        {
            float f = *(float *)&value;
            data[events * measurement + event] = (_float)f;
        };

    private:
        _float *data;
    };

    template <typename _float>
    class PointerSample : public Sample
    {
    public:
        PointerSample(const int measurements,
                      const long events,
                      _float **data) : Sample(measurements, events), data(data){};
        PointerSample(const int measurements,
                      const long events,
                      _float *data,
                      std::string key) : Sample(measurements, events, key), data(data){};

    protected:
        epp_word get_word(int measurement, long event)
        {
            float f = data[measurement][event];
            return *(epp_word *)&f;
        };

        void put_word(int measurement, long event, epp_word value)
        {
            float f = *(float *)&value;
            data[measurement][event] = (_float)f;
        };

    private:
        _float **data;
    };

    class SampleStream : public Aws::IOStream
    {
    protected:
        class sample_buffer : public std::streambuf
        {

        public:
            explicit sample_buffer(Sample &sample);
            virtual ~sample_buffer();

        protected:
            virtual std::streambuf::int_type underflow();
            virtual std::streambuf::int_type overflow(std::streambuf::int_type value);
            virtual std::streambuf::int_type sync();

        private:
            Sample *sample;
            epp_word *buffer;
            long next_event;
        };

    public:
        explicit SampleStream(Sample &sample);
    };

    class Subset : public std::vector<bool>
    {
    public:
        explicit Subset(Sample &sample);
        Subset(Sample &sample, std::string key);
        Sample *sample;
        std::string get_key();

    private:
        std::string key;
        friend class SubsetStream;
    };

    class SubsetStream : public std::iostream
    {
    protected:
        class subset_buffer : public std::streambuf
        {
        public:
            explicit subset_buffer(Subset &subset);
            virtual ~subset_buffer();
            virtual std::streambuf::int_type underflow();
            virtual std::streambuf::int_type overflow(std::streambuf::int_type value);
            virtual std::streambuf::int_type sync();

        private:
            Subset *subset;
            uint8_t *buffer;
            long next_event;
            friend class SubsetStream;
        };

    public:
        explicit SubsetStream(Subset &subset);
    };

    class Client
    {
    public:
        Client();
        ~Client();
        json ajax(const std::string &endpoint, const json &request);
        bool stage(Sample &sample);
        bool fetch(Sample &sample);
        bool stage(Subset &subset);
        bool fetch(Subset &subset);

    private:
        CURL *curl = nullptr;
        struct curl_slist *slist = nullptr;
        Aws::S3::S3Client *s3_client;
        Aws::S3::S3Client &s3();
    };

    void Init();
    void Finish();
}
#endif  /* client.h */