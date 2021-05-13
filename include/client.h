#include <ios>
#include <sstream>
#include <fstream>
#include <memory>
#include <curl/curl.h>
#include <nlohmann/json.hpp>
#include <aws/core/Aws.h>
#include <aws/core/utils/logging/LogLevel.h>
#include <aws/s3/S3Client.h>
#include <openssl/sha.h>
#include <fftw3.h>

using json = nlohmann::json;

namespace EPP
{
    typedef uint32_t epp_word;
    typedef unsigned char hash_t[SHA256_DIGEST_LENGTH];
    static fftwf_plan DCT;
    static fftwf_plan IDCT;

    // resolution of the density estimator, best if it's a power of 2 for FFTW
    const int N = 1 << 9;

    class Sample
    {
    public:
        const int measurments;
        const long events;
        std::string get_key();

    protected:
        Sample(const int measurments,
               const long events);
        Sample(const int measurments,
               const long events,
               std::string key);
        std::string key;
        hash_t hash;

    private:
        virtual epp_word get_word(int measurment, long event) = 0;
        virtual void put_word(int measurment, long event, epp_word data) = 0;
        friend class SampleStream;
    };

    template <typename _float>
    class DefaultSample : public Sample
    {
    public:
        DefaultSample(const int measurments,
                      const long events,
                      _float *data) : Sample(measurments, events), data(data){};
        DefaultSample(const int measurments,
                      const long events,
                      _float *data,
                      std::string key) : Sample(measurments, events, key), data(data){};

    protected:
        epp_word get_word(int measurment, long event)
        {
            float f = data[measurments * event + measurment];
            return *(epp_word *)&f;
        };

        void put_word(int measurment, long event, epp_word value)
        {
            float f = *(float *)&value;
            data[measurments * event + measurment] = (_float)f;
        };

    private:
        _float *data;
    };

    template <typename _float>
    class TransposeSample : public Sample
    {
    public:
        TransposeSample(const int measurments,
                        const long events,
                        _float *data) : Sample(measurments, events), data(data){};
        TransposeSample(const int measurments,
                        const long events,
                        _float *data,
                        std::string key) : Sample(measurments, events, key), data(data){};

    protected:
        epp_word get_word(int measurment, long event)
        {
            float f = data[events * measurment + event];
            return *(epp_word *)&f;
        };

        void put_word(int measurment, long event, epp_word value)
        {
            float f = *(float *)&value;
            data[events * measurment + event] = (_float)f;
        };

    private:
        _float *data;
    };

    template <typename _float>
    class PointerSample : public Sample
    {
    public:
        PointerSample(const int measurments,
                      const long events,
                      _float **data) : Sample(measurments, events), data(data){};
        PointerSample(const int measurments,
                      const long events,
                      _float *data,
                      std::string key) : Sample(measurments, events, key), data(data){};

    protected:
        epp_word get_word(int measurment, long event)
        {
            float f = data[measurment][event];
            return *(epp_word *)&f;
        };

        void put_word(int measurment, long event, epp_word value)
        {
            float f = *(float *)&value;
            data[measurment][event] = (_float)f;
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
            sample_buffer(Sample &sample);
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
        SampleStream(Sample &sample);
    };

    class Subset : public std::vector<bool>
    {
    public:
        Subset(Sample &sample);
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
            subset_buffer(Subset &subset);
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
        SubsetStream(Subset &subset);
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
        CURL *curl = NULL;
        struct curl_slist *slist = NULL;
        Aws::S3::S3Client *s3_client;
        Aws::S3::S3Client &s3();
    };

    void Init();
    void Finish();
}
