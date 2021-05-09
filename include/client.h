#include <openssl/sha.h>
#include <curl/curl.h>
#include <nlohmann/json.hpp>
#include <aws/core/Aws.h>
#include <aws/core/utils/logging/LogLevel.h>
#include <aws/s3/S3Client.h>

using json = nlohmann::json;

namespace EPP
{
    typedef uint32_t epp_word;
    typedef unsigned char hash_t[SHA256_DIGEST_LENGTH];

    class Sample
    {
    public:
        const int measurments;
        const int events;
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
                      std::string key) : Sample(measurments, events, key){};

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

    class SampleStream : public std::iostream
    {
    protected:
        class sample_buffer : public std::streambuf
        {

        public:
            sample_buffer(Sample &sample);
            virtual ~sample_buffer();
            virtual std::streambuf::int_type underflow();
            virtual std::streambuf::int_type overflow(std::streambuf::int_type value);

        private:
            Sample *sample;
            epp_word *buffer;
            long next_event;
        };

    public:
        SampleStream(Sample &sample) : std::iostream(new sample_buffer(sample)){};
    };

    class Subset
    {
    };

    class Client
    {
    public:
        Client();
        ~Client();
        json ajax(const std::string &endpoint, const json &request);
        bool stage(Sample &sample);
        void fetch(Sample &sample);

    private:
        CURL *curl = NULL;
        struct curl_slist *slist = NULL;
        Aws::SDKOptions aws_options;
        Aws::String aws_region;
        Aws::Client::ClientConfiguration aws_config;
        Aws::S3::S3Client s3_client;
    };
}
