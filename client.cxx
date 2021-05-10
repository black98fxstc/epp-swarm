#include <client.h>
#include <string>
#include <cstring>
#include <sstream>

#include <aws/core/auth/AWSCredentialsProvider.h>
#include <aws/s3/model/HeadObjectRequest.h>
#include <aws/s3/model/GetObjectRequest.h>
#include <aws/s3/model/PutObjectRequest.h>

#include <credentials.h>

namespace EPP
{
    Aws::SDKOptions aws_options;

    void Init()
    {
        aws_options.loggingOptions.logLevel = Aws::Utils::Logging::LogLevel::Debug;
        Aws::InitAPI(aws_options);
    };

    void Finish()
    {
        Aws::ShutdownAPI(aws_options);
    };

    using json = nlohmann::json;

    struct ajax_action
    {
        std::string request;
        std::string response;
        size_t req_pos;
    };

    static size_t ajax_write(void *data, size_t size, size_t nmemb, void *userdata)
    {
        size_t realsize = size * nmemb;
        struct ajax_action *ajax = (struct ajax_action *)userdata;

        ajax->response.append((const char *)data, realsize);

        return realsize;
    }

    static size_t ajax_read(char *ptr, size_t size, size_t nmemb, void *userdata)
    {
        size_t realsize = size * nmemb;
        struct ajax_action *ajax = (struct ajax_action *)userdata;

        size_t remaining = ajax->request.size() - ajax->req_pos;
        if (realsize > remaining)
            realsize = remaining;
        memcpy(ptr, &(ajax->request[ajax->req_pos]), realsize);
        ajax->req_pos += realsize;

        return realsize;
    }

    Aws::S3::S3Client &Client::s3()
    {
        if (!s3_client)
        {
            Aws::Client::ClientConfiguration aws_config;
            aws_config.region = "us-west-2";
            s3_client = new Aws::S3::S3Client(aws_credentials, aws_config);
        }
        return *s3_client;
    };

    Client::Client()
    {
        curl_global_init(CURL_GLOBAL_ALL);
        curl = curl_easy_init();
        if (!curl)
        {
            curl_easy_cleanup(curl);
        }
        slist = curl_slist_append(slist, "Accept: application/json");
        slist = curl_slist_append(slist, "Content-Type: application/json");
    }

    Client::~Client()
    {
        curl_slist_free_all(slist);
        curl_easy_cleanup(curl);
        curl_global_cleanup();

        delete s3_client;
    }

    json Client::ajax(const std::string &endpoint, const json &request)
    {
        struct ajax_action ajax;
        ajax.request = request.dump();
        ajax.req_pos = 0;
        ajax.response = std::string("");

        CURLcode res;
        char error_message[CURL_ERROR_SIZE];
        long http_code;

        curl_easy_setopt(curl, CURLOPT_URL, endpoint.c_str());
        curl_easy_setopt(curl, CURLOPT_POST, 1L);
        curl_easy_setopt(curl, CURLOPT_HTTPHEADER, slist);
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, ajax_write);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, (void *)&ajax);
        curl_easy_setopt(curl, CURLOPT_READFUNCTION, ajax_read);
        curl_easy_setopt(curl, CURLOPT_READDATA, (void *)&ajax);
        curl_easy_setopt(curl, CURLOPT_ERRORBUFFER, error_message);

        res = curl_easy_perform(curl);
        if (res != CURLE_OK)
            fprintf(stderr, "curl_easy_perform() failed: %s\n", error_message);
        curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &http_code);
        if (http_code != 200)
            fprintf(stderr, "curl_easy_perform() HTTP status: %d\n", (int)http_code);

        json response = json::parse(ajax.response);
        return response;
    }

    bool Client::stage(Sample &sample)
    {
        Aws::S3::Model::HeadObjectRequest head_request;
        head_request.SetBucket("stanford-facs-epp-data");
        head_request.SetKey(sample.get_key().c_str());
        Aws::S3::Model::HeadObjectOutcome head_outcome =
            s3().HeadObject(head_request);
        if (head_outcome.IsSuccess())
            return false;

        std::shared_ptr<Aws::IOStream> input_data =
            Aws::MakeShared<SampleStream>("EPP", sample);
        Aws::S3::Model::PutObjectRequest put_request;
        put_request.SetBucket("stanford-facs-epp-data");
        put_request.SetKey(sample.get_key().c_str());
        put_request.SetContentLength(sample.measurments * sample.events * sizeof(epp_word));
        put_request.SetBody(input_data);
        Aws::S3::Model::PutObjectOutcome put_outcome =
            s3().PutObject(put_request);

        if (put_outcome.IsSuccess())
        {
            return true;
        }
        else
        {
            return false;
        }
    };

    void Client::fetch(Sample &sample)
    {
        SampleStream s(sample);
        char buf[40];
        for (int i = 0; i < 40; i++)
            s.put('x');
        s.write(buf, 40);
        s.write(buf, 40);
        s.write(buf, 40);
        s.write(buf, 40);
        s.write(buf, 40);
        s.write(buf, 40);
        s.write(buf, 40);
        s.write(buf, 40);
        s.write(buf, 40);
        s.flush();
        
        Aws::S3::Model::GetObjectRequest request;
        request.SetBucket("stanford-facs-epp-data");
        request.SetKey(sample.get_key().c_str());
        request.SetResponseStreamFactory(
            [&sample]() {
                return new SampleStream(sample);
            });

        Aws::S3::Model::GetObjectOutcome get_object_outcome = s3().GetObject(request);
        if (get_object_outcome.IsSuccess())
            return;
    };

    bool Client::stage(Subset &subset)
    {
        Aws::S3::Model::PutObjectRequest request;
        request.SetBucket("stanford-facs-epp-data");
        request.SetKey(subset.sample->get_key().c_str());

        std::shared_ptr<Aws::IOStream> input_data =
            Aws::MakeShared<SubsetStream>("SampleAllocationTag", subset);

        request.SetBody(input_data);

        Aws::S3::Model::PutObjectOutcome outcome =
            s3().PutObject(request);

        if (outcome.IsSuccess())
        {
            return true;
        }
        else
        {
            return false;
        }
    };

    bool Client::fetch(Subset &subset)
    {
        Aws::S3::Model::GetObjectRequest request;
        request.SetBucket("stanford-facs-epp-data");
        request.SetKey(subset.sample->get_key().c_str());
        request.SetResponseStreamFactory(
            [&subset]() {
                return new SubsetStream(subset);
            });

        Aws::S3::Model::GetObjectOutcome get_object_outcome = s3().GetObject(request);
        return true;
    };

}