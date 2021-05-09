#include <client.h>
#include <string>
#include <cstring>
#include <curl/curl.h>
#include <nlohmann/json.hpp>

#include <aws/s3/model/GetObjectRequest.h>
#include <aws/s3/model/PutObjectRequest.h>

namespace EPP
{
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

        aws_options.loggingOptions.logLevel = Aws::Utils::Logging::LogLevel::Debug;
        Aws::InitAPI(aws_options);
        aws_config.region = "us-east-1";
        s3_client = Aws::S3::S3Client(aws_config);
    }

    Client::~Client()
    {
        curl_slist_free_all(slist);
        curl_easy_cleanup(curl);
        curl_global_cleanup();

        ShutdownAPI(aws_options);
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
        Aws::S3::Model::PutObjectRequest request;
        request.SetBucket("stanford-facs-epp-data");
        request.SetKey(sample.get_key().c_str());

        std::shared_ptr<Aws::IOStream> input_data =
            Aws::MakeShared<SampleStream>("SampleAllocationTag", sample);

        request.SetBody(input_data);

        Aws::S3::Model::PutObjectOutcome outcome =
            s3_client.PutObject(request);

        if (outcome.IsSuccess())
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
        Aws::S3::Model::GetObjectRequest request;
        request.SetBucket("stanford-facs-epp-data");
        request.SetKey(sample.get_key().c_str());
        request.SetResponseStreamFactory(
            [&sample]() {
                return new SampleStream(sample);
            });

        Aws::S3::Model::GetObjectOutcome get_object_outcome = s3_client.GetObject(request);
    };

    bool Client::stage(Subset &subset)
    {
        Aws::S3::Model::PutObjectRequest request;
        request.SetBucket("stanford-facs-epp-data");
        // request.SetKey(subset.sample->get_key().c_str());

        std::shared_ptr<Aws::IOStream> input_data =
            Aws::MakeShared<SubsetStream>("SampleAllocationTag", subset);

        request.SetBody(input_data);

        Aws::S3::Model::PutObjectOutcome outcome =
            s3_client.PutObject(request);

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

        Aws::S3::Model::GetObjectOutcome get_object_outcome = s3_client.GetObject(request);
        return true;
    };

}