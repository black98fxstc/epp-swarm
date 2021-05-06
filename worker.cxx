#include <cstring>
#include <string>
#include <iostream>
#include <curl/curl.h>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

struct ajax_action
{
    std::string request;
    std::string response;
    size_t req_pos;
};

static size_t write_callback(void *data, size_t size, size_t nmemb, void *userdata)
{
    size_t realsize = size * nmemb;
    struct ajax_action *ajax = (struct ajax_action *)userdata;

    ajax->response.append((const char *)data, realsize);
    
    return realsize;
}

static size_t read_callback(char *ptr, size_t size, size_t nmemb, void *userdata)
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

json epp_ajax(const std::string &endpoint, const json &request)
{
    struct ajax_action ajax;
    ajax.request = request.dump();
    ajax.req_pos = 0;
    ajax.response = std::string("");

    CURL *curl;
    CURLcode res;
    struct curl_slist *slist = NULL;

    curl_global_init(CURL_GLOBAL_ALL);
    curl = curl_easy_init();
    if (curl)
    {
        curl_easy_setopt(curl, CURLOPT_URL, endpoint.c_str());
        curl_easy_setopt(curl, CURLOPT_POST, 1L);
        slist = curl_slist_append(slist, "Accept: application/json");
        slist = curl_slist_append(slist, "Content-Type: application/json");
        curl_easy_setopt(curl, CURLOPT_HTTPHEADER, slist);
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_callback);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, (void *)&ajax);
        curl_easy_setopt(curl, CURLOPT_READFUNCTION, read_callback);
        curl_easy_setopt(curl, CURLOPT_READDATA, (void *)&ajax);

        res = curl_easy_perform(curl);
        if (res != CURLE_OK)
            fprintf(stderr, "curl_easy_perform() failed: %s\n",
                    curl_easy_strerror(res));

        curl_slist_free_all(slist);
        curl_easy_cleanup(curl);
    }
    curl_global_cleanup();

    json response = json::parse(ajax.response);
    return response;
}

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " endpoint\n";
        return 1;
    }

    json request;
    request["action"] = "Do something!";
    request["argument"] = "something to work on";

    json response = epp_ajax(argv[1], request);

    std::cout << request.dump(4) << std::endl;
    std::cout << response.dump(4) << std::endl;

    return 0;
}