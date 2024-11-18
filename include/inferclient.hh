/*
 * inferclient.hh
 */


#ifdef EBONY

#ifndef INFERENCE_CLIENT_H
#define INFERENCE_CLIENT_H

#include <string>
#include <vector>
#include <curl/curl.h>
#include "json.hpp"
#include "orthoexon.hh"

using namespace std;


/**
 * @brief client for requesting clamsa predictions from inference server
 *
 * @author Hannah Thierfeldt
 */


vector<string> splitMSAData(string &msa_data, size_t chunkSize = 100 * 1024 * 1024); // default chunk size is 5 Mb

void parseResponse(json &response, vector<bit_vector> &ssbound, list<OrthoExon> &hects);



class Connection {
    public:
        Connection();
        Connection(string server_url, long connect_timeout_ms, long response_timeout_ms, int max_tries);
        ~Connection();

        // setter
        void setServerURL(string url) {server_url = url;}
        void setConnectTimeout(int timeout_ms) {connect_timeout_ms = timeout_ms;}
        void setResponseTimeout(int timeout_ms) {response_timeout_ms = timeout_ms;}
        void setMaxTries(int tries) {max_tries = tries;}

        json sendRequest(string &request_data);


    private:
	CURL* curl;
        string server_url;
        unsigned long connect_timeout_ms;
        unsigned long response_timeout_ms;
        unsigned long max_tries;  // max number of request tries

        static size_t writeCallback(void *contents, size_t size, size_t nmemb, void *userp);
};


class ConnectionHandler {
    public:
        ConnectionHandler();
        ~ConnectionHandler();

        void readConfigFile();

        void setDataHead(string data_head) {request_frame["input"] = json::array({data_head});}
        json getRequestFrame() const {return request_frame;}
        int getFlanking() const {return flanking;}

        Connection* createConnection();
        //Connection createConnection(string server_url, long connect_timeout_ms, long response_timeout_ms, int max_tries);

    private:
        string server_url;
        unsigned long connect_timeout_ms;
        unsigned long response_timeout_ms;
        unsigned long max_tries;
        json request_frame; // head of every request data
        unsigned long flanking;
};

#endif
#endif
