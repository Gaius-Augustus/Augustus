#ifdef EBONY

#include <iostream>
#include <map>
#include <thread>
#include <chrono>
#include "inferclient.hh"
#include "properties.hh"

using namespace std;


Connection::Connection() : connect_timeout_ms(10000), response_timeout_ms(10000), max_tries(5) {
    cout << "Creating cURL connection" << endl;
    curl = curl_easy_init();
    if(!curl) {
            throw runtime_error("Failed to create CURL connection.");
        }
}

Connection::Connection(string server_url, long connect_timeout_ms, long response_timeout_ms, int max_tries) : server_url(server_url), connect_timeout_ms(connect_timeout_ms), response_timeout_ms(response_timeout_ms), max_tries(max_tries) {
    cout << "Creating cURL connection" << endl;
    curl = curl_easy_init();
    if(!curl) {
        throw runtime_error("Failed to create CURL connection.");
    }
}

Connection::~Connection() {
    if (curl) {
            curl_easy_cleanup(curl);
        }
}

size_t Connection::writeCallback(void *contents, size_t size, size_t nmemb, void *userp) {
    size_t total_size = size * nmemb;
    ((string*)userp)->append((char*)contents, total_size);
    return total_size;
}

json Connection::sendRequest(string& request_data) {

    cout << "Querying inference server (ebony)" << endl;
    // set up
    curl_easy_setopt(curl, CURLOPT_URL, server_url.c_str());
    curl_easy_setopt(curl, CURLOPT_CONNECTTIMEOUT_MS, connect_timeout_ms);
    curl_easy_setopt(curl, CURLOPT_TIMEOUT_MS, response_timeout_ms);

    // set up data
    curl_easy_setopt(curl, CURLOPT_POSTFIELDS, request_data.c_str());
    curl_easy_setopt(curl, CURLOPT_POSTFIELDSIZE, request_data.length());

    struct curl_slist* headers = nullptr;
    headers = curl_slist_append(headers, "Content-Type: application/json");
    curl_easy_setopt(curl, CURLOPT_HTTPHEADER, headers);

    // send request
    int retries = 0;
    string response_data;
    json parsed_json;

    while(retries < max_tries){

        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, writeCallback);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response_data);

        CURLcode res = curl_easy_perform(curl);
        if (res == CURLE_OK) {

            try {  // convert response to JSON format
                parsed_json = json::parse(response_data);
            } catch (const exception &e) {
               cerr << "ERROR: Invalid response format from server. Could not convert to JSON: " << e.what() << endl;
               return json{};
            }

            if (parsed_json.contains("preds")) {  // if response includes predictions
                curl_slist_free_all(headers);
                return parsed_json;

            } else { // retry
                cerr << "CURL request attempt " << (retries + 1) << " failed: server response did not include any predictions." << endl;
                ++retries;

                if (retries < max_tries){
                    cerr << "Retrying..." << endl;
                    // exponential backoff: sleep and increase waiting times for curl
                    this_thread::sleep_for(chrono::milliseconds(5000 * (1 << retries)));
                    curl_easy_setopt(curl, CURLOPT_CONNECTTIMEOUT_MS, connect_timeout_ms * (1 << (retries + 1)));
                    curl_easy_setopt(curl, CURLOPT_TIMEOUT_MS, response_timeout_ms * (1 << (retries + 1)));
                }
            }

        } else if (res == CURLE_OPERATION_TIMEDOUT || res == CURLE_COULDNT_CONNECT) {  // timed out, retry
            cerr << "CURL request attempt " << (retries + 1) << " failed: " << curl_easy_strerror(res) << endl;
            ++retries;

            if (retries < max_tries){
                cerr << "Retrying..." << endl;
                // exponential backoff: sleep and increase waiting times for curl
                this_thread::sleep_for(chrono::milliseconds(5000 * (1 << retries)));
                curl_easy_setopt(curl, CURLOPT_CONNECTTIMEOUT_MS, connect_timeout_ms * (1 << (retries + 1)));
                curl_easy_setopt(curl, CURLOPT_TIMEOUT_MS, response_timeout_ms * (1 << (retries + 1)));
            }

        } else {
            cout << "Query failed." << endl;
            cerr << "ERROR: CURL request failed. Non-retriable error occurred: " << curl_easy_strerror(res) << endl;
            curl_slist_free_all(headers);
            return json{};
        }
    }
    cout << "Query did not yield ebony predictions." << endl;
    return json{};  // all request retries failed
}


ConnectionHandler::ConnectionHandler(): connect_timeout_ms(1000), response_timeout_ms(10000), max_tries(5) {
    curl_global_init(CURL_GLOBAL_DEFAULT);

    readConfigFile();

    request_frame["tuples_overlap"] = false;
    request_frame["trans_dict"] = json::array();
    request_frame["input"] = json::array({""});
}

ConnectionHandler::~ConnectionHandler() {
    curl_global_cleanup();
}

Connection* ConnectionHandler::createConnection() {
    return new Connection(server_url, connect_timeout_ms, response_timeout_ms, max_tries);
}


void ConnectionHandler::readConfigFile() {
    string filename;
    filename = Properties::getProperty( INFERCLI_KEY ); // = "inferCliCfgFile"

    // open file
    ifstream file(filename);
    if(!file) {
        string altfname = Properties::getConfigFilename("") + "infer_cli/" + filename;
        file.clear();
        file.close();
        file.open(altfname);
        if (!file) {
            throw ios_base::failure("ConnectionHandler::readConfigFile(): Could neither find inference client config file " + filename
				   + " nor " + altfname +".");
        } else
            filename = altfname;
    }

    // read contents
    map<string,string> params;
    string line;
    while (getline(file, line)) {
        if (line.empty() || line[0] == '#') {  // skip empty lines and comments
            continue;
        }
        istringstream iss(line);
        string pname, value;
        if (iss >> pname >> value) {
            params[pname] = value;
        }
    }

    // save client config

    // mandatory parameters
    if (params.find("server_url") != params.end())
         server_url = params["server_url"];
    else
        throw invalid_argument("ConnectionHandler::readConfigFile(): Could not find parameter 'server_url' in inference client config file " + filename + ".");

    if (params.find("model") != params.end())
        request_frame["model"] = params["model"];
    else
        throw invalid_argument("ConnectionHandler::readConfigFile(): Could not find parameter 'model' in inference client config file " + filename + ".");

    if (params.find("flanking") != params.end()) {
        try {
            size_t pos;
            flanking = stoul(params["flanking"], &pos);
            if (params["flanking"].size() != pos)
                throw invalid_argument("String contains non-integer characters");
        } catch (const exception &e) {
            throw invalid_argument("ConnectionHandler::readConfigFile(): Could not convert parameter 'flanking' of inference client config file " + filename + " to int (" + e.what() + "). ");
        }
    } else
         throw invalid_argument("ConnectionHandler::readConfigFile(): Could not find parameter 'flanking' in inference client config file " + filename + ".");

    // optional parameters
    if (params.find("connect_timeout_ms") != params.end()) {
        try {
            size_t pos;
            unsigned long c = stoul(params["connect_timeout_ms"], &pos);
            if (params["connect_timeout_ms"].size() ==  pos)
                connect_timeout_ms = c;
        } catch (...) {}
    }

    if (params.find("response_timeout_ms") != params.end()) {
         try {
             size_t pos;
             unsigned long c = stoul(params["response_timeout_ms"], &pos);
             if (params["response_timeout_ms"].size() ==  pos)
                 response_timeout_ms = c;
         } catch (...) {}
    }

    if (params.find("max_tries") != params.end()) {
        try {
            size_t pos;
            unsigned long c = stoul(params["max_tries"], &pos);
            if (params["max_tries"].size() == pos)
                max_tries = c;
        } catch (...) {}
    }

}


vector<string> splitMSAData(string &msa_data, size_t chunkSize){
    vector<string> chunks;
    stringstream ss(msa_data);
    string line;
    string currentChunk;
    size_t currentSize;

    if (msa_data.size() < chunkSize){ // if data small enough, return without splitting
         chunks.push_back(msa_data);
    } else {
        // split into chunkSize big chunks
        while (getline(ss, line)) {
            line += "\n";  // add stripped newline
            currentChunk += line;
            currentSize += line.size();

            if (line == "\n" && currentSize >= chunkSize){  // only split at empty line
                chunks.push_back(currentChunk);
                currentChunk.clear();
                currentSize =  0;
            }
        }
        // add last chunk
        if (!currentChunk.empty()){
            chunks.push_back(currentChunk);
        }
    }
    return chunks;
}

void parseResponse(json &response, vector<bit_vector> &ssbound, list<OrthoExon> &hects){

    if (response.empty())
        return;

    // convert json array to vectors with preds
    vector<vector<double>> preds;
    try {
        preds = response["preds"].get<vector<vector<double>>>();
    } catch(const exception &e) {
        cout << "WARNING: Could not read predictions." << endl;
        cerr << "ERROR: Could not read ebony predictions. " << e.what() << endl;
        cerr << "Received server response: " << response << endl;
        return;
    }

    // assign ebony scores to ortho exon boundaries

    if (hects.size() != ssbound.size()){
        throw length_error("Size mismatch: hects and ssbound must have the same length.");
    }

    // iterator for OrthoExons, index for ebony scores
    auto hectIt = hects.begin();
    size_t predIdx = 0;

    // loop over ssbound and corresponding OrthoExons
    for (size_t i = 0; i < ssbound.size(); ++i, ++hectIt) {
        const bit_vector& bounds = ssbound[i];

        if (bounds[0] && bounds[1]) { // set both boundaries
            double left = preds[predIdx++][1];
            double right = preds[predIdx++][1];
            hectIt->setEbony(left, right);
        }  else if (bounds[0])  // set left boundary
            hectIt->setEbony(preds[predIdx++][1], -1.0);
        else if (bounds[1])  // set right boundary
            hectIt->setEbony(-1.0, preds[predIdx++][1]);
    }
    cout << "Parsed server response " << endl;
}


#endif
