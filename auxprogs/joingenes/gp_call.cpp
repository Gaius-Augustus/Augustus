#include <iostream>	// input/output
#include <string>	// strings
#include <vector>	// vectors
#include <fstream>	// input/output on hard drive
#include <stdlib.h>	// atof (change into double)
#include <stack>	// stack
#include <list>		// list
#include <unordered_map>		// hash
#include <sstream>	// string stream; to convert int into string

#include "joingenes.h"
#include "gp_ios.h"

using namespace std;

int main(int argc, char* argv[])
{
	if (argc <= 1) {
		cerr << "Usage: " << argv[0] << " no --genesets" << endl;
		return 1;
	}
	list<string> filenames;
	list<int> priorities;
	string filename_out;
	for (int i = 1; i < argc; ++i) {
		if (strstr(argv[i], "--genesets")!=NULL && strstr(argv[i], "=")!=NULL) {
			if (filenames.size() == 0){
				char *temp_inside;
				string str_temp = argv[i];
				unsigned int n_comma = count(str_temp.begin(), str_temp.end(), ',');
				strtok(argv[i], "=");
				for (unsigned int j = 0; j <= n_comma; j++){
					temp_inside = strtok(NULL, ",");
					filenames.push_front(temp_inside);
				}
			}else{
				cerr << "to much times --genesets option" << endl;
				return 1;
			}
		}else if (strstr(argv[i], "--priorities")!=NULL && strstr(argv[i], "=")!=NULL){
			if (priorities.size() == 0){
				char *temp_inside;
				string str_temp = argv[i];
				unsigned int n_comma = count(str_temp.begin(), str_temp.end(), ',');
				strtok(argv[i], "=");
				for (unsigned int j = 0; j <= n_comma; j++){
					temp_inside = strtok(NULL, ",");
					priorities.push_front(atoi(temp_inside));
				}


			}else{
				cerr << "to much times --priorities option" << endl;
			}
		}else if (strstr(argv[i], "--out")!=NULL && strstr(argv[i], "=")!=NULL){
			strtok(argv[i], "=");
			char *temp_inside;
			temp_inside = strtok(NULL, ",");
			filename_out = temp_inside;
		}else{
			cerr << argv[i] << " is an unknown option" << endl;
			return 1;
		}
	}
	if (filenames.size() == 0){
		cerr << "no filenames added by --genesets=x1,x2 option" << endl;
		return 1;
	}
	if (filename_out == ""){
		cerr << "no output filename added by --out=x option" << endl;
		return 1;
	}
	if (priorities.size() == 0){
		if (filenames.size()>1){
			cout << "all priorities set to 1 (default); set priorities by --priorities=x1,x2" << endl;
		}
		for (list<string>::iterator it_f = filenames.begin(); it_f != filenames.end(); it_f++){
			priorities.push_front(1);
		}
	}
	unordered_map<string,Gene> gene_map;
	list<Transcript> transcript_list;

	if (priorities.size() == filenames.size()){
		list<int>::iterator it_p = priorities.begin();
		for (list<string>::iterator it_f = filenames.begin(); it_f != filenames.end(); it_f++){
			load(gene_map, transcript_list, (*it_f) ,(*it_p));
			cout << "After loading " << (*it_f) << " (Priority " << (*it_p) << ") there are " << transcript_list.size() << " transcripts in transcript_list." << endl;
			it_p++;
		}
	}else{
		cerr << "number of input files and priorities dont match" << endl;
		return 1;
	}
	//check_frame_annotation(transcript_list);
	transcript_list.sort();
	cout << transcript_list.size() << endl;
	seek_overlaps(transcript_list);
	//save(gene_map, transcript_list, filename_output);
	cout << endl;
	return 0;
}



/*int main(int argc, char* argv[])
{
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << "--destination DESTINATION SOURCE" << std::endl;
        return 1;
    }
    std::vector <std::string> sources;
    std::string destination;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--destination") {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                destination = argv[i++]; // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                  std::cerr << "--destination option requires one argument." << std::endl;
                return 1;
            }  
        } else {
            sources.push_back(argv[i]);
        }
    }
    return move(sources, destination);
}*/
