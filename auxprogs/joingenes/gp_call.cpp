#include <iostream>		// input/output
#include <string>		// strings
#include <fstream>		// input/output on hard drive
#include <list>			// list
#include <unordered_map>		// hash
#include <unistd.h>		// for getopt()
#include <getopt.h>		// for getopt() with long option

#include "joingenes.h"
#include "gp_ios.h"

using namespace std;

static const struct option longOpts[] = {
    { "genesets", required_argument, NULL, 'g' },
    { "priorities", required_argument, NULL, 'p' },
    { "output", required_argument, NULL, 'o' },
    { "help", no_argument, NULL, 'h' },
    { NULL, no_argument, NULL, 0 }
};

void display_help(void)
{
	cout << "This is the help documentation of the joingenes programm" << endl;
	cout << "Options:" << endl;
	cout << "\tRequired parameters:" << endl;
	cout << "\t\t--genesets=file1,file2,...\t-g file1,file2,...\twhere \"file1,file2,...,fileN\" is a list of GFF formatted files with gene sets" << endl;
	cout << "\t\t--output=ofile\t\t\t-o ofile\t\twhere \"ofile\" is a name for an output file in GTF format" << endl;
	cout << "\tOptional parameters:" << endl;
	cout << "\t\t--priorities=pr1,pr2,...\t-p pr1,pr2,...\t\twhere \"pr1,pr2,...,prn\" have to be integers different from 0" << endl;
	cout << "\t\t\t\t\t\t\t\t\tThere must be as many priorities as input files." << endl;
	cout << "\t\t\t\t\t\t\t\t\tLarger numbers means a higher priority." << endl;
	cout << "\t\t\t\t\t\t\t\t\tIf no priorities are added, all priorties are set to 1" << endl;
	cout << "\t\t\t\t\t\t\t\t\tThis option is only useful, if there is more then one geneset." << endl;
	cout << "\t\t--help \t\t\t\t-h\t\t\tprints the help documentation" << endl;
	cout << endl;
	exit( EXIT_FAILURE );
}

void display_error(string const &error)
{
	cerr << "Usage error: " << error << endl;
	cerr << "Use option --help or -h for usage explanations." << endl;
	exit( EXIT_FAILURE );
}

list<string> get_filenames(char* filename_string){
	list<string> filenames;
	char *temp_inside;
	string str_temp = filename_string;
	unsigned int n_comma = count(str_temp.begin(), str_temp.end(), ',');
	if (n_comma){
		temp_inside = strtok(filename_string, ",");
		filenames.push_front(temp_inside);
		for (unsigned int j = 1; j <= n_comma; j++){
			temp_inside = strtok(NULL, ",");
			filenames.push_front(temp_inside);
		}
	}else{
		filenames.push_front(filename_string);
	}
	return filenames;
}

list<int> get_priorities(char* priority_string){
	list<int> priorities;
	char *temp_inside;
	string str_temp = priority_string;
	unsigned int n_comma = count(str_temp.begin(), str_temp.end(), ',');
	if (n_comma){
		temp_inside = strtok(priority_string, ",");
		priorities.push_front(atoi(temp_inside));
		for (unsigned int j = 1; j <= n_comma; j++){
			temp_inside = strtok(NULL, ",");
			priorities.push_front(atoi(temp_inside));
		}
	}else{
		priorities.push_front(atoi(priority_string));
	}
	return priorities;
}

int main(int argc, char* argv[])
{
int opt = 0;
static const char *optString = "g:p:o:h?";
	list<string> filenames;
	list<int> priorities;
	string filename_out;
	int longIndex;

	opt = getopt_long(argc, argv, optString, longOpts, &longIndex);
	while (opt != -1) {
		switch (opt) {
			case 'g':
				filenames = get_filenames(optarg);
				break;

            case 'p':
				priorities = get_priorities(optarg);
				break;

			case 'o':
				filename_out = optarg;
				break;
			case 'h':
				display_help();
				break;
			case '0':				// long option without short form
				break;

			case '?':
				display_error("You entered an invalid option.");
				break;

			default:
				break;
        }
        opt = getopt_long(argc, argv, optString, longOpts, &longIndex);
    }
	if (filenames.size() == 0){
		display_error("Missing input filenames.");
	}
	if (filename_out == ""){
		display_error("Missing output filename.");
	}
	if (priorities.size() == 0){
		if (filenames.size()>1){
			cout << "All priorities set to 1 (default)." << endl;
		}
		for (list<string>::iterator it_f = filenames.begin(); it_f != filenames.end(); it_f++){
			priorities.push_front(1);
		}
	}else{
		for (list<int>::iterator it_p = priorities.begin(); it_p != priorities.end(); it_p++){
			if ((*it_p) == 0){
				display_error("It is forbidden that a priority is equal to 0.");
			}
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
		for (list<Transcript>::iterator it = transcript_list.begin(); it != transcript_list.end(); it++){
			(*it).initiate();
		} 
	}else{
		display_error("Number of input files and priorities is not equal.");
	}
	for (list<Transcript>::iterator it = transcript_list.begin(); it != transcript_list.end(); it++){
		if (!check_frame_annotation(*it)){
			cerr << (*it).t_id << " has wrong frame annotation." << endl;
			//display_error("Frames are not correct.");
		}
	}
	seek_overlaps(transcript_list, filename_out);
	return 0;
}
