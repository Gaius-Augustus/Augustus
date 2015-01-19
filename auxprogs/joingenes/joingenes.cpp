#include <iostream>		// input/output
#include <string>		// strings
#include <fstream>		// input/output on hard drive
#include <list>			// list
#include <unordered_map>		// hash
#include <unistd.h>		// for getopt()
#include <getopt.h>		// for getopt() with long option

#include "jg_transcript.h"
#include "jg_ios.h"

using namespace std;

static const struct option longOpts[] = {
    { "genesets", required_argument, NULL, 'g' },
    { "priorities", required_argument, NULL, 'p' },
    { "output", required_argument, NULL, 'o' },
    { "errordistance", required_argument, NULL, 'e' },
    { "genemodel", required_argument, NULL, 'm'},
    { "onlycompare", required_argument, NULL, 'c'},
    { "suppress", required_argument, NULL, 's' },
    { "join", required_argument, NULL, 'j' },
    { "help", no_argument, NULL, 'h' },
    { NULL, no_argument, NULL, 0 }
};

void display_help(void)
{
    cout << "This is the help documentation of the joingenes program." << endl;
    cout << "Options:" << endl;
    cout << "\tNecessary:" << endl;
    cout << "\t\t--genesets=file1,file2,...\t-g file1,file2,...\twhere \"file1,file2,...,filen\" have to be datafiles with genesets in gtf format." << endl;
    cout << "\t\t--output=ofile\t\t\t-o ofile\t\twhere \"ofile\" have to be a name for an output datafile (gtf)." << endl;
    cout << "\tOptional:" << endl;
    cout << "\t\t--priorities=pr1,pr2,...\t-p pr1,pr2,...\t\twhere \"pr1,pr2,...,prn\" have to be positiv integers (different from 0)." << endl;
    cout << "\t\t\t\t\t\t\t\t\tHave to be as many as filenames are added." << endl;
    cout << "\t\t\t\t\t\t\t\t\tBigger numbers means a higher priority." << endl;
    cout << "\t\t\t\t\t\t\t\t\tIf no priorities are added, the program will set all priorties to 1." << endl;
    cout << "\t\t\t\t\t\t\t\t\tThis option is only usefull, if there are more then one geneset." << endl;
    cout << "\t\t--errordistance=x\t\t-e x\t\t\twhere \"x\" have to be non-negativ integer" << endl;
    cout << "\t\t\t\t\t\t\t\t\tIf a prediction is <=x bases nexto a prediction range border, the program suppose, that there could be a mistake." << endl;
    cout << "\t\t\t\t\t\t\t\t\tDefault is, that it dont cares about prediction ranges." << endl;
    cout << "\t\t--genemodel=x\t\t\t-m x\t\t\twhere \"x\" have to be a genemodel out of \"eukaryote\"." << endl;
    cout << "\t\t\t\t\t\t\t\t\tDefault is eukaryotic." << endl;
    cout << "\t\t--onlycompare=x\t\t\t-c x\t\t\twhere \"x\" can be set to \"true\"." << endl;
    cout << "\t\t\t\t\t\t\t\t\tDefault is false." << endl;
    cout << "\t\t\t\t\t\t\t\t\tIf this value is set to true, it disables the normal function of the program." << endl;
    cout << "\t\t\t\t\t\t\t\t\tactivates a compare and seperate mode to seperate equal transcripts from non equal ones." << endl;
    cout << "\t\t--suppress=pr1,pr2,..\t\t-s pr1,pr2,...\t\twhere \"pr1,pr2,...,prm\" have to be positiv integers (different from 0)." << endl;
    cout << "\t\t\t\t\t\t\t\t\tDefault is none." << endl;
    cout << "\t\t\t\t\t\t\t\t\tif the core of a joined/non-joined transcript has one of these priorities it will not occur in the output file." << endl;
    cout << "\t\t--join=x\t\t\t-j x\t\t\twhere \"x\" is a boolean value (default is true)" << endl;
    cout << "\t\t\t\t\t\t\t\t\tIf this parameter is set to false, the program will not join/merge/shuffle;" << endl;
    cout << "\t\t\t\t\t\t\t\t\tit will only decide between the unchanged input transcipts and output them." << endl;
    cout << "\t\t--selecting=x\t\t\t-l x\t\t\twhere \"x\" is a boolean value (default is true)" << endl;
    cout << "\t\t\t\t\t\t\t\t\tIf this parameter is set to false, the program will not select at the end between \"contradictory\" transcripts." << endl;
    cout << "\t\t\t\t\t\t\t\t\t\"contradictory\" is self defined with respect to known biological terms." << endl;
    cout << "\t\t\t\t\t\t\t\t\tThe selection works with a self defined scoring function." << endl;
    cout << "\t\t--help \t\t\t\t-h\t\t\tprints the help documentation." << endl;
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
	filenames.push_back(temp_inside);
	for (unsigned int j = 1; j <= n_comma; j++){
	    temp_inside = strtok(NULL, ",");
	    filenames.push_back(temp_inside);
	}
    }else{
	filenames.push_back(filename_string);
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
	priorities.push_back(atoi(temp_inside));
	for (unsigned int j = 1; j <= n_comma; j++){
	    temp_inside = strtok(NULL, ",");
	    priorities.push_back(atoi(temp_inside));
	}
    }else{
	priorities.push_back(atoi(priority_string));
    }
    return priorities;
}

list<int> getSuppressionList(char* suppString){
    list<int> supprList;
    char *temp_inside;
    string str_temp = suppString;
    unsigned int n_comma = count(str_temp.begin(), str_temp.end(), ',');
    if (n_comma){
	temp_inside = strtok(suppString, ",");
	supprList.push_back(atoi(temp_inside));
	for (unsigned int j = 1; j <= n_comma; j++){
	    temp_inside = strtok(NULL, ",");
	    supprList.push_back(atoi(temp_inside));
	}
    }else{
	supprList.push_back(atoi(suppString));
    }
    return supprList;
}

int main(int argc, char* argv[])
{
    int opt = 0;
    static const char *optString = "g:p:o:e:m:c:s:j:l:h?";
    list<string> filenames;
    list<int> priorities;
    list<int> supprList;
    string filename_out;
    int longIndex;
    Properties properties;
    properties.errordistance = -1;
    properties.genemodel = "eukaryote";
    properties.onlyCompare = false;
    properties.join = true;
    properties.selecting = true;

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
	    properties.outFileName = optarg;
	    break;

	case 'e':
	    if (strstr(optarg, "rror") !=NULL){display_error("You have to write two \"-\" before the optionname \"errordistance\".");}
	    properties.errordistance = atoi(optarg);
	    break;

	case 'm':
	    if (strstr(optarg, "bacterium") == NULL){display_error("You used a wrong parameter for genemodel (allowed are \"bacterium\").");}
	    properties.genemodel = optarg;
	    break;

	case 'c':
	    if (strstr(optarg, "true") != NULL || strstr(optarg, "1") != NULL)
		properties.onlyCompare = true;
	    else
		cerr << "Unknown option for \"--onlycompare\", so this parameter is set to default. (possible is \"true\" or \"1\")" << endl;
	    break;

	case 's':
	    supprList = getSuppressionList(optarg);
	    break;

	case 'j':
	    if (strstr(optarg, "false") != NULL || strstr(optarg, "0") != NULL)
		properties.join = false;
	    else
		cerr << "Unknown value for \"--join\", so this parameter is set to default. (possible is \"false\" or \"0\")" << endl;
	    break;

	case 'l':
	    if (strstr(optarg, "false") != NULL || strstr(optarg, "0") != NULL)
		properties.join = false;
	    else
		cerr << "Unknown value for \"--selecting\", so this parameter is set to default. (possible is \"false\" or \"0\")" << endl;
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
    if (properties.outFileName == ""){
	display_error("Missing output filename.");
    }
    if (priorities.size() == 0){
	if (filenames.size()>1){
	    cout << "All priorities set to 1 (default)." << endl;
	}
	for (list<string>::iterator it_f = filenames.begin(); it_f != filenames.end(); it_f++){
	    priorities.push_back(1);
	}
    }else{
	for (list<int>::iterator it_p = priorities.begin(); it_p != priorities.end(); it_p++){
	    if ((*it_p) == 0){
		display_error("It is forbidden that a priority is equal to 0.");
	    }
	}
    }

    properties.filenames = filenames;
    properties.priorities = priorities;
    properties.supprList = supprList;

    if (properties.onlyCompare == true && (properties.priorities.size() != 2 || properties.priorities.begin() == properties.priorities.end())){
	display_error("CompareMode is only working with exact 2 priorities, which have to be different.");
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
    unordered_map<string,list<Transcript>> splitted_transcript_list;
    for (list<Transcript>::iterator it = transcript_list.begin(); it != transcript_list.end(); it++){
	splitted_transcript_list[(*it).getChr()].push_back(*it);
    }

    fstream outfile;
    outfile.open(properties.outFileName, ios::out);		// delete content of file filename
    outfile.close();

    if (properties.onlyCompare){
	string filenameEqual1 = "equal1.gtf", filenameEqual2 = "equal2.gtf", filenameStop1 = "same_stop1.gtf", filenameStop2 = "same_stop2.gtf", filenameUne1 = "unequal1.gtf", filenameUne2 = "unequal2.gtf";
	fstream outfile;
	outfile.open(filenameEqual1, ios::out);
	outfile.close();
	outfile.open(filenameEqual2, ios::out);
	outfile.close();
	outfile.open(filenameStop1, ios::out);
	outfile.close();
	outfile.open(filenameStop2, ios::out);
	outfile.close();
	outfile.open(filenameUne1, ios::out);
	outfile.close();
	outfile.open(filenameUne2, ios::out);
	outfile.close();
    }
    for(auto it = splitted_transcript_list.begin(); it != splitted_transcript_list.end(); it++){
        divideInOverlapsAndConquer(it->second, properties);
    }
    return 0;
}
