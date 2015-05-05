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
    { "selection", required_argument, NULL, 'l' },
    { "mergeTaxa", required_argument, NULL, 't' },
    { "alternatives", required_argument, NULL, 'a' },
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
    cout << "\t\t\t\t\t\t\t\t\tDefault is 1000." << endl;
    cout << "\t\t\t\t\t\t\t\t\tTo disable the function, set errordistance to a negative number (e.g. -1)." << endl;
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
//    cout << "\t\t--mergetaxa=x\t\t\t-t x\t\t\twhere \"x\" is a boolean value (default is false)" << endl;
//    cout << "\t\t\t\t\t\t\t\t\tIf this parameter is set to true, the program merges the genes and transcripts with same taxa from different input data." << endl;
    cout << "\t\t--alternatives=x\t\t-a\t\t\tis a flag" << endl;
    cout << "\t\t\t\t\t\t\t\t\tIf this flag is set, the program joines different genes, if the transcripts of the genes are alternative variants." << endl;
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

void check_cds_stop_combination(Transcript* tx, pair<string,string> &stopVariants){
    if (tx->strand == '+'){
	if ((tx->tl_complete.second && tx->stop_list.empty()) || (!tx->tl_complete.second && !tx->stop_list.empty())){display_error("Semantic Error: stop-complete but no stop_codon.");}
	if (!tx->stop_list.empty()){
	    if (tx->tes == tx->stop_list.back().to){
		stopVariants.first = tx->t_id;
	    }else{
		stopVariants.second = tx->t_id;
	    }
	}
    }else{
	if (!tx->stop_list.empty()){
	    if (tx->tis == tx->stop_list.front().from){
		stopVariants.first = tx->t_id;
	    }else{
		stopVariants.second = tx->t_id;
	    }
	}
    }

    if (!stopVariants.first.empty() && !stopVariants.second.empty()){
	display_error("E.g. the stop_codon of "+stopVariants.first+" is inside of the CDS and the stop_codon of "+stopVariants.second+" is outside of the CDS. You have to make it equally!");
    }
}

int main(int argc, char* argv[])
{
    int opt = 0;
    static const char *optString = "g:p:o:e:m:c:s:j:l:t:ah?";
    list<string> filenames;
    list<int> priorities;
    list<int> supprList;
    string filename_out;
    int longIndex;
    Properties properties;
    properties.errordistance = 1000;
    properties.genemodel = "eukaryote";
    properties.onlyCompare = false;
    properties.join = true;
    properties.selecting = true;
    properties.mergeTaxa = false;
    properties.alternatives = false;
    properties.nrOfPrintedGenes = 0;
    properties.unknownCount = 1;
    properties.minimumIntronLength = 20; // hard coded (not optional at the moment)

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
		properties.selecting = false;
	    else
		cerr << "Unknown value for \"--selecting\", so this parameter is set to default. (possible is \"false\" or \"0\")" << endl;
	    break;

	case 't':
	    if (strstr(optarg, "true") != NULL || strstr(optarg, "1") != NULL)
		properties.mergeTaxa = true;
	    else
		cerr << "Unknown value for \"--mergetaxa\", so this parameter is set to default. (possible is \"true\" or \"1\")" << endl;
	    break;

	case 'a':
	    properties.alternatives = true;
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
    properties.minimumIntronLength = 20;

    if (properties.onlyCompare == true && (properties.priorities.size() != 2 || properties.priorities.begin() == properties.priorities.end())){
	display_error("CompareMode is only working with exact 2 priorities, which have to be different.");
    }

    unordered_map<string,Gene*> geneMap;
    properties.geneMap = &geneMap;
    unordered_map<string,Transcript*> transcriptMap;
    properties.transcriptMap = &transcriptMap;

//    list<Transcript> transcript_list;
//    properties.txList = &transcript_list;

    if (priorities.size() == filenames.size()){
	unordered_map<string,bool> taxaMap;
	list<int>::iterator it_p = priorities.begin();
	for (list<string>::iterator it_f = filenames.begin(); it_f != filenames.end(); it_f++){
	    load(geneMap, (*it_f) ,(*it_p), taxaMap, properties);
	    cout << "After loading " << (*it_f) << " (Priority " << (*it_p) << ") there are " << (*properties.transcriptMap).size() << " transcripts in transcript list." << endl;
	    it_p++;
	}


//(*properties.transcriptMap)

//(*(pointer->second))

	for (auto pointer = (*properties.transcriptMap).begin(); pointer != (*properties.transcriptMap).end(); pointer++){
	//for (list<Transcript>::iterator it = transcript_list.begin(); it != transcript_list.end(); it++){
	    if ((*(pointer->second)).exon_list.size() == 0){
		deleteTx(pointer->second, properties);
	    }else{
		(*(pointer->second)).initiate(properties);
	    }
	}
    }else{
        display_error("Number of input files and priorities is not equal.");
    }

// just for tests:
/*    for (auto pointer = (*properties.transcriptMap).begin(); pointer != (*properties.transcriptMap).end(); pointer++){
	if (pointer->second->parent->children.size() > 0){
	    cout << pointer->second->parent->children.size() << endl;
	    for (list<Transcript*>::iterator it_f = pointer->second->parent->children.begin(); it_f != pointer->second->parent->children.end(); it_f++){
		cout << (*it_f)->t_id << ": " << (*it_f)->exon_list.size() << " " << (*it_f)->tes << endl;
	    }
	}
    }cout << "URFURFURF" << endl; sleep(10);*/


    pair<string,string> stopVariants;		// <stop_codon is in CDS, stop_codon is not in CDS>

    for (auto pointer = (*properties.transcriptMap).begin(); pointer != (*properties.transcriptMap).end(); pointer++){
	check_cds_stop_combination(pointer->second, stopVariants);
	if (!check_frame_annotation(*pointer->second)){
	    pointer->second->isNotFrameCorrect = true;
	    cerr << pointer->second->t_id << " has wrong frame annotation and will be ignored in joining step." << endl;
	    //display_error("Frames are not correct.");
	}
    }
    unordered_map<string,list<Transcript*>> splitted_transcript_list;
    for (auto pointer = (*properties.transcriptMap).begin(); pointer != (*properties.transcriptMap).end(); pointer++){
	splitted_transcript_list[(*pointer->second).getChr()].push_back(pointer->second);
    }

    fstream outfile;
    outfile.open(properties.outFileName, ios::out);		// delete content of file filename
    outfile.close();


    if (properties.onlyCompare){
	string filenameEqual = "cEqual.gtf", filenameStop1 = "cSameStop1.gtf", filenameStop2 = "cSameStop2.gtf", filenameUne1 = "cUnequal1.gtf", filenameUne2 = "cUnequal2.gtf";
	fstream outfile;
	outfile.open(filenameEqual, ios::out);
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
