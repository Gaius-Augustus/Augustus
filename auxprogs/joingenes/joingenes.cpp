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
    { "onlycompare", no_argument, NULL, 'c'},
    { "suppress", required_argument, NULL, 's' },
    { "nojoin", no_argument, NULL, 'j' },
    { "noselection", no_argument, NULL, 'l' },
    { "alternatives", no_argument, NULL, 'a' },
    { "stopincoding", no_argument, NULL, 'i' },
    { "help", no_argument, NULL, 'h' },
    { NULL, no_argument, NULL, 0 }
};

void display_help(void)
{
    cout << "joingenes - merge several gene sets into one, including the combination of transcripts to new ones not present in input" << endl;
    cout << "This program works in several steps:" << endl;
    cout << "  1. divide the set of all transcripts into smaller sets, in which all transcripts are on the same sequence and are overlapping at least with one other transcript in this set (set is called \"overlap\")" << endl;
    cout << "  2. delete all duplications of transcripts and save the variant with the highest \"score\"" << endl;
    cout << "  3. if sequence ranges are set for some transcripts, the program detects, whether the distance to that range is dangerously close" << endl;
    cout << "  4. join:" << endl;
    cout << "    4.1. if there is a transcript dangerously close to one/both end(s) of a sequence range, the program creates a copy without the corresponding terminal exon" << endl;
    cout << "    4.2. if there is a transcript with start or stop codon in a set and a second one without this codon and they are \"joinable\", than this step joines the corresponesponding terminal exons" << endl;
    cout << "  5. selection: selects the \"best\" gene structure out of all possible \"maximum\" gene structures" << endl;
    cout << "    - \"maximum\" gene structure is a set of transcripts from an overlap so that there is no other transcript in the overlap, which can be added to the set without producing a \"contradiction\"" << endl;
    cout << "    - a gene structure is \"better\" than another one, if it has the transcript with the highest \"score\", which is not present in the other gene structure" << endl;
    cout << endl;
    cout << "  Options:" << endl;
    cout << "    Necessary:" << endl;
    cout << "\t--genesets=file1,file2,...\t-g file1,file2,...\twhere \"file1,file2,...,filen\" have to be data files with genesets in GTF format." << endl;
    cout << "\t--output=ofile\t\t\t-o ofile\t\twhere \"ofile\" is the name for an output file (GTF)." << endl;
    cout << "    Optional:" << endl;
    cout << "      Parameters:" << endl;
    cout << "\t--priorities=pr1,pr2,...\t-p pr1,pr2,...\t\twhere \"pr1,pr2,...,prn\" have to be positiv integers (different from 0)." << endl;
    cout << "\t\t\t\t\t\t\t\tHave to be as many as filenames are added." << endl;
    cout << "\t\t\t\t\t\t\t\tBigger numbers means a higher priority." << endl;
    cout << "\t\t\t\t\t\t\t\tIf no priorities are added, the program will set all priorties to 1." << endl;
    cout << "\t\t\t\t\t\t\t\tThis option is only useful, if there is more than one geneset." << endl;
    cout << "\t\t\t\t\t\t\t\tIf there is a conflict between two transcripts, so that they can not be picked in the same genestructure, joingenes decides for the one with the highest priority." << endl;
    cout << "\t--errordistance=x\t\t-e x\t\t\twhere \"x\" is a non-negative integer" << endl;
    cout << "\t\t\t\t\t\t\t\tIf a prediction is <=x bases next to a prediction range border, the program supposes, that there could be a mistake." << endl;
    cout << "\t\t\t\t\t\t\t\tDefault is 1000." << endl;
    cout << "\t\t\t\t\t\t\t\tTo disable the function, set errordistance to a negative number (e.g. -1)." << endl;
    cout << "\t--genemodel=x\t\t\t-m x\t\t\twhere \"x\" is a genemodel from the set {eukaryote, bacterium}." << endl;
    cout << "\t\t\t\t\t\t\t\tDefault is eukaryotic." << endl;
    cout << "\t--alternatives\t\t\t-a\t\t\tis a flag" << endl;
    cout << "\t\t\t\t\t\t\t\tIf this flag is set, the program joines different genes, if the transcripts of the genes are alternative variants." << endl;
    cout << endl;
    cout << "      Format changes of input/output:" << endl;
    cout << "\t--suppress=pr1,pr2,..\t\t-s pr1,pr2,...\t\twhere \"pr1,pr2,...,prm\" have to be positiv integers (different from 0)." << endl;
    cout << "\t\t\t\t\t\t\t\tDefault is none." << endl;
    cout << "\t\t\t\t\t\t\t\tif the core of a joined/non-joined transcript has one of these priorities it will not occur in the output file." << endl;
    cout << "\t--stopincoding\t\t\t-i\t\t\tis a flag" << endl;
    cout << "\t\t\t\t\t\t\t\tIf this flag is set, the program joines the stop_codons to the CDS." << endl;
    cout << endl;
    cout << "      Enable/Disable program parts:" << endl;
    cout << "\t--nojoin\t\t\t-j\t\t\tis a flag (to disable joining)." << endl;
    cout << "\t\t\t\t\t\t\t\tIf this flag is set, the program will not join/merge/shuffle;" << endl;
    cout << "\t\t\t\t\t\t\t\tit will only decide between the unchanged input transcipts and output them." << endl;
    cout << "\t--noselection\t\t\t-l\t\t\tis a flag (to disable selection)." << endl;
    cout << "\t\t\t\t\t\t\t\tIf this flag is set, the program will NOT select at the end between \"contradictory\" transcripts." << endl;
    cout << "\t\t\t\t\t\t\t\t\"contradictory\" is self defined with respect to known biological terms." << endl;
    cout << "\t\t\t\t\t\t\t\tThe selection works with a self defined scoring function." << endl;
    cout << endl;
    cout << "      Secondary program functions:" << endl;
    cout << "\t--onlycompare\t\t\t-c\t\t\tis a flag." << endl;
    cout << "\t\t\t\t\t\t\t\tIf this flag is set, it disables the normal function of the program and" << endl;
    cout << "\t\t\t\t\t\t\t\tactivates a compare and seperate mode to seperate equal transcripts from non equal ones." << endl;
    cout << endl;
    cout << "      This help:" << endl;
    cout << "\t--help \t\t\t\t-h\t\t\tprints the help documentation." << endl;
    cout << endl;
    exit( EXIT_FAILURE );
}

void display_error(string const &error)
{
    cerr << "Usage error: " << error << endl;
    display_help();
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

void check_cds_stop_combination(Transcript* tx, pair<string,string> &stopVariants, Properties &properties){
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
	if ((tx->tl_complete.second && tx->stop_list.empty()) || (!tx->tl_complete.second && !tx->stop_list.empty())){display_error("Semantic Error: stop-complete but no stop_codon.");}
	if (!tx->stop_list.empty()){
	    if (tx->tes == tx->stop_list.front().from){
		stopVariants.first = tx->t_id;
	    }else{
		stopVariants.second = tx->t_id;
	    }
	}
    }

    if (!stopVariants.first.empty() && !stopVariants.second.empty()){
	display_error("E.g. the stop_codon of "+(*properties.transcriptMap)[stopVariants.first]->originalId+" from "+(*properties.transcriptMap)[stopVariants.first]->inputFile+" is inside of the CDS and the stop_codon of "+(*properties.transcriptMap)[stopVariants.second]->originalId+" from "+(*properties.transcriptMap)[stopVariants.second]->inputFile+" is outside of the CDS. You have to make it equally!");
    }
}

int main(int argc, char* argv[])
{
    int opt = 0;
    static const char *optString = "g:p:o:e:m:s:jlaich?";
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
    properties.alternatives = false;
    properties.nrOfPrintedGenes = 0;
    properties.unknownCount = 1;
    properties.minimumIntronLength = 20; // hard coded (not optional at the moment)
    properties.stopincoding = false;

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
	    properties.onlyCompare = true;
	    break;

	case 's':
	    supprList = getSuppressionList(optarg);
	    break;

	case 'j':
	    properties.join = false;
	    break;

	case 'l':
	    properties.selecting = false;
	    break;

	case 'a':
	    properties.alternatives = true;
	    break;

	case 'i':
	    properties.stopincoding = true;
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

    if (priorities.size() == filenames.size()){
	unordered_map<string,bool> taxaMap;
	list<int>::iterator it_p = priorities.begin();
	for (list<string>::iterator it_f = filenames.begin(); it_f != filenames.end(); it_f++){
	    load(geneMap, (*it_f) ,(*it_p), taxaMap, properties);
	    cout << "# After loading " << (*it_f) << " (Priority " << (*it_p) << ") there are " << (*properties.transcriptMap).size() << " transcripts in transcript list." << endl;
	    it_p++;
	}

	for (auto pointer = (*properties.transcriptMap).begin(); pointer != (*properties.transcriptMap).end(); pointer++){
	    bool hasCDS = false;
	    for (list<Exon>::iterator it = pointer->second->exon_list.begin(); it != pointer->second->exon_list.end(); it++){
		if ((*it).feature == "CDS"){
		    hasCDS = true;
		    break;
		}
	    }
	    if (!hasCDS){
		cerr << (*properties.transcriptMap)[pointer->second->t_id]->originalId << " from file " << (*properties.transcriptMap)[pointer->second->t_id]->inputFile << " has no CDS and will be deleted." << endl;
		deleteTx(pointer->second, properties);
		continue;
	    }
	    if ((*(pointer->second)).exon_list.size() == 0){
		deleteTx(pointer->second, properties);
	    }else{
		(*(pointer->second)).initiate(properties);
	    }
	}
    }else{
        display_error("Number of input files and priorities is not equal.");
    }

    pair<string,string> stopVariants;		// <stop_codon is in CDS, stop_codon is not in CDS>

    for (auto pointer = (*properties.transcriptMap).begin(); pointer != (*properties.transcriptMap).end(); pointer++){
	check_cds_stop_combination(pointer->second, stopVariants, properties);
	if (!check_frame_annotation(*pointer->second)){
	    pointer->second->isNotFrameCorrect = true;
	    displayWarning((*properties.transcriptMap)[pointer->second->t_id]->originalId+" from file "+(*properties.transcriptMap)[pointer->second->t_id]->inputFile+" has wrong frame annotation and will be ignored in joining step.", properties, "frame");
	    //display_error("Frames are not correct.");
	}
    }
    warningSummary("There are ", " transcripts with wrong frame in the input files.", properties, "frame");

    // tests, if the alternatives definition in the input is more strict than the one in this program or not
//    testInputAlternatives(properties);

    // take a look at length of introns
/*    list<int> intronLengths;
    for (auto pointer = (*properties.transcriptMap).begin(); pointer != (*properties.transcriptMap).end(); pointer++){
	for (list<Exon>::iterator it1 = pointer->second->exon_list.begin(); it1 != pointer->second->exon_list.end(); it1++){
	    list<Exon>::iterator it2 = it1;
	    it2++;
	    if (it2 != pointer->second->exon_list.end()){
		int intronLength = (*it2).from - (*it1).to - 1;
		intronLengths.push_back(intronLength);
if (intronLength == 89114){cout << pointer->first << " " << intronLength << endl;}
	    }
	}
    }
    intronLengths.sort();
    for (list<int>::iterator it = intronLengths.begin(); it != intronLengths.end(); it++){
	cout << (*it) << endl;
    }*/

    // splits the transcript list by the scaffold of the transcript
    unordered_map<string,list<Transcript*>> splitted_transcript_list;
    for (auto pointer = (*properties.transcriptMap).begin(); pointer != (*properties.transcriptMap).end(); pointer++){
	splitted_transcript_list[(*pointer->second).getChr()].push_back(pointer->second);
    }

    fstream outfile;
    outfile.open(properties.outFileName, ios::out);		// delete content of file filename
    outfile.close();


    if (properties.onlyCompare){
	string filenameEqual = "cEqual.gtf", filenameStop1 = "cAlternative1.gtf", filenameStop2 = "cAlternative2.gtf", filenameUne1 = "cUnequal1.gtf", filenameUne2 = "cUnequal2.gtf";
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
