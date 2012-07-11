/*	Prints the coverage of a BAM alignment file

	Created: 13-Sep-11
	Last modified: 11-July-2012
*/

#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <api/BamAlignment.h> 
#include "bamtools_count.h"
#include "bamtools_sort.h" 
#include <api/algorithms/Sort.h> 
#include <algorithm>
#include <iostream>
#include <vector>
#include <string> 
#include <stdio.h>
#include <stdexcept>
#include <time.h> 
#include <getopt.h> 
#include <sstream>
#include <stdlib.h>
#include <iomanip> 
#include <unordered_map>
#include "filterBam.h"

using namespace std; 
using namespace BamTools;
using namespace BamTools::Algorithms; 

struct optionalCounters_t {
  int outPaired;
  int outUniq;
  int outBest;
};


int main(int argc, char *argv[])
{
  // Variable definition
  BamReader reader;
  BamWriter writer;
  vector<BamAlignment> qali;
  unordered_map<string, int> qNameStems;
  unordered_map<string, int>::iterator it;
  time_t tStart, tEnd, tElapsed; 
  vector<CigarOp> cigar;
  char cigarType;
  uint32_t cigarLength;
  uint32_t qLength;
  uint32_t RefID;
  uint32_t editDistance;
  string rName;
  string qName;
  string qNameStem = ""; 
  string oldQnameStem = "";
  // Alignment
  string qSuffix;
  bool strand;
  int rstart, rend;
  int cigarSize; 
  int sumMandI;
  float coverage;
  float percId;
  int32_t baseInsert;  // baseInsert = qBaseInsert+tBaseInsert, in some cases
  // Counters
  int line = 0;
  int outNoMapping = 0;
  int outMinId = 0;
  int outMinCover = 0;
  int outIntrons = 0;
  int outPaired = 0;
  int outUniq = 0;
  int outBest = 0;
  optionalCounters_t optionalCounters;
  int counter=0;
  int sizeSal;
  int sizePtrSal;
  int alignmentCount = 0;

  // Initialising options
  struct globalOptions_t globalOptions;
  globalOptions = initOptions(argc, argv);
  bool best = globalOptions.best;
  bool help = globalOptions.help;
  bool noIntrons = globalOptions.noIntrons;
  bool paired = globalOptions.paired;
  bool uniq = globalOptions.uniq;
  bool verbose = globalOptions.verbose;
  int insertLimit = globalOptions.insertLimit;
  int maxIntronLen = globalOptions.maxIntronLen;
  int maxSortesTest = globalOptions.maxSortesTest;
  int minCover = globalOptions.minCover;
  int minId = globalOptions.minId;
  int minIntronLen = globalOptions.minIntronLen;
  float uniqThresh = globalOptions.uniqThresh;
  const char* commonGeneFile = globalOptions.commonGeneFile;
  const char* inputFile = globalOptions.inputFile;
  const char* outputFile = globalOptions.outputFile;	
  const char* pairbedFile = globalOptions.pairBedFile;
  // cout << "best=" << best << endl;
  // cout << "help=" << help << endl;
  // cout << "noIntrons=" << noIntrons << endl;
  // cout << "paired=" << paired << endl;
  // cout << "uniq=" << uniq << endl;
  // cout << "verbose=" << verbose << endl;
  // cout << "Input file is " << inputFile << endl;
  // cout << "Output file is " << outputFile << endl;
  // cout << "Value of paired is " << paired << endl;
  // cout << "insertLimit=" << insertLimit << endl;
  // cout << "maxIntronLen=" <<  maxIntronLen << endl;
  // cout << "minCover=" <<  minCover << endl;
  // cout << "minId=" <<  minId << endl;
  // cout << "minIntronLen=" << minIntronLen << endl;
  // cout << "uniqThresh=" << uniqThresh << endl;

  // Starting timer
  tStart = time(NULL);    


  // Opening BAM file
  if (!reader.Open(inputFile)) 
	{
	cerr << "Could not open input BAM files. " << endl;
	return 0;
  	}
 
  // Counting number of alignments, so that memory can be allocated
  BamAlignment al;
  while (reader.GetNextAlignment(al)) {++alignmentCount;} 	
  reader.Close();
  cout << "The number of alignments in this file is: " << alignmentCount << endl;
  BamAlignment *pal = new BamAlignment[alignmentCount];
  vector<SingleAlignment*> vecSingleAlignments;
  SingleAlignment *ptrSal2;
  cout << "Allocating space for (" << alignmentCount << ") single-alignments x " << sizeof(ptrSal2) 
	   << " (bytes) = " << alignmentCount*sizeof(ptrSal2)/1024 << " kbytes" << endl;

  // Re-Opening BAM file and filtering
  reader.Open(inputFile);

  // Read reference data
  const RefVector refData = reader.GetReferenceData();
  const SamHeader header = reader.GetHeader();

  // Open BamWriter, HEADER and REFERENCE will be written into the file
  if (!writer.Open(outputFile, header, refData)) 
	{
	cerr << "Could not open output BAM file" << endl;
	return 0;
	}

  // Scanning through alignments, reading CIGAR strings and printing
  // useful stats
  nextAlignment:
    while (reader.GetNextAlignment(pal[line])) 
	  {

		// Displaying only data whose Reference seq ID has a defined mapping
		if (!(pal+line)->IsMapped()) // Unmapped alignment
		  {	 
			if (verbose) cout << "alignment " << line << ": unmapped" << endl;
			outNoMapping++;
			goto nextAlignment;
		  }
 
		cigar = pal[line].CigarData;
		qName = pal[line].Name; // query name
		qNameStem = qName.substr(0, qName.find("/")); // Equal to: $qNameStem =~ s/[\.\/]([fr12])$//; 
		qSuffix = qName.substr(qName.find("/")+1, qName.length());	
		qLength = pal[line].Length; // query length
		RefID = pal[line].RefID;
		sumMandI = 0; 
		pal[line].GetTag("NM", editDistance); // edit distance

 
		// Computing the pctge of Identity
		percId = 100*(qLength-editDistance)/qLength;
		// Filtering by percentage identity
		if (percId < minId)
		  {
			outMinId++;
			if (verbose) cout << "alignment " << line << ": out by percId=" << percId << ", minId=" << minId << endl;
			goto nextAlignment;
		  }		


		// Computing the coverage
		sumMandI = sumMandIOperations(cigar, "no");
		coverage = (float)100*sumMandI/qLength;
		// Filtering by coverage
		if (coverage < minCover)
		  {
			outMinCover++;
			if (verbose) cout << "alignment " << line << ": out by coverage=" << coverage << ", minCover=" << minCover << endl;
			goto nextAlignment;
		  }		

  		// Intron gap filter
  		baseInsert = sumDandIOperations(cigar, "no");
  		// Save if complying with insertLimit
  		if (noIntrons && baseInsert > insertLimit)
  		  {
  			outIntrons++;
			if (verbose) cout << qName << " filtered out by intron criterion= " << baseInsert << " < minId=" << insertLimit << endl;
  			goto nextAlignment;
  		  }


		// Storing alignment info into location addressed by vector of pointers
		SingleAlignment *ptrSal = new SingleAlignment(pal+line, coverage, percId, &refData);
		sizePtrSal = sizeof(ptrSal);
		vecSingleAlignments.push_back(ptrSal);
		// cout << "alignment " << line << ":" << *ptrSal << endl;
		// cout << *vecSingleAlignments.at(counter) << endl;
		counter++;

		// Saving alignments that passed the filters
  		writer.SaveAlignment(pal[line]); 

 		line++;

	  } // end while	
  

	// Closing file handles
    reader.Close();
	writer.Close();


	// Displaying results to STDOUT
	cout <<  "\n\n------------------------------------- " << endl;
	cout <<  "\nSummary of filtered alignments: " << endl;
	cout <<  "------------------------------------- " << endl;
	cout <<  "unmapped        : " << outNoMapping << endl;
	cout <<  "percent identity: " << outMinId << endl;
	cout <<  "coverage        : " << outMinCover << endl;
	if (noIntrons) {cout <<  "nointrons       : " << outIntrons << endl;}
	if (paired) 
	  {
		cout << "not paired      : " << outPaired << endl;
      	cout << "quantiles of unspliced insert lengths:" << endl;
		// try { // catches: instance of 'std::out_of_range'
		//   	  for (int it=1; it<10; it++)
		// 		{cout << "q[" << 10*it << "%]=" << insertlen.at(ceil(it*insertlen.size()/10)) << ",";}
		// 	} catch (out_of_range& oor) { 
		//   		  cerr << "[Quantiles]:" << oor.what() << "; insertlen.size()=" << insertlen.size() << endl; 
		// 	}
		cout << endl;
 	  } 
	if (uniq) {cout << "unique          : " << outUniq << endl;}
	if (best) {cout << "best            : " << outBest << endl;}



	// Displaying command line
	cout <<  "------------------------------------- " << endl;
	cout << "Cmd line: " << endl;
	for (int it=0; it<argc;it++)
	  {
		cout << argv[it] << " ";
	  }
	cout << endl;
	cout <<  "------------------------------------- " << endl;


	// // Printing sizes associated to vector of: vecSingleAlignments
	// cout << "----------------------------------------------------" << endl;
	// cout << "Summary statistics: vecSingleAlignments " << endl;
	// cout << "----------------------------------------------------" << endl;
	// // cout << "Size of 1 alignment: " << sizeof(al) << " bytes" << endl;
	// cout << "No. alignments stored in vecSingleAlignments vector: " << (int)vecSingleAlignments.size() << "\n";
	// cout << "capacity of vecSingleAlignments vector: " << (int) vecSingleAlignments.capacity() << "\n";
	// cout << "max_size of vecSingleAlignments vector: " << (int) vecSingleAlignments.max_size() << "\n";
	// cout << "Size of 1 pointer to SingleAlignment is: " << sizePtrSal << " bytes " << endl;
	// cout << "Size of (" << (int) vecSingleAlignments.size() << ") pointers to single-alignments is: "  
	// 	 << sizePtrSal << " (bytes) x " << vecSingleAlignments.size() << " (pointer to alignments) = " 
	// 	 << sizePtrSal*vecSingleAlignments.size()<< " bytes \n";
	// cout << "Size of vector container: " << sizeof(vecSingleAlignments) << " bytes " << endl;
	// cout << "Size of whole vector of pointers (container+contents) = " << sizeof(vecSingleAlignments)+sizePtrSal*vecSingleAlignments.size() << " bytes" << endl;

	// cout << "Out of printQali, vecSingleAlignments is located at " << &vecSingleAlignments << endl; 
	// cout << "-----------------------------------------------------------" << endl;
	// cout << "Before sorting [Using memory-infefficient function]:" << endl;
	// printQaliOld(vecSingleAlignments);
	// cout << "-----------------------------------------------------------" << endl;
	// cout << "After sorting [Using memory-friendlier function]:" << endl;
	// std::stable_sort(vecSingleAlignments.begin(), vecSingleAlignments.end(), byScore);
	// printQali(vecSingleAlignments);
	// cout << "-----------------------------------------------------------" << endl;

 
}


