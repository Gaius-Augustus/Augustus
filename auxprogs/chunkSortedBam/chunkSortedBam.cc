/*	Prints the coverage of a BAM alignment file

	Created: 13-Sep-11
	Last modified:  7-May-2012
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


int main(int argc, char *argv[])
{
  // Variable definition
  BamReader reader;
  BamWriter writer;
  time_t tStart, tEnd, tElapsed; 
  uint32_t qLength;
  uint32_t RefID;
  string qName = "";
  string qNameStem = ""; 
  string oldQnameStem = "";
  // Alignment
  string qSuffix;
  float coverage = 0;
  float percId = 0;
  // Counters
  int line = 0;
  // Initialising options
  struct globalOptions_t globalOptions;
  globalOptions = initOptions(argc, argv);
  const char* inputFile = globalOptions.inputFile;
  const char* outputFile = globalOptions.outputFile;
  int counter=0;
  int sizeSal;
  int sizePtrSal;
  int alignmentCount = 0;


  // Starting timer
  tStart = time(NULL);    


  // Opening BAM file
  if (!reader.Open(inputFile)) 
	{
	cerr << "Could not open input BAM files. " << endl;
	return 0;
  	}
 

  cout << "Allocating memory [Counting alignments]" << endl;
  BamAlignment al;
  while (reader.GetNextAlignment(al)) {++alignmentCount;} 	
  reader.Close();
  cout << "Number of alignments is: " << alignmentCount << endl;
  BamAlignment *pal = new BamAlignment[alignmentCount];
  vector<SingleAlignment*> vecSingleAlignments;
  SingleAlignment *ptrSal2;
  cout << "Allocating space for (" << alignmentCount << ") single-alignments x " << sizeof(ptrSal2) 
	   << " (bytes) = " << alignmentCount*sizeof(ptrSal2) << " bytes" << endl;

  // Re-Opening BAM file and scanning
  reader.Open(inputFile);

  // Read reference and header data
  const RefVector refData = reader.GetReferenceData();
  const SamHeader header = reader.GetHeader();

  // Open BamWriter, HEADER and REFERENCE will be written into the file
  if (!writer.Open(outputFile, header, refData)) 
	{
	cerr << "Could not open output BAM file number: " << endl;
	return 0;
	}

  // Scanning through alignments
    while (line < 10 && oldQnameStem.compare(qNameStem) == 0)  
	  {

		reader.GetNextAlignment(pal[line]);
		qName = pal[line].Name; // query name
		qNameStem = qName.substr(0, qName.find("/")); // Equal to: $qNameStem =~ s/[\.\/]([fr12])$//; 
		cout << "qNameStem = " << qNameStem << endl;

		// Storing alignment info into location addressed by vector of pointers
		SingleAlignment *ptrSal = new SingleAlignment(pal+line, coverage, percId, &refData);
		cout << *ptrSal << endl;

		// Saving alignments that passed the filters
  		writer.SaveAlignment(pal[line]); 
 		line++;

		oldQnameStem = qNameStem;

	  } // end while	
  

	// Closing file handles
    reader.Close();
	writer.Close();


} // end main


