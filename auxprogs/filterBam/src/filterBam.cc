/*	
	NOTE: Filters BAM alignments according to:

	1. Percentage identity (percId)
	2. Coverage (coverage)
	3. Intron gaps (noIntrons, baseInsert)
	4. Uniqueness (determined by percId and coverage)

*/

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
#include <cstdio>
#include <algorithm>
#include <unordered_map> 
#include <fstream>  
#include <map>
#include <set>
#include <math.h>  
#include "filterBam.h"
#include "bamaccess.hh"

#ifdef USE_SEQLIB
#include "bamseqlibaccess.hh"
#else
#include "bamtoolsaccess.hh"
#endif


struct optionalCounters_t { 
  int outPaired;
  int outUniq;
  int outBest;
}; 

using namespace std;

void printQali(const vector<BamAlignmentRecord_> &qali, bool pairwiseAlignments);
float scoreMate(const BamAlignmentRecord_ &al1, const BamAlignmentRecord_ &al2, int &dist, const globalOptions_t &globalOptions);
void printMatedMap(const map<int,int> mated);
void processQuery(vector<BamAlignmentRecord_> &qali, const globalOptions_t &globalOptions, std::shared_ptr<BamFileWriter> &writer, string oldQnameStem, optionalCounters_t &optionalCounters, vector<PairednessCoverage> &pairCovSteps, vector<int> &insertlen, map<string, multimap<int,int>> &pairCovSteps2);
void printPairCovSteps(const vector<PairednessCoverage> &pairCovSteps);
void printChrOfPairCovSteps(vector<PairednessCoverage> &pairCovSteps, string chr);
vector<PairednessCoverage> compactifyBed(vector<PairednessCoverage> &pairCovSteps, const globalOptions_t &globalOptions);
void printSizeOfCoverInfo(const vector<PairednessCoverage> &pairCovSteps);
string bool_cast(const bool b);
vector<string> uniqueKeys(const vector<PairednessCoverage> &m);


int main(int argc, char *argv[])
{
    struct globalOptions_t globalOptions;
    globalOptions = initOptions(argc, argv);

  // Variable definition
#ifdef USE_SEQLIB
  std::shared_ptr<BamFileReader> reader = std::make_shared<BamSeqLibReader>(globalOptions);
  std::shared_ptr<BamFileWriter> writer= std::make_shared<BamSeqLibWriter>();
#else  
  std::shared_ptr<BamFileReader> reader = std::make_shared<BamToolsReader>(globalOptions);
  std::shared_ptr<BamFileWriter> writer = std::make_shared<BamToolsWriter>();
#endif

  BamAlignmentRecord_ al;
  vector<BamAlignmentRecord_> qali;
  unordered_map<string,int> qNameStems;
  time_t tStart, tEnd;
  int32_t qLength;
  int32_t editDistance;
  string qName;
  string qNameStem = ""; 
  string oldQnameStem = "";
  // Alignment
  string qSuffix;
  float coverage;
  float percId;
  uint32_t baseInsert;  // baseInsert = qBaseInsert+tBaseInsert, in some cases
  bool useEqualSigns = false;
  bool checkForEqualSigns = true;
  uint32_t numEquals;
  bool is_NM_WarningPrinted = false;
  // Counters
  int line = 0;
  int outMap = 0;
  int outMinId = 0;
  int outMinCover = 0;
  int outIntrons = 0;
  int singleton = 0, notMateMapped = 0, notPaired = 0, notOnSameTarget = 0;
  bool isSingleton; 
  optionalCounters_t optionalCounters;
  // Pairedness coverage
  map<string, multimap<int,int>> pairCovSteps2;
  vector<PairednessCoverage> pairCovSteps;
  vector<PairednessCoverage> compactPairCovSteps;  
  // Optional counters
  optionalCounters.outPaired = 0;
  optionalCounters.outUniq = 0;
  optionalCounters.outBest = 0;
  // Initializing options
  bool best = globalOptions.best;
  bool help = globalOptions.help;
  bool noIntrons = globalOptions.noIntrons;
  bool paired = globalOptions.paired;
  bool uniq = globalOptions.uniq;
  bool verbose = globalOptions.verbose;
  bool pairwiseAlignments = globalOptions.pairwiseAlignments;
  uint32_t insertLimit = (unsigned) globalOptions.insertLimit;
  int maxIntronLen = globalOptions.maxIntronLen;
  int maxSortesTest = globalOptions.maxSortesTest;
  int minCover = globalOptions.minCover;
  int minId = globalOptions.minId;
  float uniqThresh = globalOptions.uniqThresh;
  const char* commonGeneFile = globalOptions.commonGeneFile;
  const char* inputFile = globalOptions.inputFile;
  const char* outputFile = globalOptions.outputFile;	
  const char* pairBedFile = globalOptions.pairBedFile;
  vector<int> insertlen;

  if (verbose)
	{
	  cout << "------------------------------------------------" << endl;
	  cout << "Selected options are: " << endl;
	  cout << "Input file: " << inputFile << endl;
	  cout << "Output file: " << outputFile << endl;
	  cout << "uniq=" << uniq << endl;
	  cout << "uniqThresh=" << uniqThresh << endl;
	  cout << "best=" << best << endl;
	  cout << "minId=" <<  minId << endl;
	  cout << "minCover=" <<  minCover << endl;
	  cout << "noIntrons=" << noIntrons << endl;
	  cout << "maxIntronLen=" <<  maxIntronLen << endl;
	  cout << "insertLimit=" << insertLimit << endl;
	  cout << "paired=" << paired << endl;
	  cout << "pairwiseAlignments=" << pairwiseAlignments << endl;
	  cout << "pairBedFile=" << pairBedFile << endl;
	  cout << "commonGeneFile=" << commonGeneFile << endl;
	  cout << "threads=" << globalOptions.threads << endl;
	  cout << "verbose=" << verbose << endl;
	  cout << "help=" << help << endl;
	  cout << "------------------------------------------------" << endl;
	}

  // Starting timer
  tStart = time(NULL);    

  // Opening BAM file
  if (!reader->openReader(inputFile))
	{
	cerr << "Could not open input BAM files. " << endl;
	return 0;
  	}

  // Open BamWriter, HEADER and REFERENCE will be written into the file
  if (!writer->openWriter(outputFile, *reader))
	{
	cerr << "Could not open output BAM file" << endl;
	return 0;
	}

  // Sweeping through alignments
  	while (reader->getNextAlignmentRecord(al))
	  {
		line++;
	   
       if (line%100000==1)
	   { 
	   	 printf("\r processed line %d", line);
       }		

		 // Call to compactify bed
	   if ((unsigned)strlen(pairBedFile) > 0 && line%10000000 == 0)
		  {
			cout << "\nCompactifying coverage after " << line << " lines..." << endl;
			compactPairCovSteps = compactifyBed(pairCovSteps, globalOptions);
			cout << "done\n" << endl;
		  }


		// Update qnamestem with the current line's query name 
 		// if option 'paired' used then it expects .f,.r or /1,/2 
		// suffixes of mate pairs;
		qName = al->getQueryName();
		qNameStem = qName; 

		if (paired)
		  {	

			if (!pairwiseAlignments)
			  {
				if (verbose) {cout << "qName=" << qName << endl;}
				// Taking out suffix: {f,r} or {1,2}
				if (qName.find('/') != std::string::npos)
				  {
					qNameStem = qName.substr(0, qName.find("/"));
					qSuffix = qName.substr(qName.find("/")+1, qName.length());	
				  } else {
				  	if (qName.find_last_of('-') != std::string::npos) 
						{
					  		if (verbose) {cout << "Found slash as indicator of pairedness" << endl;}
					  		qNameStem = qName.substr(0, qName.find_last_of("-"));
					  		qSuffix = qName.substr(qName.find_last_of("-")+1, qName.length());	
							if (verbose) {cout << "qNameStem=" << qNameStem << " and qSuffix=" << qSuffix << endl;}
						} 
				  } 

			  } else {//No signs to separate reads
			  	  qNameStem = qName.substr(0, qName.length());
					if (al->isFirstMate()) {
							qSuffix = "1";
					} else if (al->isSecondMate()) {
							qSuffix = "2";
					} else { if(verbose){cout << "Multiple fragments" << endl;}}
				  if (verbose) {cout << "qNameStem=" << qNameStem << " and qSuffix=" << qSuffix << endl;}

			}
		  } 

		// Filter for data whose Reference seq ID is not defined; i.e. RNAME= * in SAM format;
		// i.e. unmapped fragment without coordinate 
		if (!al->isMapped())
		  {	  
			if (verbose)
			  {
				cout << qName << " filtered out because: unmapped " << endl;
			  }
			outMap++;
			continue; // go to next alignment
		  }


        // skip if alignment is not paired, mapped, nor mate is mapped
		if (pairwiseAlignments)
		  {
			isSingleton = al->isPaired() && al->isMapped() && !al->isMateMapped();
			if( isSingleton ) {singleton++;}
			if ( !al->isPaired() ) {notPaired++; continue;} // go to next alignment
			if ( !al->isMateMapped() ) {notMateMapped++; continue;} // go to next alignment
	        // skip if alignment & mate not on same reference sequence
			if ( al->getRefID() != al->getMateRefID() ) {notOnSameTarget++; continue;} // go to next alignment
		}

	
		// Verifying file is sorted by query name
		// cout << line << ": oldQnamestem=" << oldQnameStem << ", qNameStem=" << qNameStem << endl;
		if (oldQnameStem.compare(qNameStem) && oldQnameStem.compare(""))   
		  { 

			// Checking whether file is sorted by query name
			if (line <= maxSortesTest && qNameStems[qNameStem])   
			  {
				cerr << "\nInput file " << inputFile << " not sorted by query name!\n" << qNameStem.c_str() << 
				  " occurred previously. " << endl;
				cerr << "Do either of the following: " << endl;
				cerr << "1) Convert the file into SAM with e.g. the 'bamtools' software. \n" << 
				  "   Then sort lexicographically by queryname, i.e. use the command" << endl;
				cerr << "   'export LC_ALL=C' and then 'sort -k 1,1'\n" 
					 << "   Convert back again into BAM format, with e.g. 'samtools' software. " << endl;
				cerr << "2) Sort BAM file directly with your preferred software package, " << 
				  		"   e.g. 'samtools' or 'bamtools' " << endl;  	
				cerr << "3) In the case of samtools, the command is the following: \n" <<
				  "   'samtools sort [-n] file.bam' \n" <<  
				  "   [-n] should sort by query name, just as 'sort -k 10,10' would do in a PSL file." << endl;
				cerr << "   Without options, the sorting will be done by reference name and target coordinate," << endl;  	
				cerr << "   just as a 'sort -n -k 16,16 | sort -k 14,14' would do with PSL." << endl;  	
				cerr << "   For more information check the man page included in samtools distribution." << endl;
				cerr << "4) bamtools can also sort bam files 'bamtools sort -queryname -in file.bam'," << endl;
				cerr << "   but only provides the option to do it by queryname." << endl;	

				exit(1);
			  }	

			if (qali.size()>0)
			  {
				processQuery(qali, globalOptions, writer, oldQnameStem, optionalCounters, pairCovSteps, insertlen, pairCovSteps2);
			  }
		  }  // end outer if

  		/////////////////////////////////////////////////////////////////
  		// Filters for: "percId", "coverage" and intron gaps "noIntrons"   
  		/////////////////////////////////////////////////////////////////

  		// Fetching alignment information
  		qName = al->getQueryName(); // query name
  		qLength = al->getQuerySequenceLength(); // query length (TODO: consider situations where qLength=0,undefined)
		
		if (checkForEqualSigns) {
		    if ((line - outMap - notPaired - notMateMapped - notOnSameTarget) > 100) {
		        // check first 100 (valid) lines for equal signs in aligned sequence
		        // equal signs present indicates "samtools calmd -e" was applied on bam file
		        checkForEqualSigns = false;
		    }
		}

		if (useEqualSigns) {
		    // Percentage Identity filter; compute equal sign count
		    percId = (float) 100 * al->countEqualSignsInQuerySequence() / qLength;
		}
		else if (checkForEqualSigns && (numEquals = al->countEqualSignsInQuerySequence()) != 0) {
		    // cout << "BAM file seems to have been pre-processed with calmd." << endl;
		    // cout << "Computing percentage identity by counting number of (=) signs in SEQ field." << endl;
		    percId = (float) 100 * numEquals / qLength;
		    useEqualSigns = true;
		    checkForEqualSigns = false;
		}
		else if (al->getTagData("NM", editDistance)) {
			// "NM" - Edit distance tag, which records the Levenshtein distance between the read and the reference.
			percId = (float)100*(qLength-editDistance)/qLength;  
		  } else {
			percId = 100.0f;
			if (!is_NM_WarningPrinted) {
				std::cout << "WARNING: Some alignments were missing the NM tag, and could therefore "
				          << "not be filtered based on edit distance." << std::endl;
				is_NM_WarningPrinted = true;
			}
		  }

  		if (percId < minId)
  	  	{
  			outMinId++;
			if (verbose)
			  {
				cout << qName << " filtered out by percid=" << percId << " < minId=" << minId << endl;
			  }
			continue; // go to next alignment
  	  	}


  		// Coverage filter
		coverage = (float) 100 * al->sumMandIOperations() / qLength; // sumMandIOperations() - Equiv to $qEnd-$qStart in PSL
  		if (coverage < minCover)
  		  {	
  			outMinCover++;
			if (verbose)
			  {
				cout << qName << " filtered out by coverage= " << coverage << " < minCover=" <<minCover << endl;
			  }
			continue; // go to next alignment
  		  }	


  		// Intron gap filter
  		baseInsert = al->sumDandIOperations();
  		// Save if complying with insertLimit
  		if (noIntrons && baseInsert > insertLimit)
  		  {
  			outIntrons++;
			if (verbose)
			  {
				cout << qName << " filtered out by intron criterion= " << baseInsert << " > insertLimit=" << 
					insertLimit << endl;
			  }
			continue; // go to next alignment
  		  }

	
		if (verbose)
		  {
			cout << qName << " passed with parameters (percId, coverage)=(" << percId << "," << 
			  coverage << ")"; 
			if (noIntrons)
			 { cout << ", baseInsert=" << baseInsert << " > insertLimit=" << insertLimit << endl;} 
				else {cout << endl;}; 
		  }


		// Appending coverage and percId into alignment
		al->setPercId(percId);
		al->setCoverage(coverage);

  		// Push @qali, [$_, $targetname, $qsuffix, $strand, $tstart, $tend, $percId, $coverage]; 
        qali.push_back(al);

		oldQnameStem = qNameStem;
		if (line <= maxSortesTest)
		  {
			qNameStems[qNameStem] = 1;
		  }

	  } ////// end while 


	// Calling of: "processQuery() if ($qnamestem ne "");
	if (qNameStem.compare(""))
	  {
		processQuery(qali, globalOptions, writer, oldQnameStem, optionalCounters, pairCovSteps, insertlen, pairCovSteps2);

		if (verbose) {printPairCovSteps(pairCovSteps);}
	  }


	// Write pairedness coverage data into pairBedFile
	if ((unsigned)strlen(pairBedFile) > 0)
	  {
		  ofstream bedFile;
		  bedFile.open(pairBedFile);
		  bedFile << "track type=bedGraph name=\"pairedness coverage\" description=\"pairedness coverage\"";
		  bedFile << " visibility=full color=200,100,0 altColor=200,100,0\n";
		  cout << "pairBedFile selected. Calling compactifyBed..." << endl;
		  printSizeOfCoverInfo(pairCovSteps);
		  vector<PairednessCoverage> compactPairCovSteps = compactifyBed(pairCovSteps, globalOptions); 
		  printSizeOfCoverInfo(compactPairCovSteps);
 		  vector<string> chrNames = uniqueKeys(compactPairCovSteps);
		  vector<PairednessCoverage> pairCovStepsOfChr;
		  string chr;

		  cout << "Number of chromosomes is " << chrNames.size() << endl;
		  for (unsigned int chrIt = 0; chrIt < chrNames.size(); chrIt++ ) // display chromosomes in alph order
			{
			  chr = chrNames.at(chrIt);
			  if (chr.empty()) continue;

			  // Sweep through each chromosome and stuff their position values into a vector 
			  for (unsigned int it=0; it<compactPairCovSteps.size(); it++)
				{
				  if (compactPairCovSteps.at(it).chr == chr)
					pairCovStepsOfChr.push_back(compactPairCovSteps.at(it)); 
				}

			  int cov = 0;
			  int pos = 0;
			  PairednessCoverage step;
			  for (unsigned int jit=0; jit<pairCovStepsOfChr.size(); jit++)
				{
				  step = pairCovStepsOfChr.at(jit); 
				  if (pos < step.coord && cov > 0)
				  	{
					  bedFile << chr << "\t" << pos << "\t" << step.coord << "\t" << cov << endl;
					} 
				  pos = step.coord;
				  cov += step.label;
				}

			  if (cov!=0)
				{cout << "Inconsistent" << endl;}

			} // end outer for

		  // Closing bed file
		  bedFile.close();
	  }


	// Sorting insertlengths; i.e. @insertlen = sort {$a <=> $b} @insertlen;
	std::stable_sort(insertlen.begin(), insertlen.end());


	// Displaying results to STDOUT
	cout <<  "\n\n------------------------------------- " << endl;
	cout <<  "Processed alignments: " << line << endl;
	cout <<  "\nSummary of filtered alignments: " << endl;
	cout <<  "------------------------------------- " << endl;
	cout <<  "unmapped        : " << outMap << endl;
	cout <<  "percent identity: " << outMinId << endl;
	cout <<  "coverage        : " << outMinCover << endl;
	if (noIntrons) {cout <<  "nointrons       : " << outIntrons << endl;}
	if (paired) 
	  {
		cout << "not paired (our criterion)	: " << optionalCounters.outPaired << endl;
      	cout << "quantiles of unspliced insert lengths:" << endl;
		try { // catches: instance of 'std::out_of_range'
			for (unsigned int it=1; it<10; it++) {
			    unsigned int pos = ceil(it*insertlen.size()/10);
			    unsigned int quantile = 0;
			    if (pos < insertlen.size()) {
			        quantile = insertlen.at(pos);
			    }
			    cout << "q[" << 10*it << "%]=" << quantile << ",";
			}
			} catch (out_of_range& oor) { 
		  		  cerr << "[Quantiles]:" << oor.what() << "; insertlen.size()=" << insertlen.size() << endl; 
			}
		cout << endl;
 	  } 
	if (uniq) {cout << "not unique	: " << optionalCounters.outUniq << endl;}
	if (best) {cout << "not best	: " << optionalCounters.outBest << endl;}
	if (pairwiseAlignments)
 	  {
		cout << "singleton	: " << singleton << endl;
		cout << "on diff. target : " << notOnSameTarget << endl;
	  }


	// Displaying command line
	cout <<  "------------------------------------- " << endl;
	cout << "Cmd line: " << endl;
	for (int it=0; it<argc;it++)
	  {
		cout << argv[it] << " ";
	  }
	cout << endl;
	cout <<  "------------------------------------- " << endl;

	if (verbose)
	  {
		// Printing sortedness test results
		cout << "----------------------------------------" << endl;
		cout << "File is correctly sorted " << endl;
		cout << "----------------------------------------" << endl;
	  }

	// Closing file handles
	reader->close();
	writer->close();
    
	// Ending timer and printing elapsed time
	tEnd = time(NULL);    
	printElapsedTime(tEnd, tStart);


} // end main


string bool_cast(const bool b) 
{
    ostringstream ss;
    ss << boolalpha << b;
    return ss.str();
}



vector<string> uniqueKeys(const vector<PairednessCoverage> &m)
{
  vector<string> unKeys;

  if (m.empty()) {return unKeys;}

  // find the number of different keys
  size_t nKeys = 1;
  vector<PairednessCoverage>::const_iterator it = m.begin();
  string lastkey = it->chr;
  unKeys.push_back(it->chr);
  ++it;

  while (it != m.end()) 
	{
      if (lastkey < it->chr) 
		{
		  ++nKeys;
		  lastkey = it->chr;
		  unKeys.push_back(it->chr);
		}
	  ++it;
	}
  return unKeys;
}


void printQali(const vector<BamAlignmentRecord_> &qali, bool pairwiseAlignments)
{
  vector<BamAlignmentRecord_>::const_iterator it = qali.begin();
  std::stringstream ss_rstart, ss_rend;
  string rName;
  string qName;
  string qSuffix;
  bool strand;
  float percId, coverage, score;

  for (; it != qali.end(); it++)
	{
	  rName = (*it)->getReferenceName();
	  qName = (*it)->getQueryName();
	  if (!pairwiseAlignments) {
		  qSuffix = qName.substr(qName.find("/")+1, qName.length()); 
	  } else {
		if ((*it)->isFirstMate()) {
					qSuffix = "1";
		} else if ((*it)->isSecondMate()) {
					qSuffix = "2";
				} else { cout << "Multiple fragments" << endl;}
	  }
	  strand = (*it)->isReverseStrand();
	  ss_rstart << (*it)->getStartPosition();
	  ss_rend << (*it)->getEndPosition();
	  percId = (*it)->getPercId();
	  coverage = (*it)->getCoverage();
	  score = (*it)->getScore();

	  // Printing to stdout
	  cout << "\t[ ";
	  cout << qName << " " << rName << " " << qSuffix << " " << bool_cast(strand) << " " << 
			ss_rstart.str() << " " << ss_rend.str() << " " << percId << " " << 
		    coverage << " " << score; 
	  cout << " ]," << endl;
	  // clear strings 
	  ss_rstart.str(""); ss_rend.str("");
	}
}

// for comparing quality of two mate-pair read alignments (it, jit)
float scoreMate(const BamAlignmentRecord_ &al1, const BamAlignmentRecord_ &al2, int &dist, const globalOptions_t &globalOptions)
{
  int maxIntronLen = globalOptions.maxIntronLen;
  bool best = globalOptions.best;
  float percId1, percId2;
  string coverage1, coverage2;

  // Retrieving percent identity
  percId1 = al1->getPercId();
  percId2 = al2->getPercId();
  // Retrieving coverage values
  coverage1 = al1->getCoverage();
  coverage2 = al2->getCoverage();

  float coverageIt = atof(coverage1.c_str());
  float coverageJit = atof(coverage2.c_str());

  float score = (percId1 + percId2)/100 // percent identity
				 + (coverageIt + coverageJit)/100; // percent coverage
  if (!best)
  	{ // penalty for distance between mates. Do not use if option 'best' is chosen, 
	  // otherwise a one base difference may cause a difference
	  score -= dist/maxIntronLen/10; 
	}

  return score;
}



//
// checking whether two alignments (or two alignment pairs) are similar
// Purpose: Due to separate handling of spliced and unspliced alignments it can happen
// that very similar alignments are reported, e.g. an unspliced read going approximately up to an intron
// and a spliced read with a few base pairs on one exon.
// These should not be considered ambiguous when --uniq is specified.
// TODO: check whether boolean comparisons can be made with "uint32_t"
bool similar(BamAlignmentRecord_ &alR, BamAlignmentRecord_ &alS, const globalOptions_t &globalOptions) 
{
  // Extracting options
  bool verbose = globalOptions.verbose;
	
  // First alignment
  int32_t rStart = alR->getStartPosition();
  int32_t rEnd = alR->getEndPosition();
  string rName = alR->getQueryName();
  // Second alignment
  int32_t sStart = alS->getStartPosition();
  int32_t sEnd = alS->getEndPosition();
  string sName = alS->getQueryName();

  if (verbose)
	{
	  cout << "[SIMILAR]: checking whether " << rName << " and " << sName << " are approx. the same:"<< endl;
	  cout << "[SIMILAR]: Overlap: (" << rName << "):" << rEnd << " <= (" << sName 
		   << "):" << sStart << " ?" << endl; 
	  cout << "[SIMILAR]: Overlap: (" << rName << "):" << rStart << " >= (" << sName 
		   << "):" << sEnd << " ?" << endl; 
	}

  if (rEnd <= sStart || rStart >= sEnd) // here: similar = overlapping target range
	{
	  if (verbose) {cout << "[SIMILAR]: Alignments are NOT SIMILAR" << endl;}
	  return false;
	} else {
		if (verbose) {cout << "[SIMILAR]: Ans: Alignments are SIMILAR" << endl;}
		return true;
  		   }
}

void printMatedPairsInfo(vector<BamAlignmentRecord_> &qali, vector<MatePairs> &matepairs)
{
  vector<MatePairs>::iterator itMp = matepairs.begin();

  string qName1, qName2;
  cout << "Printing Mate-pairs info: [matepairs.size()=" << matepairs.size() << "]" << endl;
  for (; itMp != matepairs.end(); itMp++)
	{
	  int it = (*itMp).alIt;
	  int jit = (*itMp).alJit;
	  float score = (*itMp).score;
	  qName1 = qali.at(it)->getQueryName();
	  qName2 = qali.at(jit)->getQueryName();
	  cout << "(" << it << "," << jit << ")=" << qName1 << ", " << qName2 << ", scoreMate=" << score << endl;
 	}
}


void printMatedMap(const map<int,int> mated)
{
  map<int,int>::const_iterator itMated = mated.begin();
  cout << "Printing mated summary:" << endl;
  for (; itMated!=mated.end(); itMated++)
	{
	  cout << "mate:" << (*itMated).first << ", mated:" << (*itMated).second << " times" << endl;
	}
}


void printPairCovSteps(const vector<PairednessCoverage> &pairCovSteps)
{
  vector<string> chrNames = uniqueKeys(pairCovSteps);
  string chr;

  cout << "Printing pairCovSteps contents: key=>[value(s)]" << endl;
  // Sweep through each chromosome and print only values associated to a chromosome name
    for (unsigned int chrIt = 0; chrIt < chrNames.size(); chrIt++ ) 
  	  {

		chr = chrNames.at(chrIt);
		cout << chr << "=>" << endl;
		for (unsigned int it = 0; it<pairCovSteps.size(); it++)
		  {
			if (pairCovSteps.at(it).chr == chr)
			  cout << "\t[" << pairCovSteps.at(it).coord << "," << pairCovSteps.at(it).label << "]" << endl; 
		  }
	  }
}


void printChrOfPairCovSteps(vector<PairednessCoverage> &pairCovSteps, string chr)
{

  cout << "Printing coverage of chr=" << chr << endl;
  // Print only values of pairCovSteps that coincide with given chr name
  cout << chr << "=>" << endl;
  for (unsigned int it=0; it<pairCovSteps.size(); it++)
	{
	  if (pairCovSteps.at(it).chr == chr)
		cout << "\t[" << pairCovSteps.at(it).coord << "," << pairCovSteps.at(it).label << "]" << endl; 
	}
}


void printSizeOfCoverInfo(const vector<PairednessCoverage> &pairCovSteps)
{
  vector<string> chrNames = uniqueKeys(pairCovSteps);
  string chr;

  // Sweep through each chromosome and print only values associated to a chromosome name
    for (unsigned int chrIt = 0; chrIt < chrNames.size(); chrIt++ ) 
  	  {
		chr = chrNames.at(chrIt);
		int keySize = 0;
		for (unsigned int it=0; it<pairCovSteps.size(); it++)
		  {
			if (pairCovSteps.at(it).chr == chr)
			  keySize++;
		  }
		cout << "Size of cover. info of chr=" << chr << " is " << keySize << endl;
	  }
}


// If several steps coincide then summarize them equivalently by one step in order to 
// 1) Save memory or 
// 2) output a bed file
vector<PairednessCoverage> compactifyBed(vector<PairednessCoverage> &pairCovSteps, const globalOptions_t &globalOptions)
{
  bool verbose = globalOptions.verbose;
  vector<string> chrNames = uniqueKeys(pairCovSteps);
  std::stable_sort(chrNames.begin(), chrNames.end());
  string chr;
  int32_t before=0, after=0;
  vector<PairednessCoverage> pairCovStepsOfChr;
  vector<PairednessCoverage> compactPairCovSteps;

  // Now sweep through each chromosome and stuff into a vector all the position values corresponding to 
  // that vector
    for (unsigned int chrIt = 0; chrIt < chrNames.size(); chrIt++ ) // display chromosomes in alph order
  	  {
		chr = chrNames.at(chrIt);
		if (chr.empty()) continue; 

 		// Sweep through each chromosome and stuff their position values into a vector 
		for (unsigned int it=0; it<pairCovSteps.size(); it++)
		  {
			if (pairCovSteps.at(it).chr == chr)
			  pairCovStepsOfChr.push_back(pairCovSteps.at(it)); 
		  }


		if (verbose)
		  {
			cout << "------------------------------------------------------\n";
			cout << "Pairedness coverage BEFORE sorting:\n";
			printChrOfPairCovSteps(pairCovSteps, chr);
		  }
		// Foreach "$chr" in pairCovSteps, sort the contents of the "coord" field in ascending order
        // {-1} corresponds to pEnd, and {1} corresponds to pStart
		std::stable_sort(pairCovStepsOfChr.begin(), pairCovStepsOfChr.end());

		if (verbose)
		  {
			cout << "------------------------------------------------------\n";
			cout << "Pairedness coverage AFTER sorting:\n";
			printChrOfPairCovSteps(pairCovSteps, chr);
			cout << "------------------------------------------------------\n";
			cout << "-------------------------------------------------------" << endl;
			cout << "Before removal" << endl;
			printChrOfPairCovSteps(pairCovStepsOfChr, chr);
		  }
 
		before += pairCovStepsOfChr.size();



		unsigned int jit=0;
		while(jit < pairCovStepsOfChr.size()-1)
		  {
			if (pairCovStepsOfChr.at(jit).coord == pairCovStepsOfChr.at(jit+1).coord)
			  {
				pairCovStepsOfChr.at(jit).label += pairCovStepsOfChr.at(jit+1).label;
				pairCovStepsOfChr.erase(pairCovStepsOfChr.begin()+jit+1); // Remove element jit+1
			  } else {
			  	  jit++;
			  }
		  }

		if (verbose)
		  {
			cout << "After removal" << endl;
			printChrOfPairCovSteps(pairCovStepsOfChr, chr);
			cout << "-------------------------------------------------------" << endl;
		  }
		after += pairCovStepsOfChr.size();

		// Insert compact pairCovStepsOfChr into a new vector
		compactPairCovSteps.insert(compactPairCovSteps.end(), pairCovStepsOfChr.begin(), pairCovStepsOfChr.end());
		pairCovStepsOfChr.clear();
	  }

		cout << "\nPairedness coverage before compactifying: " << before << ", after: " << after << endl;

	return compactPairCovSteps;
}



void processQuery(vector<BamAlignmentRecord_> &qali, const globalOptions_t &globalOptions, std::shared_ptr<BamFileWriter> &writer, string oldQnameStem, optionalCounters_t &optionalCounters, vector<PairednessCoverage> &pairCovSteps, vector<int> &insertlen, map<string, multimap<int,int>> &pairCovSteps2)
{
  // Expanding options
  bool best = globalOptions.best;
  bool paired = globalOptions.paired;
  bool uniq = globalOptions.uniq;
  bool verbose = globalOptions.verbose;
  bool pairwiseAlignments = globalOptions.pairwiseAlignments;
  int maxIntronLen = globalOptions.maxIntronLen;
  float uniqThresh = globalOptions.uniqThresh;
  const char* commonGeneFile = globalOptions.commonGeneFile;
  const char* pairBedFile = globalOptions.pairBedFile;
  // Pairedness 
  multimap<int,int> pairCoord;

  // Mate variables
  string itQname, jitQname;
  string itQsuffix, jitQsuffix;
  bool itStrand, jitStrand;
  // Matepairs
  vector<MatePairs> matepairs;
  MatePairs mp;
  map<int,int> mated;
  int32_t inslen, dist;
  // unsigned int inslen, dist;
  unsigned int it, jit;
  uint32_t jitTstart;
  uint32_t jitTend;
  uint32_t itTstart;
  uint32_t itTend;

  if (verbose)
	{
	  if (qali.size() > 2){
	  cout << "----------------------------" << endl;
	  cout << "Options: (paired, uniq, best) = (" << paired << ", " << uniq << ", " << best << ")" << endl;
		cout << "Alignments passed for processing (qali.size())= " << qali.size() << endl;
}
	  cout << "----------------------------" << endl;
	}

  // Defines whether to treat reads as single- or pair-ended ones
  if (paired && !qali.empty())
	{
 
	  // cout << "The size of qali is: " << qali.size() << endl;

	  if (verbose && qali.size()>1)
		{
		  // Printing Qali before sorting
		  cout << "------------------------------------------------------" << endl;
		  cout << "qali BEFORE sorting by ref. Name and ref. start position:" << endl;
		  printQali(qali, pairwiseAlignments);
		}

	  if (qali.size() > 1)
		{
		  // Sorting by $tname and then by $tstart (check BAM format specification within SAM specs)

		  // sort by QueryName in ascending order
		  std::stable_sort(qali.begin(), qali.end(), [](const BamAlignmentRecord_ lhs, const BamAlignmentRecord_ rhs) {
		    std::less<std::string> comp;
		    return comp(lhs->getQueryName(), rhs->getQueryName());
		  });

		  //sort by StartPosition in ascending order
		  std::stable_sort(qali.begin(), qali.end(), [](const BamAlignmentRecord_ lhs, const BamAlignmentRecord_ rhs) {
		      // force unmapped aligmnents to end
		      if (lhs->getRefID() == -1) return false;
		      if (rhs->getRefID() == -1) return true;
		      // if on same reference, sort on position
		      if (lhs->getRefID() == rhs->getRefID()) {
		          return lhs->getStartPosition() < rhs->getStartPosition();
		      }
		      // otherwise sort on reference ID
		      return lhs->getRefID() < rhs->getRefID();
		  });
		}

	  if (verbose && qali.size()>1)
		{
		  // Printing Qali after sorting
		  cout << "------------------------------------------------------" << endl;
		  cout << "qali AFTER sorting by ref. Name and ref. start position:" << endl;
		  printQali(qali, pairwiseAlignments);
		  cout << "------------------------------------------------------" << endl;
		}

	  // Evaluating candidate mate-pairs.
	  // Sweeping through all possible pairings of alignments within Qali, looking 
	  // for candidates that belong to different: mate (qSuffix), strand (+,-) 
	  // and whose distance < maxInsertLength and insert length ...
	  for (it=0; it<qali.size()-1; it++)
		{
		  BamAlignmentRecord_& bar = qali.at(it);
		  itQname = bar->getQueryName(); 
		  if (!pairwiseAlignments) {
		  	itQsuffix = itQname.substr(itQname.find("/")+1, itQname.length()); 
		  } else {
				  	if (bar->isFirstMate()) {
							itQsuffix = "1";
					} else if (bar->isSecondMate()) {
							itQsuffix = "2";
					}
		  }
		  itStrand = bar->isReverseStrand();
		  itTstart = bar->getStartPosition();
		  itTend = bar->getEndPosition();

		  // Only loop until chromosome is different
		  for (jit=it+1;
		       jit < qali.size() && bar->getReferenceName() == qali.at(jit)->getReferenceName();
		       jit++)
			{
			  BamAlignmentRecord_& barj = qali.at(jit);
			  jitQname = barj->getQueryName();
			  if (!pairwiseAlignments) {
				jitQsuffix = jitQname.substr(jitQname.find("/")+1, jitQname.length()); } else {	
				  	if (barj->isFirstMate()) {
							jitQsuffix = "1";
					} else if (barj->isSecondMate()) {
							jitQsuffix = "2";
					}
			  }
			  jitStrand = barj->isReverseStrand();

			  if (verbose)
				{
				  cout << "comparing pair: " << it << "," << jit << ": [" 
						<< itQname << ", " << jitQname  << "]=[" << bar->getReferenceName() << ", mate " << itQsuffix 
					   << ", strand " << bool_cast(itStrand) << "; " << barj->getReferenceName() 
					   << ", mate " << jitQsuffix << ", strand " << bool_cast(jitStrand) << "]" << endl; 
				}

			  if (itQsuffix!=jitQsuffix) // different mates: (f,r) or (1,2)
  		  	  	{
  		  	  	  if (itStrand!=jitStrand) //different strands: (false, true)=(+,-)
  		  	  	  	{
		  			  jitTstart = barj->getStartPosition();
		  			  jitTend = barj->getEndPosition();
  		  	  		  dist = jitTstart - itTend - 1;
  		  	  		  if (itTstart>jitTstart) {dist = itTstart - jitTend - 1;}

					  if (verbose)
						{
						  cout << "Alignments " << it << "," << jit << " above point to common target, "  
							   << "with different mates, with different strands, dist=" << dist 
							   << ", maxIntronLen=" << maxIntronLen << endl;
						}

					  // If alignments not too far apart but not overlapping either, then they're mated
  		  			  if (dist < maxIntronLen) // Steffi commented this out: && (dist>=0 || (dist<0 && dist>(int32_t)(-1)*qali.at(it).Length)))
  		  			  	{

						  if ((verbose && (dist<0 && dist>(int32_t)(-1)*bar->getQuerySequenceLength())))
							{
							  cout << "dist=" << dist << " with a query length of " << bar->getQuerySequenceLength() << " indicates overlapping reads" << endl;
							}

  		  			  	  //push @matepairs, [$i,$j,scoreMate($i,$j,$dist)];
  		  			  	  mp.setValues(it, jit, scoreMate(bar, barj, dist, globalOptions));
						  matepairs.push_back(mp);

  		  			  	  if (!mated[it]) {mated[it]=0;}
  		  			  	  if (!mated[jit]) {mated[jit]=0;} 
  		  			  	  mated[it]++;
  		  			  	  mated[jit]++;
  		  			  	  inslen = jitTend - itTstart - 1;
  		  			  	  if (inslen<0) {inslen = itTend - jitTstart - 1;} 

						  if (verbose)
							{
							 cout << ">>>found mate pair:"<< it <<"," <<jit <<" (above alignment)!!!" << endl;
							  // cout << ": [" 
							  // 	   << itQname << ", " << jitQname  << "]=[" << bar->getReferenceName() << ", mate " 
							  // 	   <<itQsuffix << ", strand " << bool_cast(itStrand) << "; " << barj->getReferenceName() 
							  // 	   << ", mate " << jitQsuffix << ", strand " << bool_cast(jitStrand) 
							  // 	   << "]" << endl; 
							  cout << "mated[" << it << "]=" << mated[it] << ",mated[" << jit <<"]="
								   << mated[jit] << ", inslen=" << inslen << ",dist=" << dist << endl;
							}

  		  			  	  // push @insertlen, $inslen;
  		  			  	  insertlen.push_back(inslen); 

  		  			  	} else {
							if (verbose)
							  {
								cout <<"dist=" << dist << " not right between " << it << " and " 
									<< jit << endl;
							  }
  		  			  } // end if(dist < maxIntronLen && dist>=0)

  		  			} // end middle if
  		  		} // end outer if

			} // end inner for
		} // end middle for


	  if (verbose)
		{
		  cout << "Summary of evaluation of candidate mate-pairs" << endl;
		  cout << "found " << matepairs.size() << " mate pairs, involving "  
			   << mated.size() << " mates." << endl;
		  cout << "------------------------------------------------------" << endl;
		  printMatedPairsInfo(qali, matepairs);
		  printMatedMap(mated);
		  cout << "------------------------------------------------------" << endl;
		  cout << "(paired,uniq,best)=[" << paired << "," << uniq << "," << best << "]\n";

		}

	  // Counting all alignments that were not paired
	  optionalCounters.outPaired += (int)qali.size() - (int)mated.size();

	  if ((!uniq && !best) || matepairs.size()<2)
	  	{ // (uniq,best= (0,0) or size of mated pairs=1, let pass all mate-paired alignments
		  if (verbose)
			{
			  cout << "------------------------------------------------\n";
			  cout << "Letting pass all mated-paired alignments= " << mated.size() << ", listed below:\n";
			  cout << "Size of matepairs=" << matepairs.size() << "\n";
			  cout << "------------------------------------------------\n";
			}

		  if (matepairs.size()>0)
			{
			  map<int,int>::iterator m_it;
			  for (m_it=mated.begin(); m_it!=mated.end(); m_it++)
				{ 
				  if (verbose)
					{cout << "Letting pass paired-alignments: (" << (*m_it).first << ") : " 
					      << qali.at((*m_it).first)->getQueryName() << ", "
					      << qali.at((*m_it).first)->getReferenceName() << endl;
					}
				  writer->saveAlignment(qali.at((*m_it).first));
				}
			}

	  	} else { // (uniq or best) selected and matepair.size()>2

   			if (verbose)
			  {
				cout << "------------------------------------------------------" << endl;
				cout << "Sort (descending) mated pairs by scoreMate" << endl;
				cout << "BEFORE sorting [vector size]: " << matepairs.size() << endl;
				printMatePairs(matepairs, qali);
			  }

			// Sort matepairs by score (in descending order) 	
			std::stable_sort(matepairs.begin(), matepairs.end()); 

			if (verbose)
			  {
				cout << "AFTER sorting" << endl;
				printMatePairs(matepairs, qali);
		        cout << "------------------------------------------------------" << endl;
			  }  
 
			// Selection of Uniq-ue mate-pairs takes precedence over selecting Best mate-pairs
	  	  	if (uniq)
	  		  {// let pass only one mate pair, and only if second is significantly worse
			  	unsigned int second; float score_topPair, score_2ndPair, ratio;
				string mate_1_topPair, mate_1_2ndPair, mate_2_topPair, mate_2_2ndPair;
				second = 1;

				if (verbose)
				  {
					cout << "------------------------------------------------------" << endl;
					cout << "Comparing similarity between pairs:" << endl;
					cout << matepairs.at(0).alIt << " and " << matepairs.at(second).alIt << endl;
					cout << matepairs.at(0).alJit << " and " << matepairs.at(second).alJit << endl;
					cout << "------------------------------------------------------" << endl;
				  }

	  			while(second < matepairs.size() && 
					  similar(qali.at(matepairs.at(0).alIt),qali.at(matepairs.at(second).alIt),globalOptions) 
					  							&& 
					  similar(qali.at(matepairs.at(0).alJit),qali.at(matepairs.at(second).alJit),globalOptions))
	  			  {
	  				second++;
	  			  }
				if (verbose) {cout << "second=" << second << endl;}

				// "Second" mate-pair is the one with lowest score, but that is still similar
				if (second < matepairs.size())
				  {
					// Fetching data from top-scored mate-pair:
					mate_1_topPair = qali.at(matepairs.at(0).alIt)->getQueryName();
					mate_2_topPair = qali.at(matepairs.at(0).alJit)->getQueryName();
					score_topPair = matepairs.at(0).score;
					// Fetching data from second-worst mate-pair:
					mate_1_2ndPair = qali.at(matepairs.at(second).alIt)->getQueryName();
					mate_2_2ndPair = qali.at(matepairs.at(second).alJit)->getQueryName();
					score_2ndPair = matepairs.at(second).score;

					// Computing ratio between both sets of mate-pairs:
					ratio = score_2ndPair/score_topPair;

					if (verbose)
					  {
						cout << "Selecting a unique mate-pair in terms of its score" << endl;
						cout << "------------------------------------------------------------------------" 
							<< endl;
						cout << "Position of last similar mate-pair (indexed by second)=" << second << endl;
						cout << "Comparing scores between optimal and second mate-pairs: " << endl;
						cout << "[" << mate_1_topPair << " paired with " << mate_2_topPair << "; score=" 
							 << score_topPair << "]" << endl; 
						cout << "[" << mate_1_topPair << " paired with " << mate_2_2ndPair << "; score=" 
							 << score_2ndPair << "]" << endl;
						cout << "Ratio between these two mate-pairs: " << ratio << endl;
						cout << "------------------------------------------------------------------------\n"; 
					  }

					//... significantly worse in terms of "scoreMate"
					if (ratio < uniqThresh)
					  {
						if (verbose)
						  {
							cout << "Letting pass unique mate-pair: (" << mate_1_topPair << " and " 
								 << mate_2_topPair << ") because:" << endl;
							cout << "'lowest-scored' but similar mate-pair: (" << mate_1_2ndPair << " and " 
								 << mate_2_2ndPair << ")" << endl;
							cout << "is significantly worse, ratio=" << ratio << "<uniqThresh=" 
							 	<< uniqThresh << endl;
						  }
						writer->saveAlignment(qali.at(matepairs.at(0).alIt));
						writer->saveAlignment(qali.at(matepairs.at(0).alJit));
						optionalCounters.outUniq += mated.size()-2; // Drop all mated-pairs except the optimal two
					  } else {// dropping all mate-pairs belonging to this query

					  	  if (verbose)
						    {
							  cout << "(" << matepairs.size() <<") mate-pairs filtered out by uniqueness "; 
							  cout << "because:";
							  cout << "\nThey are similar and have score ratio=" << ratio 
								   << ">uniqThresh=" << uniqThresh << endl;
							  for (MatePairs matePair : matepairs)
								{
								  cout << qali.at(matePair.alIt)->getQueryName() << ", "
									   << qali.at(matePair.alJit)->getQueryName() << ") filtered out by ";
								  cout << "uniqueness criterion." << endl;	
								}	
							  cout << "Clearing contents of matepairs" << endl;				   	
							}

						  optionalCounters.outUniq += mated.size(); // Drop all mated alignments
						  matepairs.clear();
					  }

				  } else {// '(second == matepairs.size()' => all mate-pairs in "matepairs" are similar
					  if (verbose)
						{
						  cout << "------------------------------------------------------------------------" 
								<< endl; 
				  		  cout << "Suboptimal mate-pairs are all similar" << endl;
				  	 	  cout << "Letting pass only top-scored mate-pair: " 
							   << qali.at(matepairs.at(0).alIt)->getQueryName() << " and " 
							   << qali.at(matepairs.at(0).alJit)->getQueryName() << ", score=" 
							   << matepairs.at(0).score << endl;
						  cout << "(" << matepairs.size()-1<< ") mate-pairs filtered out by uniqueness" << endl;
						}

					  writer->saveAlignment(qali.at(matepairs.at(0).alIt));
					  writer->saveAlignment(qali.at(matepairs.at(0).alJit));
					  optionalCounters.outUniq += mated.size()-2; // Drop all alignments except the first pair of alignments
				}

				// keep only the best pair (if any)
				matepairs.resize(1);

	    	  } else { // (uniq, best) = (0, 1); let pass only "best" mate-pairs that share maximum score

			  	  if (verbose)
					{
					  cout << "------------------------------------------------------------------------" 
						   << endl;
					  cout << "Options: (paired, uniq, best) = (" << paired << ", " << uniq << ", " << best 
						   << ")" << endl;
					}

				  float optScore, tempScore; 
				  vector<string> bestTnames; 
				  optScore = matepairs.at(0).score;
				  tempScore = optScore;
				  unsigned int numBest = 0;

				  // Taking provisions for repeated indices
				  map<int,int> writtenIndices;
				  int mateIt, mateJit;

				  while (numBest < matepairs.size() && matepairs.at(numBest).score == optScore)
					{

					  // Verifying whether mate-pair(s) have been written before
					  mateIt = matepairs.at(numBest).alIt;
					  if (!writtenIndices.count(mateIt))
					  	{
						  if (verbose)
							{
							  cout << "Letting pass alignment: " << qali.at(matepairs.at(0).alIt)->getQueryName()
								   << ", i.e. pair: " << mateIt << ", score=" << tempScore 
								   << ", optScore=" << optScore << endl;
							}
						  writer->saveAlignment(qali.at(matepairs.at(numBest).alIt));
						  // Updating written indices with newly saved alignments
						  writtenIndices[mateIt] = 1;
						}

					  // Verifying whether mate-pair(s) have been written before
					  mateJit = matepairs.at(numBest).alJit;
					  if (!writtenIndices.count(mateJit))
					  	{
					  	  if (verbose)
					  		{
					  		  cout << "Letting pass alignment: " << qali.at(matepairs.at(0).alJit)->getQueryName()  
					  			   << ", i.e. pair: " << mateJit << ", score=" << tempScore 
					  			   << ", optScore=" << optScore << endl;
					  		}
					  	  writer->saveAlignment(qali.at(matepairs.at(numBest).alJit));
					  	  // Updating written indices with newly saved alignemnts
					  	  writtenIndices[mateJit] = 1;
					  	}

					  // Keeping name of reference for: geneFile
					  if (verbose)
						{
						  cout << "Save at bestTnames:" 
						       << qali.at(mateIt)->getReferenceName() << endl;
						}
					  bestTnames.push_back(qali.at(mateIt)->getReferenceName());
					  numBest++;

					}

				  optionalCounters.outBest += mated.size()-writtenIndices.size(); 
				  // Retaining only mate-pairs with optimal score
				  matepairs.resize(numBest);

				if (verbose)
				  {
					cout << "------------------------------------------------------------------------" << endl;
					cout << "Filtered out alignments (best criterion): " << qali.size() << endl;
					printQali(qali, pairwiseAlignments);
					cout << "------------------------------------------------------------------------" << endl;
					cout << "Summary [Common target names]: " << endl;
					cout << "bestTnames.size()=" << bestTnames.size() << endl;
				  }

				if (bestTnames.size() > 1 && commonGeneFile)
					{
					  std::set<std::string> geneNames;

					  for (std::string& tName : bestTnames)
						{
						// Replace suffixes of the type chrX.t123
						  tName = tName.substr(0, tName.find(".t")); 
						  geneNames.insert(tName);
						  if (verbose) {cout << "tName= " << tName << endl;}
						}
					
					if (geneNames.size() > 1)
						{
						  ofstream geneFile;
						  geneFile.open(commonGeneFile);
						  for (const std::string& tName : geneNames)
							{
							  if (verbose)
								{
								  cout << oldQnameStem << "\t"<< tName << endl;
								}
							  geneFile << oldQnameStem <<"\t"<< tName;
							}
						  geneFile << endl;
						  geneFile.close();  
						}
					} 

			} // end if (unique)

	  	} // (!uniq && !best) || matepairs.size()<2)

	  // Write pairedness coverage info into the pairBedFile	
	  // output pairedbed info: go through list of all mate pairs and store start and end position

	  if (strlen(pairBedFile)>0)
		{

		  if (verbose)
			{
			  cout << "---------------------Start: Pair bed section ------------------------\n";
			  cout << "The size of matepairs=" << matepairs.size() << "\n";
			  cout << "---------------------End: Pair bed section ------------------------\n";
			}

		  int pEnd, pStart;
		  string chr;
		  PairednessCoverage pc;

		  while (matepairs.size()>0)
		  	{
		  	  chr = qali.at(matepairs.at(0).alIt)->getReferenceName();
		  	  pEnd = qali.at(matepairs.at(0).alJit)->getEndPosition();
		  	  pStart = qali.at(matepairs.at(0).alIt)->getStartPosition();

			  // Storing info using map structures
			  pairCoord.insert(pair<int,int>(pStart, 1));
			  pairCovSteps2.insert(pair<string,multimap<int,int>>(chr, pairCoord));
			  pairCoord.insert(pair<int,int>(pEnd, -1));
			  pairCovSteps2.insert(pair<string,multimap<int,int>>(chr, pairCoord));

			  // Storing info using PairednessCoverage class
			  pc.setValues(pStart, 1, chr);
			  pairCovSteps.push_back(pc);
			  pc.setValues(pEnd, -1, chr);
			  pairCovSteps.push_back(pc);
			  matepairs.erase(matepairs.begin()); // Delete first member of matepairs
		  	}
		}

	} else {// IF NOT PAIRED, single alignments

		if ((uniq || best) && qali.size()>1)
		  {

			// Computing scores for each alignment, score=percId+coverage;
			for (BamAlignmentRecord_& bar : qali)
			{
			    bar->setScore(bar->getCoverage() + bar->getPercId());
			}

			if (verbose)
			  {
				cout << "Scoring alignments and sorting them according to such score." << endl;
				cout << "------------------------------------------------------------------------" << endl;
				cout << "qali BEFORE sorting by score:" << endl;
				printQali(qali, pairwiseAlignments);
				cout << "------------------------------------------------------------------------" << endl;
			  }

			// Sorting alignments by score in descending order
			std::stable_sort(qali.begin(), qali.end(), [](const BamAlignmentRecord_& lhs, const BamAlignmentRecord_& rhs) {
			    return lhs->getScore() > rhs->getScore(); // desc
			});

			if (verbose)
			  {
				cout << "------------------------------------------------------------------------" << endl;
				cout << "qali AFTER sorting by score:" << endl;
				printQali(qali, pairwiseAlignments);
				cout << "------------------------------------------------------------------------" << endl;
			  }

			BamAlignmentRecord_& firstAl = qali.at(0);
			// Selection of Uniq-ue alignments takes precedence over selecting the Best alignments
			if (uniq) 
			  { // let pass only one alignment, and only if second is significantly worse
			  	unsigned int second; float scoreFirst, scoreSecond, ratio;
				second = 1;

				// Sweep through all alignments and stop until a di-similar one is found
				while (second < qali.size() && similar(firstAl, qali.at(second), globalOptions))
				  {
					second++;
				  }


				// "Second" alignment is the one with lowest score, but that is still similar
				if (second < qali.size())
				  {
					BamAlignmentRecord_& secondAl = qali.at(second);
					// Comparing scores between best-scored alignment and one indexed by "second"
					
					scoreFirst = firstAl->getScore();
					scoreSecond = secondAl->getScore();
					// Computing ratio between first and second alignments
					ratio = scoreSecond/scoreFirst;

					if (verbose)
					  {
						cout << "Selecting a unique alignment in terms of its score" << endl;
						cout << "------------------------------------------------------------------------" 
							<< endl;
						cout << "Position of second worst alignment (index)=" << second << endl;
						cout << "Comparing scores between optimal and second: " << endl;
						cout << "[" << firstAl->getReferenceName() << "]<-"
							 << firstAl->getQueryName() << "; score=" << scoreFirst << " and " << endl;
						cout << "[" << secondAl->getReferenceName() << "]<-"
							<< secondAl->getQueryName() << "; score=" << scoreSecond << endl;
						cout << "Ratio between these two alignments: " << ratio << endl;
						cout << "------------------------------------------------------------------------\n"; 
					  }
	
					// ... significantly worse in terms of score
					if (ratio < uniqThresh) 
					  {		
						if (verbose)
						  {
							cout << "Letting pass unique alignment " << firstAl->getQueryName() << " because:" << endl;
							cout << "'lowest-scored' but similar alignment " << secondAl->getQueryName() << endl;
							cout << "is significantly worse, ratio=" << ratio << "<uniqThresh=" 
							 	<< uniqThresh << endl;
						  }
						writer->saveAlignment(firstAl);
						optionalCounters.outUniq += qali.size()-1;
						if (verbose)
						  {	cout << "(" << optionalCounters.outUniq << "/" << qali.size() << ") alignments filtered out " 
								 << " by uniqueness " << endl;
							for (unsigned int it=1; it<qali.size(); it++)
							  {cout << qali.at(it)->getQueryName() << " filtered out by uniqueness criterion." << endl;}
						  }
					  } else { // dropping all alignments belonging to the same query				  	

						if (verbose)
						  {
							cout << "(" << qali.size() << ") alignments filtered out by uniqueness because:";
							cout << "\nratio=" << ratio << ">uniqThresh=" << uniqThresh 
								 << " between top and lowest-scored similar alignments." << endl;
							for (unsigned int it=0; it<qali.size(); it++)
							  {cout << qali.at(it)->getQueryName() << " filtered out by uniqueness criterion." << endl;}
						  }					  	
					  	  optionalCounters.outUniq += qali.size();//Filtered out alignments by uniq. increases by size of Qali

					  }
				  } else {// '(second == qali.size()' => all alignments in "qali" are similar
				  	  if (verbose)
						{	
						  scoreFirst = firstAl->getScore();
						  cout << "------------------------------------------------------------------------" 
							<< endl;
				  	 	  cout << "Suboptimal alignments are all similar " << endl;
				  	 	  cout << "Letting pass only top-scored alignment: " << firstAl->getQueryName()
							   << ", score=" << scoreFirst << endl;
						  cout << "(" << qali.size()-1 << ") alignments filtered out by uniqueness" << endl;
						  for (unsigned int it=1; it<qali.size(); it++)
						 	{cout << qali.at(it)->getQueryName() << " filtered out by uniqueness criterion." << endl;}
						  cout << "------------------------------------------------------------------------" 
							<< endl;
						}
					  // Removing percId and coverage Tags before writing into file
					  writer->saveAlignment(firstAl);
				  	  optionalCounters.outUniq += qali.size()-1;
				  }

			  } else { // (uniq, best) = (0, 1); let pass only (best) alignments that share maximum score

			  if (verbose)
				{
				  cout << "------------------------------------------------------------------------" << endl;
				  cout << "Options: (paired, uniq, best) = (" << paired << ", " << uniq << ", " << best 
					<< ")" << endl;
				}

			  	float optScore = firstAl->getScore();
				vector<string> bestTnames; 

				while (qali.size() > 0 && qali.at(0)->getScore() == optScore)
				  {
					if (verbose)
					  {
						cout << "Letting pass alignment (best): " << firstAl->getQueryName() << ", score="
							 << qali.at(0)->getScore() << ", optScore=" << optScore << endl;
						cout << "Storing at bestTnames=" << firstAl->getReferenceName() << endl;
					  }

					// Removing percId and coverage Tags before writing into file
					writer->saveAlignment(qali.at(0));
					// Keeping name of reference for: geneFile
					bestTnames.push_back(qali.at(0)->getReferenceName());
					qali.erase(qali.begin()); // Deletes first member of qali
					if (verbose)
					  {cout << "After deleting alignment, qali.size()=" << qali.size() << endl;}
				  } 

				optionalCounters.outBest += qali.size();
				if (verbose)
				  {
					cout << "------------------------------------------------------------------------" << endl;
					cout << "Filtered out alignments (best criterion): " << qali.size() << endl;
					printQali(qali, pairwiseAlignments);
					cout << "------------------------------------------------------------------------" << endl;
					cout << "Summary [Common target names]: " << endl;
					cout << "bestTnames.size()=" << bestTnames.size() << endl;
				  }

				if (bestTnames.size() > 1 && commonGeneFile)
				  {
					std::set<std::string> geneNames;

					for (std::string& tName : bestTnames)
					  {
						// Replace suffixes of the type chrX.t123
						tName = tName.substr(0, tName.find(".t")); 
						geneNames.insert(tName);
						if (verbose) {cout << "tName= " << tName << endl;}
					  }
					
					if (verbose) {cout << "Size of geneNames=" << geneNames.size() << endl;}
					if (geneNames.size() > 1)
					  {
						ofstream geneFile;
						geneFile.open(commonGeneFile);
						
						// Writing name of query first..
						geneFile << oldQnameStem;

						// Then writing the common genes for such query
						for (const std::string& tName : geneNames)
						  {
							if (verbose) 
							  {
								cout << "commonGeneFile:" << oldQnameStem << "\t"<< tName<< endl;
							  }
							 geneFile << "\t"<< tName;
						  }
						geneFile << endl;
						geneFile.close();
					  }
				  } 

			  } // end if (unique)	

		} else { // (uniq, best) = (0, 0) && qali.size()>1)"

		  if (verbose)
			{
		  		cout << "Letting pass alignments, i.e. a total of: " << qali.size() << endl;
			}

			for(BamAlignmentRecord_& bar : qali)
				{
				  // Removing percId and coverage Tags before writing into file
				  writer->saveAlignment(bar); // Prints alignment line into file
				}

		} // end if !(uniq || paired)


  } // end if !(paired)


  // Clearing contents of qali
  qali.clear();

} // end processQuery()
