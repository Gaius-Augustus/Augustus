/*	list:sort
	
	NOTE: 
	This script only works with previously SORTED files!!!

	Created: 4-November-2011
	Last modified: 4-January-2012
*/ 
 
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <api/BamAlignment.h>
#include <api/algorithms/Sort.h> 
#include "bamtools_sort.h" 
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
#include <boost/lexical_cast.hpp>
#include <cstdio>
#include <algorithm>
#include <unordered_map>
#include "header.h"


using namespace BamTools;
using namespace BamTools::Algorithms; 
using namespace std;
using namespace boost;

string bool_cast(const bool b) 
{
    ostringstream ss;
    ss << boolalpha << b;
    return ss.str();
}

struct optionalCounters_t {
  int outPaired;
  int outUniq;
  int outBest;
};

//
// Generate a random number between 0 and 1
// return a uniform number in [0,1].
double unifRand()
{
    return rand() / double(RAND_MAX);
}

//
// Generate a random number in a real interval.
// param a one end point of the interval
// param b the other end of the interval
// return a inform rand numberin [a,b].
double unifRand(double a, double b)
{
    return (b-a)*unifRand() + a;
}

//
// Reset the random number generator with the system clock.
void seed()
{
    srand(time(0));
}


void sortQali(vector<BamAlignment>& qali);
void printQali(vector<BamAlignment> &qali, const RefVector refData);
float scoreMate(vector<BamAlignment> qali, int it, int jit, int dist, globalOptions_t globalOptions);
void printMatedPairsInfo(list<MatePairs> matepairs, list<int> insertlen);
void printMatedMap(map<int,int> mated);
optionalCounters_t processQuery(vector<BamAlignment> &qali, const RefVector refData, globalOptions_t globalOptions, BamWriter* ptrWriter);


int main(int argc, char *argv[])
{
  seed(); // For random number generation
  // Variable definition
  BamReader reader;
  BamWriter writer;
  BamAlignment al;
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
  const char* pairbedFile = globalOptions.pairbedFile;
  cout << "------------------------------------------------" << endl;
  cout << "Selected options are: " << endl;
  cout << "best=" << best << endl;
  cout << "help=" << help << endl;
  cout << "noIntrons=" << noIntrons << endl;
  cout << "paired=" << paired << endl;
  cout << "uniq=" << uniq << endl;
  cout << "verbose=" << verbose << endl;
  cout << "Input file is " << inputFile << endl;
  cout << "Output file is " << outputFile << endl;
  cout << "Value of paired is " << paired << endl;
  cout << "insertLimit=" << insertLimit << endl;
  cout << "maxIntronLen=" <<  maxIntronLen << endl;
  cout << "minCover=" <<  minCover << endl;
  cout << "minId=" <<  minId << endl;
  cout << "minIntronLen=" << minIntronLen << endl;
  cout << "uniqThresh=" << uniqThresh << endl;
  cout << "------------------------------------------------" << endl;

  // Starting timer
  tStart = time(NULL);    

  // Opening BAM file
  if (!reader.Open(inputFile)) 
	{
	cerr << "Could not open input BAM files. " << endl;
	return 0;
  	}

  // Retrieve 'metadata' from BAM files, these are required by BamWriter
  const RefVector refData = reader.GetReferenceData();
  const SamHeader header = reader.GetHeader();

  // Open BamWriter, HEADER and REFERENCE will be written into the file
  if (!writer.Open(outputFile, header, refData)) 
	{
	cerr << "Could not open output BAM file" << endl;
	return 0;
	}

  // Sweeping through alignments
  nextAlignment:
  	while (reader.GetNextAlignment(al)) 
	  {
		line++;
		qName = al.Name;
		qNameStem = qName; 
		RefID = al.RefID;
 		
		if (paired)
		  {	
			qNameStem = qName.substr(0, qName.find("/"));
			qSuffix = qName.substr(qName.find("/")+1, qName.length());	
		  } 

		// Displaying only data whose Reference seq ID has a defined mapping
		rName = getReferenceName(refData, RefID);
		if (rName.find("printReferenceName")!=-1)
		  {	  
			// cout << qName << " filtered out because it did not map with: " << rName << endl;
			outNoMapping++;
			goto nextAlignment;
		  }

	
		// Verifying file is sorted in ascending order
		if (oldQnameStem.compare("") && oldQnameStem.compare(qNameStem))   
		  { 
			// Checking whether 10th field is sorted in ascending order
			if (line <= maxSortesTest && qNameStems[qNameStem])   
			  {
				cerr << "Input file not sorted by query name!\n" << qNameStem.c_str() << 
				  " occurred previously. Set LC_ALL=C and sort -k 1,1" << endl; 	
				exit(1);
			  }	
			/////////////////////////////////////////////////////////////////	  	
			// 	Calling of: "processQuery() if (@qali);" goes here
			/////////////////////////////////////////////////////////////////
		  }  // end outer if

  		/////////////////////////////////////////////////////////////////
  		//
  		//     Filters by "$percId", "$coverage" and intron gaps "$nointrons" go here  
  		//     (see checkFileSortedness.pl)
  		//
  		/////////////////////////////////////////////////////////////////

  		// Fetching alignment information
  		cigar = al.CigarData;
  		cigarSize = cigar.size();
  		qName = al.Name; // query name
  		qLength = al.Length; // query length
		// if (qLength == 0) {qLength = 1;} //////////////////////////////////////////// TAKE OUT!!!!!!!!!!
  		RefID = al.RefID; // ID of reference seq. (later used)
  		sumMandI = 0; // Equiv to $qEnd-$qStart in PSL
  		baseInsert = 0;
  		al.GetTag("NM", editDistance); // edit distance (see SAM spec)

		// Percentage Identity filter
  		percId = 100*(qLength-editDistance)/qLength;  
		// percId = unifRand(75,100);
  		if (percId < minId)
  	  	{
  			outMinId++;
			cout << qName << " filtered out by percid=" << percId << " < minId=" << minId << endl;
  			goto nextAlignment;
  	  	}

  		// Coverage filter
  		sumMandI = sumMandIOperations(cigar, "no");
   		coverage = (float)100*sumMandI/qLength; 
		// coverage = unifRand(75,100);
  		// Save if complying with minCover
  		if (coverage < minCover)
  		  {	
  			outMinCover++;
			cout << qName << " filtered out by coverage= " << coverage << " < minId=" << minCover << endl;
  			goto nextAlignment;
  		  }	

  		// Intron gap filter
  		baseInsert = sumDandIOperations(cigar, "no");
  		// Save if complying with insertLimit
  		if (noIntrons && baseInsert > insertLimit)
  		  {
  			outIntrons++;
			cout << qName << " filtered out by intron criterion= " << baseInsert << " < minId=" << 
				insertLimit << endl;
  			goto nextAlignment;
  		  }


	
		cout << qName << " passed with parameters (percId, coverage)=(" << percId << "," << 
				coverage << ")" << endl; 

		// Appending coverage and percId into alignment
		std::stringstream field;
		field << percId;
		al.AddTag("pi", "Z", field.str());
 		field.str("");
		field << coverage;
		al.AddTag("co", "Z", field.str());
		field.str("");

  		/* Push @qali, [$_, $targetname, $qsuffix, $strand, $tstart, $tend, $percId, $coverage]; 
  		   Generating array of strings to do equivalent to: */
        qali.push_back(al);


		oldQnameStem = qNameStem;
		if (line <= maxSortesTest)
		  {
			qNameStems[qNameStem] = 1;
		  }

	  } ////// end while 


	// cout << "----------------------------" << endl;
	// cout << "qali.size()=" << qali.size() << endl;
	// cout << "----------------------------" << endl;

	// /* 	Calling of: "processQuery() if ($qnamestem ne "");" goes here */
	// try {
	  optionalCounters = processQuery(qali, refData, globalOptions, &writer);
	  outPaired = optionalCounters.outPaired;
	  outUniq = optionalCounters.outUniq;
	  outBest = optionalCounters.outBest;
	// } catch  (out_of_range& oor) {
	//   cout << "There was an error here: " << oor.what() << endl;
	//   cout << "Size of qali is: " << qali.size() << endl;
	// }

	// Displaying results to STDOUT
	cout <<  "\n  filtered out: " << endl;
	cout <<  "----------------: " << endl;
	cout <<  "unmapped        : " << outNoMapping << endl;
	cout <<  "percent identity: " << outMinId << endl;
	cout <<  "coverage        : " << outMinCover << endl;
	if (noIntrons)
	  {cout <<  "nointrons       : " << outIntrons << endl;}
	if (paired) 
	  {
		cout << "not paired      : " << outPaired << endl;
      	cout << "quantiles of unspliced insert lengths: " << endl;
	  	/// Print quantiles
 	  } 
	if (uniq) 
	  {cout << "unique          : " << outUniq << endl;}
	if (best)
	  {cout << "best            : " << outBest << endl;}



	// Displaying command line
	cout << "Cmd line: " << endl;
	for (int it=0; it<argc;it++)
	  {
		cout << argv[it] << " ";
	  }
	cout << endl;

	// Printing sortedness test results
	cout << "----------------------------------------" << endl;
	cout << "File is correctly sorted " << endl;
	cout << "----------------------------------------" << endl;


	// Closing file handles
	reader.Close();
	writer.Close();

	// Ending timer and printing elapsed time
	tEnd = time(NULL);    
	tElapsed = printElapsedTime(tEnd, tStart);     

} // end main



void sortQali(vector<BamAlignment>& qali) 
{
        sort( qali.begin(), qali.end(), Sort::ByName() );
        sort( qali.begin(), qali.end(), Sort::ByPosition() );
}

void printQali(vector<BamAlignment> &qali, const RefVector refData)
{
  vector<BamAlignment>::iterator it = qali.begin();
  std::stringstream ss_rstart, ss_rend, ss_percId, ss_coverage, ss_score;
  string rName;
  string qName;
  string qSuffix;
  bool strand;
  string percId, coverage, score;

  for (it; it != qali.end(); it++)
	{
	  rName = getReferenceName(refData, (*it).RefID); 
	  qName = (*it).Name; 
	  qSuffix = qName.substr(qName.find("/")+1, qName.length());	
	  strand = (*it).IsReverseStrand();
	  ss_rstart << (*it).Position; 
	  ss_rend << (*it).GetEndPosition();
  	  (*it).GetTag("pi", percId);
  	  (*it).GetTag("co", coverage);
  	  (*it).GetTag("sc", score);
 	  ss_percId << percId;  
	  ss_coverage << coverage;  
	  ss_score << score;  

	  // Printing to stdout
	  cout << "\t[ ";
	  cout << qName << " " << rName << " " << qSuffix << " " << bool_cast(strand) << " " << 
			ss_rstart.str() << " " << ss_rend.str() << " " << ss_percId.str() << " " << 
		    ss_coverage.str() << " " << ss_score.str(); 
	  cout << " ]," << endl;
	  // clear strings
	  ss_rstart.str(""); ss_rend.str(""); ss_percId.str(""); ss_coverage.str(""); ss_score.str("");
	}
}

// parameter later used for comparing quality between two alignments
vector<BamAlignment> scoreAli(vector<BamAlignment>& qali)
{
  string s_percId, s_coverage;
  float percId, coverage, score;
  std::stringstream ss_score;
  vector<BamAlignment>::iterator it = qali.begin();

  for (it; it!=qali.end(); it++)
	{
	  	(*it).GetTag("pi", s_percId);
  		(*it).GetTag("co", s_coverage);
 		ss_score << atof(s_coverage.c_str())+atof(s_percId.c_str()); 
		// cout << "Adding score=" << ss_score.str() << " to " << (*it).Name << endl;
		(*it).AddTag("sc", "Z", ss_score.str());
		ss_score.str("");
	}

  return qali;

}

// for comparing quality of two mate-pair read alignments (it, jit)
float scoreMate(vector<BamAlignment> qali, int it, int jit, int dist, globalOptions_t globalOptions)
{
  int maxIntronLen = globalOptions.maxIntronLen;
  bool best = globalOptions.best;

  string percId1, coverage1, percId2, coverage2;
  string rf1, rf2;

  qali.at(it).GetTag("pi", percId1);
  qali.at(jit).GetTag("pi", percId2);
  qali.at(it).GetTag("co", coverage1);
  qali.at(jit).GetTag("co", coverage2);

  float percIdIt = atof(percId1.c_str());
  float percIdJit = atof(percId2.c_str());
  float coverageIt = atof(coverage1.c_str());
  float coverageJit = atof(coverage2.c_str());

  float score = (percIdIt + percIdJit)/100 // percent identity
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
bool similar(BamAlignment al1, BamAlignment al2)
{
  uint32_t jitTstart = al1.Position;
  uint32_t jitTend = al2.GetEndPosition();  
  uint32_t itTstart = al1.Position; 
  uint32_t itTend = al2.GetEndPosition();

  if (jitTend <= itTstart || jitTstart >= itTend) // here: similar = overlapping target range
	{
	  cout << "they are not similar" << endl;
	  return false;
	} else {
			cout << "alignments " << al1.Name << " and " << al2.Name << " are similar" << endl;
			return true;
  		   }
}


void printMatedPairsInfo(list<MatePairs> matepairs, list<int> insertlen)
{
  list<MatePairs>::iterator it = matepairs.begin();
  list<int>::iterator itInsertlen = insertlen.begin();
  cout << "i j scoreMate insertLen" << endl;
  for (it, itInsertlen; it != matepairs.end(); it++, itInsertlen++)
	{
	  // cout << *it << "," << *itInsertlen;
	  cout << *it << " " << *itInsertlen << endl;
	}
}


void printMatedMap(map<int,int> mated)
{
  map<int,int>::iterator itMated = mated.begin();
  cout << "key=>value" << endl;
  for (itMated; itMated!=mated.end(); itMated++)
	{
	  cout << (*itMated).first << "=>" << (*itMated).second << endl;
	}
}


optionalCounters_t processQuery(vector<BamAlignment> &qali, const RefVector refData, globalOptions_t globalOptions, BamWriter* ptrWriter)
{
  struct optionalCounters_t optionalCounters;
  // Optional counters
  int outPaired = 0;
  int outUniq = 0;
  int outBest = 0;

  // Expanding options
  bool best = globalOptions.best;
  bool paired = globalOptions.paired;
  bool uniq = globalOptions.uniq;
  bool verbose = globalOptions.verbose;
  int insertLimit = globalOptions.insertLimit;
  int maxIntronLen = globalOptions.maxIntronLen;
  int maxSortesTest = globalOptions.maxSortesTest;
  int minCover = globalOptions.minCover;
  int minId = globalOptions.minId;
  int minIntronLen = globalOptions.minIntronLen;
  int uniqThresh = globalOptions.uniqThresh;

  // Mate variables
  string itRname, jitRname;
  string itQname, jitQname;
  string itQsuffix, jitQsuffix;
  bool itStrand, jitStrand;
  // Matepairs
  list<MatePairs> matepairs;
  MatePairs mp;
  list<int> insertlen;
  map<int, int> mated;
  int32_t inslen, dist;
  // unsigned int inslen, dist;
  float scoremates, correction;
  int it, jit;
  float scoreArray[SIZE];
  uint32_t jitTstart;
  uint32_t jitTend;
  uint32_t itTstart;
  uint32_t itTend;
  string reservedField;

  if (paired)
	{
	  // //// Printing Qali before and after sorting
	  // cout << "------------------------------------------------------" << endl;
	  // cout << "qali BEFORE sorting by Name and position:" << endl;
	  // printQali(qali, refData);
	  sortQali(qali); 	  // Sorting by $tName and then by $tstart 
	  // cout << "------------------------------------------------------" << endl;
	  // cout << "qali AFTER sorting by Name and position:" << endl;
	  // printQali(qali, refData);


	  for (it=0; it<qali.size()-1; it++)
		{
		  itRname = getReferenceName(refData, qali.at(it).RefID); 
		  itQname = qali.at(it).Name; 
		  itQsuffix = itQname.substr(itQname.find("/")+1, itQname.length());	
		  itStrand = qali.at(it).IsReverseStrand();

		  for (jit=it+1; jit<qali.size() && getReferenceName(refData, qali.at(it).RefID)==getReferenceName(refData, qali.at(jit).RefID); jit++) 
			{
			  jitRname = getReferenceName(refData, qali.at(jit).RefID);
			  jitQname = qali.at(jit).Name;
			  jitQsuffix = jitQname.substr(jitQname.find("/")+1, jitQname.length());	
			  jitStrand = qali.at(jit).IsReverseStrand();
			  // cout << "comparing " << it << "," <<  getReferenceName(refData, qali.at(it).RefID) << 
			  // 		", with " << jit << "," << getReferenceName(refData, qali.at(jit).RefID) << endl;
			  if (itQsuffix!=jitQsuffix) //qSuffix.qali.at(it)==qSuffix.qali.at(jit)
  		  	  	{
  		  	  	  if (itStrand!=jitStrand) //strand.qali.at(it)==strand.qali.at(jit)
  		  	  	  	{
		  			  jitTstart = qali.at(jit).Position; 
		  			  jitTend = qali.at(jit).GetEndPosition(); ///////////+1
		  			  itTstart = qali.at(it).Position; 
    	  			  itTend = qali.at(it).GetEndPosition(); //////////+1
  		  	  		  dist = jitTstart - itTend - 1;
  		  	  		  if (itTstart>jitTstart)
  		  	  		  	{
  		  	  		  	  dist = itTstart - jitTend - 1;
  		  	  		  	}
  		  			  if (dist < maxIntronLen && dist>=0)
  		  			  	{
  		  			  	  // printf("found mate pair %d,%d\n", it, jit);
  		  			  	  ///////////////////////////////////////////////////////////////////////////////
  		  			  	  //push @matepairs, [$i,$j,scoreMate($i,$j,$dist)];
  		  			  	  scoreArray = {(float)it, (float)jit, scoreMate(qali, it, jit, dist, globalOptions)};
  		  			  	  mp.push(scoreArray); matepairs.push_back(mp);
  		  			  	  ///////////////////////////////////////////////////////////////////////////////
  		  			  	  if (!mated[it]) {mated[it]=0;}
  		  			  	  if (!mated[jit]) {mated[jit]=0;} 
  		  			  	  mated[it]++;
  		  			  	  mated[jit]++;
  		  			  	  inslen = jitTend - itTstart - 1;
  		  			  	  if (inslen<0)
  		  			  	  	{
  		  			  	  	  inslen = itTend - jitTstart - 1;
  		  			  	  	} 
  		  			  	  scoremates = scoreMate(qali, it, jit, dist, globalOptions);
  		  			  	  correction = dist/maxIntronLen/10;
  	  	  			  	  // printf("mated[%d]=%d,mated[%d]=%d,inslen=%I32d,dist=%I32df,score=%2.5,
  		  			  	  // 		  correction=%2.5f\n", it, mated[it], jit, mated[jit], inslen, dist, 
  		  			  	  // 		  scoremates, correction);
  		  			  	  printf("mated[%d]=%d,mated[%d]=%d,inslen=%I32d,dist=%I32d\n", it, mated[it], 
								 jit, mated[jit], inslen, dist);
  		  			  	  ///////////////////////////////////////////////////////////////////////////////
  		  			  	  // push @insertlen, $inslen;
  		  			  	  insertlen.push_back(inslen); 
  		  			  	  ///////////////////////////////////////////////////////////////////////////////
  		  			  	} else {
  		  					// printf("dist=%I32d not right between %d and %d\n", dist, it, jit);
  		  			  } // end if(dist < maxIntronLen && dist>=0)

  		  			} // end middle if
  		  		} // end outer if

			} // end inner for
		} // end middle for

 
	  printf("found %d mate pairs, involving %d mates.\n", (int)matepairs.size(), (int)mated.size());

	  /////////////////////////////////////////////
	  // Printing by-products of computeMatePairs()
	  printMatedPairsInfo(matepairs, insertlen);
	  printMatedMap(mated);
	  // printQali(qali, refData);
	  /////////////////////////////////////////////


	  //////////////////////////////
	  // Insert commented code here
	  // $outPaired += @qali - scalar(keys %mated);
	  int outPaired=0;
	  outPaired += (int)qali.size() - (int)mated.size();
	  cout << "scalarQali=" << (int)qali.size() << endl;
	  cout << "scalarMated=" << (int)mated.size() << endl;
	  cout << "outPaired=" << (int)outPaired << endl;
	  if (!uniq && !best || matepairs.size()<2)
		{ // let pass all read alignments that are involved in mate pairs
		  map<int,int>::iterator m_it;
		  for (m_it=mated.begin(); m_it!=mated.end(); m_it++)
			{ 
			  // try {

				(*ptrWriter).SaveAlignment(qali.at((*m_it).first)); // Saving alignment
			  // } catch  (out_of_range& oor) {
			  // 	cout << "There was an error at line " << line << ", " << oor.what() << endl;
			  // }

			}
		} else {
			// Sort matepairs by score, i.e. @matepairs = sort {$b->[2] <=> $a->[2]} @matepairs;
			matepairs.sort(); 
			printList(matepairs);
	  }
	  // 	  if (uniq)
	  // 		{// let pass only best mate pair, and only if second is 
	  // 		  // significantly worse
	  // int second = 1;
	  // list <MatePairs>::iterator matesIter = matepairs.begin();
	  // while(second < matepairs.size()) // && similar())
	  // 	{
	  // 	  second++;
	  // 	  matesIter++;
	  // 	}
	  // cout << "second=" << second << endl;

	  // 		  if (second < (int)matepairs.size())
	  // 			  {
	  // matesIter = matepairs.begin();
	  // float ratio = (*matesIter).values[2]/((*matepairs.begin()).values[2]);
	  // cout << "ratio=" << ratio << endl;
	  // if (verbose)
	  //   {
	  // printf("\nbest two mates\n");
	  // matesIter = matepairs.begin();
	  // cout << qali.at((int)(*matesIter).values[0]) << "paired with\n" << 
	  //   		qali.at((int)(*matesIter).values[1]) << "score=" << 
	  //   		(*matesIter).values[2] << "\n";
	  // matesIter++;
	  // cout << qali.at((int)(*matesIter).values[0]) << "paired with\n" << 
	  //   		qali.at((int)(*matesIter).values[1]) << "score=" << 
	  //   		(*matesIter).values[2] << "\n";
	  // } // end if(verbose)

	  // 				if (ratio < uniqthresh)
	  // 				  {
	  // // print the two alignments for best mate pair only
	  // printf("%d", matepairs); // $qali[$matepairs[0]->[0]]->[0]);
	  // printf("%d", matepairs); // $qali[$matepairs[0]->[1]]->[0]);
	  // 				  } else {
	  // 				  		matepairs.clear();
	  // 				  }
	  // 			  } else {
	  // 			  if (verbose)
	  // 				{
	  // 				  printf( "suboptimal mate pairs are similar\n");
	  // 				}
	  // 			  // printf($qali[$matepairs[0]->[0]]->[0]);
	  // 			  // printf($qali[$matepairs[0]->[1]]->[0]);			
	  // 			  }// end if(second<@matepairs.size())
	  // 			///// splice @matepairs, 1; # keep only the best pair (if any)

	  // 		} else { // best: take all best alignment pairs
	  // 		int optscore;
	  // 		list<string> bestTnames;
	  // 		int numbest = 0;
	  // 		while (numbest < matepairs.size()) 
	  // 		  {
	  // 			// printf($qali[$matepairs[$numbest]->[0]]->[0]);
	  // 			// printf($qali[$matepairs[$numbest]->[1]]->[0]);
	  // 			// push @bestTnames, $qali[$matepairs[$numbest]->[0]]->[1];
	  // 			numbest++;
	  // 		  }
	  // 		// $outBest += @matepairs - $numbest;
	  // 		// splice @matepairs, $numbest; # keep only the first $numbest pairs
	  // 		if (bestTnames.size()>1)
	  // 		  {
	  // 			map<string,int> genenames; 
	  // 			int itTname = 0; 
	  //           for (itTname; itTname<genenames.size(); itTname++) //my $Tname (@bestTnames) 
	  // 				{ 
	  // 				  // if (genenames[itTname]==1)
	  // 				  // 	{
	  // 				  // 	//$Tname =~ s/\.t\d+//; $genenames{$Tname}=1;
	  // 				  // 	}
	  // 				}
	  // 			// print COMMON $oldqnamestem . "\t" . join(" ", keys %genenames) . "\n" if (%genenames > 1 && defined($commongenefile));

	  // 		  }		
	  // 	  } // end if(uniq)

	  // } // end if (!uniq && !best || matepairs.size()<2)
	
	  // // output pairedbed info: go through list of all mate pairs and store start and end position
	  // if (pairbedfile)
	  // 	{
	  // 	  while (matepairs.size() > 0)
	  // 		{
	  // 	  	// DO SOMETHING ELSE
		
	  // 		}
	  // 	}
	  //////////////////////////////

	} else {// IF NOT PAIRED, single read
		if ((uniq || best) && qali.size()>1)
		  {
			// cout << "------------qali BEFORE sorting by score:---------------------------" << endl;
			// printQali(qali, refData);
			// Computing scores for each alignment and sorting them in descending order
			qali = scoreAli(qali);
			std::stable_sort( qali.begin(), qali.end(), Sort::ByTag<std::string>("sc", Sort::DescendingOrder) );
			// cout << "------------qali AFTER sorting by score:----------------------------" << endl;
			// printQali(qali, refData);
			if (uniq)
			  {
				int second; int score0, score2; float ratio;
				second = 1;
				// Sweep through all alignments and stop until a di-similar one is found
				while (second < qali.size() && similar(qali.at(0), qali.at(second)))
				  {
					second++;
				  }
				cout << "The position of second=" << second << endl;
				// The $second alignment is the one with lowest $score, but that is still similar
				if (second < qali.size())
				  {
					qali.at(0).GetTag("sc", score0);
					qali.at(second-1).GetTag("sc", score2);
					cout << "score0 = " << score0 << "; score2 = " << score2 << endl;
					ratio = (float)score2/score0;
					// if (verbose)
					//   {
					cout << "best two alignments " << endl;
					cout << getReferenceName(refData, qali.at(0).RefID) << "\nscore=" << score0 << endl;
					cout << getReferenceName(refData, qali.at(0).RefID) << "\nscore=" << score2 << endl;
					cout << "ratio = " << ratio << endl;
					// }
	
					if (ratio < uniqThresh) // ... significantly worse in terms of "scoreAli"
					  {		
						cout << "Letting pass only one alignment alignment: " << qali.at(0).Name << 
						  " because second alignment " << qali.at(second).Name << 
						  " is significantly worse." << endl;
						(*ptrWriter).SaveAlignment(qali.at(0)); // Prints alignment line into file
						outUniq += qali.size()-1;
					  } else {
					  outUniq += qali.size(); // dropping all alignments belonging to the same query
					  cout << "All queries (" << qali.size() << " alignments): " << 
						qali.at(0).Name << " filtered out by uniqueness ratio=" << 
						ratio << "> uniqThresh=" << uniqThresh << endl;
					}
				  } else {// Implies: (second == qali.size()) => all alginments in "qali" are similar	
				  cout << "Suboptimal alignments are all similar " << endl;
				  cout << "Letting pass only one alignment: " << qali.at(0).Name << endl;
				  (*ptrWriter).SaveAlignment(qali.at(0)); // Letting pass only best
				  outUniq += qali.size()-1;
				}

			} else { // if !(unique) but (best)], take all best alignments that share maximum score

			    cout << "----------------------------" << endl;
				cout << "qali.size()=" << qali.size() << endl;
				cout << "----------------------------" << endl;
				string s_optScore, s_tempScore;
				float optScore, tempScore; 
				qali.at(0).GetTag("sc", s_optScore);
				optScore = atof(s_optScore.c_str());
				tempScore = optScore;
				vector<string> bestTnames;
				while(tempScore == optScore && qali.size()>0)
				  {
					(*ptrWriter).SaveAlignment(qali.at(0)); 
				  // cout << "Wrote to file " << qali.at(0).Name << ", with tempScore=" << tempScore << endl;
					bestTnames.push_back(getReferenceName(refData, qali.at(0).RefID));
					qali.erase(qali.begin()); // Delete first member of qali
					qali.at(0).GetTag("sc", s_tempScore);
					tempScore = atof(s_tempScore.c_str());
				  }

				outBest += qali.size();
				cout << "-------------------------------------" << endl;
				cout << "Filtered out alignments by best criterion: " << qali.size() << endl;
				// printQali(qali, refData);		
				cout << "-------------------------------------" << endl;
				if (bestTnames.size()>1)
				  {
					string geneNames;
					for (int it=0; it<bestTnames.size(); it++)
					  {
						geneNames.append(bestTnames.at(it));
						geneNames.append(" ");
					  }
					cout << "Records for commongenefile: " << qali.at(0).Name <<  " " << 
					  bestTnames.size() << " genes " << endl;
				  }

			  } // end if (unique)
	
		} else { //This else corresponds to: "if !((uniq || best) && qali.size()>1)"
		  	cout << "Letting pass all alignments, that is a total of: " << qali.size() << endl;
			for(int it=0; it<qali.size(); it++)
				{
				  (*ptrWriter).SaveAlignment(qali.at(it)); // Prints alignment line into file
				}

		} // end if !(uniq || paired)

	// Clearing contents of qali
	qali.clear();


  } // end if !(paired)


  // Updating filter counters	
  optionalCounters.outPaired = outPaired;
  optionalCounters.outUniq = outUniq;
  optionalCounters.outBest = outBest;

  return optionalCounters;
} // end processQuery()


