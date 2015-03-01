/*****************************************************************************\
 * Filename : extrinsicinfo.hh
 * Author   : Mario Stanke
 *
 *
 * Description: EST information or homology information or information from inter species
 *              comparisons. Hints in general.
 *
 * Date       |   Author              |  Changes
 *------------|-----------------------|---------------------------------
 * 26.11.02   | Mario Stanke          | creation of the file
 * 27.08.03   | Mario Stanke          | grades of hints
 * 23.06.06   | Mario Stanke          | introduce CDS, CDSpart, UTR, UTRpart hints
 * 26.07.06   | Mario Stanke          | introduce nonexonpart hints
 * 28.08.06   | Mario Stanke          | HintGroup
 * 01.07.08   | Mario Stanke          | added new hint type nonirpartF
 \******************************************************************************/

#ifndef __EXTRINSICINFO_HH
#define __EXTRINSICINFO_HH

// project includes
#include "hints.hh"
#include "properties.hh"
#include "gene.hh"

// standard C/C++ includes
#include <set>
#include <fstream>
#include <vector>
#include <iostream>


#define forwDSS    0
#define revDSS     1
#define forwASS    2
#define revASS     3
#define maxStoreMalus 3000

/*
 * Range of prediction and a set of HintGroups turned off
 */
class PredictionRun {
public:
    PredictionRun(){
	begin = end = -1;
	omittedGroups = NULL;
	allHints = false;
    } 
    PredictionRun(int begin, int end, list<HintGroup*> *omittedGroups, bool allHints = false){
	this->begin = begin;
	this->end = end;
	this->omittedGroups = omittedGroups;
	this->allHints = allHints;
    }
    void print(int beginPos = -1);
    int begin;
    int end;
    bool allHints; // this is true for the (first) prediction run that includes all hints, this is not quite the same as omittedGroups=NULL
    list<HintGroup*> *omittedGroups;
};

/*
 * Plan of individiual prediction steps with each a range and a set of hint groups turned on/off
 */
class PredictionScheme {
public:
    PredictionScheme(int n){
	this->n = n;
    }
    void addRun(PredictionRun run);
    void print(int beginPos = -1);
    int n;
    list<PredictionRun> predictionRuns; 
};


class FeatureCollection;

/*
 * SequenceFeatureCollection
 * holds all extrinsic feature information for one sequence
 * 
 */
class SequenceFeatureCollection {
public:
    SequenceFeatureCollection(FeatureCollection *collection) {
	this->collection = collection;
	featureLists = new list<Feature>[NUM_FEATURE_TYPES];
	groupList = NULL;
	groupGaps = NULL;
	predictionScheme = NULL;
	sorted = true;
	groupsSorted = false;
	hintedSites = NULL;
	hasLocalSSmalus = NULL;
	seqlen = 0;
	K = BLOCKSIZE;
	firstEnd = NULL;
	lastStart = NULL;
    }
    SequenceFeatureCollection(SequenceFeatureCollection& other, int from, int to, bool rc = false);
    ~SequenceFeatureCollection(){
	if (featureLists)
	    delete[] featureLists;
	if (predictionScheme)
	    delete predictionScheme;
	if (groupList)
	    delete groupList;
	if (groupGaps)
	    delete groupGaps;
	if (hintedSites)
	    delete [] hintedSites;
	if (hasLocalSSmalus)
	  delete [] hasLocalSSmalus;
	if (firstEnd) {
	    for (int type=0; type < NUM_FEATURE_TYPES; type++)
		delete [] firstEnd[type];
	    delete [] firstEnd;
	}
	if (lastStart) {
	    for (int type=0; type < NUM_FEATURE_TYPES; type++)
		delete [] lastStart[type];
	    delete [] lastStart;
	}
	
    }

    void addFeature(Feature f);
    void printFeatures(ostream& out);
    void sortFeatureLists();
    void checkGroupConsistency(AnnoSequence *seq);
    void warnInconsistentHints(AnnoSequence *seq);
    static void deleteEqualElements(list<Feature> &list);
    list<Feature> getFeatureList(FeatureType type) {
	return featureLists[(int) type];
    }
    Feature *getFeatureAt(FeatureType type, int endPosition, Strand strand);
    Feature *getFeatureListAt(FeatureType type, int endPosition, Strand strand);
    Feature *getAllActiveFeatures(FeatureType type);
    Feature *getFeatureListInRange(FeatureType type, int startPosition, int endPosition, 
				   Strand strand, int seqRelFrame=-1);
    Feature *getFeatureListBeginningInRange(FeatureType type, int startPosition, int endPosition, 
					    Strand strand, int seqRelFrame=-1);
    Feature *getFeatureListOvlpingRange(FeatureType type, int startPosition, int endPosition, Strand strand);
    Feature *getFeatureListOvlpingRange(Bitmask featuretypes, int startPosition, int endPosition, Strand strand);
    Feature *getFeatureListContaining(Bitmask featuretypes, int position, Strand strand);
    Feature *getExonListInRange(int startPosition, int endPosition, Strand strand, int seqRelFrame=-1);
    Feature *getExonListOvlpingRange(int startPosition, int endPosition, Strand strand, int seqRelFrame=-1);
    double localSSMalus(FeatureType type, int pos, Strand strand);
    void shift(int offset);
    void computeHintedSites(const char* dna);
    bool validDSSPattern(const char* dna) const {
	try {
	    return validDSS[Seq2Int(2)(dna)];
	} catch (...) {
	    return false;
	}
    }
    bool validRDSSPattern(const char* dna) const {
	try {
	    return validDSS[Seq2Int(2).rc(dna)];
	} catch (...) {
	    return false;
	}
    }
    bool validASSPattern(const char* dna) const {
	try {
	    return validASS[Seq2Int(2)(dna)];
	} catch (...) {
	    return false;
	}
    }
    bool validRASSPattern(const char* dna) const {
	try {
	    return validASS[Seq2Int(2).rc(dna)];
	} catch (...) {
	    return false;
	}
    }
    bool validSplicePattern(string s) const {
	return validHintedSites.count(s) > 0;
    }
    bool validRSplicePattern(string s) const {
	string reverseComplement(4,0);
	putReverseComplement(reverseComplement.begin(), s.c_str(), 4);
	return validSplicePattern(reverseComplement);
    }
    bool isHintedDSS(int pos, Strand strand) const {
	return hintedSites && 
	    hintedSites[pos][strand == plusstrand ? forwDSS : revDSS];
    }
    bool isHintedASS(int pos, Strand strand) const {
	return hintedSites &&
	    hintedSites[pos][strand == plusstrand ? forwASS : revASS];
    }
    void deleteFeatureAt(FeatureType type, int endPosition, Strand strand);
    void cleanRedundantFeatures();
    void setSeqLen(int len){seqlen = len;}
    void makeGroups();
    void printGroups();
    void findGenicGaps(); // find intra-group gaps: nonirpart hints
    void findGroupGaps(); // find inter-group gaps: good positions to break up sequence
    void determineInterGroupRelations();
    void resetConformance();
    void emptyTrash();
    void computeIndices();
    list<Feature>::iterator getPosFirstEndAtOrAfter(int type, int e);
    list<Feature>::iterator getPosStartAfter(int type, int s);
    void rescaleBoniByConformance();
    void createPredictionScheme(list<AltGene> *genes);
    void setActiveFlag(list<HintGroup*> *groups, bool flag);
    list<AltGene> *joinGenesFromPredRuns(list<list<AltGene> *> *genesOfRuns, int maxtracks, bool uniqueCDS);
    void sortIncompGroupsOfGroups(); // sort by pointers, for '==' operator of lists
    set<HintGroup*> *getCausingGroups(PredictionRun &pr);
    void prepare(AnnoSequence *annoseq, bool print, bool withEvidence=true);
    void prepareLocalMalus(const char* dna); // precompute tables for computing a local exonpart/UTRpart/CDSpart malus later
    int numZeroCov(int start, int end, int type, Strand strand);
    static void initHintedSplicesites(string ssList);
    FeatureCollection *collection;

    static Bitmask     validDSS;
    static Bitmask     validASS;
    static set<string> validHintedSites;
 
    list<Feature>    *featureLists;
    list<HintGroup>  *groupList;
    list<Feature>    *groupGaps; // sorted list of gaps between groups (irpart)
    bool             groupsSorted;
    PredictionScheme *predictionScheme;
    bool             sorted;
    Bitmask          *hintedSites;
    Bitmask          *hasLocalSSmalus;
    int              seqlen;
    int              K; // block size for firstEnd and lastStart
    // firstEnd[t][k] holds an iterator to featureList[t] to the first element f in the list such that
    // f->end >= k*K, where k>=0 and k ist such that k*K <= seqlen
    list<Feature>::iterator **firstEnd;
    // lastStart[t][k] holds an iterator to featureList[t] to the list such that all following list elements f have
    // f->start > k*K, where k>=0 and k ist such that k*K <= seqlen 
    list<Feature>::iterator **lastStart;
    vector<int> cumCovUTRpartPlus; // cumulative number of positions not covered by UTRpart hints on the plus strand
    vector<int> cumCovUTRpartMinus; // cumulative number of positions not covered by UTRpart hints on the minus strand
    vector<int> cumCovCDSpartPlus; // cumulative number of positions not covered by CDSpart hints on the plus strand
    vector<int> cumCovCDSpartMinus; // cumulative number of positions not covered by CDSpart hints on the minus strand
    vector<int> cumCovExonpartPlus; // cumulative number of positions not covered by exonpart hints on the plus strand
    vector<int> cumCovExonpartMinus; // cumulative number of positions not covered by exonpart hints on the minus strand

private:
    void addCumCov(vector<bool> &cov, const list<Feature>& flist, Strand strand);
};

/*
 * FeatureTypeInfo
 * holds all type-specific extrinsic feature information
 *
 */
struct FeatureTypeInfo {
    FeatureTypeInfo(int numSources=1, double b=-1.0, double m=1.0, double lm=1.0) :
	bonus(b), // -1 means not initialised
	malus(m), // don't change anything unless we have at least searched 
	localMalus(lm),
	gradeclassbounds(numSources, vector<double>()),
	gradequots(numSources, vector<double>(1, 1.0))
    {}

    /*
     * gradeclass returns the grade grade associated with the score (e-value).
     */ 
    int gradeclass(int source, double score) {
	int klasse = 0;
	while (klasse < gradeclassnums(source)-1
	       && score >= gradeclassbounds[source][klasse])
	    klasse++;
	return klasse;
    }

    int gradeclassnums(int source) {  // for each source the number of classes   
	return gradequots[source].size();
    }
    void read(istream& datei, int source) {
	int numclasses;
	datei >> numclasses;
	if (numclasses < 0 || numclasses > 10){
	    cerr << "Error: number of classes=" << numclasses << endl; 
	    throw ProjectError("Error: number of classes out of range.");
	}
	gradeclassbounds[source].resize(numclasses-1);
	gradequots[source].resize(numclasses);
	for (int i=0; i<numclasses-1; i++) {
	    datei >> gradeclassbounds[source][i];
	    if (i>0 && gradeclassbounds[source][i] < gradeclassbounds[source][i-1]) {
		cerr << "Error: class bounds not increasing!" << endl;
		throw ProjectError("Error reading class bounds");
	    }
	}
	for (int i=0; i < numclasses; i++) 
	    datei >> gradequots[source][i];
    }
    double bonus;
    double malus;
    double localMalus; // (stronger) additional malus given the gene has evidence at other regions

    vector<vector<double> > gradeclassbounds; // for source and class the boundary

    /*
     * gradequot equals p+(g)/p-(g) for a given type
     * where g is the grade associated with the score (e-value).
     */ 
    vector<vector<double> > gradequots; // for each source and class the grade quotient
};

/*
 * FeatureCollection
 * holds all extrinsic feature information for a set of sequences
 * 
 */
class FeatureCollection {
public:
    FeatureCollection() :
	hasHintsFile(false),
	numSeqsWithInfo(0),
	numSources(1),
	malustable(NULL),
	localmalustable(NULL)
    {
	try {
	    offset = Properties::getIntProperty( "predictionStart" ) - 1;
	} catch (...) {
	    offset = 0;
	}
	string ssList = "gtag,gcag,";
	try {
	    ssList += Properties::getProperty("allow_hinted_splicesites");
	    ssList += ",";
	} catch (KeyNotFoundError) {}
	SequenceFeatureCollection::initHintedSplicesites(ssList);

	if (offset < 0)
	    offset = 0;
	Feature::offset = offset;
    }
    FeatureCollection(const FeatureCollection &other);
  // TODO deep assignment operator
  //FeatureCollection & operator = (const FeatureCollection & other){
  //    return *this;}

    ~FeatureCollection(){
	for (map<string, SequenceFeatureCollection*>::iterator it = collections.begin();
	     it != collections.end();
	     it++) {
	    delete it->second;
	    it->second = NULL;
	}
	if (malustable){
	  for (int type=0; type < NUM_FEATURE_TYPES; type++)
	    if (malustable[type])
	      delete [] malustable[type];
	  delete [] malustable;
	}
	if (localmalustable){
	  for (int type=0; type < NUM_FEATURE_TYPES; type++)
	    if (localmalustable[type])
	      delete [] localmalustable[type];
	  delete [] localmalustable;
	}
	/*if(sourceKey)
	    delete [] sourceKey;
	if(individual_liability)
	    delete [] individual_liability;
	if(oneGroupOneGene)
	    delete [] oneGroupOneGene;
	*/
    }

    SequenceFeatureCollection& getSequenceFeatureCollection(const char *seqname){
	SequenceFeatureCollection*& psfc = collections[seqname];
	if (psfc == NULL)
	    psfc = new SequenceFeatureCollection(this);
	return *psfc;
    }
    SequenceFeatureCollection* getSequenceFeatureCollection(string seqname){
	return collections[seqname];
    }
    bool isInCollections(string seqname){return collections.count(seqname)>0;}
    void readGFFFile(const char *filename);
    void setBonusMalus(Feature& f);
    void readExtrinsicCFGFile();
    // reading the extrinsic config file is split into two parts:
    void readSourceRelatedCFG(istream& datei); // reading in all source related parameters
    void readTypeInfo(istream& datei); // reading in the bonus/malus table
    int getNumSeqsWithInfo() { return numSeqsWithInfo;}
    int getNumCommonSeqs(AnnoSequence *annoseq);
    void printAccuracyForSequenceSet(const AnnoSequence* annoseqs, bool cleanRedundancies=true);
    void printAccuracyForSequenceSetOld(const AnnoSequence* annoseqs, bool cleanRedundancies=true);
    bool skeyExists(string skey);
    int esource(string skey);
    bool getIndividualLiability(string skey){return individual_liability[esource(skey)];}
    bool get1group1gene(string skey){return oneGroupOneGene[esource(skey)];}
    double malus(FeatureType type){
	return typeInfo[type].malus;
    }
    double localMalus(FeatureType type){
     return typeInfo[type].localMalus;
    }

    double partMalus(FeatureType type, int len);
    double localPartMalus(FeatureType type, int len, Double bonus, int nindep);
    bool hasHintsFile;
    int offset;
private:
    int numSeqsWithInfo;
    ifstream datei;

    int numSources;              // number of different sources of hints
    string *sourceKey;           // for each source the key string
    bool *individual_liability;  // for each source the flag individual_liability
    bool *oneGroupOneGene;       // for each source the flag 1group1gene
    FeatureTypeInfo typeInfo[NUM_FEATURE_TYPES]; // the type-specific information
    map<string, SequenceFeatureCollection*> collections;
    double** malustable;
    double** localmalustable;
};

#endif    //__EXTRINSICINFO_HH
