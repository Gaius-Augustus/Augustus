/*****************************************************************************\
 * Filename : hints.hh
 * Author   : Mario Stanke
 *
 *
 * Description: Hints on the gene structure
 *
 * Date       |   Author              |  Changes
 *------------|-----------------------|---------------------------------
 * 20.10.06   | Mario Stanke          | creation of the file by splitting the source file extrinsicinfo.hh
 * 01.07.08   | Mario Stanke          | added nonirpart hint type
 \******************************************************************************/

#ifndef __HINTS_HH
#define __HINTS_HH

// project includes
#include "types.hh"

// standard C/C++ includes
#include <cmath>  // for pow
#include <list>

#define NUM_FEATURE_TYPES 17
#define BONUS_FACTOR 1
#define QUOT_PSEUDOCOUNT 1
#define BLOCKSIZE 1000

using namespace std;

// Note: the order of FeatureTypes is important for the function compatibleWith
enum FeatureType {startF = 0, stopF, assF, dssF, tssF, ttsF, exonpartF, exonF, intronpartF, intronF, irpartF, CDSF, CDSpartF, UTRF, UTRpartF, nonexonpartF, nonirpartF};
extern const char* featureTypeNames[NUM_FEATURE_TYPES];

// start:       begin and end position of a start codon
// start:       begin and end position of a stop codon
// ass:         begin=end=the first intron base upstream of an exon
// dss:         begin=end=the first intron base downstream of an exon
// exonpart:    interval contained in a coding exon
// exon:        exact coding exon
// intronpart:  interval contained in an intron 
// intron:      exact intron 
// tss:         short region that contains the transcription initiation site
// tts:         short region that contains the transcription termination site
// irpart:      interval that is part of an intergenic region
// CDS:         exact coding sequence
// CDSpart:     interval contained in a coding sequence
// UTR:         exact utr exon (only the untranslated part of the exon)
// UTRpart:     part of utr exon
// nonexonpart: part of intergenic region or part of intron
// nonirpart:   part of genic region

bool isSignalType(FeatureType type);
bool isGFF(ifstream &istrm);

class Feature {
public:
    Feature(){
	bonus = malus = 1.0;
	esource = '?';
	active = true;
	discard = false;
	numContradicting = 0.0;
	numSupporting = 0;
	mult = 1;
    }

    Feature(int anfang, int ende, FeatureType typ, Strand strang, int leserahmen, string equelle) {
	start = anfang;
	end = ende;
	type = typ;
	strand = strang;
	frame = leserahmen;
	esource = equelle;
	score = 0.0;
	active = true;
	discard = false;
	numContradicting = 0.0;
	numSupporting = 0;
	mult = 1;
    }
  
    ~Feature(){
    }
   
    static FeatureType getFeatureType(string typestring);
    static FeatureType getFeatureType(int typeint);
    double exonpartMalus(int len){
	return pow(malus, len);
    }
    double distance_faded_bonus(int pos);
    bool compatibleWith(Feature &other);
    bool weakerThan(Feature &other, bool &strictly);
    double conformance();
    int length() {return end - start + 1;}
    void shiftCoordinates(int start,int end,bool rc = false);
    void setFrame(string f);
    void setStrand(string s);
    // fields of the GFF-format
    string seqname;
    string source;
    string feature;
    string groupname;
    int priority; // >=0 higher priority -> more important. -1 reserved for not specified
    int start, end;
    double score;
    double bonus;
    double malus;
    Strand strand;
    /* frame definition gff: One of '0', '1', '2' or '.'. '0' indicates that the specified region is in frame, 
       i.e. that its first base corresponds to the first base of a codon. '1' indicates 
       that there is one extra base, i.e. that the second base of the region corresponds 
       to the first base of a codon, and '2' means that the third base of the region is 
       the first base of a codon. If the strand is '-', then the first base of the region 
       is value of <end>, because the corresponding coding region will run from <end> to 
       <start> on the reverse strand. */
    int frame; 
    string attributes;
    int gradeclass;
    string esource; // 'annotrain' is reserved for annotation in the training
    FeatureType type;
    bool active;
    bool discard;
    Feature *next;  // used for making a partial list in SequenceFeatureCollections
    float numContradicting; // fractional number of other hints that contradict this one
    int numSupporting;
    static int offset;
    int mult; // multiplicity for summarizing several identical hints
};

ostream& operator<<(ostream&out, Feature& feature);
istream& operator>>( istream& in, Feature& feature );
 
bool operator<(const Feature& f1, const Feature& f2);
bool operator==(const Feature& f1, const Feature& f2);

/*
 * HintGroup
 * Hints that are known to belong to the same gene, 
 * for example, because they come from the same mRNA, form a group.
 */

class HintGroup{
public:
    HintGroup(){
	hints = NULL;
	name = "";
	incompGroups = strongerGroups = NULL;
	begin = end = -1;
	geneBegin = geneEnd = -1;
	priority = -1;
	copynumber = 1;
	trashy = false;
    }
    ~HintGroup(){
	if (hints)
	    delete hints;
	if (incompGroups)
	    delete incompGroups;
	if (strongerGroups)
	    delete strongerGroups;
    }
    friend bool operator<(const HintGroup& g1, const HintGroup& g2);
    friend bool operator==(const HintGroup& g1, const HintGroup& g2);
    string getName() const {return name;}
    int getPriority() const {return priority;}
    int getBegin() const {return begin;}
    int getEnd() const {return end;}
    int getGeneBegin() const {return geneBegin;}
    int getGeneEnd() const {return geneEnd;}
    int getCopyNumber() const {return copynumber;}
    void addCopyNumber(int n) {copynumber += n;}
    int getSize() const {if (hints) return hints->size(); else return 0;}
    string getSource() const {if (hints){ return hints->front()->esource;} else return "";}
    list<HintGroup*> *getIncompGroups(){return incompGroups;}
    list<HintGroup*> *getStrongerGroups(){return strongerGroups;}
    list<Feature*> *getHints() {return hints;}
    void print(ostream& out, bool withHints=false);
    void sortFeatures();
    void addFeature(Feature *hint);
    bool compatibleWith(HintGroup &other, Feature *&rascal1, Feature *&rascal2, bool &weakerThan);
    void updateFeatureConformance(HintGroup &other);
    bool nestedGenePossible(HintGroup &other);
    bool isTrashy();
    bool canCauseAltSplice();
    void setActiveFlag(bool active);
    void setDiscardFlag(bool discard);
    void addIncompGroup(HintGroup *otherGroup);
    void addStrongerGroup(HintGroup *otherGroup);
    void sortIncompGroup(){if (incompGroups) {incompGroups->sort();}}
private:
    list<Feature*> *hints;
    list<HintGroup*> *incompGroups; // incompatible HintGroups
    list<HintGroup*> *strongerGroups; // groups that are properly stronger
    string name;
    int priority;
    int begin;
    int end;
    int geneBegin;
    int geneEnd;
    int copynumber;
    bool trashy;
};

void printSrcGroupEvidence(list<HintGroup*> *groupList);

#endif    //__HINTS_HH
