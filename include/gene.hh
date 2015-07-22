/**********************************************************************
 * file:    gene.hh
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  
 * authors: Mario Stanke, mario.stanke@uni-greifswald.de
 *
 **********************************************************************/

#ifndef _GENE_HH
#define _GENE_HH

// project includes
#include "motif.hh"
#include "pp_scoring.hh"
#include "hints.hh"


// Forward declarations
class State;
class Transcript;
class Gene;
class AltGene;
class AnnoSequence;

template<class T>
struct ptr_comparison
{
    bool operator()(T* a, T* b) { return *a < *b; } 
};

class SrcEvidence {
public:
    SrcEvidence(string srcname){
	this->srcname = srcname;
	freq = 0;
    }
    string srcname;
    int freq;
    list<string> groupnames;
};

bool operator<(const SrcEvidence& e1, const SrcEvidence& e2);


/*
 * class SrcEvidence
 * A summary of extrinsic evidence by source of evidence.
 * Example:
 * E : 7 (e1,e2,e3,e4,e5,e6,e7)
 * P : 2 (p1,p2)
 * T : 1 (t1)
 * numEvidence = 10. numSources = 3. srcnames = [E,P,T], freq= [7,2,1], groupnames[0]= {e1,e2,...,e7}
 */
class Evidence {
public:
    Evidence(bool withNames) {
	numEvidence = 0;
	this->withNames = withNames;
    }
    ~Evidence() {
    }
    void add(string source, string name);
    void add(string source);
    void print();
    list<SrcEvidence> sourceEvidence;
    int numEvidence;
    bool withNames;
};

Evidence *stateEvidenceSummary(State *state1, State *state2 = NULL);


/**
 * class State
 * 
 */
class State {
public:
  /**
   * The beginning position of the current state.
   */
  int     begin;
  /**
   * The ending position of the current state.
   */
  int     end;
  /**
   * The next state.
   */
  State *next;
  StateType type;
  Double prob;
  bool hasScore;
  float apostprob; // score = posterior probability 
  int sampleCount;
  Evidence *evidence;
  char truncated; // truncated left, truncated right?
  int framemod; // reading frame modifier, frame = frame given by type + framemod (for truncation)
    
  State(int first, int last, StateType t) :
      begin(first), end(last),
      next(NULL),
      type(t),
      prob(0),
      hasScore(false),
      apostprob(0),
      sampleCount(1),
      evidence(NULL),
      truncated(0),
      framemod(0)
  {}

  State() : 
      begin(0), 
      end(0), 
      next(NULL), 
      type(TYPE_UNKNOWN),
      prob(0),
      hasScore(false),
      apostprob(0),
      sampleCount(1),
      evidence(NULL),
      truncated(0),
      framemod(0) {}
  ~State(){
    if (evidence)
	delete evidence;
  }
  State *cloneStateSequence(){
    State *erg = new State(*this);
    if (next)
      erg->next = next->cloneStateSequence();
    else 
      erg->next = NULL;
    return erg;
  }
  // copy constructor
    State(const State& other) :
	begin(other.begin),
	end(other.end),
	next(other.next),
	type(other.type),
	prob(other.prob),
	hasScore(other.hasScore),
	apostprob(other.apostprob),
	sampleCount(other.sampleCount),
	evidence(NULL),
	truncated(other.truncated),
	framemod(other.framemod)
    {
	if (other.evidence)
	    evidence = new Evidence(*other.evidence);
    }
  int frame(); // reading frame, in case it is a (coding) exon, usually depends only on the type
  bool frame_compatible(const Feature *hint);
  void addEvidence(string srcname) {if (!evidence) evidence = new Evidence(false); evidence->add(srcname);}
  bool operator< (const State &other) const;
  bool operator== (const State &other) const;
  int length() {
    return end - begin + 1;
  }
  Strand strand() {return isOnFStrand(type)? plusstrand : minusstrand;}
  State *getBiologicalState();
  void setTruncFlag(int end, int predEnd, int dnalen);
  void includeFrameModIntoType();
};

bool frame_compatible(State *ex1, State* ex2);

/**
 * A path through the Hidden-Markov-Model states
 * author Mario Stanke
 */
class StatePath {
public:
    StatePath(){
	first = NULL;
	pathemiProb = 0.0;
	intron_d = 0;
	seqname = "";
    }
    ~StatePath(){
	State *tmp;
	while( first ){
	    tmp = first->next;
	    delete first;
	    first = tmp;
	}
    }

    void push(State *st){
	st->next = first;
	first = st;
    }

    int size() {
	int n = 0;
	for (State* st = first; st != NULL; st = st->next) {
	    n++;
	}
	return n;
    }

    void print();
    static StatePath* condenseStatePath(StatePath *oldpath);
    Transcript* projectOntoGeneSequence(const char *genenames);
    static StatePath* getInducedStatePath(Transcript *genelist, int dnalen, bool printErrors=true);
    void reverse();
    bool operator== (const StatePath &other) const;
    bool operator< (const StatePath &other) const;
private:
    void pushIntron(int begin, int end, int frame, bool onFStrand);

public:
    string seqname;
    State* first;
    Double pathemiProb;
    int intron_d;
    list<PP::Match> proteinMatches;
};


int lenStateList(State *head);

/**
 * class Transcript
 * a single splice version, not neccessarily coding
 */
class Transcript {
public:
    Transcript() {
	exons = introns = (State*) NULL;
	next = NULL;
	transstart = transend = -1;
	seqname = id = source = "";
	strand = plusstrand;
	complete = true;
	apostprob = 0.0;
	hasProbs = false;
	throwaway = false;
	viterbi = true;
    }
    Transcript (const Transcript& other);
    Transcript& operator= (const Transcript& other);
    virtual ~Transcript(){
	State *st, *tmp;
	list<State*> sl = getExInInHeads();
	for (list<State*>::iterator it = sl.begin(); it != sl.end(); ++it){
	    st = *it;
	    while( st ){
		tmp = st->next;
		delete st;
		st = tmp;
	    }
	}
    }
    virtual Transcript* clone() {return new Transcript(*this);}
    Transcript *cloneGeneSequence(){
	Transcript *res = clone(); // returns a Transcript or a Gene as appropriate
	if (next)
	    res->next = next->cloneGeneSequence();
	else 
	    res->next = NULL;
	return res;
    }
    static void destroyGeneSequence(Transcript *head) {
	Transcript* nextHead = head;
	while (nextHead) {
	    head = nextHead;
	    nextHead = nextHead->next;
	    delete head;
	}
    }
    void addStatePostProbs(float p);
    void setStatePostProbs(float p);
    void addSampleCount(int k);
    void setSampleCount(int k);
    virtual int geneBegin() const { return transstart;}
    virtual int geneEnd() const { return transend;}
    virtual bool isCoding() const { return false; }
    bool operator< (const Transcript &other) const;
    bool operator== (const Transcript &other) const;
    void normPostProb(float n);
    void updatePostProb(Transcript* other);
    virtual list<State*> getExInHeads() const {
	list<State*> L;
	L.push_back(exons);
	L.push_back(introns);
	return L;
    }
    virtual list<State*> getExInInHeads() const { return getExInHeads();}
    double meanStateProb();
    char* getExonicSequence(AnnoSequence *annoseq = NULL,
			    bool noOffset = false) const; // CDS or whole RNA, respectively
    virtual void shiftCoordinates(int d);
    virtual bool almostIdenticalTo(Transcript *other);
    virtual void printCodingSeq(AnnoSequence *annoseq) const {}; // print nothing
    virtual void printProteinSeq(AnnoSequence *annoseq) const {};// for noncoding
    virtual void printBlockSequences(AnnoSequence *annoseq) const {}; // genes
    virtual void printGFF() const;
    virtual void printEvidence() const {}; // implemented only for coding genes
    void setStateHasScore(bool has);
    static Transcript* getGenesOnStrand(Transcript* genes, Strand strand);
    static void filterTranscriptsByMaxTracks(list<Transcript*> &gl, int maxTracks);
    virtual double supportingFraction(HintGroup *group) {return 0.0;}
public:
    State*      exons; // the exons (not UTR)
    State*      introns; // the introns between 'exons'
    int         transstart, transend; // transcription boundaries, -1 if not known
    Transcript* next;
    string  id;       // transcript id
    string  seqname;
    string  source;
    Strand  strand;
    Boolean complete; //coding region is complete
    string  geneid;   // id of the AltGene
    float   apostprob;
    bool    hasProbs;
    bool    throwaway;
    bool    viterbi;

    static bool gff3;
    static bool print_tss;
    static bool print_tts;
};


void filterGenePrediction(list<Transcript*> &gl, list<Transcript*> &filteredTranscripts, const char *seq, Strand strand, bool noInFrameStop, double minmeanexonintronprob=0.0, double minexonintronprob=0.0);

/**
 * class Gene
 * a single splice version (transcript) of a coding gene
 */
class Gene : public Transcript {
public:
    Gene() {
	utr5exons = utr3exons = utr5introns = utr3introns = (State*) 0;
	length   = clength = 0;
	weight   = 1;
	strand   = plusstrand;
	frame    = 0;
	complete5utr = complete3utr = true;
	codingstart = codingend = -1;
	supportingEvidence = incompatibleEvidence = CDSexonEvidence = CDSintronEvidence = UTR5stateEvidence = UTR3stateEvidence = NULL;
    }
    Gene(const Gene& other);
    virtual Gene* clone() { return new Gene(*this); }
    ~Gene(){ // this calls ~Transcript implicitly, which deletes exons and introns
	if (supportingEvidence)
	    delete supportingEvidence;
	if (incompatibleEvidence)
	    delete incompatibleEvidence;
	if (CDSintronEvidence)
	    delete CDSintronEvidence;
	if (CDSexonEvidence)
	    delete CDSexonEvidence;
	if (UTR5stateEvidence)
	    delete UTR5stateEvidence;
	if (UTR3stateEvidence)
	    delete UTR3stateEvidence;
    }
    list<State*> getExInHeads() const { list<State*> L; L.push_back(exons); L.push_back(introns); L.push_back(utr5exons); L.push_back(utr3exons); return L;}
    list<State*> getExInInHeads() const { list<State*> L = getExInHeads(); L.push_back(utr5introns); L.push_back(utr3introns); return L;}
    bool hasInFrameStop(AnnoSequence *annoseq) const;
    // void computeBC(char *seq);
    int numExons() const;
    State *lastExon() const;
    bool identicalCDS(Gene *other);
    using Transcript::almostIdenticalTo; // not really, but to prevent the compiler warning on MAC
    virtual bool almostIdenticalTo(Gene *other);
    void shiftCoordinates(int d);
    int geneBegin() const { return (transstart>=0)? transstart : codingstart;}
    int geneEnd() const { return (transend>=0)? transend : codingend;}
    virtual bool isCoding() const { return true; }
    void addUTR(State *mrnaRanges, bool complete_l=true, bool complete_r=true);
    void compileExtrinsicEvidence(list<HintGroup>  *groupList);
    double supportingFraction(HintGroup *group);
    void addSupportedStates(HintGroup *group);
    double getPercentSupported() const;
    int getCDSCoord(int loc, bool comp) const;  
    bool completeCDS() const;
    void print();
    void printGFF() const;
    void printCodingSeq(AnnoSequence *annoseq) const;
    void printProteinSeq(AnnoSequence *annoseq) const;
    void printBlockSequences(AnnoSequence *annoseq) const;
    virtual void printEvidence() const;
    void truncateMaskedUTR(AnnoSequence *annoseq);

    static void init();

    /// members for UnTranslated Region
    State*  utr5exons;
    State*  utr3exons;
    State*  utr5introns;
    State*  utr3introns;
    int     codingstart, codingend; // transstart <= codingstart <= codingend <= transend, if not -1
    bool    complete5utr;
    bool    complete3utr;
  
    /// The length of the span of the coding part (with introns)
    int     length;
    /// The coding length of the gene
    int     clength;
    /// The reading frame position of the first base (usually 0)
    int     frame;
    BaseCount bc;
    int weight;
    Evidence *supportingEvidence;
    Evidence *incompatibleEvidence;
    Evidence *CDSexonEvidence;
    Evidence *CDSintronEvidence;
    Evidence *UTR5stateEvidence;
    Evidence *UTR3stateEvidence;

    list<PP::Match> proteinMatches;  // true if gene matches protein profile

    /// output options
    static bool print_start;
    static bool print_stop;
    static bool print_introns;
    static bool print_cds;
    static bool print_exonnames;
    static bool stopCodonExcludedFromCDS;
    static bool print_utr;
    static bool print_blocks;
};


/*
 * AltGene is a gene in the original sense. It can contain several transcripts (coding or not).
 * 
 */
class AltGene {
public:
    list<Transcript*> transcripts;
    int mincodstart; // leftmost position of any start codon (if coding)
    int maxcodend; //
    Strand strand;
    string id;
    string seqname;
    float apostprob;
    bool hasProbs;

    AltGene() {
	mincodstart = maxcodend = -1;
	strand = plusstrand;
	apostprob = 0.0;
	hasProbs = 0;
    }
    bool operator< (const AltGene &other) const;
    void addGene(Transcript* tx);
    bool overlaps(Transcript *tx);
    void shiftCoordinates(int d);
    void sortTranscripts(int numkeep=-1);
    void deleteSuboptimalTranscripts(bool uniqueCDS);
    int minTransBegin();
    int maxTransEnd();
    bool isCoding(){ return transcripts.empty() || dynamic_cast<Gene*> (transcripts.front());}
};

Transcript* getPtr(list<AltGene> *gl);

void printGeneList(list<AltGene> *genelist, AnnoSequence *annoseq, bool withCS, bool withAA, bool withEvidence);
void printGeneList(Transcript* seq, AnnoSequence *annoseq, bool withCS, bool withAA);
void printGeneSequence(Transcript* seq, AnnoSequence *annoseq = NULL, bool withCS=false, bool withAA=true);
list<Gene*>* sortGenePtrList(list<Gene*>);
list<AltGene> *reverseGeneList(list<AltGene> *altGeneList, int endpos);
list<AltGene>* groupTranscriptsToGenes(list<Transcript*> &transcripts);
 
void reverseGeneSequence(Transcript* &seq, int endpos);
void postProcessGenes(list<AltGene> *genes, AnnoSequence *annoseq);
Gene* promoteToCoding(Transcript* tx, AnnoSequence *as); // make a coding gene from a non-coding if it contains a good ORF

class Annotation {
public:
  Annotation() {
    genes = lastGene = (Gene*) 0;
    forwardGenes =  backwardGenes  = lastForwardGene = lastBackwardGene = (Gene*) 0;
    forwardPath = condensedForwardPath = backwardPath = condensedBackwardPath = (StatePath*) 0;
    path = condensedPath = (StatePath*) 0;
    emiProb = 1.0;
    forwardEmiProb = 1.0;
    backwardEmiProb = 1.0;
  }
  ~Annotation() {
    Transcript::destroyGeneSequence(genes);
    if (path)
      delete path;
    if (condensedPath)
      delete condensedPath;

    Transcript::destroyGeneSequence(forwardGenes);
    Transcript::destroyGeneSequence(backwardGenes);
    if (forwardPath)
      delete forwardPath;
    if (backwardPath)
      delete backwardPath;
    if (condensedForwardPath)
      delete condensedForwardPath;
    if (condensedBackwardPath)
      delete condensedBackwardPath;
  }

  void appendGene(Transcript* gene);
  void appendForwardGene(Transcript* gene);
  void appendBackwardGene(Transcript* gene);
  const void printGFF() const;

  StatePath *path;
  StatePath *condensedPath;

  StatePath *forwardPath;            
  StatePath *backwardPath;           
  StatePath *condensedForwardPath;   
  StatePath *condensedBackwardPath;  

  Transcript *genes;
  Transcript *forwardGenes;  
  Transcript *backwardGenes; 

  Double emiProb;
  Double forwardEmiProb;  
  Double backwardEmiProb; 

private:
  Transcript *lastGene;
  Transcript *lastForwardGene, *lastBackwardGene; 
};


class AnnoSequence {
public:
    AnnoSequence(){
	length = 0;
	sequence = 0;
	seqname = 0;
	next = (AnnoSequence*) 0;
	anno = (Annotation*) 0;
	weight = 1;
	offset = 0;
    }
    ~AnnoSequence(){
	if (sequence)
	    delete[] sequence;
	if (seqname)
	    delete [] seqname;
	if (anno)
	    delete anno;
    }
    static void deleteSequence(AnnoSequence *head){
	AnnoSequence* nextHead = head;
	while(nextHead) {
	    head = nextHead;
	    nextHead = nextHead->next;
	    delete head;
	}
    }
    void setWeight(int w) {
	weight = w;
	if (!anno)
	    return;
	
	for (Transcript* t = anno->genes; t != NULL; t = t->next){
	    Gene *g = dynamic_cast<Gene*> (t);
	    if (g)
		g->weight = w;
	}
    }
    /*
     * reverse complement the sequence and all genes
     */
    AnnoSequence* getReverseComplement();
    void printGFF();

    char*        seqname;
    int          length;
    int          offset; // if >0 this is the number of bases that have been cut out in the front of the original input sequence
    char*        sequence;
    BaseCount    bc;
    AnnoSequence *next;
    Annotation*  anno;
    int          weight;
};

class AnnoSeqGeneIterator {
public:
  AnnoSeqGeneIterator(const AnnoSequence* annoseqHead){
    annoseq = annoseqHead;
    while (annoseq && !(annoseq->anno && annoseq->anno->genes))
      annoseq = annoseq->next;
    if (annoseq)
      gene = annoseq->anno->genes;
    else 
      gene = NULL;
  };
  friend const AnnoSeqGeneIterator& operator++(AnnoSeqGeneIterator& gi);
    
  bool hasMoreElements(){
    return (annoseq != NULL && gene != NULL);
  }
  const AnnoSequence *annoseq;
  const Transcript *gene;
};


class CompareStatePathPtr{
public: 
  CompareStatePathPtr(){ }
  ~CompareStatePathPtr(){ }
  bool operator()(StatePath *first, StatePath *second){
    return (first->pathemiProb > second->pathemiProb);
  }
};

class StatePathCollection {
public:
  StatePathCollection(){
    sorted = false;
  };
  ~StatePathCollection(){
    StatePath *p;
    while (!pathlist.empty()){
      p = pathlist.front();
      pathlist.pop_front();
      delete p;
    }
  }
  void addPath(StatePath *p){
    if (!containsPath(p))
      pathlist.push_back(p);
    sorted = false;
  }
  int size(){
    return pathlist.size();
  }
  void sort(){
    pathlist.sort(CompareStatePathPtr());
    sorted = true;
  }
  void printAPosterioriProbs(ostream& out) {
    for(list<StatePath*>::iterator it=pathlist.begin(); it!=pathlist.end(); it++)
      out << (*it)->pathemiProb << " , ";
    out << endl;
  }
  bool containsPath(StatePath *p);
  int positionInCollection(StatePath *p);

private:
  list<StatePath*> pathlist;
  bool sorted;
};

void examineBaseCountOfGeneSeq(AnnoSequence  *as);

/* class CodonUsage: not used right now
class CodonUsage {
public:
  CodonUsage(){
    init();
  };
  CodonUsage(char *seq, int len, int frame=0){
    init();
    addCUofSeq(seq, len, frame);
  };
  CodonUsage(char *seq, int frame=0){
    init();
    addCUofSeq(seq, strlen(seq), frame=0);
  };
  CodonUsage(TrainingData *genedata) {	
    init();
    while (genedata) {
      addCUofSeq(genedata->seq, genedata->seqLen, 0);
      genedata = genedata->next;
    }
    computeUsage();
  };

  double meanLogProb(char *seq, int len, int frame=0);
  void addCUofSeq(char *seq, int len, int frame=0);
  void init();
  void computeUsage();
  void print(ostream &out);
private:
  //GeneticCode code;
  int codoncount[64];
  double codonusage[64];
  double logcodonusage[64];
  double aaFrequencies[20];
};
*/ // class CodonUsage


class FreqSegment{
public:
    FreqSegment(int s, int f){
	start = s;
	freq = f;
    }
    int start;
    int freq;
};

#endif   //  _GENE_HH
