 /*****************************************************************************\
 * Filename : utrmodel.hh
 * Author   : Mario Stanke
 * Description: Untranslated Region Model Header File
 *
 *
 * Date       |   Author        |  Changes
 *------------|-----------------|----------------------------------------
 * 21.09.2005 | Mario Stanke    | creation of the file
 * 27.03.2006 | Mario Stanke    | introduced UTR intron
 * 05.04.2007 | Mario Stanke    | distance distribution of tata box to tss
 \******************************************************************************/

#ifndef _UTRMODEL_HH
#define _UTRMODEL_HH

#include "statemodel.hh"


/**
 * The UTR model class.
 *
 * author Mario Stanke
 */
 
class UtrModel : public StateModel {
public:
    UtrModel();
    ~UtrModel();

    StateType getStateType( ) const {
	return utype;
    }

  /**
   * Build all needed probabilities from the given
   * gene set.
   *
   * @param   annoseq A single linked list annotated sequences
   */
  void buildModel( const AnnoSequence* annoseq, int parIndex );
  void registerPars( Parameters* parameters);
  void processStates( const Gene* gene );
  void process5SingleExon( const State* exon, bool withLen=true );
  void process5InitialExon( const State* exon, bool withLen=true );
  void process5InternalExon( const State* exon);
  void process5TerminalExon( const State* exon);
  void process5Intron( int begin, int end);
  void process3SingleExon( const State* exon, bool withLen=true );
  void process3InitialExon( const State* exon, bool withLen=true );
  void process3InternalExon( const State* exon);
  void process3TerminalExon( const State* exon, bool withLen);
  void process3Intron( int begin, int end);
  
  
    /**
     * Print the wanted probabilities in files.
     * The names of the files should be given with
     * the Properties object /UtrModel/outfile
     */
    void printProbabilities   ( int zusNumber, BaseCount *bc, const char* suffix = NULL );
    void initAlgorithms       ( Matrix<Double>&, int);
    void viterbiForwardAndSampling(ViterbiMatrixType&, ViterbiMatrixType&, int, int, 
				   AlgorithmVariant, OptionListItem&);
    Double emiProbUnderModel  (int begin, int end) const;
    Double endPartEmiProb     (int begin, int end, int endOfBioExon) const;
    Double notEndPartEmiProb  (int begin, int end, int endOfBioExon, Feature *exonparts) const;
    void getEndPositions      ( int end, int &beginOfEndPart, int &endOfBioExon) const;
    Double tssupSeqProb       ( int left, int right, bool reverse) const;  
    Double tssProb            ( int left) const;
    static void computeTtsProbs    (int from, int to);
    static void init();
    static void resetPars(){
	if (utrcount == 0)
	    return;
	initSnippetProbs();
	initAlgorithmsCalled = false;
    }
    static void updateToLocalGC(int from = -1, int to = -1);
    static void clearSegProbs();
    static void readProbabilities(int zusNumber);
    static void readAllParameters();
    static void storeGCPars(int idx);
    static void resetModelCount(){utrcount = 0;};
    static void setTtsSpacing(int spacing){ ttsSpacing = spacing; };
    
private:
  Double seqProb            ( int left, int right, bool reverse, int type) const; // deprecated
  static void computeLengthDistributions( );
  static void fillTailsOfLengthDistributions( );
  void process5InitSequence( const char* start, const char* end);
  void process5Sequence( const char* start, const char* end);
  void process3Sequence( const char* start, const char* end);
  void processTssupSequence( const char* start, const char* end);
  void buildProbabilities ( const AnnoSequence* annoseq );
  void buildTSSModel( const AnnoSequence* annoseq );
  void buildTTSModel( const AnnoSequence* annoseq );
  int findTATA(const char* seq, int maxpos, bool reverseComplement=false) const;
  void processTSS(const char* start);
  void initCountVars      ( );
  Double longIntronProb(int internalBegin, int internalEnd) const;
  static void initSnippetProbs();
  void decrementEndOfPred( int &endOfPred, list<int>::iterator &eopit, bool inCache);
  void updatePossibleEOPs(list<int>::iterator &eopit, int endOfBioExon, bool &inCache);
private:
  StateType              utype;
  Integer                gweight;
  EOPList                eop;
  static Integer         utrcount;
  static vector<Integer> utr5_emicount;
  static vector<Integer> utr5init_emicount;
  static vector<Integer> utr3_emicount;
  static Double          utr_patpseudo;
  static PatMMGroup      utr5_emiprobs;
  static PatMMGroup      *GCutr5_emiprobs;
  static PatMMGroup      utr5init_emiprobs;
  static PatMMGroup      *GCutr5init_emiprobs;
  static PatMMGroup      utr3_emiprobs;
  static PatMMGroup      *GCutr3_emiprobs;
  static Integer         utr5init_gesbasen;
  static Integer         utr5_gesbasen;
  static Integer         utr3_gesbasen;
  static Integer         k;
  static double          utr5patternweight; // old way: this is applied AFTER reading from the parameters file
  static double          utr3patternweight; // old way: this is applied AFTER reading from the parameters file
  static double          utr5prepatternweight; // for computing a mixture directly after HMM training and BEFORE writing down to parameter file
  static double          utr3prepatternweight; // for computing a mixture directly after HMM training and BEFORE writing down to parameter file
  static vector<Integer> tssup_emicount;
  static Double          tssup_patpseudo;
  static vector<Double>  tssup_emiprobs;
  static vector<Double>  *GCtssup_emiprobs;
  static Integer         tssup_gesbasen;
  static Integer         tssup_k;
  static vector<Integer> lenCount5Single;       // Length count of single exons
  static vector<Integer> lenCount5Initial;      // Length count of initial exons
  static vector<Integer> lenCount5Internal;     // Length count of internal exons
  static vector<Integer> lenCount5Terminal;     // Length count of terminal exons
  static vector<Double>  lenDist5Single;        // Length distribution of single exons
  static vector<Double>  lenDist5Initial;       // Length distribution of initial exons
  static vector<Double>  lenDist5Internal;      // Length distribution of internal exons
  static vector<Double>  lenDist5Terminal;      // Length distribution of terminal exons
  static vector<Double>  tailLenDist5Single;    // Tail probabilities of the length distribution of single exons
  static vector<Integer> lenCount3Single;       // Length count of single exons
  static vector<Integer> lenCount3Initial;      // Length count of initial exons
  static vector<Integer> lenCount3Internal;     // Length count of internal exons
  static vector<Integer> lenCount3Terminal;     // Length count of terminal exons
  static vector<Double>  lenDist3Single;        // Length distribution of single exons
  static vector<Double>  lenDist3Initial;       // Length distribution of initial exons
  static vector<Double>  lenDist3Internal;      // Length distribution of internal exons
  static vector<Double>  lenDist3Terminal;      // Length distribution of terminal exons
  static vector<Double>  tailLenDist3Single;    // Tail probabilities of the length distribution of single exons
  static vector<Double>  tssProbsPlus;          // to store tss probabilities
  static vector<Double>  tssProbsMinus;         // to store tss probabilities
  static Integer         num5Single, num5Initial, num5Internal, num5Terminal, num5Introns;
  static Integer         numHuge5Single, numHuge5Initial, numHuge5Internal, numHuge5Terminal; 
  static Integer         num3Single, num3Initial, num3Internal, num3Terminal, num3Introns;
  static Integer         numHuge3Single, numHuge3Initial, numHuge3Internal, numHuge3Terminal; 
  static Integer         exonLenD;            // use detailed length distribution up to this number
  static Integer         max_exon_length;
  static Integer         max3singlelength;
  static Integer         max3termlength;
  static double          slope_of_bandwidth;  // for smoothing
  static Integer         minwindowcount;      // see class Smooth in commontrain.hh
  static Boolean         hasLenDist;
  static Integer         tss_start;
  static Integer         tss_end;
  static Integer         tata_start;
  static Integer         tata_end;
  static Integer         tata_pseudocount;
  static Integer         d_tss_tata_min;
  static Integer         d_tss_tata_max;
  static Motif           *tssMotif;           // motif of the transcription start site of tata-less promotors
  static Motif           *GCtssMotif;
  static Motif           *ttsMotif;           // motif of the transcription termination site (downstream of polyA signal)
  static Motif           *GCttsMotif;
  static Motif           *tssMotifTATA;       // motif of the transcription start site of tata promotors
  static Motif           *GCtssMotifTATA;
  static Motif           *tataMotif;          // motif of the tata box (if existent)
  static Motif           *GCtataMotif;
  // UTR intron related member variables
  static vector<Integer> intron_emicount;
  //static vector<Double>  intron_emiprobs;
  //static Integer         intron_k;            // order of the markov chain
  //  static SnippetProbs    *rInitSnippetProbs5, *rSnippetProbs3, *intronSnippetProbs;
  static SegProbs        *initSegProbs5, *segProbs5, *rInitSegProbs5, *rSegProbs5, *rSegProbs3, *segProbs3, *intronSegProbs;
  static bool            initAlgorithmsCalled, haveSnippetProbs;
  static vector<Integer> aataaa_count;
  static vector<Double>  aataaa_probs;
  static int             aataaa_boxlen;
  static string          polyasig_consensus;
  static int             d_polya_cleavage_min;
  static int             d_polya_cleavage_max;
  static double          prob_polya;
  static int             tts_motif_memory;
  static double pUtr5Intron, pUtr3Intron, prUtr5Intron, prUtr3Intron;
  static Double          *ttsProbPlus, *ttsProbMinus;
  static vector<Integer> distCountTata;
  static int             lastParIndex;
  static int             verbosity;
  static int             ttsSpacing; // without hints allow 3' end only every ttsSpacing bases for speed
};

class UtrModelError : public ProjectError {
public:
    UtrModelError(string msg) : ProjectError(msg) {}
};

#endif    //  _UTRMODEL_HH
