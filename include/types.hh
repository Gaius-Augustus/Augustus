/*****************************************************************************\
 * Filename : types.hh
 * Authors  : Emmanouil Stafilarakis, Mario Stanke
 *
 * Copyright: ©Stafilarakis, Stanke
 *
 * Description: Several typedefs and definitions.
 *
 *
 * Date       |   Author              |  Changes
 *------------|-----------------------|----------------------------------------
 * 26.09.2001 | Stafilarakis Emm.     | Creation of the file
 * 07.05.2002 | Mario Stanke          | Included the defs of dss and ass sizes
 * 24.03.2003 | Mario Stanke          | introducing shadow states 
 * 19.05.2003 | Stanke                | recursively include other property files
 * 10.05.2006 | Mario Stanke          | added 10 new species in the usage info
 * 15.05.2006 | Mario Stanke          | added 3' UTR states
 * 29.02.2012 | Mario Stanke          | added quantile function
 * 25.02.2014 | Stefanie Koenig       | added ftoa function (double to string)
 * 28.09.2014 | Mario Stanke          | added noncoding (nc) states
\******************************************************************************/

#ifndef _TYPES_HH
#define _TYPES_HH

// project includes
#include "lldouble.hh"

// standard C/C++ includes
#include <string>
#include <exception>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <limits>
#include <algorithm> // for std:sort
#include <map>

using namespace std;


enum Strand {STRAND_UNKNOWN=-1, plusstrand, minusstrand, bothstrands};
char strandChar (Strand s);
ostream& operator<< (ostream& strm, const Strand s);

/**
 * ASS       = acceptor splice site, between intron and exon
 * DSS       = donor splice site , between exon and intron
 * stop codons depend on genetic code, prefer to use Geneticcode::isStopcodon
 */
#define ASS_MIDDLE      2
#define DSS_MIDDLE      2
#define STARTCODON_LEN  3
#define STOPCODON_LEN   3
#define OCHRECODON      "taa"
#define OPALCODON       "tga"
#define AMBERCODON      "tag"
#define RCOCHRECODON    "tta"
#define RCOPALCODON     "tca"
#define RCAMBERCODON    "cta"
#define TRUNC_LEFT      1
#define TRUNC_RIGHT     2 

#define SPECIES_SUBDIR "species/"
#define MODEL_SUBDIR "model/"
#define EXTRINSIC_SUBDIR "extrinsic/"

#define VERSION "3.2.1"

#define PREAMBLE "# This output was generated with AUGUSTUS (version " VERSION ").\n\
# AUGUSTUS is a gene prediction tool written by M. Stanke (mario.stanke@uni-greifswald.de),\n\
# O. Keller, S. König, L. Gerischer and L. Romoth.\n\
# Please cite: Mario Stanke, Mark Diekhans, Robert Baertsch, David Haussler (2008),\n\
# Using native and syntenically mapped cDNA alignments to improve de novo gene finding\n\
# Bioinformatics 24: 637-644, doi 10.1093/bioinformatics/btn013"

#define GREETING "AUGUSTUS (" VERSION ") is a gene prediction tool\n\
written by M. Stanke, O. Keller, S. König, L. Gerischer and L. Romoth."

#define SPECIES_LIST "usage:\n\
augustus [parameters] --species=SPECIES queryfilename\n\
\n\
where SPECIES is one of the following identifiers\n\
\n\
identifier                               | species\n\
-----------------------------------------|----------------------\n\
human                                    | Homo sapiens\n\
fly                                      | Drosophila melanogaster\n\
arabidopsis                              | Arabidopsis thaliana\n\
brugia                                   | Brugia malayi\n\
aedes                                    | Aedes aegypti\n\
tribolium                                | Tribolium castaneum\n\
schistosoma                              | Schistosoma mansoni\n\
tetrahymena                              | Tetrahymena thermophila\n\
galdieria                                | Galdieria sulphuraria\n\
maize                                    | Zea mays\n\
toxoplasma                               | Toxoplasma gondii\n\
caenorhabditis                           | Caenorhabditis elegans\n\
(elegans)                                | Caenorhabditis elegans \n\
aspergillus_fumigatus                    | Aspergillus fumigatus\n\
aspergillus_nidulans                     | Aspergillus nidulans\n\
(anidulans)                              | Aspergillus nidulans\n\
aspergillus_oryzae                       | Aspergillus oryzae\n\
aspergillus_terreus                      | Aspergillus terreus\n\
botrytis_cinerea                         | Botrytis cinerea\n\
candida_albicans                         | Candida albicans\n\
candida_guilliermondii                   | Candida guilliermondii\n\
candida_tropicalis                       | Candida tropicalis\n\
chaetomium_globosum                      | Chaetomium globosum\n\
coccidioides_immitis                     | Coccidioides immitis\n\
coprinus                                 | Coprinus cinereus\n\
coprinus_cinereus                        | Coprinus cinereus\n\
cryptococcus_neoformans_gattii           | Cryptococcus neoformans gattii\n\
cryptococcus_neoformans_neoformans_B     | Cryptococcus neoformans neoformans\n\
cryptococcus_neoformans_neoformans_JEC21 | Cryptococcus neoformans neoformans\n\
(cryptococcus)                           | Cryptococcus neoformans\n\
debaryomyces_hansenii                    | Debaryomyces hansenii\n\
encephalitozoon_cuniculi_GB              | Encephalitozoon cuniculi\n\
eremothecium_gossypii                    | Eremothecium gossypii\n\
fusarium_graminearum                     | Fusarium graminearum\n\
(fusarium)                               | Fusarium graminearium\n\
histoplasma_capsulatum                   | Histoplasma capsulatum\n\
(histoplasma)                            | Histoplasma capsulatum\n\
kluyveromyces_lactis                     | Kluyveromyces lactis\n\
laccaria_bicolor                         | Laccaria bicolor\n\
lodderomyces_elongisporus                | Lodderomyces elongisporus\n\
magnaporthe_grisea                       | Magnaporthe grisea\n\
neurospora_crassa                        | Neurospora crassa\n\
(neurospora)                             | Neurospora crassa\n\
phanerochaete_chrysosporium              | Phanerochaete chrysosporium\n\
(pchrysosporium)                         | Phanerochaete chrysosporium\n\
pichia_stipitis                          | Pichia stipitis\n\
rhizopus_oryzae                          | Rhizopus oryzae\n\
saccharomyces_cerevisiae_S288C           | Saccharomyces cerevisiae\n\
saccharomyces_cerevisiae_rm11-1a_1       | Saccharomyces cerevisiae\n\
(saccharomyces)                          | Saccharomyces cerevisiae\n\
schizosaccharomyces_pombe                | Schizosaccharomyces pombe\n\
ustilago_maydis                          | Ustilago maydis\n\
(ustilago)                               | Ustilago maydis\n\
yarrowia_lipolytica                      | Yarrowia lipolytica\n\
nasonia                                  | Nasonia vitripennis\n\
tomato                                   | Solanum lycopersicum\n\
chlamydomonas                            | Chlamydomonas reinhardtii\n\
amphimedon                               | Amphimedon queenslandica\n\
pea_aphid                                | Acyrthosiphon pisum\n\
\n"

#define HELP_USAGE "usage:\n\
augustus [parameters] --species=SPECIES queryfilename\n\
\n\
'queryfilename' is the filename (including relative path) to the file containing the query sequence(s)\n\
in fasta format.\n\
\n\
SPECIES is an identifier for the species. Use --species=help to see a list.\n\
\n\
parameters:\n\
--strand=both, --strand=forward or --strand=backward\n\
--genemodel=partial, --genemodel=intronless, --genemodel=complete, --genemodel=atleastone or --genemodel=exactlyone\n\
  partial      : allow prediction of incomplete genes at the sequence boundaries (default)\n\
  intronless   : only predict single-exon genes like in prokaryotes and some eukaryotes\n\
  complete     : only predict complete genes\n\
  atleastone   : predict at least one complete gene\n\
  exactlyone   : predict exactly one complete gene\n\
--singlestrand=true\n\
  predict genes independently on each strand, allow overlapping genes on opposite strands\n\
  This option is turned off by default.\n\
--hintsfile=hintsfilename\n\
  When this option is used the prediction considering hints (extrinsic information) is turned on.\n\
  hintsfilename contains the hints in gff format.\n\
--AUGUSTUS_CONFIG_PATH=path\n\
  path to config directory (if not specified as environment variable)\n\
--alternatives-from-evidence=true/false\n\
  report alternative transcripts when they are suggested by hints\n\
--alternatives-from-sampling=true/false\n\
  report alternative transcripts generated through probabilistic sampling\n\
--sample=n\n\
--minexonintronprob=p\n\
--minmeanexonintronprob=p\n\
--maxtracks=n\n\
  For a description of these parameters see section 4 of README.TXT.\n\
--proteinprofile=filename\n\
  When this option is used the prediction will consider the protein profile provided as parameter.\n\
  The protein profile extension is described in section 7 of README.TXT.\n\
--progress=true\n\
  show a progressmeter\n\
--gff3=on/off\n\
  output in gff3 format\n\
--predictionStart=A, --predictionEnd=B\n\
  A and B define the range of the sequence for which predictions should be found.\n\
--UTR=on/off\n\
  predict the untranslated regions in addition to the coding sequence. This currently works only for a subset of species.\n\
--noInFrameStop=true/false\n\
  Do not report transcripts with in-frame stop codons. Otherwise, intron-spanning stop codons could occur. Default: false\n\
--noprediction=true/false\n\
  If true and input is in genbank format, no prediction is made. Useful for getting the annotated protein sequences.\n\
--uniqueGeneId=true/false\n\
  If true, output gene identifyers like this: seqname.gN\n\
\n\
For a complete list of parameters, type \"augustus --paramlist\".\nAn exhaustive description can be found in the file README.TXT.\n"

#define HELP_USAGE_ETRAINING "usage:\n\
etraining --species=SPECIES trainfilename\n\
\n\
'trainfilename' is the filename (including relative path) to the file in genbank format containing the training sequences. These can be multi-gene sequences.\n\
SPECIES is a name for your species.\n\
\n\
parameters:\n\
--/genbank/verbosity=n\n\
  Choose one of 0,1,2 or 3. The larger the verbosity, the more (error) messages you get.\n"


#define POWER4TOTHE(x) ((int) (pow((double) 4,(double) (x)) + 0.1))

/*
 * bonus for (almost sure) hint with negative score
 * the user sets an (almost sure) anchor
 */

#define BONUSSUREHINT 1e100

/**
 * 
 */
typedef long int    Integer;
/**
 * 
 */
typedef int         Short;
/**
 * 
 */
typedef bool        Boolean;
/**
 * 
 */
typedef LLDouble Double;
/**
 * 
 */
typedef double      Float;


class Constant {
public:
    static void init();
    static string fullSpeciesPath() {return configPath + speciesDir;}
    static string modelPath() {return configPath + MODEL_SUBDIR;}
    static string extrinsicPath() {return configPath + EXTRINSIC_SUBDIR;}

    // size of splicesite patterns
    static int ass_size() { return ass_start + ass_end; }           // excluding "AG"
    static int ass_whole_size() { return ass_size() + ASS_MIDDLE; } // including "AG"
    static int dss_size() { return dss_start + dss_end; }           // excluding "GT"
    static int dss_whole_size() { return dss_size() + DSS_MIDDLE; } // including "GT"
    
    static string configPath;
    static string speciesDir;
    static int decomp_num_at;
    static int decomp_num_gc;
    static int decomp_num_steps;
    static int trans_init_window;
    static int tis_maxbinsize;
    static int ass_upwindow_size;
    static int init_coding_len;
    static int et_coding_len;
    static int ass_start;    //  size of ASS pattern before "AG" 
    static int ass_end;      //                      after  "AG"
    static int dss_start;    //  size of DSS pattern before "GT" 
    static int dss_end;      //                      after  "GT"
    static int dss_maxbinsize;
    static int ass_maxbinsize;
    static int tss_upwindow_size;   // win from start of model up to tss
    static int tss_start;
    static int min_coding_len;
    static int max_exon_len;
    static Integer d_polyasig_cleavage;
    static bool keep_viterbi;
    static double gc_range_min;
    static double gc_range_max;
    static double probNinCoding;
    static double opalprob;
    static double amberprob;
    static double ochreprob;
    static bool utr_option_on;
    static bool nc_option_on;
    static Integer augustus_verbosity;
    static bool alternatives_from_evidence;
    static double subopt_transcript_threshold;
    static Integer almost_identical_maxdiff;  // maximum allowed difference for the ends of the transcript
    static bool uniqueGeneId;
    static double max_contra_supp_ratio;
    static bool reportUtrOnlyGenes;
    static bool useCRFtraining;
    static bool CRFtrainTIS;
    static bool CRFtrainSS;
    static bool CRFtrainIntron;
    static bool CRFtrainIgenic;
    static bool CRFtrainCDS;
    static bool CRFtrainUTR;
    static bool dss_gc_allowed;
    static Boolean tieIgenicIntron; // whether to tie igenic model parameters to intron model parameters, i.e. use just one content model, that of the intron
    static Boolean exoncands;
    static Boolean proteinOutput;
    static Boolean codSeqOutput;
    static Boolean contentmodels; // whether to use content models, default: true
    static Integer min_intron_len; // minimal intron length (for hints)
    static bool MultSpeciesMode; // whether we do comparative gene prediction in multiple species
    static string  treefile; // file name in which a tree is specified in Newick format
    static string speciesfilenames; // file name to file which contains the names of species and the corresponding file names
    static string dbaccess; // comma separated string with database access (hostname, database name, table name, user, passwd
    static string alnfile; // name of file that contains MSA of genomes
    static bool overlapmode; // whether overlapping exons are allowed in Viterbi algorithm
    static string orthoexons; //name of file that contains list of orthologous exons
    static Integer maxOvlp; // parameters for overlapping coding regions in bacteria
    static vector<Double> head2tail_ovlp;
    static vector<Double> head2head_ovlp;
    static vector<Double> tail2tail_ovlp;
    static unsigned temperature; // heating the distribution for sampling, 0=cold, 7=hottest, take probs to the power of (8-temperature)/8
    static bool softmasking; // if true, lower-case character regions give rise to nonexonpart hints
    static bool dbhints;
    // scores from logistic regression  
    static bool logreg;
    static vector<double> ex_sc;
    static vector<double> in_sc;
    static vector<double> lg_es;
};


extern bool inCRFTraining;

extern const int power2id[31];

#define A_SET_FLAG(x) power2id[x]

struct Bitmask {
    Bitmask (int n=0) : value(n) {}
    Bitmask (const Bitmask& other) : value(other.value) {}
    bool operator[] (int n) const {
	return value & A_SET_FLAG(n);
    }
    void set(int n) {
	value |= A_SET_FLAG(n);
    }
    void unset(int n) {
	value &= ~A_SET_FLAG(n);
    }
    static Bitmask any() {
	return Bitmask(-1);
    }
    int value;
};
    

/*
 * Base class for all exception classes in the project.
 *
 * @author Emmanouil Stafilarakis
 */
class ProjectError : public exception {
public:
    /**
     * Constructor
     * @param msg The describing error message
     */
    ProjectError( string msg ) throw () : exception(), message( msg ) {
    }
    /**
     * Constructor
     */
    ProjectError( ) throw () : exception() { }

    /**
     * Destructor
     */
    ~ProjectError() throw() {}

    /**
     * Get the message describing the error.
     * @return A string message.
     */
    string getMessage( ) { return message; }
private:
    string message;
};


class InvalidNucleotideError : public ProjectError {
public:
    InvalidNucleotideError( char t ) : ProjectError( "Invalid nucleotide '" + string(1,t) + "' encountered." ) { }
};


#define NUM_TYPES 86

enum StateType{TYPE_UNKNOWN = -1, igenic, 
	       // forward strand
	       singleG, initial0, initial1, initial2, internal0, internal1, internal2, terminal,
	       lessD0, longdss0, equalD0, geometric0, longass0,  // The five intron states for frame 0
	       lessD1, longdss1, equalD1, geometric1, longass1,  // The five intron states for frame 1
	       lessD2, longdss2, equalD2, geometric2, longass2,  // The five intron states for frame 2
	       utr5single, utr5init, utr5intron, utr5intronvar, utr5internal, utr5term, // 5'UTR states
	       utr3single, utr3init, utr3intron, utr3intronvar, utr3internal, utr3term, // 3'UTR states
	       // reverse strand
	       rsingleG, rinitial, rinternal0, rinternal1, rinternal2, rterminal0, rterminal1, rterminal2,
	       rlessD0, rlongdss0, requalD0, rgeometric0, rlongass0,  // The five intron states for frame 0
	       rlessD1, rlongdss1, requalD1, rgeometric1, rlongass1,  // The five intron states for frame 1
	       rlessD2, rlongdss2, requalD2, rgeometric2, rlongass2,  // The five intron states for frame 2
	       rutr5single, rutr5init, rutr5intron, rutr5intronvar, rutr5internal, rutr5term, // 5'UTR states
	       rutr3single, rutr3init, rutr3intron, rutr3intronvar, rutr3internal, rutr3term, // 3'UTR states
	       
	       intron_type, rintron_type, exon_type,
	       // non-protein-coding states
	       ncsingle, ncinit, ncintron, ncintronvar, ncinternal, ncterm, // forward strand
	       rncsingle, rncinit, rncintron, rncintronvar, rncinternal, rncterm // reverse strand
};

extern const char* stateTypeNames[NUM_TYPES];
extern const char* stateTypeIdentifiers[NUM_TYPES];
extern const int stateReadingFrames[NUM_TYPES];



/*
 * mod3, returns the number in {0,1,2} being congruent to the argument modulo 3
 */
inline int mod3 (int k) {
    return (k>=0)? k%3 : (k%3+3)%3;
}

inline int modm(int k, int m) {
    return (k>=0)? k%m : (k%m+m)%m;
}


/* 
 * state type functions
 */

inline StateType toStateType( const char* str ){
    int i;
    for (i=0; i<NUM_TYPES; i++)
	if (strcmp(str, stateTypeIdentifiers[i]) == 0)
	    return (StateType) i;
    return TYPE_UNKNOWN;
}

inline Boolean isInitialExon(StateType type){
    return (type==initial0 || type==initial1 || type==initial2);
}

inline Boolean isInternalExon(StateType type){
    return (type==internal0 || type==internal1 || type==internal2);
}

inline Boolean isRTerminalExon(StateType type){
    return (type==rterminal0 || type==rterminal1 || type==rterminal2);
}

inline Boolean isRInternalExon(StateType type){
    return (type==rinternal0 || type==rinternal1 || type==rinternal2);
}

inline Boolean isFirstExon(StateType type){
    return (isInitialExon(type) || isRTerminalExon(type) || type==singleG || type==rsingleG);
}

inline Boolean isLastExon(StateType type){
    return (type==terminal || type==rinitial || type==singleG || type == rsingleG);
}

inline Boolean is5UTR(StateType type){
    return (utr5single <= type && type <= utr5term)
	|| (rutr5single <= type && type <= rutr5term);
}

inline Boolean is5UTRIntron(StateType type){
    return (type == utr5intron || type == utr5intronvar 
	    || type == rutr5intron || type == rutr5intronvar);
}

inline Boolean is5UTRExon(StateType type){
    return is5UTR(type) && !is5UTRIntron(type);
}

inline Boolean is3UTR(StateType type){
    return  (utr3single <= type && type <= utr3term)
	|| (rutr3single <= type && type <= rutr3term);
}

inline Boolean is3UTRIntron(StateType type){
    return (type == utr3intron || type == utr3intronvar 
	    || type == rutr3intron || type == rutr3intronvar);
}

inline Boolean is3UTRExon(StateType type){
    return is3UTR(type) && !is3UTRIntron(type);
}

inline Boolean isFirstUTRExon(StateType type){
    return (type == utr5single || type == utr5init || type == rutr3single || type == rutr3term);
}

inline Boolean isLastUTRExon(StateType type){
    return (type == utr3single || type == utr3term || type == rutr5single || type == rutr5init);
}

inline Boolean isLongUTRIntron(StateType type){
    return (type == utr5intronvar || type == utr3intronvar || type == rutr5intronvar || type == rutr3intronvar);
}

inline Boolean isCodingExon(StateType type) {
    return  (singleG <= type && type <= terminal) 
	|| (rsingleG <= type && type <= rterminal2);
}

inline Boolean isNcIntron(StateType type) {
    return  (type == ncintron || type == ncintronvar || type == rncintron || type == rncintronvar);
}

inline Boolean isNcExon(StateType type) {
    return  (ncsingle <= type && !isNcIntron(type));
}

inline Boolean isNc(StateType type) {
    return  (ncsingle <= type);
}

inline Boolean isExon(StateType type) {
    return (isCodingExon(type) || is5UTRExon(type) || is3UTRExon(type) || isNcExon(type));
}

inline Boolean isCodingIntron(StateType type) {
    return  (lessD0 <= type && type <= longass2)
	|| (rlessD0 <= type && type <= rlongass2);
}

inline Boolean isIntron(StateType type) {
    return (isCodingIntron(type) || is5UTRIntron(type) || is3UTRIntron(type) 
	    || type == intron_type || type == rintron_type || isNcIntron(type));
}

inline Boolean isGeometricIntron(StateType type){
    return (type==geometric0 || type==geometric1 || type==geometric2);
}

inline Boolean isRGeometricIntron(StateType type){
    return (type==rgeometric0 || type==rgeometric1 || type==rgeometric2);
}

inline Boolean isLongAssIntron(StateType type){
    return (type==longass0 || type== longass1 || type==longass2);
}

inline Boolean isRLongDssIntron(StateType type){
    return (type==rlongdss0 || type== rlongdss1 || type==rlongdss2);
}

inline Boolean isLongDssIntron(StateType type){
    return (type==longdss0 || type== longdss1 || type==longdss2);
}

inline Boolean isOnFStrand(StateType type){
    return ((type >= singleG && type < rsingleG) || (type >= ncsingle && type < rncsingle));
}

inline StateType initialExon(int frame){
    frame = mod3(frame);
    if (frame == 0)
	return initial0;
    if (frame == 1)
	return initial1;
    if (frame == 2)
	return initial2;
    return TYPE_UNKNOWN;
}

inline StateType internalExon(int frame){
    frame = mod3(frame);
    if (frame == 0)
	return internal0;
    if (frame == 1)
	return internal1;
    if (frame == 2)
	return internal2;
    return TYPE_UNKNOWN;
}

inline StateType lessDIntron(int frame){
    frame = mod3(frame);
   if (frame == 0)
	return lessD0;
    if (frame == 1)
	return lessD1;
    if (frame == 2)
	return lessD2;
    return TYPE_UNKNOWN;
}


inline StateType geometricIntron(int frame){
    frame = mod3(frame);
    if (frame == 0)
	return geometric0;
    if (frame == 1)
	return geometric1;
    if (frame == 2)
	return geometric2;
    return TYPE_UNKNOWN;
}

inline StateType longdssIntron(int frame){
    frame = mod3(frame);
    if (frame == 0)
	return longdss0;
    if (frame == 1)
	return longdss1;
    if (frame == 2)
	return longdss2;
    return TYPE_UNKNOWN;
}

inline StateType longassIntron(int frame){
    frame = mod3(frame);
    if (frame == 0)
	return longass0;
    if (frame == 1)
	return longass1;
    if (frame == 2)
	return longass2;
    return TYPE_UNKNOWN;
}

inline StateType equalDIntron(int frame){
    frame = mod3(frame);
    if (frame == 0)
	return equalD0;
    if (frame == 1)
	return equalD1;
    if (frame == 2)
	return equalD2;
    return TYPE_UNKNOWN;
}

inline StateType rterminalExon(int frame){
    frame = mod3(frame);
    if (frame == 0)
	return rterminal0;
    if (frame == 1)
	return rterminal1;
    if (frame == 2)
	return rterminal2;
    return TYPE_UNKNOWN;
}

inline StateType rinternalExon(int frame){
    frame = mod3(frame);
    if (frame == 0)
	return rinternal0;
    if (frame == 1)
	return rinternal1;
    if (frame == 2)
	return rinternal2;
    return TYPE_UNKNOWN;
}
inline StateType rlongdssIntron(int frame){
    frame = mod3(frame);
    if (frame == 0)
	return rlongdss0;
    if (frame == 1)
	return rlongdss1;
    if (frame == 2)
	return rlongdss2;
    return TYPE_UNKNOWN;
}

inline StateType rlongassIntron(int frame){
    frame = mod3(frame);
    if (frame == 0)
	return rlongass0;
    if (frame == 1)
	return rlongass1;
    if (frame == 2)
	return rlongass2;
    return TYPE_UNKNOWN;
}

inline StateType rlessDIntron(int frame){
    frame = mod3(frame);
   if (frame == 0)
	return rlessD0;
    if (frame == 1)
	return rlessD1;
    if (frame == 2)
	return rlessD2;
    return TYPE_UNKNOWN;
}

inline StateType rgeometricIntron(int frame){
    frame = mod3(frame);
    if (frame == 0)
	return rgeometric0;
    if (frame == 1)
	return rgeometric1;
    if (frame == 2)
	return rgeometric2;
    return TYPE_UNKNOWN;
}

inline StateType requalDIntron(int frame){
    frame = mod3(frame);
    if (frame == 0)
	return requalD0;
    if (frame == 1)
	return requalD1;
    if (frame == 2)
	return requalD2;
    return TYPE_UNKNOWN;
}

int howOftenOccursIt(const char* haystack, const char* needle, const char* endhaystack=NULL);
bool containsJustNonNucs(const char *dna, int dnalen);
bool isNuc(const char *dna);

/*
 * expand the ~ to the $Home directory
 */
string expandHome(string filename);

/*
 * convert int to string
 */
inline string itoa(int n) {
    ostringstream strm;
    strm << n;
    return strm.str();
}

/*
 * convert double to string
 */
inline string ftoa(double n) {
    ostringstream strm;
    strm << n;
    return strm.str();
}

/*
 * copy string to newly allocated character array
 */
inline char* newstrcpy(const char* s, int len) {
    char* result = new char[len + 1];
    result[len] = '\0';
    return strncpy(result, s, len);
}
inline char* newstrcpy(const char* s) {
    int len = strlen(s);
    char* result = new char[len + 1];
    return strcpy(result, s);
}

inline char* newstrcpy(const string& s, int len) {
    char* result = new char[len + 1];
    s.copy(result, len);
    result[len] = '\0';
    return result;
}
inline char* newstrcpy(const string& s) {
    return newstrcpy(s, s.length());
}

inline void strip_newlines(char* p) {
    char* q = p;
    while (p != 0 && *p != '\0') {
        if (*p == '\n') {
            p++;
            *q = *p;
        } 
        else {
            *q++ = *p++;
        }
    }
    *q = '\0';
}

char *getRandomDNA(int len);

/* for 0 <= q <=1 get the q-th quantile of the values store in
 * vector v. 
 */
Double quantile(const vector<Double> &v, float q);
int quantile(const vector<int> &v, float q);

/* obtain a "hashtable" from a list of unique names, e.g.
 * [human, mouse, dog] becomes
 * human => 0, mouse => 1, dog => 2
 * Mario: This could become an unordered_map (expected constant time) once we commit to C++11,
 * rather than a logarithmic tree data structure
 */
map<string, size_t> *getMap (vector<string> names) throw(ProjectError);

/*
 * functions used in earlier versions

 *
 * sort an array increasingly p from index a to index b (included).
 */
// void QuickSort(Double *p, int a, int b);

/*
 * compute the unbiased empirical standard deviation of a distribution given by an array of doubles of size n
 * a distribution on the numbers {0,1, ..., n-1}
 */
// double variance(double *p, int n);

/*
 * compute the chi-square statistics for 
 * a multinomial distribution on the numbers {0,1, ..., r-1} given by p
 * reference: uniform distribution
 */
// double chiSquareUniform(int *p, int n);

// double chisquare(int a[4][4]);
// int branchPointPosition(char *seq, int asspos);

#endif   //  _TYPES_HH
