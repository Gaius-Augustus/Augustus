/**********************************************************************
 * file:    types.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  
 * authors: Mario Stanke (mario.stanke@uni-greifswald.de)
 * 
 * date    |   author      |  changes 
 * --------|---------------|------------------------------------------ 
 * 16.06.02| Mario Stanke  | Creation of the file
 * 02.01.03| Mario Stanke  | changing intron to model with 4 states
 * 24.03.03| Mario Stanke  | introducing shadow states
 * 06.08.03| Mario Stanke  | configPath
 * 04.10.05| Mario Stanke  | UTR adjustments
 * 15.05.06| Mario Stanke  | added 3' UTR
 **********************************************************************/

#include "types.hh"

// project includes
#include "properties.hh"  // for Constant::init

// standard C/C++ includes
#include <iostream>
#include <vector>
#include <cstdlib>

// declaration of "global" variables
//---------------------------------------
string Constant::configPath;
string Constant::speciesDir;
int Constant::trans_init_window = 12;
int Constant::ass_upwindow_size = 20;
int Constant::init_coding_len = 16;
int Constant::et_coding_len = 5;
int Constant::ass_start = 2;
int Constant::ass_end = 2;
int Constant::dss_start = 2;
int Constant::dss_end = 5;
int Constant::dss_maxbinsize = 0;
int Constant::ass_maxbinsize = 0;
int Constant::tis_maxbinsize = 0;
int Constant::tss_upwindow_size = 0;
int Constant::tss_start = 5;
int Constant::decomp_num_at = 1;
int Constant::decomp_num_gc = 1;
int Constant::decomp_num_steps = 1;
int Constant::min_coding_len = 102;
int Constant::max_exon_len = 12000;
Integer Constant::d_polyasig_cleavage = 20;
bool Constant::keep_viterbi = false;
double Constant::gc_range_min = 0.32;
double Constant::gc_range_max = 0.73;
double Constant::probNinCoding = 0.23;
double Constant::opalprob = 0.333;  // default value, was 0.24
double Constant::amberprob = 0.333; // default value, was 0.48
double Constant::ochreprob = 0.333; // default value, was 0.28
bool Constant::utr_option_on = false;
bool Constant::nc_option_on = false;
Integer Constant::augustus_verbosity = 1;
bool Constant::alternatives_from_evidence = true;
double Constant::subopt_transcript_threshold = 0.94;
Integer Constant::almost_identical_maxdiff = 10;
bool Constant::uniqueGeneId = false;
// class GeneticCode is taking care of these now
// bool OpenReadingFrame::ochre = true;
// bool OpenReadingFrame::opal = true;
// bool OpenReadingFrame::amber = true;
double Constant::max_contra_supp_ratio = 9.0;
bool Constant::reportUtrOnlyGenes = false;
bool Constant::useCRFtraining = false;
bool Constant::CRFtrainCDS = true;
bool Constant::CRFtrainIgenic = true;
bool Constant::CRFtrainIntron = true;
bool Constant::CRFtrainSS = true;
bool Constant::CRFtrainTIS = true;
bool Constant::CRFtrainUTR = false;
bool Constant::dss_gc_allowed = false;
Boolean Constant::tieIgenicIntron = true;
Boolean Constant::proteinOutput = true;
Boolean Constant::codSeqOutput = false;
Boolean Constant::contentmodels = true;
Boolean Constant::exoncands = false;
Integer Constant::min_intron_len = 39;
bool Constant::MultSpeciesMode = 0; // whether we do comparative gene prediction in multiple species
string Constant::treefile; // file name in which a tree is specified in Newick format
string Constant::speciesfilenames; // file name to file which contains the names of species and the corresponding file names
string Constant::dbaccess; // comma separated string with database access (database name,host name,user,passwd, table name"
string Constant::alnfile; // name of file that contains MSA of genomes
bool Constant::overlapmode = false;
string Constant::orthoexons; //name of file that contains list of orthologous exons
Integer Constant::maxOvlp = 60; // maximum overlap of coding regions for bacteria
vector<Double> Constant::head2tail_ovlp;
vector<Double> Constant::head2head_ovlp;
vector<Double> Constant::tail2tail_ovlp;
unsigned Constant::temperature = 0; // heating the distribution for sampling, 0=cold, 7=hottest
bool Constant::softmasking = false;
bool inCRFTraining = false;
bool Constant::dbhints = false;
// scores from logistic regression                                                                                                  
bool Constant::logreg;
// features are explaned in the default config file config/cgp/log_reg_parameters_default.cfg
vector<double>Constant::ex_sc;
vector<double>Constant::in_sc;
vector<double>Constant::lg_es;

// moved here from hints.cc
const int power2id[31] = {1,2,4,8,16,32,64,128,
			   256,512,1024,2048,4096,8192,16384,32768,
			   65536,131072,262144,524288,
			   1048576,2097152,4194304,8388608,
			   16777216,33554432,67108864,134217728,
			   268435456,536870912,1073741824};

const char* stateTypeNames[NUM_TYPES]= 
{"intergenic region", 
 "single exon", 
 "initial exon 0", "initial exon 1", "initial exon 2", 
 "internal exon 0", "internal exon 1", "internal exon 2", 
 "terminal exon",
 "intron lessD0", "intron dss 0", "intron equalD0", "intron geometric 0", "intron ass 0",
 "intron lessD1", "intron dss 1", "intron equalD1", "intron geometric 1", "intron ass 1",
 "intron lessD2", "intron dss 2", "intron equalD2", "intron geometric 2", "intron ass 2",
 "5'UTR single exon", "5'UTR initial exon", "5'UTR intron", "5'UTR intron var", "5'UTR internal exon", "5'UTR terminal exon",
 "3'UTR single exon", "3'UTR initial exon", "3'UTR intron", "3'UTR intron var", "3'UTR internal exon", "3'UTR terminal exon",
 "reverse single exon", 
 "reverse initial exon", 
 "reverse internal exon 0", "reverse internal exon 1", "reverse internal exon 2", 
 "reverse terminal exon 0","reverse terminal exon 1","reverse terminal exon 2",
 "reverse intron lessD0", "reverse long intron dss 0", "reverse intron equalD0", 
"reverse intron geometric 0", "reverse long intron ass 0",
 "reverse intron lessD1", "reverse long intron dss 1", "reverse intron equalD1", 
"reverse intron geometric 1", "reverse long intron ass 1",
 "reverse intron lessD2", "reverse long intron dss 2", "reverse intron equalD2", 
"reverse intron geometric 2", "reverse long intron ass 2",
"reverse 5'UTR single exon", "reverse 5'UTR initial exon", "reverse 5'UTR intron", "reverse 5'UTR intron var", "reverse 5'UTR internal exon", "reverse 5'UTR terminal exon", 
"reverse 3'UTR single exon", "reverse 3'UTR initial exon", "reverse 3'UTR intron", "reverse 3'UTR intron var", "reverse 3'UTR internal exon", "reverse 3'UTR terminal exon", 
 "intron", "reverse intron", "exon",
 "noncoding single exon", "noncoding initial exon", "noncoding intron", "noncoding intron var", "noncoding internal exon", "noncoding terminal exon",
 "reverse noncoding single exon", "reverse noncoding initial exon", "reverse noncoding intron", "reverse noncoding intron var", "reverse noncoding internal exon", "reverse noncoding terminal exon"
};

const char* stateTypeIdentifiers[NUM_TYPES]= 
{"igenic",
 "single", "initial0", "initial1", "initial2", "internal0", "internal1", "internal2", "terminal",
 "lessD0", "longdss0", "equalD0", "geometric0", "longass0",
 "lessD1", "longdss1", "equalD1", "geometric1", "longass1",
 "lessD2", "longdss2", "equalD2", "geometric2", "longass2",
 "utr5single", "utr5init", "utr5intron", "utr5intronvar", "utr5internal", "utr5term",
 "utr3single", "utr3init", "utr3intron", "utr3intronvar", "utr3internal", "utr3term",
 "rsingle", "rinitial", "rinternal0", "rinternal1", "rinternal2", "rterminal0", "rterminal1", "rterminal2",
 "rlessD0", "rlongdss0", "requalD0", "rgeometric0", "rlongass0",
 "rlessD1", "rlongdss1", "requalD1", "rgeometric1", "rlongass1",
 "rlessD2", "rlongdss2", "requalD2", "rgeometric2", "rlongass2",
 "rutr5single", "rutr5init", "rutr5intron", "rutr5intronvar", "rutr5internal", "rutr5term",
 "rutr3single", "rutr3init", "rutr3intron", "rutr3intronvar", "rutr3internal", "rutr3term",
 "intron", "rintron", "exon",
 "ncsingle", "ncinit", "ncintron", "ncintronvar", "ncinternal", "ncterm",
 "rncsingle", "rncinit", "rncintron", "rncintronvar", "rncinternal", "rncterm"
};

// don't forget to change shift the reading frames, when you introduce new states!!
const int stateReadingFrames[NUM_TYPES]=
{0,                        // intergenic region
 0,0,1,2,0,1,2,0,          // forward exons
 0,0,0,0,0,                // --  
 1,1,1,1,1,                //   |- forward introns
 2,2,2,2,2,                // -- 
 0,0,0,0,0,0,0,0,0,0,0,0,  // forward utr
 2,2,0,1,2,0,1,2,          // reverse exons
 0,0,0,0,0,                // --
 1,1,1,1,1,                //   |- reverse introns
 2,2,2,2,2,                // --
 0,0,0,0,0,0,0,0,0,0,0,0,  // reverse utr
 0,0,0,                    // other
 0,0,0,0,0,0,0,0,0,0,0,0   // noncoding
};    

char strandChar (Strand s){
    switch (s){
    case plusstrand:
	return '+';
    case minusstrand:
	return '-';
    case bothstrands:
	return '.';
    default:
	return '?';
    }
}

ostream& operator<< (ostream& strm, const Strand s){
    strm << strandChar(s);
    return strm;
}

void Constant::init(){
    configPath = Properties::getProperty(CFGPATH_KEY);
    speciesDir = Properties::getProperty(SPECIESDIR_KEY);
    try {
	trans_init_window = Properties::getIntProperty("/Constant/trans_init_window");
    } catch (ProjectError e) {
	cerr << e.getMessage();
    }
    try {
	ass_upwindow_size = Properties::getIntProperty("/Constant/ass_upwindow_size");
    } catch (ProjectError e) {
	cerr << e.getMessage();
    }
    try {
	ass_start = Properties::getIntProperty("/Constant/ass_start");
    } catch (ProjectError e) {
	cerr << e.getMessage();
    }
    try {
	ass_end = Properties::getIntProperty("/Constant/ass_end");
    } catch (ProjectError e) {
	cerr << e.getMessage();
    }
    try {
	dss_start = Properties::getIntProperty("/Constant/dss_start");
    } catch (ProjectError e) {
	cerr << e.getMessage();
    }
    try {
	dss_end = Properties::getIntProperty("/Constant/dss_end");
    } catch (ProjectError e) {
	cerr << e.getMessage();
    } 
    try {
	useCRFtraining = Properties::getBoolProperty("CRF");
    } catch (ProjectError e) {}
    try {
	dss_maxbinsize = Properties::getIntProperty("/Constant/dss_maxbinsize");
    } catch (ProjectError e) {
	if (useCRFtraining)
	    dss_maxbinsize = 100; // default, when not set and CRF training
    }
    try {
      CRFtrainSS = Properties::getBoolProperty("CRFtrainSS");
    } catch (ProjectError e) {
      CRFtrainSS = true;
    }
    try {
      CRFtrainIntron = Properties::getBoolProperty("CRFtrainIntron");
    } catch (ProjectError e) {
      CRFtrainIntron = true;
    }
    try {
      CRFtrainTIS = Properties::getBoolProperty("CRFtrainTIS");
    } catch (ProjectError e) {
      CRFtrainTIS = true;
    }
    try {
      CRFtrainCDS = Properties::getBoolProperty("CRFtrainCDS");
    } catch (ProjectError e) {
      CRFtrainCDS = true;
    }
    try {
      CRFtrainUTR = Properties::getBoolProperty("CRFtrainUTR");
    } catch (ProjectError e) {
      CRFtrainUTR = false;
    }
    try {
      CRFtrainIgenic = Properties::getBoolProperty("CRFtrainIgenic");
    } catch (ProjectError e) {
      CRFtrainIgenic = true;
    }
    try {
        ass_maxbinsize = Properties::getIntProperty("/Constant/ass_maxbinsize");
    } catch (ProjectError e) {
	if (useCRFtraining)
	    ass_maxbinsize = 100; // default, when not set and CRF training
    }
    try {
        tis_maxbinsize = Properties::getIntProperty("/Constant/tis_maxbinsize");
    } catch (ProjectError e) {
	if (useCRFtraining)
	    tis_maxbinsize = 30; // default, when not set and CRF training
    }
    try {
	init_coding_len = Properties::getIntProperty("/Constant/init_coding_len");
    } catch (ProjectError e) {
	cerr << e.getMessage();
    }
    try {
	et_coding_len = Properties::getIntProperty("/Constant/intterm_coding_len");
    } catch (ProjectError e) {
	cerr << e.getMessage();
    }
    try {
	tss_upwindow_size = Properties::getIntProperty("/Constant/tss_upwindow_size");
    } catch (ProjectError e) {
	//cerr << e.getMessage();
    }
    try {
	tss_start = Properties::getIntProperty("/UtrModel/tss_start");
    } catch (ProjectError e) {
      //cerr << e.getMessage();  TODO: error message only if --UTR=5
    }
    try {
	decomp_num_at = Properties::getIntProperty("/Constant/decomp_num_at");
    } catch (ProjectError e) {
	cerr << e.getMessage();
    }
    try {
	decomp_num_gc = Properties::getIntProperty("/Constant/decomp_num_gc");
    } catch (ProjectError e) {
	cerr << e.getMessage();
    }
    try {
	decomp_num_steps = Properties::getIntProperty("/Constant/decomp_num_steps");
    } catch (ProjectError e) {
	cerr << e.getMessage();
    } 
    try {
	min_coding_len = Properties::getIntProperty("/Constant/min_coding_len");
    } catch (ProjectError e) {
	cerr << e.getMessage();
    }
    try {
	max_exon_len = Properties::getIntProperty("/ExonModel/maxexonlength");
    } catch (ProjectError e) {
	cerr << e.getMessage();
    }
    
    // scores for logistic regression
    try {
	logreg = Properties::getBoolProperty("/CompPred/logreg");
    } catch (...) {
	logreg = true;
    }
    for(int i=0; i<16; i++){
	try {
	    ex_sc.push_back(Properties::getdoubleProperty("/CompPred/exon_score" + itoa(i) ));
	} catch (...) {
	    ex_sc.push_back(0);
	}
    }
    for(int i=0; i<4; i++){
	try {
	    in_sc.push_back(Properties::getdoubleProperty("/CompPred/intron_score" + itoa(i) ));
	} catch (...) {
	    in_sc.push_back(0);
	}
    }
    for(int i=0; i<4; i++){
	try {
	    lg_es.push_back(Properties::getdoubleProperty("lg_exon_score" + itoa(i) ));
	} catch (...) {
	    lg_es.push_back(0);
	}
    }
    
    Properties::assignProperty("/UtrModel/d_polyasig_cleavage", d_polyasig_cleavage);
    Properties::assignProperty("keep_viterbi", 	keep_viterbi);
    Properties::assignProperty("/Constant/gc_range_min", gc_range_min);
    Properties::assignProperty("/Constant/gc_range_max", gc_range_max);
    Properties::assignProperty("/Constant/probNinCoding", probNinCoding);

    // for default values, see above
    Properties::assignProperty("/Constant/opalprob", opalprob);
    Properties::assignProperty("/Constant/amberprob", amberprob);
    Properties::assignProperty("/Constant/ochreprob", ochreprob);

    Properties::assignProperty("/augustus/verbosity", augustus_verbosity);
    Properties::assignProperty("alternatives-from-evidence", alternatives_from_evidence);
    Properties::assignProperty("/Constant/subopt_transcript_threshold", subopt_transcript_threshold);
    Properties::assignProperty("/Constant/almost_identical_maxdiff", almost_identical_maxdiff);
    Properties::assignProperty("uniqueGeneId", uniqueGeneId);
    Properties::assignProperty("/Constant/max_contra_supp_ratio", max_contra_supp_ratio);
    Properties::assignProperty("reportUtrOnlyGenes", reportUtrOnlyGenes);
    Properties::assignProperty(UTR_KEY, utr_option_on);
    Properties::assignProperty(NONCODING_KEY, nc_option_on);
    Properties::assignProperty("/IntronModel/allow_dss_consensus_gc", dss_gc_allowed);
    Properties::assignProperty("tieIgenicIntron", tieIgenicIntron);
    Properties::assignProperty("protein", proteinOutput);
    Properties::assignProperty("codingseq", codSeqOutput);
    Properties::assignProperty("contentmodels", contentmodels);
    Properties::assignProperty("exoncands", exoncands);
    Properties::assignProperty("min_intron_len", min_intron_len);
    Properties::assignProperty("orthoexons", orthoexons);
    Properties::assignProperty("maxOvlp", maxOvlp);
    Properties::assignProperty("temperature", temperature);
    if (temperature > 7){
	cerr << "No temperature >7 allowed. temperature must be one of 0 1 2 3 4 5 6 7. Will use temperature=7." << endl;
        temperature = 7;
    }     	
    LLDouble::setTemperature(temperature);
 
}

int howOftenOccursIt(const char* haystack, const char* needle, const char *endhaystack){
    if (!haystack)
	return 0;
    int n=-1;
    if (endhaystack == NULL)
	endhaystack = haystack + strlen(haystack)-1;
    const char* pos = haystack;
    while (pos && pos <= endhaystack) {
	pos = strstr(pos, needle);
	if (pos)
	    pos++;
	n++;
    }
    return n;
}

bool containsJustNonNucs(const char *dna, int dnalen){
    bool haveNuc=false;
    int i=0;
    while (!haveNuc && i<dnalen) {
	char c = toupper(dna[i]);
	if(c == 'A' || c == 'C' || c == 'G' || c == 'T')
	    haveNuc = true;
	i++;
    }
    return !haveNuc;
}

bool isNuc(const char *dna){
    char c = toupper(*dna);
    return (c == 'A' || c == 'C' || c == 'G' || c == 'T');
}

/*
 * expand the ~ to the $Home directory
 */
string expandHome(string filename){
    if (filename.length()>0 && filename[0]=='~') {
	filename.replace(0,1,getenv("HOME"));
    }
    return filename;
}


char *getRandomDNA(int len) {
  char dnamap[5] = "acgt"; // 0->a, 1->c, 2->g, 3->t
  char *seq = new char[len+1];
  for (int i=0; i<len; i++) {
    seq[i] = dnamap[(int) (4.0 * rand() / (1.0 + RAND_MAX))];
  }
  return seq;
}

Double quantile(const vector<Double> &v, float q){
  if (v.size() == 0)
    return -numeric_limits<double>::infinity();
  if (q<0)
    q=0;
  if (q>1)
    q=1;
  vector<Double> w(v);
  std::sort(w.begin(), w.end()); // O(n log n) but finding a quantile would be possible also in linear time
  int threshindex = (int) (q * w.size());
  // cout << "quantile: size=" << w.size() << " index=" << threshindex << endl;
  return w[threshindex];
}

int quantile(const vector<int> &v, float q){
    if (v.size() == 0)
	return -numeric_limits<int>::infinity();
    if (q<0)
	q=0;
    if (q>1)
	q=1;
    vector<int> w(v);
    std::sort(w.begin(), w.end()); // O(n log n) but finding a quantile would be possible also in linear time                                                                                                                                                                            
    int threshindex = (int) (q * w.size());
    return w[threshindex];
}                                                                                                                                                                                                                                                                                        
 

map<string, size_t> *getMap (vector<string> names) throw(ProjectError) {
    map<string, size_t> *hashtable = new map<string, size_t>;
    for (size_t i=0; i<names.size(); i++){
	if (hashtable->find(names[i]) != hashtable->end())
	    throw ProjectError(string("getMap: list of names contains multiple entries: ") + names[i]);
	hashtable->insert(pair<string, size_t>(names[i], i));
    }
    return hashtable;
}

/*
 * functions used in earlier versions of AUGUSTUS
 */
// void QuickSort(Double *p, int a, int b){
//     if (a >= b)
// 	return;
//     Double pivot = p[a], temp;
//     bool partitioned = false;
//     int i=a-1, j=b+1;
//     while(!partitioned) {
// 	do {j--;} while(p[j] > pivot);
// 	do {i++;} while(p[i] < pivot);
// 	if (i<j) {      // exchange p[i] and p[j]
// 	    temp = p[i];
// 	    p[i] = p[j];
// 	    p[j] = temp;
// 	} else 
// 	    partitioned = true;
//     }
//     QuickSort(p, a, j);
//     QuickSort(p, j+1, b);
// }

// double variance(double *p, int n){
//     double sumsq=0, sum=0;
//     for (int i=0; i<n; i++) {
// 	sum += p[i] * i;
// 	sumsq += p[i] * i * i;
//     }
//     return sqrt((sumsq - sum*sum/n)/(n-1));
// }

// double chiSquareUniform(int *p, int r){
//     int n =0;
//     double sum=0;
//     for (int i=0; i<r; i++)
// 	n += p[i];
//     for (int i=0; i<r; i++){
// 	double expect = (double) n/r;
// 	sum += (expect-p[i])*(expect-p[i])/expect;
//     }    
//     return sum;
// }


// double chisquare(int a[4][4]) {
//     double chisq(0.0);
//     int margin1[4], margin2[4];
//     int n=0;
//     double diff, expect;
//     for (int i = 0; i<4 ; i++) {
// 	margin1[i] = a[i][0] + a[i][1] + a[i][2] + a[i][3];
// 	margin2[i] = a[0][i] + a[1][i] + a[2][i] + a[3][i];
// 	/*cout << "margin1[" << i << "]=" << margin1[i] 
// 	  << "  margin2[" << i << "]=" << margin2[i] << endl; */
// 	n += margin1[i];
//     }
//     //cerr << "n=" << n << endl;
//     for (int i=0; i<4; i++)
// 	for (int j=0; j<4; j++){
// 	    expect = (double) margin1[i]*margin2[j]/n;
// 	    diff = expect - a[i][j];
// 	    if (a[i][j] >= 5) {
// 		//cout << "beitrag[" << i << "][" << j <<"] = " << diff*diff/expect << endl;
// 		chisq += diff*diff/expect;
// 	    }
// 	}
//     return chisq;
// }
