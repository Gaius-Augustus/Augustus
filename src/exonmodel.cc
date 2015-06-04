/**********************************************************************
 * file:    exonmodel.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  coding exons only, UTR exons are modelled in utrmodel.cc
 * authors: Mario Stanke, mario@gobics.de
 *
 * date    |   author      |  changes
 * --------|---------------|------------------------------------------
 * 06.05.02| Mario Stanke  | debugging and testing
 * 18.05.02| Mario Stanke  | added length distributions
 * 31.05.02| Mario Stanke  | implementation of viterbi algorithm
 * 06.10.02| Mario Stanke  | making length distr independent of content
 * 18.11.02| Mario Stanke  | internal exon terminal emission part
 * 30.12.02| Mario Stanke  | moving splice site models to intron
 * 25.03.03| Mario Stanke  | shadow exons on reverse strand
 * 02.07.04| Mario Stanke  | malus for all exons, not only for those without a bonus
 * 19.08.04| Mario Stanke  | corrected rare bug with exon state at end of dna
 * 23.12.04| Mario Stanke  | bug fixed: wrong reading frame in call to getExonListInRange
 * 08.02.05| Mario Stanke  | bug fixed: exon hints were missed for reverse terminal or reverse single exons.
 * 22.08.05| Mario Stanke  | training on genes with incomplete 3' end, incomplete exon not used for lendist
 * 09.09.08| Mario Stanke  | introduce fuzzyness in dss and ass hints
 **********************************************************************/

#include "exonmodel.hh"

// project includes
#include "motif.hh"
#include "pp_scoring.hh"
#include "properties.hh"
#include "projectio.hh"
#include "extrinsicinfo.hh"

// standard C/C++ includes
#include <fstream>

/*
 * Initialisation of static data members
 */

vector<Integer> ExonModel::patterncount[3];      // {0,1,2}x{acgt}^(k+1), the reading frame is 
                                                 // the positiion of the emitted (last) nucleotide
vector<Integer> ExonModel::initpatterncount[3];
vector<Integer> ExonModel::etpatterncount[3];
Integer ExonModel::k = 4;
Integer ExonModel::etorder = 1;
Integer ExonModel::etpseudocount = 3;
Integer ExonModel::min_exon_length = 1;
Integer ExonModel::trans_init_window = 12;
FramedPatMMGroup ExonModel::emiprobs("exon emiprob old");
FramedPatMMGroup* ExonModel::GCemiprobs = NULL;
vector<Double> ExonModel::initemiprobs[3];
vector<Double> **ExonModel::GCinitemiprobs = NULL;
vector<Double> ExonModel::etemiprobs[3];
vector<Double> **ExonModel::GCetemiprobs = NULL;
vector<Integer> ExonModel::numExonsOfType;
vector<Integer> ExonModel::numHugeExonsOfType;  // number of exons exceeding the maximal length 
                                                    // modelled by the length distribution
vector<Double> ExonModel::lenDistSingle;   // Length distribution of Single exons 
                                               // (length of biol. ExonModel::exon)
vector<Double> ExonModel::lenDistInitial;  // Length distribution of Initial exons 
                                               // (length of biol. exon)
vector<Double> ExonModel::lenDistInternal; // Length distribution of Internal exons 
                                               // (length of biol. exon)
vector<Double> ExonModel::lenDistTerminal; // Length distribution of Terminal exons 
                                               // (length of biol. exon)
Motif*   ExonModel::transInitMotif = NULL;      // weight matrix before the translation initiation
Motif*   ExonModel::GCtransInitMotif = NULL;    // weight matrix before the translation initiation
BinnedMMGroup  ExonModel::transInitBinProbs("tis", 1); // CRF-features based on transInitMotif, monotonic increasing
BinnedMMGroup* ExonModel::GCtransInitBinProbs = NULL;  // for all GC content classes
Motif**  ExonModel::etMotif = NULL;             // weight matrices before the donor splice site (3 frames)
Motif*** ExonModel::GCetMotif = NULL;           // array of above
Integer  ExonModel::numSingle=0, ExonModel::numInitial=0, ExonModel::numInternal=0, 
         ExonModel::numTerminal=0;
Integer  ExonModel::numHugeSingle=0, ExonModel::numHugeInitial=0, 
         ExonModel::numHugeInternal=0, ExonModel::numHugeTerminal=0; 
Matrix<vector<Double> > ExonModel::Pls;
Matrix<vector<Double> >* ExonModel::GCPls = NULL;
Integer ExonModel::exoncount = 0;
Boolean ExonModel::hasLenDist = false;
//Boolean ExonModel::hasAAdep = false;      // not in use right now
Integer ExonModel::gesbasen[3]  = { 0, 0, 0 };
Double  ExonModel::patpseudo = 1;         // pseudocount for patterns in sequence
Integer ExonModel::exonLenD = 4000;       // use detailled length distribution up to this number
Integer ExonModel::minPatSum = 0;         // for the decision to shorten the emission pattern
vector<Integer> ExonModel::lenCountSingle;   // length count of Single exons (length of biol. exon)
vector<Integer> ExonModel::lenCountInitial;  // length count of Initial exons (length of biol. exon)
vector<Integer> ExonModel::lenCountInternal; // length count of Internal exons (length of biol. exon)
vector<Integer> ExonModel::lenCountTerminal; // length count of Terminal exons (length of biol. exon)
double ExonModel::slope_of_bandwidth = 0.1;  // for smoothing
Integer ExonModel::minwindowcount = 5;       // see class Smooth in commontrain.hh
Integer ExonModel::tis_motif_memory = 3;     // see class Smooth in commontrain.hh
Integer ExonModel::tis_motif_radius = 2;     // see class Smooth in commontrain.hh
list<string>* ExonModel::tiswins = NULL;     // holds translation initiation windows (for CRF training)

int             ExonModel::numModels = 1;
//BaumWelch*      ExonModel::baumwelch = NULL;
Double*         ExonModel::modelStartProbs = NULL;
//AADependency    ExonModel::aadep = 0;
int             ExonModel::ilend = 550;
OpenReadingFrame* ExonModel::orf = NULL;
int             ExonModel::ochrecount = 0; // frequencies of the stop codons
int             ExonModel::ambercount = 0;
int             ExonModel::opalcount  = 0;
bool            ExonModel::initAlgorithmsCalled = false;
bool            ExonModel::haveORF = false;
int             ExonModel::lastParIndex = -1; // GC-index of current parameter set
int             ExonModel::verbosity;
int             ExonModel::startcounts[64] = {0};
int             ExonModel::lenboostL = INT_MAX; // parameters for boosting length dist of single and initial exons
double          ExonModel::lenboostE = 0; // lengths above L are improved, the more the larger E is

/* --- OpenReadingFrame methods ------------------------------------ */

/*
 * constructor
 */
OpenReadingFrame::OpenReadingFrame(const char *dna, int _max_exon_length, int _n) :
    n(_n), max_exon_length(_max_exon_length)
{
    nearestStopForward.resize(n);
    nearestStopReverse.resize(n);
    int stopcodpos, i;
    stopcodpos=-1;
    for (i=0; i<=n - STOPCODON_LEN; i+=3) {
	if (isStopcodon(dna+i)) 
	    stopcodpos = i;
	nearestStopForward[i] = stopcodpos;
    }
    stopcodpos=-1;
    for (i=1; i<=n - STOPCODON_LEN; i+=3) {
	if (isStopcodon(dna+i))
	    stopcodpos = i;
	nearestStopForward[i] = stopcodpos;
    }
    stopcodpos=-1;
    for (i=2; i<=n - STOPCODON_LEN; i+=3) {
	if (isStopcodon(dna+i))
	    stopcodpos = i;
	nearestStopForward[i] = stopcodpos;
    }
    stopcodpos=-1;
    for (i=0; i<=n - STOPCODON_LEN; i+=3) {
	if (isRCStopcodon(dna+i))
	    stopcodpos = i;
	nearestStopReverse[i] = stopcodpos;
    }
    stopcodpos=-1;
    for (i=1; i<=n - STOPCODON_LEN; i+=3) {
	if (isRCStopcodon(dna+i))
	    stopcodpos = i;
	nearestStopReverse[i] = stopcodpos;
    }
    stopcodpos=-1;
    for (i=2; i<=n - STOPCODON_LEN; i+=3) {
	if (isRCStopcodon(dna+i))
	    stopcodpos = i;
	nearestStopReverse[i] = stopcodpos;
    }
    if (n>5) {
	// nearestStopForward[n - STOPCODON_LEN] = nearestStopForward[n - STOPCODON_LEN - 3];
	nearestStopForward[n - STOPCODON_LEN +1] = nearestStopForward[n - STOPCODON_LEN - 2];
	nearestStopForward[n - STOPCODON_LEN +2] = nearestStopForward[n - STOPCODON_LEN - 1];
	// nearestStopReverse[n - STOPCODON_LEN] = nearestStopReverse[n - STOPCODON_LEN - 3];
	nearestStopReverse[n - STOPCODON_LEN +1] = nearestStopReverse[n - STOPCODON_LEN - 2];
	nearestStopReverse[n - STOPCODON_LEN +2] = nearestStopReverse[n - STOPCODON_LEN - 1];
    } else {// TODO (unwichtig)
    }
}


/*
 * Get the position furthest to the left of base such that with respect to 
 * frame at position base there is no stop codon in the window from 
 * leftmostExonBegin to base. If forward is false the reverse complement is taken.
 * Initialize the tables of stopcodons.
 */
int OpenReadingFrame::leftmostExonBegin(int frame, int base, bool forward){
    int pos;
    if (forward) {
	if (frame==0 || frame==1)
	    pos=base-frame-3;
	else 
	    pos=base-frame;
    } else {
	if (frame==1 || frame==2)
	    pos=base+frame-5;
	else 
	    pos=base-2;
    }
    
    Integer leftmostbegin;
    if (pos >= n)
	pos -= 3*((pos-n+3)/3);
    // alternativ: pos = n-3 + (pos-n)%3

    if (pos >= 0) {
	leftmostbegin = forward? nearestStopForward[pos]+1 : nearestStopReverse[pos]+1;
    }
    else 
	leftmostbegin = 0;
   

    // this has the effect that no exon can be longer than max_exon_length
    // which is actually not a necessary condition but it is easier this way
    // The 10 is a security distance, because the biological exon can be longer than the inner sequence part 
    int max_allowed_len = max_exon_length - Constant::ass_upwindow_size - Constant::ass_start - ASS_MIDDLE - DSS_MIDDLE - Constant::dss_start;
    if (leftmostbegin < base-max_allowed_len)
      leftmostbegin = base-max_allowed_len;
    return leftmostbegin;
}

/*
 * containsInFrameStopcodon: not needed, can be implemented by GeneticCode
 */

bool OpenReadingFrame::isStopcodon( const char* dna ){
    if (*dna != tolower(*dna))
	throw ProjectError("Internal Error: upper case dna"); // TEMP: just for now
    // Oliver, 2009/06/24: commented this out to give translation table priority
    // over {ochre,amber,opal}prob values; this basically causes
    // stopcodons with prob=0 to appear neither prematurely nor as stopcodons
    // if( *dna == 't' ){
    // 	if( *(dna+1) == 'a' )
    // 	    return ((*(dna+2) == 'a' && ochre) || (*(dna+2) == 'g' && amber));
    // 	if( *(dna+1) == 'g' )
    // 	    return (*(dna+2) == 'a' && opal);
    // }
    // return false;

    // give the chosen translation table priority over {ochre,amber,opal}prob
    // 
    return GeneticCode::isStopcodon(dna);
}

bool OpenReadingFrame::isRCStopcodon( const char* dna ){
    return GeneticCode::isRCStopcodon(dna);
}


/* --- ExonModel methods ------------------------------------------- */

/*
 * constructor
 */
ExonModel::ExonModel() : gweight(1) {
    etype = toStateType( Properties::getProperty("/ExonModel/type", exoncount++) );
    win = stateReadingFrames[etype];
    switch (etype) {
	case singleG: case initial0: case initial1: case initial2:
	    // start codon at left end of exon
	    beginPartLen = STARTCODON_LEN + trans_init_window;
	    innerPartOffset = STARTCODON_LEN;
	    break;
	case rsingleG: case rterminal0: case rterminal1:case rterminal2:  
	    // reverse stop codon at left end of exon
	    beginPartLen = innerPartOffset = STOPCODON_LEN;
	    break;
	default:   // splice site at left end of exon
	    beginPartLen = 0;
	    innerPartOffset = isOnFStrand(etype) ? Constant::ass_end : Constant::dss_start;
    }
    baseOffset = getBaseOffset(etype);
    innerPartEndOffset = getInnerPartEndOffset(etype);
}


/*
 * getBaseOffset
 */
int ExonModel::getBaseOffset(StateType type){
    if (type == singleG || type == terminal) // stop codon at right end of exon
	return 0;
    
    if (type == rsingleG || type == rinitial) // reverse start codon at right end of exon
	return -trans_init_window; 
    
    // splice site at right end of exon
    return isOnFStrand(type) ? Constant::dss_start : Constant::ass_end;
}

/*
 * getInnerPartEndOffset
 */
int ExonModel::getInnerPartEndOffset(StateType type){
    if (type == singleG || type == terminal) // stop codon at right end of exon
	return STOPCODON_LEN;
    
    if (type == rsingleG || type == rinitial) // reverse start codon at right end of exon
	return STARTCODON_LEN; 
    
    // splice site at right end of exon
    return isOnFStrand(type) ? Constant::dss_start : Constant::ass_end;
}

/*
 * destructor
 */
ExonModel::~ExonModel( ){
    // TODO: delete what is neccessary
    if( --exoncount == 0 )
	lastParIndex = -1;
}

/*
 * ===[ ExonModel initialisation of class variables ]======================
 */
void ExonModel::init() {
    /*
     * Default values for the parameters.
     */
    k = 4;
    patpseudo  = Double(1);   
    slope_of_bandwidth = 0.1;   // for smoothing
    minwindowcount=3;           // see class Smooth in commontrain.hh	
    try{
	verbosity = Properties::getIntProperty("/ExonModel/verbosity");
    }catch( ProjectError e) {    
	cerr << e.getMessage();
    }
    try{
	k = Properties::getIntProperty( "/ExonModel/k" );
    }catch( ProjectError e) {    
	cerr << e.getMessage();
    }
    try{
	patpseudo = Properties::getDoubleProperty( "/ExonModel/patpseudocount" );
    }catch( ProjectError e) { 
	cerr << e.getMessage();
    }
    try{
	exonLenD = Properties::getIntProperty( "/ExonModel/exonlengthD" );
    }catch( ProjectError e) { 
	cerr << e.getMessage();
    }
    try{
	slope_of_bandwidth = Properties::getdoubleProperty( "/ExonModel/slope_of_bandwidth");
    }catch( ProjectError e) { 
	cerr << e.getMessage();
    }
    try{
	minwindowcount = Properties::getIntProperty( "/ExonModel/minwindowcount");
    }catch( ProjectError e) { 
	cerr << e.getMessage();
    }
    try{
	min_exon_length = Properties::getIntProperty( "/ExonModel/minexonlength");
    } catch( ProjectError e) { 
	// optional parameter
    }
    try{
	minPatSum = Properties::getIntProperty( "/ExonModel/minPatSum");
    } catch( ProjectError e) { 
	cerr << e.getMessage();
    }
    try{
	etorder = Properties::getIntProperty( "/ExonModel/etorder");
    } catch( ProjectError e) { 
	cerr << e.getMessage();
    }
    try{
	etpseudocount = Properties::getIntProperty( "/ExonModel/etpseudocount");
    } catch( ProjectError e) { 
	cerr << e.getMessage();
    }
    try{
	tis_motif_memory = Properties::getIntProperty( "/ExonModel/tis_motif_memory");
    } catch( ProjectError e) { 
	cerr << e.getMessage();
    }
    try{
	tis_motif_radius = Properties::getIntProperty( "/ExonModel/tis_motif_radius");
    } catch( ProjectError e) { 
	cerr << e.getMessage();
    }
    try{
	lenboostE = Properties::getdoubleProperty( "/ExonModel/lenboostE");
	if (lenboostE < 0.0){
	    cout << "Warning: /ExonModel/lenboostE < 0 does not make sense. Will use 0 instead." << endl;
	    lenboostE = 0.0;
	}
    } catch( ProjectError e) {}
    try{
	lenboostL = Properties::getIntProperty( "/ExonModel/lenboostL");
    } catch( ProjectError e) {}
	
    trans_init_window = Constant::trans_init_window;

    // reserve space for GC content dependent arrays IF NOT DONE BEFORE
    if (!GCPls)
      GCPls = new Matrix<vector<Double> >[Constant::decomp_num_steps];
    if (!GCtransInitMotif)
      GCtransInitMotif = new Motif[Constant::decomp_num_steps];
    if (!GCtransInitBinProbs){
      GCtransInitBinProbs = new BinnedMMGroup[Constant::decomp_num_steps];
      for (int idx = 0; idx < Constant::decomp_num_steps; idx++){
	GCtransInitBinProbs[idx].setName("tis");
	GCtransInitBinProbs[idx].monotonic = 1; // increasing
      }
    }
    if (!GCetMotif){
      GCetMotif = new Motif**[Constant::decomp_num_steps];
      for (int idx = 0; idx < Constant::decomp_num_steps; idx++){
	GCetMotif[idx] = new Motif*[3];
	GCetMotif[idx][0] = new Motif();
	GCetMotif[idx][1] = new Motif();
	GCetMotif[idx][2] = new Motif();
      }
    }
    if (!GCemiprobs)
      GCemiprobs = new FramedPatMMGroup[Constant::decomp_num_steps];
    if (!GCinitemiprobs){
      GCinitemiprobs = new vector<Double>*[Constant::decomp_num_steps];
      for (int idx = 0; idx < Constant::decomp_num_steps; idx++)
	GCinitemiprobs[idx] = new vector<Double>[3];
    }
    if (!GCetemiprobs){
      GCetemiprobs = new vector<Double>*[Constant::decomp_num_steps];
      for (int idx = 0; idx < Constant::decomp_num_steps; idx++)
	GCetemiprobs[idx] = new vector<Double>[3];
    }
}

/*
 * ===[ ExonModel::readProbabilities ]====================================
 */
void ExonModel::readProbabilities(int parIndex) {
    exit(1);
    if( parIndex == lastParIndex )
	return;
    
    string filename = Constant::fullSpeciesPath() + Properties::getProperty("/ExonModel/infile");
    ifstream istrm; 
    istrm.open(filename.c_str(), ifstream::in);
    if (istrm) {
	int size;
	Seq2Int s2i_e(k+1);
	int dummyl, dummyi, dummyk;
	Double dbl;

	// for backward compatibility, check whether STARTCODON parameters exist at all
	streampos spos = istrm.tellg();
	istrm >> goto_line_after( "[STARTCODONS]" );
	if (!istrm){
	    istrm.clear();
	    istrm.seekg(spos); // go back to where you were and use default start codons
	} else {
	    GeneticCode::readStart(istrm);
	}

	if (!hasLenDist) {
	    // read length distributions
	    istrm >> goto_line_after( "[LENGTH]");
	    istrm >> comment >> exonLenD;
	    istrm >> comment >> slope_of_bandwidth;
	    istrm >> comment >> minwindowcount;
	    istrm >> comment >> numSingle >> numInitial >> numInternal >> numTerminal;
	    istrm >> comment >> numHugeSingle >> numHugeInitial >> numHugeInternal >> numHugeTerminal;
	    istrm >> comment;
	    lenDistSingle.assign(Constant::max_exon_len+1, 0);
	    lenDistInitial.assign(Constant::max_exon_len+1, 0);
	    lenDistInternal.assign(Constant::max_exon_len+1, 0);
	    lenDistTerminal.assign(Constant::max_exon_len+1, 0);
	    double boostfactor = 1.0;
	    for( int i = 0; i <= exonLenD; i++ ){
		istrm >> dummyi;
		istrm >> dbl;
		lenDistSingle[i] = dbl / 1000;
		istrm >> dbl; 
		lenDistInitial[i]= dbl / 1000; 
		istrm >> dbl;
		if (i > lenboostL){
		    // if requested (e.g. on bacteria), boost the probabilities of lengths > L by (1+E)^{i-L}
		    boostfactor *= (1 + lenboostE);
		    lenDistInitial[i] *= boostfactor;
		    lenDistSingle[i] *= boostfactor;
		}
		lenDistInternal[i] = dbl / 1000;
		istrm >> dbl;
		lenDistTerminal[i] = dbl / 1000;
	    }
	    // single exon can't be shorter than min_coding_len
	    for ( int i = 0; i < Constant::min_coding_len; i++) 
		lenDistSingle[i] = 0;
	    fillTailsOfLengthDistributions();
	    hasLenDist = true;
	}
	    
	/*
	 * begin of content dependent part
	 */

	char zusString[6];
	sprintf(zusString, "[%d]", parIndex);
	istrm >> goto_line_after(zusString);

	/*
	 * content model
	 */
	
	istrm >> goto_line_after( "[P_ls]" );
	istrm >> comment;      // dummy k 
	Pls.assign(k+1, 3);
	for( int l = 0; l <= k; l++ ){
	    string checkBase;
	    istrm >> comment >> dummyl;
	    int size = POWER4TOTHE(l+1);
	    Seq2Int s2i(l+1);
	    Pls[l][0].resize( size );
	    Pls[l][2] = Pls[l][1] = Pls[l][0];
	    for( int j = 0; j < size; j++ ){
		istrm >> comment;
		int pn = s2i.read(istrm);
		if (pn != j)
		    throw ProjectError("ExonModel::readProbabilities: Error reading file " + filename +
				       " at P_ls, pattern " + s2i.INV(pn));
		istrm >> Pls[l][0][j]
		      >> Pls[l][1][j]
		      >> Pls[l][2][j];
	    }
	}

	if (transInitMotif)
	    delete transInitMotif;
	transInitMotif = new Motif();
	istrm >> goto_line_after( "[TRANSINIT]" );
	transInitMotif->read(istrm);

	if (etMotif){
	    delete etMotif[0];
	    delete etMotif[1];
	    delete etMotif[2];
	}

	etMotif = new Motif*[3];
	etMotif[0] = new Motif();
	istrm >> goto_line_after( "[ETMOTIF0]" );
	etMotif[0]->read(istrm);
	etMotif[1] = new Motif();
	istrm >> goto_line_after( "[ETMOTIF1]" );
	etMotif[1]->read(istrm);
	etMotif[2] = new Motif();
	istrm >> goto_line_after( "[ETMOTIF2]" );
	etMotif[2]->read(istrm);

	// EMISSION
	// make the emission probabilities
	for (int f=0; f<3; f++) {
	    emiprobs.probs[f].resize( Pls[k][f].size() );
	    emiprobs.order = k;
	}
	// for backward compatibility, check whether EMISSION parameters exist at all, if not compute them from the Pls (old version)
	spos = istrm.tellg();
	istrm >> goto_line_after( "[EMISSION]" );
	if (!istrm){
	    istrm.clear();
	    istrm.seekg(spos); // go back to where you were
	    for (int f=0; f<3; f++)
		computeEmiFromPat(Pls[k][f], emiprobs.probs[f], k);
	} else {
	    istrm >> comment >> size;
	    istrm >> comment >> dummyk;
	    istrm >> comment >> dbl;
	    for( int i = 0; i < emiprobs.probs[0].size(); i++ ){
		istrm >> comment;
		int pn = s2i_e.read(istrm);
		if (pn != i)
		    throw ProjectError("ExonModel::readProbabilities: Error reading file " + filename +
				       " at EMISSION, pattern " + s2i_e.INV(pn));
		istrm >> emiprobs.probs[0][i]
		      >> emiprobs.probs[1][i]
		      >> emiprobs.probs[2][i];
	    }
	}
	
	istrm >> goto_line_after( "[INITEMISSION]" );
	istrm >> comment >> size;
	istrm >> comment >> dummyk;
	if (dummyk != k)
	    throw ProjectError("ExonModel::readProbabilities: Mismatch in order of exon INITEMISSION Markov chain.");
	istrm >> comment >> patpseudo;
	initemiprobs[0].resize(size);
	initemiprobs[1].resize(size);
	initemiprobs[2].resize(size);
	while( istrm >> comment >> ws, istrm && istrm.peek() != '[' ){
	    int pn = s2i_e.read(istrm);
	    istrm >> initemiprobs[0][pn]
		  >> initemiprobs[1][pn]
		  >> initemiprobs[2][pn];
	}


	istrm >> goto_line_after( "[ETEMISSION]" );
	istrm >> comment >> size;
	istrm >> comment >> k;   // TODO: k nicht von hier bestimmen lassen
	istrm >> comment >> patpseudo;
	etemiprobs[0].resize(size);
	etemiprobs[1].resize(size);
	etemiprobs[2].resize(size);
	while( istrm >> comment >> ws, istrm && istrm.peek() != '[' ){
	    int pn = s2i_e.read(istrm);
	    istrm >> etemiprobs[0][pn]
		  >> etemiprobs[1][pn]
		  >> etemiprobs[2][pn];
	}
	lastParIndex = parIndex;
    } else {
	string msg("ExonModel: Couldn't open file ");
	msg += filename;
	throw ProjectError(msg);

    }
}

/*
 * readAllParameters
 * read in species_exon_probs.pbl parameter files and store parameters for all gc content classes
 */
void ExonModel::readAllParameters(){
  string filename = Constant::fullSpeciesPath() + Properties::getProperty("/ExonModel/infile");
  ifstream istrm; 
  istrm.open(filename.c_str(), ifstream::in);
  if (istrm) {
    int size;
    Seq2Int s2i_e(k+1);
    int dummyl, dummyi, dummyk;
    Double dbl;
    
    // for backward compatibility, check whether STARTCODON parameters exist at all
    streampos spos = istrm.tellg();
    istrm >> goto_line_after( "[STARTCODONS]" );
    if (!istrm){
	istrm.clear();
	istrm.seekg(spos); // go back to where you were and use default start codons
    } else {
	GeneticCode::readStart(istrm);
    }

    if (!hasLenDist) {
      // read length distributions
      istrm >> goto_line_after( "[LENGTH]");
      istrm >> comment >> exonLenD;
      istrm >> comment >> slope_of_bandwidth;
      istrm >> comment >> minwindowcount;
      istrm >> comment >> numSingle >> numInitial >> numInternal >> numTerminal;
      istrm >> comment >> numHugeSingle >> numHugeInitial >> numHugeInternal >> numHugeTerminal;
      istrm >> comment;
      lenDistSingle.assign(Constant::max_exon_len+1, 0);
      lenDistInitial.assign(Constant::max_exon_len+1, 0);
      lenDistInternal.assign(Constant::max_exon_len+1, 0);
      lenDistTerminal.assign(Constant::max_exon_len+1, 0);
      
      for( int i = 0; i <= exonLenD; i++ ){
	istrm >> dummyi;
	istrm >> dbl;
	lenDistSingle[i] = dbl / 1000;
	istrm >> dbl; 
	lenDistInitial[i]= dbl / 1000; 
	istrm >> dbl; 
	lenDistInternal[i] = dbl / 1000;
	istrm >> dbl;
	lenDistTerminal[i] = dbl / 1000;
      }
      // single exon can't be shorter than min_coding_len
      for ( int i = 0; i < Constant::min_coding_len; i++) 
	  lenDistSingle[i] = 0;
      fillTailsOfLengthDistributions();
      hasLenDist = true;
    }

    // if requested (e.g. on bacteria), boost the probabilities of lengths > L by (1+E)^{i-L}
    Double boostfactor = 1;
    for (int i = lenboostL+1; i< lenDistSingle.size(); i++){
	boostfactor *= 1 + lenboostE;
	lenDistSingle[i] *= boostfactor;
    }

    /*
     * begin of GC content dependent part
     */
    char zusString[6];
   
    // loop over GC content classes and read in all remaining data
    for (int idx = 0; idx < Constant::decomp_num_steps; idx++) {
	  
      sprintf(zusString, "[%d]", idx+1);
      istrm >> goto_line_after(zusString);

      /*
       * content model
       */
	  
      istrm >> goto_line_after( "[P_ls]" );
      istrm >> comment;      // dummy k 
      GCPls[idx].assign(k+1, 3);
      for( int l = 0; l <= k; l++ ){
	string checkBase;
	istrm >> comment >> dummyl;
	int size = POWER4TOTHE(l+1);
	Seq2Int s2i(l+1);
	GCPls[idx][l][0].resize( size );
	GCPls[idx][l][2] = GCPls[idx][l][1] = GCPls[idx][l][0];
	for( int j = 0; j < size; j++ ){
	  istrm >> comment;
	  int pn = s2i.read(istrm);
	  if (pn != j)
	    throw ProjectError("ExonModel::readProbabilities: Error reading file " + filename +
			       " at P_ls, pattern " + s2i.INV(pn));
	  istrm >> GCPls[idx][l][0][j]
		>> GCPls[idx][l][1][j]
		>> GCPls[idx][l][2][j];
	  if (!Constant::contentmodels)
	      GCPls[idx][l][0][j] = GCPls[idx][l][1][j] = GCPls[idx][l][2][j] = 1.0/size; // use uniform distribution
	}
      }

      istrm >> goto_line_after( "[TRANSINIT]" );
      GCtransInitMotif[idx].read(istrm);
      spos = istrm.tellg();
      istrm >> goto_line_after( "[TRANSINITBIN]" );
      if (!istrm) {
	istrm.clear();
	istrm.seekg(spos); // go back to where you were
	Constant::tis_maxbinsize = 0; // no binning at all
      } else {
	GCtransInitBinProbs[idx].read(istrm);
	GCtransInitBinProbs[idx].setName(string("tis bin gc") + itoa(idx+1));
	//cout << "tis number of bins: " << GCtransInitBinProbs[idx].nbins << endl;
	//GCtransInitBinProbs[idx].write(cout);
      }
      istrm >> goto_line_after( "[ETMOTIF0]" );
      GCetMotif[idx][0]->read(istrm);
      istrm >> goto_line_after( "[ETMOTIF1]" );
      GCetMotif[idx][1]->read(istrm);
      istrm >> goto_line_after( "[ETMOTIF2]" );
      GCetMotif[idx][2]->read(istrm);

      // EMISSION
      // make the emission probabilities
      GCemiprobs[idx].setName(string("exon emiprob gc") + itoa(idx+1));
      for (int f=0; f<3; f++) {
	GCemiprobs[idx].probs[f].resize( GCPls[idx][k][f].size() );
	GCemiprobs[idx].order = k;
      }
      // for backward compatibility, check whether EMISSION parameters exist at all, if not compute them from the Psl (old version)
      spos = istrm.tellg();
      istrm >> goto_line_after( "[EMISSION]" );
      if (!istrm){
	istrm.clear();
	istrm.seekg(spos); // go back to where you were
	for (int f=0; f<3; f++)
	  computeEmiFromPat(GCPls[idx][k][f], GCemiprobs[idx].probs[f], k);
      } else {
	istrm >> comment >> size;
	istrm >> comment >> dummyk;
	istrm >> comment >> dbl;
	for( int i = 0; i < GCemiprobs[idx].probs[0].size(); i++ ){
	  istrm >> comment;
	  int pn = s2i_e.read(istrm);
	  if (pn != i)
	    throw ProjectError("ExonModel::readProbabilities: Error reading file " + filename +
			       " at EMISSION, pattern " + s2i_e.INV(pn));
	  istrm >> GCemiprobs[idx].probs[0][i]
		>> GCemiprobs[idx].probs[1][i]
		>> GCemiprobs[idx].probs[2][i];
	  if (!Constant::contentmodels)
	      GCemiprobs[idx].probs[0][i] = GCemiprobs[idx].probs[1][i] = GCemiprobs[idx].probs[2][i] = .25;
	}
      }
	
      istrm >> goto_line_after( "[INITEMISSION]" );
      istrm >> comment >> size;
      istrm >> comment >> dummyk;
      if (dummyk != k)
	throw ProjectError("ExonModel::readProbabilities: Mismatch in order of exon INITEMISSION Markov chain.");
      istrm >> comment >> patpseudo;
      GCinitemiprobs[idx][0].resize(size);
      GCinitemiprobs[idx][1].resize(size);
      GCinitemiprobs[idx][2].resize(size);
      while( istrm >> comment >> ws, istrm && istrm.peek() != '[' ){
	int pn = s2i_e.read(istrm);
	istrm >> GCinitemiprobs[idx][0][pn]
	      >> GCinitemiprobs[idx][1][pn]
	      >> GCinitemiprobs[idx][2][pn];
      }


      istrm >> goto_line_after( "[ETEMISSION]" );
      istrm >> comment >> size;
      istrm >> comment >> k;
      istrm >> comment >> patpseudo;
      GCetemiprobs[idx][0].resize(size);
      GCetemiprobs[idx][1].resize(size);
      GCetemiprobs[idx][2].resize(size);
      while( istrm >> comment >> ws, istrm && istrm.peek() != '[' ){
	int pn = s2i_e.read(istrm);
	istrm >> GCetemiprobs[idx][0][pn]
	      >> GCetemiprobs[idx][1][pn]
	      >> GCetemiprobs[idx][2][pn];
      }
    } // loop over gc content classes ([1],[2], etc)
    istrm.close();
  } else {
    string msg("ExonModel: Couldn't open file ");
    msg += filename;
    throw ProjectError(msg);

  }
}


/*
 * ===[ ExonModel::getCodonUsage]====================================
 * Compute a probability vector pi(0),... pi(63) that represents the
 * codon usage.
 * For now, codon usage is approximated by emission probabilities
 * in future on could alternatively measure it during training and output in species_exon_probs.pbl.
 */
double *ExonModel::getCodonUsage(){
    int numCodons = 20000; // number of sampled codons
    double *pi = new double[64];
    for (int c=0; c<64; c++){
      pi[c] = 0;
    }
    
    if (!GCemiprobs)
    	throw ProjectError("ExonModel::getCodonUsage emission probabilities not initialized?");
    FramedPatMMGroup &e = GCemiprobs[Constant::decomp_num_steps/2]; // assume average GC content
    char *sampledCDS = getSampledCDS(e.probs, e.order, numCodons);
    //cout << "sampled cds " << sampledCDS << endl;
    Seq2Int s2i(3);
    for (int i=0; i<numCodons; i++)
	pi[s2i(sampledCDS+e.order+i*3)] += 1;
    delete [] sampledCDS;
    // normalize pi
    int s = 0;
    for (int c=0; c<64; c++)
	s += (int) pi[c]; // here, the elements of pi are still integers
    
    //    cout << "codon usage: " << endl;
    for (int c=0; c<64; c++){
	pi[c] /= s;
	// cout << c << "\t" << s2i.inv(c) << "\t" << pi[c] << endl;
    }
    
    return pi;
}

/*
 * ===[ ExonModel::fillTailsOfLengthDistributions]====================================
 * Define the tail of the length distributions: from exonLenD up to max_exon_length
 * Make it so that there is no jump at exonLenD and so that the probability of huge
 * exons (>= d) is correct.
 * Probabilities for lengths greater than max_exon_length are 0.
 */
void ExonModel::fillTailsOfLengthDistributions( ){
    Double a,p; 
    int k;
    
    a = lenDistSingle[exonLenD];
    p = Double(1) - a/(Double(numHugeSingle+1)/Double(numSingle+1));
    for (k = exonLenD+1; k <= Constant::max_exon_len; k++){
	lenDistSingle[k] = p * lenDistSingle[k-1];
    }

    a = lenDistInitial[exonLenD];
    p = Double(1) - a/((Double(numHugeInitial)+1)/(numInitial+1));
    for (k = exonLenD+1; k <= Constant::max_exon_len; k++){
	lenDistInitial[k] = p * lenDistInitial[k-1];
    }

    a = lenDistInternal[exonLenD];
    p = Double(1) - a/((Double(numHugeInternal)+1)/(numInternal+1));
    for (k = exonLenD+1; k <= Constant::max_exon_len; k++){
	lenDistInternal[k] = p * lenDistInternal[k-1];
    }

    a = lenDistTerminal[exonLenD];
    p = Double(1) - a/((Double(numHugeTerminal)+1)/(numTerminal+1));
    for (k = exonLenD+1; k <= Constant::max_exon_len; k++){
	lenDistTerminal[k] = p * lenDistTerminal[k-1];
    }
}


/*
 * ===[ ExonModel::initAlgorithms ]=======================================
 */
void ExonModel::initAlgorithms(Matrix<Double>& trans, int cur) {
    initAlgorithmsCalled = true;
}

// stuff that needs to be called once for all exon states
// set these parameters to the one of the GC content index gcIdx    
void ExonModel::updateToLocalGC(int from, int to){
    Pls = GCPls[gcIdx];
    emiprobs = GCemiprobs[gcIdx];
    initemiprobs[0] =  GCinitemiprobs[gcIdx][0];
    initemiprobs[1] =  GCinitemiprobs[gcIdx][1];
    initemiprobs[2] =  GCinitemiprobs[gcIdx][2];
    etemiprobs[0] =  GCetemiprobs[gcIdx][0];
    etemiprobs[1] =  GCetemiprobs[gcIdx][1];
    etemiprobs[2] =  GCetemiprobs[gcIdx][2];
    transInitMotif = &GCtransInitMotif[gcIdx];
    transInitBinProbs = GCtransInitBinProbs[gcIdx];
    etMotif = GCetMotif[gcIdx];
}

/*
 * ===[ ExonModel::viterbiForwardAndSampling ]=====================================
 *
 * Find the predecessor state "predState" and the ending position "predBase" of 
 * this predecessor state in a most probable path ending in state "state" at position "base". 
 *
 */
void ExonModel::viterbiForwardAndSampling(ViterbiMatrixType& viterbi, // matrix of viterbi variables
					  ViterbiMatrixType& forward, // matrix of forward variables
					  int state,
					  int base,                   // end of the exon state
					  AlgorithmVariant algovar,   // doViterbiOnly, doViterbiAndForward, doSampling or doBacktracking
					  OptionListItem& oli) {
/* 
 *  an example:
 *  pred.State |ATG               | <= (k+1)     | *** ... ***   | STP        |
 *             |start codon       | pattern prob | emission prob | stop codon |
 *             |                  |              | >= 0          |            |
 *                                 <---- inner sequence  ------->
 *  --------------------------------------------------------------------------
 *  predProb   |beginPartProb     | startSeqProb | restSeqProb   | endPartProb|
 *              <--------------------- lenPartProb -------------------------->
 * 
 */ 

    Double fwdsum = 0;
    Double maxProb = 0;
    OptionsList *optionslist = NULL;
    PP::ExonScorer* exonScorer = 0;
    if (algovar==doSampling)
	optionslist = new OptionsList();
    bool checkSubstates = (profileModel != 0); // && isOnFStrand(etype) && etype != singleG

    /*
     * endPartProb
     * probability of the last - fixed length - part of the exon:
     * translation termination if the exon is a terminal exon
     */
    Double endPartProb = endPartEmiProb(base);
    ViterbiSubmapType substates(state);
    
    /*
     * right is the rightmost position of the inner sequence part
     * just in front of the endPart
     * endOfBioExon is the position most downstream (wrt forward strand) of the biological exon
     */
    
    int endOfBioExon = base + baseOffset;
    int right = endOfBioExon - innerPartEndOffset;
    
    /*
     * If it is impossible that a (complete) exon ends here, return immediately.
     */
    
    if (endPartProb <= 0 || right < 0) {
	// viterbi[base].erase(state);
	// if (needForwardTable(algovar))
	//     forward[base].erase(state);
	if (algovar == doSampling)
	    throw ProjectError("ExonModel: Trying to sample from empty options list. Internal error.");
	return;
    }

    /*
     * Reading frame of the position "right"
     * win=0,1,2 is the reading frame of __ endOfBioExon+1 __
     * on the reverse strand, count backwards 
     */
    int frameOfRight = isOnFStrand(etype) ? 
	mod3(win - (endOfBioExon + 1) + right) : 
	mod3(win +  endOfBioExon + 1  - right);
    
    // Need this for protein pattern evaluation
    // int endOfLastCodon = endOfBioExon - (isOnFStrand(etype) ? win : 2 - win);
    // if (endOfLastCodon >= dnalen)
    // 	checkSubstates=false;
    if (endOfBioExon >= dnalen)
	checkSubstates = false;

     /*
     * Determine the leftmost possible boundary of the inner sequence part.  
     * That is, determine the longest possible open reading frame.
     */
    int endOfNonStopcodon = (etype == terminal || etype == singleG) ?
	endOfBioExon - STOPCODON_LEN : endOfBioExon;
    if (endOfNonStopcodon > dnalen-1) // to be safe in those rare cases
	endOfNonStopcodon = dnalen-1;

    int frameOfEndOfNonStopcodon = isOnFStrand(etype) ?
	mod3(win - 1 - endOfBioExon + endOfNonStopcodon) : 
	mod3(win + 1 + endOfBioExon - endOfNonStopcodon);

    int ORFleft = orf->leftmostExonBegin(frameOfEndOfNonStopcodon, endOfNonStopcodon, isOnFStrand(etype));

    /*
     * get the extrinsic exonpart information about parts falling in this range
     * and with fitting reading frame
     */
    int seqRelFrame = isOnFStrand(etype)? mod3(right - frameOfRight) : mod3(right + frameOfRight);
    Feature* extrinsicexons = seqFeatColl->getExonListOvlpingRange(ORFleft-1, endOfBioExon, bothstrands, seqRelFrame);

    if (checkSubstates) {
	if (algovar == doBacktracking) {
	    // in backtracking, state is a FullStateId if etype!=terminal
	    SubstateId substate;
	    getStatePair(state, state, substate); // extract the main state index
#ifdef DEBUG
#ifndef DEBUG_BACK
	    if (!substate.undef())
		cerr << "substate bonus at base " << base << ": "
		     << viterbi[base].get(state, substate) / viterbi[base][state]
		     << endl;
#endif
#endif		    
	    // if it's a terminal or rinitial exon, we backtrack from igenic which is
	    // unaware of substates, so the given state won't refer to a substate; 
	    // in order to find out, we check viterbi: if the probability comes from the
	    // profile, there will be substates, otherwise not
	    if (isLastExon(etype))  {
#ifdef DEBUG
		if (!substate.undef())
		    throw ProjectError("viterbiForwardAndSampling called with"
				       " bad state/substate combination"
			               " (last exon but substate=" + itoa(substate) + ").");
#endif
		checkSubstates = viterbi[base].hasSubstates(state);
	    } else {
		checkSubstates = !substate.undef();
	    }
	    if (checkSubstates) 
		exonScorer = new PP::SingleTargetExonScorer(*profileModel, etype, substate, 
							    endOfBioExon, ORFleft);
	} else {
	    // profileModel->advanceHSColl(endOfLastCodon);
	    if (isLastExon(etype))
		exonScorer = new PP::SingleTargetExonScorer(*profileModel, etype,
							    endOfBioExon);
	    else {
		substates.activate();
		exonScorer = new PP::MultiTargetExonScorer(*profileModel, substates, etype,
							   endOfBioExon);
	    }
	}
    }
#ifdef DEBUG
    if (algovar==doBacktracking)
	oli.base = oli.state = -1;
#endif

    // determine bounds for exon start
    int startMax = endOfBioExon + innerPartOffset - min_exon_length + 1;
    int startMin;
    if (isRTerminalExon(etype) || etype == rsingleG ) 
	startMin = startMax = ORFleft+2; // single gene start candidate position
    else {
	if (ORFleft <= 0)
	    startMin = 0; // left truncated exon
	else
	    // innerPartOffset is beginOfStart - beginOfBioExon
	    startMin = ORFleft + innerPartOffset; 
	if (startMax > base + beginPartLen) // ensure that base>endOfPred
	    startMax = base + beginPartLen;
    }
    /*
     * loop over the length of the inner sequence
     */
    seqProb(-1,0,0);      // initialize the static variables
    for (int beginOfStart = startMax; beginOfStart >= startMin; beginOfStart--) {
	// main work done in this funcion call
	int endOfPred = beginOfStart - beginPartLen -1;
	Double notEndPartProb = notEndPartEmiProb(beginOfStart, right, frameOfRight, extrinsicexons);
	if (notEndPartProb <= 0 || endOfPred >= dnalen) 
	    continue;
	// - initial or single exon starts right at pos 1
	// - left truncated internal or terminal exon
	const ViterbiColumnType& predForw =  forward[endOfPred >= 0 ? endOfPred : 0];
	ViterbiColumnType& predVit = viterbi[endOfPred >= 0 ? endOfPred : 0];
	// compute the maximum over the predecessor states probs times transition probability
	int predReadingFrame;
	int beginOfBioExon  = beginOfStart - innerPartOffset;
	int exonLength = endOfBioExon - beginOfBioExon +1;
	if (checkSubstates)
	    exonScorer->newFirstCodon(beginOfBioExon);
	
	vector<Ancestor>::const_iterator it;
	for( it = ancestor.begin(); it != ancestor.end(); ++it ){
	    int predState = it->pos;
	    if (algovar == doSampling) {
		if (predForw[predState] == 0)
		    continue;
	    } else
		if (predVit.get(predState)==0)
		    continue;
	    StateType predStateType = (*stateMap)[predState];
	    Double transEmiProb = it->val * endPartProb * notEndPartProb;
	    predReadingFrame = stateReadingFrames[predStateType];
	    if (etype == singleG || etype == rsingleG || etype == rterminal0 || etype == rterminal1 || 
		etype == rterminal2 || etype == initial0 || etype == initial1 || etype == initial2 ||
		// has been checked in notEndPartProb
		win == mod3((isOnFStrand(etype))? predReadingFrame + exonLength: predReadingFrame - exonLength)) 
	    {
		// transition and exon length fit with respect to reading frame
		if (needForwardTable(algovar)) {
		    Double fwdsummand = predForw[predState] * transEmiProb.heated();
		    if (algovar == doSampling) {
			if (fwdsummand > 0) 
			    optionslist->add(predState, endOfPred, fwdsummand);
			continue; // next for
		    } else 
			fwdsum += fwdsummand;
		}
		if (checkSubstates) {
		    if (exonScorer->validSize())
			try {
			    exonScorer->score(predVit, predState, transEmiProb);
			    if (algovar==doBacktracking && exonScorer->hasNewMax()) {
				oli.base = endOfPred;
				oli.resetPredEnd();
				oli.state = getFullStateId(predState, exonScorer->getPredSubstate());
			    }
			} catch(NoSubmapFoundError e) {}
		} 
		Double predProb = predVit.get(predState) * transEmiProb;
		if (predProb > maxProb) {
		    maxProb = predProb;
		    if (algovar==doBacktracking && !checkSubstates) {
			oli.base = endOfPred;
			oli.resetPredEnd();
			oli.state = predState;
		    }
		}
	    }
	} // end of loop over predecessor state
	if (Constant::overlapmode){ // potentially increase maxProb and change oli if overlap is better
	    processOvlpOption(viterbi, forward, algovar, state, endOfPred, beginOfBioExon, maxProb,
			      endPartProb * notEndPartProb, fwdsum, optionslist, oli, base);
	    if (oli.predEnd >= base){
		cerr << "possible infinite loop " << oli.predEnd << " " << base << " " << state << " " << oli.state << " " << endOfPred << endl;
	    }
	}
    } // end of loop over the exon length
    
    switch (algovar) {
	case doSampling:
	    optionslist->prepareSampling();
	    try {
		oli = optionslist->sample();
	    } catch (ProjectError e) {
		cerr << "Sampling error in exon model. state=" << state << " base=" << base << endl;
		throw e;
	    }
	    delete optionslist;
	    break;
	case doViterbiAndForward:
	    if (fwdsum > 0)
		forward[base][state] = fwdsum;
	case doViterbiOnly:
	    if (checkSubstates) {
		exonScorer->postProcessing(maxProb);
		if (isLastExon(etype)) {
		    substates.activate();
		    exonScorer->exportSubstates(substates);
		    substates.useMaximum(maxProb);
		}
		if (!substates.empty()) {
		    if (!isFirstExon(etype))
			substates.registerPredecessors();
		    viterbi[base].addSubstates(substates);
		}
	    }
	    if (maxProb > 0)
		viterbi[base][state] = maxProb;
	    break;
	case doBacktracking:
	    if (checkSubstates) {
		exonScorer->addMatches(oli.base + beginPartLen + 1 - innerPartOffset, endOfBioExon);
#ifdef DEBUG
		static_cast<PP::SingleTargetExonScorer*>(exonScorer)->
		    checkResults(oli, viterbi, state, base, endPartProb, 
				 oli.state==-1 ? 0 : notEndPartEmiProb(oli.base + beginPartLen + 1, right, 
								       frameOfRight, extrinsicexons));
#endif			
	    }
	default:
	    break;
    }
    delete exonScorer;
} // ExonModel::viterbiForwardAndSampling

/*
 * ===[ ExonModel::processOvlpOption ]=====================================
 *
 * Check in the inner loop of the Viterbi algorithm whether introducing an overlap of exons
 * improves the probability maxProb of the best non-overlapping option. This introduces another nested loop.
 *        beginOfBioExon
 *           |
 *           >------------------------------------------------------------------->
 * ----------------------------------------------->
 *   |<----------- ovlp --------------->|         |
 * endOfPred                         endOfPred2   |
 *                                              endOfBioExon (previous)
 *           |-----------   bioOvlp  -------------|
*/
void ExonModel::processOvlpOption(ViterbiMatrixType& viterbi, ViterbiMatrixType& forward, AlgorithmVariant& algovar, 
				  int state, int endOfPred, int beginOfBioExon, Double &maxProb,
				  Double emiProb, Double &fwdsum, OptionsList *optionslist, OptionListItem &oli, int base) const {
  vector<Ancestor>::const_iterator it;
  // distinguish between maxOvlp (biological) and maximum for ovlp (states)
  // The actual biological overlap (that depends on the state combination) may not be larger than Constant::maxOvlp
  int maxStateOvlp = Constant::maxOvlp + 2 * trans_init_window;
  for (int ovlp = 0; ovlp <= maxStateOvlp && endOfPred + ovlp >= 0 && endOfPred + ovlp < dnalen; ovlp++){
    int endOfPred2 = endOfPred + ovlp;
    int endOfBioExon; // of left gene, depends on endOfPred2 and on type of exon
    int bioOvlp; // actual number of bases shared in the two exons
    Double lenCorrection = Double(4).pow(ovlp); // heuristic, because bases in overlap are evaluated twice
    Double lenProb;
    const ViterbiColumnType& predForw =  forward[endOfPred2 >= 0 ? endOfPred2 : 0];
    ViterbiColumnType& predVit = viterbi[endOfPred2 >= 0 ? endOfPred2 : 0]; 
    
    //cout << "state=" << state << " endOfPred=" << endOfPred << " endOfPred2=" << endOfPred2 << " maxProb=" << maxProb << endl;
    for (it = ancestor.begin(); it != ancestor.end(); ++it ){
      int predState = it->pos;
      if (algovar == doSampling) {
	if (predForw[predState] == 0)
	  continue;
      } else
	if (predVit.get(predState)==0)
	  continue; // predecessor state cannot end at this position
      StateType predStateType = (*stateMap)[predState];
      if (predStateType == igenic)
	  continue; // makes no sense to overlap with intergenic region
      endOfBioExon = endOfPred2  + getBaseOffset(predStateType);
      bioOvlp = endOfBioExon - beginOfBioExon + 1;
      if (bioOvlp < 0)
	  lenProb = 1;
      else if (bioOvlp <= Constant::maxOvlp) {
	  if (isOnFStrand(etype) == isOnFStrand(predStateType))
	      lenProb = Constant::head2tail_ovlp[bioOvlp];
	  else if (isOnFStrand(etype))
	      lenProb = Constant::head2head_ovlp[bioOvlp];
	  else 
	      lenProb = Constant::tail2tail_ovlp[bioOvlp];
      } else
	  lenProb = 0;
      if (lenProb == 0)
	  continue;
      Double transEmiProb = it->val * emiProb * lenProb * lenCorrection;
      if (needForwardTable(algovar)) {
	  Double fwdsummand = predForw[predState] * transEmiProb.heated();
	  if (algovar == doSampling) {
	      if (fwdsummand > 0) 
		  optionslist->add(predState, endOfPred, fwdsummand, endOfPred2);
	      continue;
	  } else 
	      fwdsum += fwdsummand;
      }
      Double predProb = predVit.get(predState) * transEmiProb;
      //cout << "possible predState=" << predState << " predProb=" << predProb << endl;
      if (predProb > maxProb && endOfPred2 < base) {
	  // cout << "# overlap improves at endOfPred=" << endOfPred << " endOfPred2=" << endOfPred2 
	  //     << " state=" << state << " predState=" << predState << " ovlp=" << ovlp
	  //     << " bioOvlp=" << bioOvlp << " lenProb=" << lenProb << endl;
	  maxProb = predProb;
	  if (algovar == doBacktracking) {
	      oli.base = endOfPred;
	      oli.predEnd = endOfPred2;
	      oli.state = predState;
	  }
      }
    }
  }
}

/*
 * ===[ ExonModel::endPartEmiProb ]=====================================
 *
 * Compute the probability that the end part ends at "end" in the test sequence.
 * Return 0 if it is impossible.   
 */

Double ExonModel::endPartEmiProb(int end) const {
    static Double endPartProb(0);
    Feature *feature;
    Double extrinsicEmiQuot(1);
    switch( etype ){
        case singleG: case terminal:
  	{
	    int stppos = end - STOPCODON_LEN + 1;
	    if( stppos < 0 || stppos > dnalen-3 || !GeneticCode::isStopcodon(sequence + stppos) )
	      endPartProb = 0;
	    else { 
		// assign probabilities to the stop codons
		// human: taa 28%, tga 48%, tag 24%
		if (strncmp(sequence + stppos, "taa", 3)==0)
		    endPartProb = Constant::ochreprob;
		else if (strncmp(sequence + stppos, "tag", 3)==0)
		    endPartProb = Constant::amberprob;
		else if (strncmp(sequence + stppos, "tga", 3)==0)
		    endPartProb = Constant::opalprob;
		else 
		    throw ProjectError("ExonModel::endPartEmiProb: internal error, unknown stop codon");
		// check if we have extrinsic information about a stop codon
		feature = seqFeatColl->getFeatureListOvlpingRange(stopF, end-2, end , plusstrand);
		if (feature) {
		  while (feature) {
		    if (feature->start <= end-2 && feature->end >= end)
		      extrinsicEmiQuot *= feature->distance_faded_bonus(end-1);
		    feature = feature->next;
		  }
		} else if (seqFeatColl->collection->hasHintsFile)
		  extrinsicEmiQuot = seqFeatColl->collection->malus(stopF);
		/*
		feature = seqFeatColl->getFeatureAt(stopF, end, plusstrand); 
		if (feature)
		    extrinsicEmiQuot = feature->bonus;
		else if (seqFeatColl->collection->hasHintsFile) {
		    extrinsicEmiQuot = seqFeatColl->collection->malus(stopF);
		    }*/
	    }
	}
	break;
	case rsingleG: case rinitial:
	{
	    int startpos = end - trans_init_window - STARTCODON_LEN + 1;
	    if (startpos >= 0 && GeneticCode::isStartcodon(sequence + startpos, true)) {
		endPartProb = GeneticCode::startCodonProb(sequence + startpos, true);
		if (endPartProb > 0){
		    if (startpos + STARTCODON_LEN + trans_init_window - 1 + tis_motif_memory < dnalen){
			endPartProb *= transInitMotif->seqProb(sequence + startpos + STARTCODON_LEN, true, true);// HMM
			if (transInitBinProbs.nbins >= 1) {
			    int idx = GCtransInitBinProbs[gcIdx].getIndex(endPartProb);// map prob to CRF score
			    if (inCRFTraining && (countEnd < 0 || (startpos >= countStart && startpos <= countEnd)))
				GCtransInitBinProbs[gcIdx].addCount(idx);
			    endPartProb = transInitBinProbs.avprobs[idx];
			}
		    }
		    else
			endPartProb = pow(0.25, (double)(dnalen-(startpos + STARTCODON_LEN)));
		}
	    } else 
		endPartProb = 0;
	    // check if we have extrinsic information about a reverse start codon
	    feature = seqFeatColl->getFeatureListOvlpingRange(startF, startpos, startpos + STARTCODON_LEN - 1 , minusstrand); 
	    if (feature) {
	      while (feature) {
		if (feature->start <= startpos && feature->end>= startpos + STARTCODON_LEN - 1)
		  extrinsicEmiQuot *= feature->distance_faded_bonus(startpos + 1);
		feature = feature->next;
	      }
	    } else if (seqFeatColl->collection->hasHintsFile)
	      extrinsicEmiQuot = seqFeatColl->collection->malus(startF);
	      
	    //feature = seqFeatColl->getFeatureAt(startF, startpos + STARTCODON_LEN - 1 , minusstrand); 
	    //if (feature)
	    //	extrinsicEmiQuot = feature->bonus;
	    //else if (seqFeatColl->collection->hasHintsFile)
	    //	extrinsicEmiQuot = seqFeatColl->collection->malus(startF);
	}
	break;
        case initial0: case initial1: case initial2: case internal0: case internal1: case internal2:
	{
	    int dsspos = end + Constant::dss_start + 1;
	    if (end == dnalen-1) {// exon is longer than dna, right truncated exon
	      endPartProb = 1;  // in this case allow that there is no splice site consensus
	    } else if ((dsspos + DSS_MIDDLE - 1 < dnalen && !isPossibleDSS(dsspos)) || 
		       end + Constant::dss_start >= dnalen || 
		       orf->leftmostExonBegin(win - 1, end + Constant::dss_start, true) >= end)
		endPartProb = 0;
	    else 
		endPartProb = 1;
	    feature = seqFeatColl->getFeatureListContaining(A_SET_FLAG(dssF), dsspos, plusstrand);
	    if (feature) {
		while (feature) {
		    extrinsicEmiQuot *= feature->distance_faded_bonus(dsspos);
		    feature = feature->next;
		}
	    } else if (seqFeatColl->collection->hasHintsFile)
	        extrinsicEmiQuot = seqFeatColl->collection->malus(dssF) * seqFeatColl->localSSMalus(dssF, dsspos, plusstrand);
	}
	break;
        case rterminal0: case rterminal1: case rterminal2: case rinternal0: case rinternal1: case rinternal2:
	{
	    int asspos = end + Constant::ass_end + 1;
	    if (end == dnalen-1){ // exon is longer than dna, right truncated exon
	      endPartProb = 1;  // in this case allow that there is no splice site consensus
	    } else if (end + Constant::ass_end + ASS_MIDDLE < dnalen && isPossibleRASS(asspos)){
#ifdef DEBUG
		if (!seqFeatColl->validRASSPattern(sequence + asspos))
		    throw ProjectError("pattern found is not valid");
#endif
		endPartProb = 1;
	    } else 
		endPartProb = 0;
	    feature = seqFeatColl->getFeatureListContaining(A_SET_FLAG(assF), asspos, minusstrand);
	    if (feature)
		while (feature) {
		    extrinsicEmiQuot *= feature->distance_faded_bonus(asspos);
		    feature = feature->next;
		}
	    else if (seqFeatColl->collection->hasHintsFile)
	        extrinsicEmiQuot = seqFeatColl->collection->malus(assF) * seqFeatColl->localSSMalus(assF, asspos, minusstrand);
	}
	    break;
	default:
	    cerr << "ExonModel::viterbiAlgorithm: unknown alternative." << endl;
	    endPartProb =  0;
    }
    return endPartProb * extrinsicEmiQuot;
}


/*
 * ===[ ExonModel::notEndPartEmiProb ]=====================================
 *
 * Probability of the emission of the exon exluding the fixed length end part.
 * This includes 
 * - the beginPart     (translation initiation, ass or reverse dss)
 * - the startSeq part (first k or less emissions after the begin part)
 * - the restSeq part  (rest of the variable length sequence emission part)
 * - the length part   (probability of the length of the biological exon)
 *
 * "right" is the last position of the inner sequence part as described in viterbiAlgorithm.
 * frameOfRight is the reading frame of the position "right".
 */

Double ExonModel::notEndPartEmiProb(int beginOfStart, int right, int frameOfRight, Feature *exonparts) const {
    Double beginPartProb;
    Feature *feature;
    Double extrinsicQuot = 1;

    /*
     * probability of the begin part
     */
    int beginOfBioExon = beginOfStart - innerPartOffset;
    switch( etype ){
	case singleG: case initial0: case initial1: case initial2:
	    // start codon at the beginning?
	    if (beginOfBioExon >= 0 && beginOfBioExon < dnalen - 2 &&
		GeneticCode::isStartcodon(sequence + beginOfBioExon)){
		beginPartProb = GeneticCode::startCodonProb(sequence + beginOfStart - STARTCODON_LEN);
		if (beginPartProb > 0){
		    // two cases ... . the normal one with enough sequence space before the gene
		    int transInitStart = beginOfBioExon - trans_init_window;
		    if (transInitStart > transInitMotif->k){
			beginPartProb *= transInitMotif->seqProb(sequence+transInitStart);
			if (transInitBinProbs.nbins >= 1) {
			    int idx = GCtransInitBinProbs[gcIdx].getIndex(beginPartProb);// map prob to CRF score
			    if (inCRFTraining && (countEnd < 0 || (transInitStart >= countStart && transInitStart <= countEnd)))
				GCtransInitBinProbs[gcIdx].addCount(idx);
			    beginPartProb = transInitBinProbs.avprobs[idx];
			}
		    } else {
			/* ... and the case where there is no place for the transInitMotif
			 * take emission probs of 1/4 for the rest up to the beginning of the seq
			 * endOfPred is negative in this case!
			 * Need this if the gene starts right after the sequence.
			 */
			beginPartProb *= pow(0.25, (double)(beginOfStart - STARTCODON_LEN));
		    }
		    feature = seqFeatColl->getFeatureListOvlpingRange(startF, beginOfStart-3, beginOfStart-1 , plusstrand); 
		    if (feature) {
			while (feature) {
			    if (feature->start <= beginOfStart-3 && feature->end >= beginOfStart-1)
				extrinsicQuot *= feature->distance_faded_bonus(beginOfStart-2);
			    feature = feature->next;
			}
		    } else if (seqFeatColl->collection->hasHintsFile)
			extrinsicQuot = seqFeatColl->collection->malus(startF);
		}
	    } else 
		beginPartProb = 0; 
	    break;
	case terminal: case internal0: case internal1: case internal2:
	    if (beginOfStart > 0) {
		// only a shortcut if there is no possible ass at the right position:
		if ((beginOfBioExon < 0) || 
		    ((beginOfBioExon - ASS_MIDDLE >=0 ) && !isPossibleASS(beginOfBioExon - 1)))
		    beginPartProb = 0;
		else 
		    beginPartProb = 1; // splice site is evaluated in other state
		if (beginPartProb > 0){
		    feature = seqFeatColl->getFeatureListContaining(A_SET_FLAG(assF), beginOfBioExon - 1, plusstrand);
		    if (feature) {
			while (feature) {
			    extrinsicQuot *= feature->distance_faded_bonus(beginOfBioExon - 1);
			    feature = feature->next;
			}
		    } else if (seqFeatColl->collection->hasHintsFile)
		        extrinsicQuot = seqFeatColl->collection->malus(assF) * seqFeatColl->localSSMalus(assF, beginOfBioExon - 1, plusstrand);
		}
	    } else if (beginOfStart == 0) // left truncated
	        beginPartProb = 1;
	    else 
	        beginPartProb = 0;
	    break;
	case rsingleG: case rterminal0: case rterminal1: case rterminal2:
	    // reverse stop codon
	    if (beginOfBioExon < 0)
		beginPartProb = 0;
	    else if (strncmp(sequence + beginOfBioExon, "tta", 3)==0)
		beginPartProb = Constant::ochreprob;
	    else if (strncmp(sequence + beginOfBioExon, "cta", 3)==0)
		beginPartProb = Constant::amberprob;
	    else if (strncmp(sequence + beginOfBioExon, "tca", 3)==0)
		    beginPartProb = Constant::opalprob;
	    else 
		beginPartProb = 0;
	    // extrinsic info about a reverse stop codon?
	    if (beginPartProb > 0){
	       feature = seqFeatColl->getFeatureListOvlpingRange(stopF, beginOfStart-3, beginOfStart-1 , minusstrand); 
	       if (feature) {
		 while (feature) {
		   if (feature->start <= beginOfStart-3 && feature->end >= beginOfStart-1)
		     extrinsicQuot *= feature->distance_faded_bonus(beginOfStart-2);
		   feature = feature->next;
		 }
	       } else if (seqFeatColl->collection->hasHintsFile)
		 extrinsicQuot = seqFeatColl->collection->malus(stopF);
	       /*
	       feature = seqFeatColl->getFeatureAt(stopF, beginOfStart-1 , minusstrand);
		if (feature)
		    extrinsicQuot *= feature->bonus;
		else  if (seqFeatColl->collection->hasHintsFile)
		extrinsicQuot *= seqFeatColl->collection->malus(stopF);*/
	    }
	    break;
	case rinitial: case rinternal0: case rinternal1: case rinternal2:
	    // reverse donor splice site
	    if (beginOfStart == 0)
	      beginPartProb = 1; // left truncated
	    else if (beginOfBioExon < 0 || (beginOfBioExon - DSS_MIDDLE > 0 && !isPossibleRDSS(beginOfBioExon - 1)))
		beginPartProb = 0;
	    else 
		beginPartProb = 1; // splice site is evaluated in other state
	    if (beginPartProb > 0){
		feature = seqFeatColl->getFeatureListContaining(A_SET_FLAG(dssF), beginOfBioExon - 1, minusstrand);
		if (feature) {
		    while (feature) {
			extrinsicQuot *= feature->distance_faded_bonus(beginOfBioExon - 1);
			feature = feature->next;
		    }
		} else if (seqFeatColl->collection->hasHintsFile && beginOfBioExon > 0)
		    extrinsicQuot = seqFeatColl->collection->malus(dssF) * seqFeatColl->localSSMalus(dssF, beginOfBioExon - 1, minusstrand);
	    }
	    break;
	default:
	    throw ExonModelError("notEndPartEmiProb: Wrong ExonType");
	    break;
    }

    if (!(beginPartProb > 0)) 
	return 0;
	
    /*
     * sequence emission probability: restSeqProb
     */
    Double restSeqProb;
    if (beginOfStart > right) {
	/*
	 * inner sequence empty or negative length = overlapping begin and end part
	 * this is for dealing with very short exons, which would be prohibited otherwise
	 * by the length of the fixed length emission of the beginning and end part.
	 * The method here is: Use the beginning and end part even if not really applicable
	 * and correct for the fact that the emissions of the overlapping part were taken into
	 * account twice by dividing by a standard emission probability of 1/4.
	 */
	restSeqProb = POWER4TOTHE( beginOfStart - right -1);
    } else if (right-beginOfStart <= k){
	/*
	 *  inner sequence has length <= k+1, a short exon but with inner sequence of positive length
	 */
	restSeqProb = 1;
	Seq2Int s2i(right-beginOfStart+1);
	try {
	    if (isOnFStrand(etype))
		restSeqProb = Pls[right-beginOfStart][frameOfRight][s2i(sequence + beginOfStart)];
	    else
		restSeqProb = Pls[right-beginOfStart][mod3(frameOfRight + right - beginOfStart)]
		    [s2i.rc(sequence + beginOfStart)];		
	} catch (InvalidNucleotideError e) {
	    // we dont assume anything in this case, take iid uniform distribution on {a,c,g,t}
	    restSeqProb = pow(Constant::probNinCoding, right-beginOfStart + 1); // 0.25
	}
    } else {
	/*
	 * inner sequence has length > k+1, this is the normal case
	 */
	
	/*
	 * inner sequence part from endOfStart+1 to right
	 * | initial pattern | initial content model| exon content model| exon terminating model (dss)|
	 * Depending on the exon type, some of the models may not apply
	 * For short exons the precedence is inital -> terminal -> inner
	 * but terminal applies only completely or not at all.
	 */
	int endOfInitial, beginOfTerm,  endOfStart; // forward
	int endOfTerm, beginOfInitP, beginOfInitial;// reverse

	endOfStart = beginOfStart + k-1;
	beginOfInitP = right - (k-1);
	try {
	    if (k==0)
		restSeqProb = 1;
	    else
		if (isOnFStrand(etype)) // init pattern at left side
		    restSeqProb = Pls[k-1][mod3(frameOfRight-right+endOfStart)][Seq2Int(k)(sequence + beginOfStart)];
		else                  // init pattern at right side
		    restSeqProb = Pls[k-1][mod3(frameOfRight+right-beginOfInitP)]
			[Seq2Int(k).rc(sequence + beginOfInitP)];
	} catch (InvalidNucleotideError e) {
	    restSeqProb = pow(Constant::probNinCoding, (int) k ); // 0.25
	}

	switch( etype ){
	    case singleG:
		/*
		 * | initial pattern | initial content model | content model |
		 *   <---   k   --->
		 */
		endOfInitial = endOfStart + Constant::init_coding_len;
		if (endOfInitial > right)
		    endOfInitial = right;
		restSeqProb *=
		    initialSeqProb(endOfStart+1, endOfInitial, mod3(frameOfRight-right+endOfInitial)) *
		    seqProb(endOfInitial+1, right, frameOfRight);
		break;
	    case initial0: case initial1: case initial2:
		/*
		 * | initial pattern | initial content model | content model | eterminal model |
		 *   <---   k   --->   <- init_coding_len ->                  < et_coding_len >  
		 */
		endOfInitial = endOfStart + Constant::init_coding_len;
		if (endOfInitial > right) {
		    endOfInitial = right;
		    beginOfTerm = right + 1;
		} else  {
		    beginOfTerm = right - Constant::et_coding_len + 1;
		    if (beginOfTerm  <= endOfInitial )
			beginOfTerm = right + 1; // no terminal part 
		}
		restSeqProb *= 
		    initialSeqProb(endOfStart+1, endOfInitial, mod3(frameOfRight-right+endOfInitial)) *
		    seqProb(endOfInitial+1, beginOfTerm-1, mod3(frameOfRight-right+(beginOfTerm-1))) *
		    eTermSeqProb(beginOfTerm, right, frameOfRight);
		break;	
	    case internal0: case internal1: case internal2: 
		/*
		 * | initial pattern | content model | eterminal model |
		 *   <---   k   --->                  < et_coding_len >  
		 */
		beginOfTerm = right - Constant::et_coding_len + 1;
		if (beginOfTerm <= endOfStart)
		    beginOfTerm = right + 1;   // exon not long enough, no terminal part
		restSeqProb *= 
		    seqProb(endOfStart+1, beginOfTerm-1, mod3(frameOfRight-right+(beginOfTerm-1))) *
		    eTermSeqProb(beginOfTerm, right, frameOfRight);
		break;
	    case terminal:
		/*
		 * | initial pattern | content model |
		 *   <---   k   --->
		 */
		restSeqProb *= seqProb(endOfStart+1, right, frameOfRight);
		break;
	    case rsingleG:
		/*
		 * | content model | initial content model | initial pattern |
		 *                   <- init_coding_len ->   <---   k   --->
		 */
		beginOfInitial = beginOfInitP - Constant::init_coding_len;
		if (beginOfInitial < beginOfStart)
		    beginOfInitial = beginOfStart;
		restSeqProb *=
		    initialSeqProb(beginOfInitial, beginOfInitP-1, mod3(frameOfRight+right-(beginOfInitP-1))) *
		    seqProb(beginOfStart, beginOfInitial-1, mod3(frameOfRight+right-(beginOfInitial-1)));
		break;
	    case rinitial:
		/*
		 * | eterminal model | content model | initial content model | initial pattern |
		 *   < et_coding_len >                 <- init_coding_len ->   <---   k   --->
		 */
		beginOfInitial = beginOfInitP - Constant::init_coding_len;
		if (beginOfInitial < beginOfStart) {
		    beginOfInitial = beginOfStart;
		    endOfTerm = beginOfStart - 1;
		} else {	    
		    endOfTerm = beginOfStart + Constant::et_coding_len - 1;
		    if (endOfTerm >= beginOfInitial)
			endOfTerm = beginOfStart - 1;   // exon not long enough, no terminal part
		}
		restSeqProb *= 
		    initialSeqProb(beginOfInitial, beginOfInitP-1, mod3(frameOfRight+right-(beginOfInitP-1))) *
		    seqProb(endOfTerm + 1, beginOfInitial-1, mod3(frameOfRight+right-(beginOfInitial-1))) *
		    eTermSeqProb(beginOfStart, endOfTerm, mod3(frameOfRight+right-endOfTerm));
		break;	
	    case rinternal0: case rinternal1: case rinternal2:
		/*
		 * | eterminal model | content model | initial pattern |
		 *   < et_coding_len >                 <---   k   ---> 
		 */
		endOfTerm = beginOfStart + Constant::et_coding_len - 1;
		if (endOfTerm >= beginOfInitP) 
		    endOfTerm = beginOfStart - 1;   // exon not long enough, no terminal part
		restSeqProb *=
		    seqProb(endOfTerm + 1, beginOfInitP-1, mod3(frameOfRight+right-(beginOfInitP-1))) *
		    eTermSeqProb(beginOfStart, endOfTerm, mod3(frameOfRight+right-endOfTerm));
		break;
	    case rterminal0: case rterminal1: case rterminal2:
		/*
		 * | content model | initial pattern | 
		 *                   <---   k   --->
		 */
		restSeqProb *=
		    seqProb(beginOfStart, beginOfInitP-1,mod3(frameOfRight+right-(beginOfInitP-1)));
		break;
	    default:
		throw ExonModelError("notEndPartEmiProb, inner sequence part: Wrong ExonType");
		break;
	}
    }
    
    /*
     * compute the probability of the length
     */
    Double lenPartProb;
    int endOfBioExon = right + innerPartEndOffset;
    int exonLength = endOfBioExon - beginOfBioExon + 1;

    if (exonLength < 1) {
	lenPartProb = 0;
    } else{
	switch( etype ){
	    case singleG: case rsingleG:
		if (exonLength%3 == 0)
		    lenPartProb = 3*lenDistSingle[exonLength];
		else
		    lenPartProb = 0;
		break;
	    case initial0: case initial1: case initial2:
		if (exonLength%3 == win && exonLength > 2) 
		    lenPartProb = 3*lenDistInitial[exonLength];
		else
		    lenPartProb = 0;
		break;
	    case rinitial:
		if (exonLength > 2)
		    lenPartProb = 3*lenDistInitial[exonLength];
		else 
		    lenPartProb = 0;
		break;
	    case internal0: case internal1: case internal2: case rinternal0: case rinternal1: case rinternal2:
		/*
		 * for internal and terminal exons the length and the 
		 * predecessor state must fit together, see viterbi algorithm
		 */
		lenPartProb = 3*lenDistInternal[exonLength];
		break;
	    case terminal:
		lenPartProb = 3*lenDistTerminal[exonLength];
		break;
	    case rterminal0: case rterminal1: case rterminal2:
		if (mod3(2-exonLength) == win) 
		    lenPartProb = 3*lenDistTerminal[exonLength];
		else 
		    lenPartProb = 0;
		break;
	    default:
		throw ExonModelError("Wrong ExonType");
		break;
	}
    }
    /*
     * Multiply a bonus/malus to extrinsicQuot for every exonpart hint that is
     * covered by this biological exon
     */
    bool exonSupport = false; // used for malus
    bool CDSSupport = false;  // used for malus
    int numEPendingInExon=0, numCPendingInExon=0, nep=0; // just used for malus. 
    bool strandOK;
    Double partBonus = 1;
    for (Feature *part = exonparts; part!= NULL; part = part->next){
	strandOK = part->strand == bothstrands || part->strand == STRAND_UNKNOWN || isOnFStrand(etype) == (part->strand == plusstrand);
	if (part->type == exonpartF || part->type == CDSpartF){
	    if (part->type == exonpartF && part->end >= beginOfBioExon && part->end <= endOfBioExon)
		numEPendingInExon++;
	    if (part->type == CDSpartF && part->end >= beginOfBioExon && part->end <= endOfBioExon)
		numCPendingInExon++;
	    if (strandOK){
	        if (part->start >= beginOfBioExon && part->end <= endOfBioExon){
		    partBonus *= part->bonus;
		    nep += 1;
  	        } else { // exonpart not completely contained in exon
		  if (part->type == exonpartF) {
		    if ((etype == singleG || etype == initial0 || etype == initial1 || etype == initial2 ||
			 etype == rsingleG || etype == rterminal0 || etype == rterminal1 || etype == rterminal2) &&
			part->end <= endOfBioExon && part->end >= beginOfBioExon){
		      partBonus *= sqrt(part->bonus);
		      nep += 1;
		    }
		    if ((etype == singleG || etype == terminal || etype == rsingleG || etype == rinitial) &&
			part->start >= beginOfBioExon && part->start <= endOfBioExon){
		      partBonus *= sqrt(part->bonus); 
		      nep += 1;
		    }
		    // could both happen if single exon gene has overlapping exonpart at both end points
		  }
		}
	    }
	}
	/* 
	 * Multiply a bonus to extrinsicQuot for every exon/CDS hint that exactly matches these
	 * exon boundaries. Malus for unsupported exons.
	 */
	if (part->type == CDSF && part->start == beginOfBioExon && part->end == endOfBioExon && strandOK){
	    CDSSupport = true;
	    extrinsicQuot *= part->bonus;
	}
	if (part->type == exonF && strandOK) {
	    if (etype == singleG || etype == rsingleG) {
		// do nothing, since both exons ends are handled in UTR states (if there is UTR)
	    } else if (etype == internal0 || etype == internal1 || etype == internal2 || etype == rinternal0 || etype == rinternal1 || etype == rinternal2) {
		if (part->start == beginOfBioExon && part->end == endOfBioExon) {
		    exonSupport = true;
		    extrinsicQuot *= part->bonus; 
		}
	    } else if (etype == terminal || etype == rinitial) {
		if (part->start == beginOfBioExon && part->end > endOfBioExon) {
		    exonSupport = true;
		    extrinsicQuot *= sqrt(part->bonus); 
		}
	    } else {
		if (part->start < beginOfBioExon && part->end == endOfBioExon) {
		    exonSupport = true;
		    extrinsicQuot *= sqrt(part->bonus); 
		}
	    }
	}
    }
    extrinsicQuot *= partBonus;
    if (seqFeatColl && seqFeatColl->collection->hasHintsFile) {
	/* We have searched for extrinsic features.
	 * Then multiply the malus for each position of the exon.
	 * We should exclude those (few) positions where an exonpart ends.
	 * Exons longer than their exonpart hint have an incentive to become shorter.
	 */
      if (nep >= 5){// local malus for partially and unevenly supported CDS
	int zeroCov = seqFeatColl->numZeroCov(beginOfBioExon, endOfBioExon, CDSpartF, isOnFStrand(etype)? plusstrand : minusstrand);
	Double localPartMalus = seqFeatColl->collection->localPartMalus(CDSpartF, zeroCov, partBonus, nep);
	if (localPartMalus < 1 / partBonus) // at least have ab initio probabilities
	  localPartMalus = 1 / partBonus;
	// cout << "partBonus[" << beginOfBioExon << ", " << endOfBioExon << "]= " << partBonus << " zeroCov = " << zeroCov << " localPartMalus = " << localPartMalus << endl;
	extrinsicQuot *= localPartMalus;
      }
	if (exonLength-numEPendingInExon > 0)
	    extrinsicQuot *= seqFeatColl->collection->partMalus(exonpartF, exonLength-numEPendingInExon);
	if (exonLength-numCPendingInExon > 0)
	    extrinsicQuot *= seqFeatColl->collection->partMalus(CDSpartF, exonLength-numCPendingInExon);
    }

    if (seqFeatColl && seqFeatColl->collection->hasHintsFile) {
	// We have searched for extrinsic features but have not found an exon hint.
	if (!exonSupport)
	    extrinsicQuot *= seqFeatColl->collection->malus(exonF);
	if (!CDSSupport)
	    extrinsicQuot *= seqFeatColl->collection->malus(CDSF);
    }
    return beginPartProb * restSeqProb * lenPartProb * extrinsicQuot;
}


/*
 * ===[ ExonModel::emiProbUnderModel ]=====================================
 *
 * Probability of emitting dna[begin]...dna[end] in this exon state.
 * Includes sequence emission and state length.
 */

Double ExonModel::emiProbUnderModel(int begin, int end) const {
    if (begin > end)
	return 1;

    /*
     * endPartProb
     */
    Double endPartProb = endPartEmiProb(end);
    
    /*
     * right is the last position of the inner sequence part
     * endOfBioExon is the last position of the biological exon state
     * beginOfStart, see viterbiAlgorithm
     */
    
    int endOfBioExon = end + baseOffset;
    int right = endOfBioExon - innerPartEndOffset;
    int beginOfStart = begin + beginPartLen;
    int beginOfBioExon = beginOfStart - innerPartOffset;
    
    if( !(endPartProb > 0) || (right < 0))
	return 0;
    
     /*
     * Reading frame of the position "right"
     * win=0,1,2 is the reading frame of __ endOfBioExon+1 __
     * on the reverse strand, count backwards 
     */
    int frameOfRight = isOnFStrand(etype) ? 
	mod3(win - (endOfBioExon + 1) + right) : 
	mod3(win +  endOfBioExon + 1  - right);
    
    // Determine the leftmost possible boundary of the inner sequence part.  
    int left = orf->leftmostExonBegin(frameOfRight, right, isOnFStrand(etype));
    
    if (beginOfStart < left) // e.g. a stop codon in the reading frame
	return 0;
    /*
     * get the extrinsic exonpart information about parts falling in this range
     * and with fitting reading frame
     */
    Feature *extrinsicexons = NULL;
    if (seqFeatColl)
	extrinsicexons = seqFeatColl->getExonListOvlpingRange(beginOfBioExon, endOfBioExon, bothstrands, mod3(right - frameOfRight));
    seqProb(-1,0,0);
    Double notEndPartProb = notEndPartEmiProb(beginOfStart, right, frameOfRight, extrinsicexons); 
    
    return endPartProb * notEndPartProb;
}


/*
 * computes the probability of the emission of the sequence from left to right
 * left and right included
 */

Double ExonModel::seqProb(int left, int right, int frameOfRight) const {
    static Double seqProb = 1;
    static int oldleft = -1, oldright = -1, oldframe = -1;
    static StateType oldtype = TYPE_UNKNOWN;

    bool reverse = !isOnFStrand(etype);
    Seq2Int s2i(k+1);
    if (left < 0) {   // new initialisation
	seqProb = 1;
        oldleft = oldright = oldframe = -1;
        oldtype = TYPE_UNKNOWN;
	return 1;
    }
    if (left > right) 
        return 1;

    if (right == oldright && frameOfRight == oldframe && left <= oldleft && etype == oldtype) {
        for (int curpos = oldleft-1; curpos >= left; curpos--){
            try {
		int f = reverse? mod3(frameOfRight+right-curpos) : mod3(frameOfRight-right+curpos);
		int pn = reverse? s2i.rc(sequence+curpos) : s2i(sequence+curpos-k);
		seqProb *= emiprobs.probs[f][pn];
            } catch (InvalidNucleotideError e) {
                seqProb *= Constant::probNinCoding; //  0.25, 1/4
            }
        }
        oldleft = left;
        return seqProb;
    }
    
    // compute everything new
    seqProb = 1.0;
    for (int curpos = right; curpos >= left; curpos--) {
	try {
	    int f = reverse? mod3(frameOfRight+right-curpos) : mod3(frameOfRight-right+curpos);
	    int pn = reverse? s2i.rc(sequence+curpos) : s2i(sequence+curpos-k);
	    seqProb  *= emiprobs.probs[f][pn];
	    if (inCRFTraining && (countEnd < 0 || (curpos >= countStart && curpos <= countEnd)))
		GCemiprobs[gcIdx].addCount(GCemiprobs[gcIdx].getIndex(f,pn));
	} catch (InvalidNucleotideError e) {
	    seqProb  *= Constant::probNinCoding; // 0.25 1/4
	}
    }
    oldleft = left;
    oldright = right;
    oldframe = frameOfRight;
    oldtype = etype;
    return seqProb;
}

/*
 * eTermSeqProb
 * terminal part of exons before donor splice site
 */ 
Double ExonModel::eTermSeqProb(int left, int right, int frameOfRight) const {
//     static Double seqProb = 1; 
//     static int oldleft = -1, oldright = -1, oldframe =-1;
//     static StateType oldtype = TYPE_UNKNOWN;
    if (left > right) 
        return 1;
    bool reverse = !isOnFStrand(etype);
    if (!(right-left+1 == Constant::et_coding_len)) 
	throw ProjectError("Assertion failed in ExonModel::eTermSeqProb");
//     the first part was never called
//     if (left == oldleft && right == oldright && frameOfRight == oldframe && etype == oldtype)
// 	return seqProb;
//     else {
	//seqProb = etMotif[mod3(frameOfRight-Constant::et_coding_len + 1)]->seqProb(
	//    sequence + left);
// old method: Markov Model
	Double seqProb = 1; 
	Seq2Int s2i(k+1);
	for (int curpos = right; curpos >= left; curpos--) {
	    try {
		int f = reverse? mod3(frameOfRight+right-curpos) : mod3(frameOfRight-right+curpos);
		int pn = reverse? s2i.rc(sequence+curpos) : s2i(sequence+curpos-k);
		seqProb *= etemiprobs[f][pn];
	    } catch (InvalidNucleotideError e) {
		seqProb *= Constant::probNinCoding; // 0.25, 1/4
	    }
	}
// 	oldleft = left;
// 	oldright = right;
// 	oldframe = frame;
// 	oldtype = etype;
	return seqProb;
//     }
}

/*
 * initialSeqProb
 * first coding part of an initial or single exon (after the initial pattern)
 */ 
Double ExonModel::initialSeqProb(int left, int right, int frameOfRight) const {
    if (left > right) 
        return 1;
    Double seqProb = 1;
    Seq2Int s2i(k+1);
    bool reverse = !isOnFStrand(etype);
    for (int curpos = right; curpos >= left; curpos--) {
	try {
	    int f = reverse? mod3(frameOfRight+right-curpos) : mod3(frameOfRight-right+curpos);
	    int pn = reverse? s2i.rc(sequence+curpos) : s2i(sequence+curpos-k);
	    seqProb *= initemiprobs[f][pn];
	} catch (InvalidNucleotideError e) {
	    seqProb *= Constant::probNinCoding; // 0.25, 1/4
	}
    }
    return seqProb;
}
