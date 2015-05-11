/**********************************************************************
 * file:    igenicmodel.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  intergenic region
 * authors: Mario Stanke, mario@gobics.de
 *
 * date    |   author      |  changes 
 * --------|---------------|------------------------------------------
 * 1.05.02 | Mario Stanke  |  Viterbi Algorithm
 * 1.11.02 | Mario Stanke  |  emiProbUnderModel
 **********************************************************************/

#include "igenicmodel.hh"

// project includes
#include "properties.hh"
#include "merkmal.hh"
#include "projectio.hh"
#include "geneticcode.hh"
#include "extrinsicinfo.hh"
#include "intronmodel.hh" // so splice sites models can be reused here

// standard C/C++ includes
#include <fstream>

/*
 * Initialisation of static data members
 */

Integer IGenicModel::k = 4;
Double IGenicModel::patpseudocount = 1.0;
Integer IGenicModel::gesbasen = 0;
vector<Integer> IGenicModel::emicount;
PatMMGroup      IGenicModel::emiprobs("igenic emiprob");
PatMMGroup*     IGenicModel::GCemiprobs = NULL;
int IGenicModel::lastParIndex = -1;
vector<vector<Double> > IGenicModel::Pls;
vector<vector<Double> >* IGenicModel::GCPls = NULL;
int IGenicModel::verbosity;
double IGenicModel::geoProb = 0.9999;

/*
 * Initialisation of class variables
 */
void IGenicModel::init() {
    try{
	k = Properties::getIntProperty( "/IGenicModel/k" );
	patpseudocount = Properties::getDoubleProperty( "/IGenicModel/patpseudocount" );
	verbosity = Properties::getIntProperty("/IGenicModel/verbosity"); 
    } catch( ProjectError ) { }
    
    // create GC content array data structures IF NOT THERE YET
    if (!GCemiprobs)
      GCemiprobs = new PatMMGroup[Constant::decomp_num_steps];
    if (!GCPls)
      GCPls = new vector<vector<Double> >[Constant::decomp_num_steps];
}

/*
 * initAlgorithms
 * this is called before running the inference algorithms
 */
void IGenicModel::initAlgorithms (Matrix<Double>& trans, int cur) {
    geoProb = trans[cur][cur].doubleValue();
}

/*
 * updateToLocalGC
 * this is called when the GC index changes
 */
void IGenicModel::updateToLocalGC(int from, int to){
  // set these parameters to the one of the GC content index gcIdx
  if (!Constant::tieIgenicIntron 
      || IntronModel::GCemiprobs == NULL 
      || IntronModel::GCemiprobs[gcIdx].probs.size()==0 // this happens for intronless species
      || IntronModel::k != k) {
      emiprobs = GCemiprobs[gcIdx];
  } else {// use the intron content model
      emiprobs = IntronModel::GCemiprobs[gcIdx];
  }
  
  Pls = GCPls[gcIdx];
}

/*
 * ===[ IGenicModel::readProbabilities ]==================================
 */
void IGenicModel::readProbabilities( int parIndex ) {
    if( parIndex == lastParIndex )
	return;

    string filename = Constant::fullSpeciesPath() + Properties::getProperty("/IGenicModel/infile");
    ifstream istrm(filename.c_str());

    if( istrm ){
        int l;

	char zusString[6];
	sprintf(zusString, "[%d]", parIndex);
	istrm >> goto_line_after(zusString);
	istrm >> comment >> k;
	if (istrm.eof())
	    throw ProjectError("IgenicModel::readProbabilities: Error reading file " + filename + ". Truncated?");

        Pls.resize( k+1 );
	istrm >> goto_line_after( "[P_ls]" );
        for( int i = 0; i <= k; i++ ){
            istrm >> comment >> l >> comment;
            Pls[i].assign( POWER4TOTHE(l+1), 0.0 );
	    Seq2Int s2i(i+1);
	    for( int j = 0; j < Pls[i].size(); j++ ) {
                istrm >> comment;
		int pn = s2i.read(istrm);
		if (pn != j)
		    throw ProjectError("IgenicModel::readProbabilities: Error reading file " + filename +
				       " at P_ls, pattern " + s2i.INV(pn));
		istrm >> Pls[i][j];
            }
	}
	
	emiprobs.probs.resize(Pls[k].size() );
	emiprobs.order = k;
	Seq2Int s2i(k+1);
	istrm >> goto_line_after( "[EMISSION]" );
	if (!istrm){ // for backward compatibility with old HMM-parameter files
	    computeEmiFromPat(Pls[k], emiprobs.probs, k);
	} else {
	    int size;
	    istrm >> comment >> size;
	    for( int j = 0; j < emiprobs.probs.size(); j++ ) {
		istrm >> comment;
		int pn = s2i.read(istrm);
		if (pn != j)
		    throw ProjectError("IgenicModel::readProbabilities: Error reading file " + filename +
				       " at EMISSION, pattern " + s2i.INV(pn));
		istrm >> emiprobs.probs[j];
	    }
	}
	istrm.close();
	lastParIndex = parIndex;
    } else {
	string msg("IGenicModel: Couldn't open file ");
	msg += filename;
	throw ProjectError(msg);
    }
    //cerr << "exiting IGenicModel::readProbabilities" << endl;
}

/*
 * readAllParameters
 */
void IGenicModel::readAllParameters(){
    string filename = Constant::fullSpeciesPath() + Properties::getProperty("/IGenicModel/infile");
    ifstream istrm(filename.c_str());

    if( istrm ){
        // all GC content dependent parameters
        int l;
	char zusString[6];
	
	// loop over GC content classes
	for (int idx = 0; idx < Constant::decomp_num_steps; idx++) {
	  GCemiprobs[idx].setName(string("igenic emiprob gc") + itoa(idx+1));
	  sprintf(zusString, "[%d]", idx+1);
	  istrm >> goto_line_after(zusString);
	  istrm >> comment >> k;
	  if (istrm.eof())
	    throw ProjectError("IgenicModel::readAllParameters: Error reading file " + filename + ". Truncated?");
	  
	  GCPls[idx].resize( k+1 );
	  istrm >> goto_line_after( "[P_ls]" );
	  for( int i = 0; i <= k; i++ ){
            istrm >> comment >> l >> comment;
	    int size = POWER4TOTHE(l+1);
            GCPls[idx][i].assign( size, 0.0 );
	    Seq2Int s2i(i+1);
	    for( int j = 0; j < size; j++ ) {
	      istrm >> comment;
	      int pn = s2i.read(istrm);
	      if (pn != j)
		throw ProjectError("IgenicModel::readProbabilities: Error reading file " + filename +
				   " at P_ls, pattern " + s2i.INV(pn));
	      istrm >> GCPls[idx][i][j];
	      if (!Constant::contentmodels)
		  GCPls[idx][i][j] = 1.0/size; // use uniform distribution
            }
	  }
	
	  GCemiprobs[idx].probs.resize(GCPls[idx][k].size() );
	  GCemiprobs[idx].order = k;
	  Seq2Int s2i(k+1);
	  streampos spos = istrm.tellg();
	  istrm >> goto_line_after( "[EMISSION]" );
	  if (!istrm){ // for backward compatibility with old HMM-parameter files
	    istrm.clear();
	    istrm.seekg(spos); // go back to where you were
	    computeEmiFromPat(GCPls[idx][k], GCemiprobs[idx].probs, k);
	  } else {
	    int size;
	    istrm >> comment >> size;
	    for( int j = 0; j < GCemiprobs[idx].probs.size(); j++ ) {
		istrm >> comment;
		int pn = s2i.read(istrm);
		if (pn != j)
		    throw ProjectError("IgenicModel::readProbabilities: Error reading file " + filename +
				       " at EMISSION, pattern " + s2i.INV(pn));
		istrm >> GCemiprobs[idx].probs[j];
		if (!Constant::contentmodels)
		    GCemiprobs[idx].probs[j] = 0.25; // use uniform distribution
	    }
	  }
	}// end for idx loop
	istrm.close();
    } else {
	string msg("IGenicModel: Couldn't open file ");
	msg += filename;
	throw ProjectError(msg);
    }
}


/*
 * ===[ IGenicModel::viterbiForwardAndSampling ]===================================
 */
void IGenicModel::viterbiForwardAndSampling(ViterbiMatrixType& viterbi, // viterbi variables
					    ViterbiMatrixType& forward, // forward variables
					    int state,               // current state
					    int base,                // position in the dna
					    AlgorithmVariant algovar,
					    OptionListItem& oli) {
    vector<Ancestor>::const_iterator it;
    Double max(-1.0), curmax; // for viterbi table
    Double fwdsum(0.0), fwdsummand;  // for forward table
    Double emiProb = emiProbUnderModel(base, base);
    Double transEmiProb;
    OptionsList *optionslist = NULL;
    
    oli.base = base - 1;    
    if (algovar == doSampling)
	optionslist = new OptionsList();
    for( it = ancestor.begin() ; it != ancestor.end(); ++it ){
	transEmiProb = it->val * emiProb;
        curmax  = viterbi[base-1].get(it->pos) * transEmiProb;
	if (needForwardTable(algovar))
	    fwdsummand = forward[base-1].get(it->pos) * transEmiProb.heated();
        if( curmax > max ){
            max = curmax;
	    oli.state = it->pos;
	}
	if (needForwardTable(algovar)) {
	    fwdsum += fwdsummand;
	    if (algovar == doSampling && fwdsummand > 0) 
		optionslist->add(it->pos, oli.base, fwdsummand);
	}
    }	
    if (doViterbi(algovar) && state != oli.state) {
	// increment the successor counts for the predecessor state
	try {
	    ViterbiSubmapType& predSubmap = viterbi[base-1].getSubstates(oli.state);
	    if (!predSubmap.empty()) 
		predSubmap.succCount()++;
	} catch(...) {}
    }

    if (algovar == doSampling) {
	optionslist->prepareSampling();
	try {
	    oli = optionslist->sample();
	} catch (ProjectError e) {
	    cerr << "Sampling error in model for intergenic region. state=" << state << " base=" << base << endl;
	    throw e;
	}
	delete optionslist;
	return;
    } else {
	if (max > 0)
	    viterbi[base][state] = max;
	if (algovar == doViterbiAndForward && fwdsum > 0)
	    forward[base][state] = fwdsum;
    }
}


/*
 * IGenicModel::emiProbUnderModel
 *
 * return the emission probability of sequence[begin]...sequence[end]
 * excluding the transition probabilities
 * from Viterbi and Forward this is only called with end = begin
 */


Double IGenicModel::emiProbUnderModel(int begin, int end) const {
    static Seq2Int s2i( k+1 );    
    Double p = 1.0;
    Double extrinsicProb = 1.0;
    bool have_irparthint=false, have_nonexonpartF=false, have_nonirparthint=false;
    
    for (; begin <= end; begin++){
	Feature * hints = seqFeatColl->getFeatureListContaining(A_SET_FLAG(irpartF) | A_SET_FLAG(nonexonpartF) | A_SET_FLAG(nonirpartF), begin, bothstrands);
	if (hints) {
	    while (hints) {
		have_irparthint |= (hints->type == irpartF);
		have_nonexonpartF |= (hints->type == nonexonpartF);
		have_nonirparthint |= (hints->type == nonirpartF);
		if (hints->type != nonirpartF)
		    extrinsicProb  *= hints->bonus;
		else 
		    extrinsicProb  /= hints->bonus; // nonirpartF bonus, divide igenic prob by this factor instead of multiplying the emiprobs of all other states
		hints = hints->next;
	    }
	} 
	if(seqFeatColl->collection->hasHintsFile){
	    if (!have_irparthint)
		extrinsicProb *= seqFeatColl->collection->malus(irpartF);
	    if (!have_nonexonpartF)
		extrinsicProb *= seqFeatColl->collection->malus(nonexonpartF);
	    if (!have_nonirparthint)
		extrinsicProb /= seqFeatColl->collection->malus(nonirpartF); // have no nonirpart=genic hint. Reward igenic base.
	}

	if( begin > k ){
	    try {
		int pn = s2i(sequence+begin-k);
		p *= emiprobs.probs[pn];
		if (inCRFTraining){
		  if (Constant::tieIgenicIntron && IntronModel::GCemiprobs != NULL 
		      && IntronModel::GCemiprobs[gcIdx].probs.size() > 0)
		    IntronModel::GCemiprobs[gcIdx].addCount(pn);
		  else 
		    GCemiprobs[gcIdx].addCount(pn);
		}
	    } catch (InvalidNucleotideError e) {
		p *= 0.25;
	    }
	} else{
	    /* 
	     * at the beginning of the sequence, too short for the order k of the Markov model 
	     * compute the emission probability based on the probs for shorter patterns 
	     */
	    Seq2Int s2it( begin+1 );
	    try {
		int basek = s2it(sequence);
		p *=  Pls[begin][basek]/ (Pls[begin][basek/4]+Pls[begin][basek/4+1]+Pls[begin][basek/4+2]+Pls[begin][basek/4+3]);
	    } catch  (InvalidNucleotideError e) {
		p *= 0.25;
	    }
	}
    }
    return p * extrinsicProb;
}
