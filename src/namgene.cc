/*******************************************************************************************************
 * file:    namgene.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  
 * authors: Mario Stanke (mario@gobics.de), Stafilarakis
 *
 * date    |   author      |  changes 
 * --------|---------------|----------------------------------------------------------------------------- 
 * 26.09.01| Stafilarakis  | creation of the class
 * ??.11.02| Stanke        | getPathEmiProb
 * 09.03.04| Stefanie Dahms| getNextCutEndPoint
 *         | Stanke        | forward algorithm, sampling, sfc, ...
 * 09.09.05| Mario Stanke  | alternative transcripts and posterior probs from sampling
 * 13.07.06| Mario Stanke  | fixed minor bug: mixed up geneBegin and codingstart for sorting transcripts
 * 31.08.06| Mario Stanke  | rewrote getNextCutEndPoint
 * 21.09.06| Mario Stanke  | getStepGenes, etc.
 * 17.01.07| Mario Stanke  | computeReachableStates: save time by not going into irrelevant states
 *******************************************************************************************************/

#include "namgene.hh"

// project includes
#include "pp_scoring.hh" 
#include "statemodel.hh"
#include "motif.hh"
#include "gene.hh"
#include "projectio.hh"  // for comment, goto_line_after
#include "mea.hh"
#include "exoncand.hh"

// standard C/C++ includes
#include <iomanip>  // for setprecision
#include <iostream>

int synchstate = 0;

NAMGene::NAMGene() {
  /*
   * Prepare the following:
   *  - determine number of states
   *  - create states
   *  - read in transition probabilities
   */
  statecount = Properties::getIntProperty( "/NAMGene/statecount" );
  try {
      noInFrameStop = Properties::getBoolProperty("noInFrameStop");
  } catch (...) {
      noInFrameStop = false;
  }
  try {
    minmeanexonintronprob = Properties::getdoubleProperty("minmeanexonintronprob");
  } catch (...) {
    minmeanexonintronprob = 0.0;
  }
  try {
    minexonintronprob = Properties::getdoubleProperty("minexonintronprob");
  } catch (...) {
    minexonintronprob = 0.0;
  }
  try {
    maxtracks = Properties::getIntProperty("maxtracks");
  } catch (...) {
    maxtracks = -1; // default, no maximum
  }
  try {
    sampleiterations = Properties::getIntProperty("sample");
  } catch (...) {
    sampleiterations = 0;
  }
  if (sampleiterations > 0) {
      if (sampleiterations < 10){
	  cerr << "Error: Number of sample iterations is too low. (sample=" << sampleiterations << ")" << endl
	       << "I will not sample (sample=0) and will not estimate posterior probabilities." << endl;
	  sampleiterations = 0;
      } else if (sampleiterations < 20){
	  cerr << "Warning: Number of sample iterations 'sample'=" << sampleiterations << " is low." << endl
	       << "Posterior probabilities will be only rough estimates." << endl;
      }
  }
  try {
      alternatives_from_sampling = Properties::getBoolProperty("alternatives-from-sampling");
  } catch (...) {
      alternatives_from_sampling = false;
  }
  try {
      show_progress = Properties::getBoolProperty("progress");
  } catch (...) {
    show_progress = false;
  }

  try {
      mea_prediction = Properties::getBoolProperty("mea");
  } catch (...) {
    mea_prediction = false;
  }

  if(mea_prediction || Constant::MultSpeciesMode){
    sampleiterations = 100;
    alternatives_from_sampling = true;
    minexonintronprob = 0;
    minmeanexonintronprob = 0;
  }

  try {
    mea_eval = Properties::getBoolProperty("mea_evaluation");
  } catch (...) {
    mea_eval = false;
  }

  // Read in protein profile
  try {
      PP::Profile* prfl;
      try {
	  prfl = new PP::Profile(Properties::getProperty("proteinprofile"));
      } catch (PP::ProfileNotFoundError e) {
	  try {
	      prfl = new PP::Profile(Properties::getConfigFilename("proteinprofile"));
	  } catch (PP::ProfileNotFoundError) {
	      throw e;
	  }
      }
      profileModel = prfl ? new PP::SubstateModel(*prfl) : NULL;
      if (prfl) {
	  cout << "# Using protein profile " << prfl->getName() << "\n";
	  cout << "# " << prfl->getInfoLine() << "\n";
      }
  } catch (KeyNotFoundError) {
      profileModel = NULL;
  } catch (PP::ProfileInsigError e) {
      cerr << e.getMessage() << ". Not using protein profile.\n";
      cerr << "Choose different block thresholds to enforce using this profile.\n";
      cout << "# Not using protein profile.";
      profileModel = NULL;
  }

  needForwardTable = (sampleiterations > 0);
  readTransAndInitProbs();
  readOvlpLenDist();
  checkProbsConsistency();
  computeReachableStates();
  createStateModels();
  
  /*
   * Create map of state types
   */
  stateMap.resize(statecount);
  for (int i=0; i< statecount; i++) {
    stateMap[i] = getStateType(i);
  }

  if (Constant::temperature && Constant::augustus_verbosity)
      cout << "# setting temperature to " << Constant::temperature << " (for sampling)" << endl;
  sampledTxs = NULL;
}

void NAMGene::setStatesInitialProbs(){  
    for( int i = 0; i < statecount; i++ ) {
	viterbi[0][i] = initProbs[i];
	if (needForwardTable)
	    forward[0][i] = initProbs[i];
    }
}

StateType NAMGene::getStateType(int i){
    if (i >= 0 && i < statecount)
	return states[i]->getStateType();
    else
	return TYPE_UNKNOWN;
}

int NAMGene::getStateIndex(StateType type){
    for( int i = 0; i < statecount; i++ ) 
	if (getStateType(i) == type)
	    return i;
    string errmsg = "Unknown state type ";
    errmsg += type;
    throw NAMGeneError( errmsg);
}

void NAMGene::viterbiAndForward( const char* dna, bool useProfile){
 
  StateModel::setPP(useProfile ? profileModel : NULL);  
  int progress, oldprogress=0;

  /*
   * Determine the length of the DNA sequence under consideration
   * and form the forward matrix
   */
  int dnalen = strlen( dna );
  viterbi.assign(dnalen, statecount);
  if (needForwardTable) 
      forward.assign(dnalen, statecount);

  /*
   * Initialize the first Viterbi column with InitialProbs
   */
  setStatesInitialProbs();

  /*
   * The heart of the algorithm.
   * The actual work is done in the states.
   */
  OptionListItem oli;
  int deletionEnd = 0;
  cs.computeStairs(dna);
  StateModel::setContentStairs(&cs);

  curGCIdx = -1; // initialize with invalid GC content class
  prepareModels(dna, dnalen);
  

  /*
   * Check if all the sequence contains of non-nucleotides (e.g. N)
   * If yes, do not execute Viterbi/Forward: everything is considered intergenic
   */

  if (containsJustNonNucs(dna, dnalen)) {
#ifdef DEBUG
    cout << "#warning: sequence contains just letters different from a,c,g,t: " << dnalen << " of them." << endl;
#endif
    curGCIdx = cs.idx[0];
    initAlgorithms(); // this indirectly calls initPredecessors required for backtracking later
    updateToLocalGCEach(curGCIdx);
    for( int j = 1; j < dnalen; j++ ) {
      for( int i = 0; i < statecount; i++ ){
	if (i == synchstate) {
	    viterbi[j][i] = viterbi[j-1][i]/4;
	  if (needForwardTable)
	      forward[j][i] = forward[j-1][i]/4;
	} else {
	    viterbi[j].erase(i);
	    if (needForwardTable)
		forward[j].erase(i);
	}
      }
    }
    return;
  }

  if (show_progress){
      cerr << "viterbi algorithm progress:\n[%]: ";
      oldprogress = 0;
  }
#ifdef DEBUG
  double capacity=0;
  double fcapacity=0;
  cerr << "debug activated!" << endl;
  vector<string> statevec;
  for (int i=0; i<statecount; i++)
      statevec.push_back(stateTypeNames[getStateType(i)]);
  double beforeCount=0.0, afterCount=0.0;
  int mincount = 1000;
#endif

  initAlgorithms(); // update GC content dependent parameters
  for( int j = 1; j < dnalen; j++ ) { // TODO: this ignores the first nucleotide
      if (cs.idx[j] != curGCIdx) {// check whether GC content has changed, this is in particular the case at the very start
          curGCIdx = cs.idx[j];
	  updateToLocalGCEach(curGCIdx, j, cs.getNextStep(j) - 1); // update GC content dependent parameters
      }
      if (show_progress) {
	  progress = 1+100*j/dnalen;
	  if (progress > oldprogress){
	      if (progress%10<9 && progress%10>0)
		  cerr << ".";
	      else if (progress%10==0)
		  cerr << progress;
	  }
	  oldprogress = progress;
      }
      if (useProfile) 
	  profileModel->advanceScores(j);
      for( int i = 0; i < statecount; i++ ){
	  if (stateReachable[i]) {
	      states[i]->viterbiForwardAndSampling(viterbi, forward, i, j, doViterbi(needForwardTable), oli);
	  }
      }
      if (j % 1000 == 0) {
#ifdef DEBUG 
	  cerr << "[" << j;
#endif
	  int deletedSubmaps=0;
	  int newDeletionEnd = StateModel::getActiveWindowStart(j);
	  while (deletionEnd < newDeletionEnd) {
	      deletedSubmaps += viterbi[++deletionEnd].eraseUnneededSubstates();
#ifdef DEBUG
	      // viterbi[deletionEnd].checkAfterWindow();
#endif
	  }

#ifdef DEBUG
	  cerr << "]";
#endif 
      }
      
#ifdef DEBUG
      try {
	  viterbi[j].testIntegrity();
      } catch (...) {
	  cerr << "j=" << j << "\n";
	  throw;
      }
      capacity += viterbi[j].capacity();
      if (needForwardTable) 
	  fcapacity += forward[j].capacity();
#endif
      if (useProfile) {
#ifdef DEBUG
	  beforeCount += viterbi[j].totalSubstateCount();
#endif
	  profileModel->clearLowScoringSubstates(viterbi[j]);
#ifdef DEBUG
	  int newcount = viterbi[j].totalSubstateCount();
	  afterCount += newcount;
	  if (newcount > mincount)  {
//	      viterbi[j].showSubstateInfo(cerr);
	      mincount = newcount;
	  }
#endif
      }
#ifdef DEBUG
# ifdef DEBUG_BACK
      if (useProfile) 
	  for ( int i=0; i<statecount; i++) 
	      if (viterbi[j].hasSubstates(i)) { 
		  if (i==8 || i==25 || i==1 || i==24) {
		      cerr << "hit at " << j << "\n";
		      states[i]->viterbiForwardAndSampling(viterbi, forward, i, j,
							   doBacktracking, oli);
		  } else {
		      const ViterbiSubmapType& submap = viterbi[j].getSubstates(i);
		      for (ViterbiSubmapType::const_iterator it = submap.begin();
			   it != submap.end();
			   ++it) {
			  states[i]->viterbiForwardAndSampling(viterbi, forward, getFullStateId(i, it->first), j,
							       doBacktracking, oli);
		      }
		  }
		  try {
		      viterbi[j].getSubstates(i).checkCounts();
		  } catch (ProjectError) {
		      throw ProjectError("checkCounts() failed at base=" + itoa(j) + ", state=" + itoa(i) + " (" + stateTypeIdentifiers[stateMap[i]] + ")");
		  }
	      }
# endif
#endif
  }
#ifdef DEBUG
  if (useProfile) {
# ifdef DEBUG_DELALL
      while (deletionEnd < dnalen-1) 
      viterbi[++deletionEnd].eraseUnneededSubstates();
      // for (int j=dnalen-1; j>0; j--) 
      //     viterbi[j].eraseUnneededSubstates();
      viterbi.showSubstateMemoryUsage();
      exit(0);
# endif
      fprintf(stderr, "Number of substates reduced to %2.0f%%, average number was %.0f, is now %.0f.\n", afterCount/beforeCount*100, beforeCount/(dnalen-1), afterCount/(dnalen-1));
      profileModel->currentResult.clear();
      profileModel->outputFullBlocks();
        cerr << "\nMemory usage of Viterbi matrix (without substates): " << 
	  (dnalen*(sizeof(ViterbiColumnType) + statecount) + capacity*sizeof(Double))/1024 << "KB.\n";
      viterbi.showSubstateMemoryUsage(0, 2*deletionEnd - dnalen);
      viterbi.showSubstateMemoryUsage(2*deletionEnd - dnalen, deletionEnd);
      viterbi.showSubstateMemoryUsage(deletionEnd);
  }
#endif
  if (show_progress)
      cerr << endl;

  /*
   * The multiplication of the terminal probabilities is not done here, because some algorithms
   * call viterbiForwardAndSampling and thus change the viterbi and forward table again.
   * This probability needs to be multiplied at the appropriate places.
   */

} // end ViterbiAndForward

StatePath* NAMGene::getSampledPath(const char *dna, const char* seqname){
  StatePath *sampledPath = new StatePath();
  if (!needForwardTable)
      return sampledPath;
  if (seqname)
      sampledPath->seqname = seqname;

  /* determine the starting state of the sampling, i.e. the last state 
   * of the sampled state path
   */
  int dnalen = forward.size(); 
  int base, state;
  Double pathEmiProb = 1.0;
  if (containsJustNonNucs(dna, dnalen)) {
      sampledPath->push(new State(0, dnalen-1, getStateType(synchstate)));
    for (int i=0; i<dnalen; i++) 
	pathEmiProb *= (float) .25;
  } else { // normal case
      const ViterbiColumnType& lastCol = forward[dnalen-1];
      OptionsList ol;
      OptionListItem oli;
      for (int i=0; i < statecount; i++) 
	  if (lastCol[i] * termProbs[i]>0.0)
	      ol.add(i, dnalen-1, lastCol[i] * termProbs[i]);
      ol.prepareSampling();
      oli = ol.sample();
    
      base = oli.base;
      state = oli.state;
      pathEmiProb *= oli.probability;
      oli.reset();

      while (base>0) {
	  if (cs.idx[base] != curGCIdx) { // check whether GC content has changed
	    curGCIdx = cs.idx[base];
	    // update GC content dependent parameters,
	    // no new computation of SegProbs (empty interval [2,1]
	    updateToLocalGCEach(curGCIdx, 2,1);
	  }
	  states[state]->viterbiForwardAndSampling( viterbi, forward, state, base, 
						    doSampling, oli);
	  State *temp = new State(oli.base+1, base, getStateType(state));
	  temp->setTruncFlag(base, oli.base, dnalen);
	  sampledPath->push(temp);
	  if (Constant::overlapmode && oli.predEnd > -INT_MAX)
	      base = oli.predEnd;
	  else
	      base = oli.base;
	  state = oli.state;
	  pathEmiProb *= oli.probability;
	  oli.reset();
      }
  }

  /*
   * this pathEmiProb is (should be equal) to what getPathEmiProb computes on the sampledPath
   */
  sampledPath->pathemiProb = pathEmiProb;
  return sampledPath;
}

/*
 * NAMGene::getViterbiPath
 * Create the Viterbi-Path
 */
StatePath* NAMGene::getViterbiPath(const char *dna, const char* seqname){
  OptionListItem oli;
  int dnalen = viterbi.size();
  StatePath *viterbiPath = new StatePath();
  int state=-1;
  Double maxValue(0.0), value;
  if (seqname != NULL) 
    viterbiPath->seqname = seqname;

  // allow viterbi path to end in any state
  const ViterbiColumnType& lastcol = viterbi[dnalen-1];
  for (int i=0; i < statecount; i++) {
      SubstateId substate;
      if (profileModel && profileModel->allow_truncated) 
	  lastcol.getMaxSubstate(i, substate, value); // allow viterbi path to end in substates too
      else 
	  value = lastcol[i];
      value *= termProbs[i];
      if (value > maxValue) {
	  maxValue = value;
	  state = getFullStateId(i, substate);
      }
  }
  if (maxValue <= 0) {
      throw ProjectError("No feasible path found in HMM");
  }

  /*
   * Set the common probability of the viterbi path and the emission
   */
  viterbiPath->pathemiProb = maxValue;
#ifdef DEBUG
  cerr << "Viterbi Score: " << maxValue << "\n";
  cerr << "ln ""        : " << maxValue.log() << "\n";
#endif
    int base = dnalen-1;
    oli.reset();
    
    while (base>0) {
	int stateidx;
	SubstateId substate;  
	getStatePair(state, stateidx, substate);
#ifdef DEBUG_STATES
	if (stateMap[stateidx] != igenic && !isGeometricIntron(stateMap[stateidx]) &&
	    !isRGeometricIntron(stateMap[stateidx])) 
	{
	    cerr << "state=" << stateidx << ", stateType=" << stateTypeIdentifiers[stateMap[stateidx]]
		 << ", base=" << base << ", ln vit=" <<  viterbi[base].get(stateidx, substate).log() << endl;
	    if (isCodingExon(stateMap[stateidx])) 
		cerr << endl;
	}
#endif
	// need to test on long example  seqs with index changes (bug fixed on Mar 25th, 2015)
	if (cs.idx[base] != curGCIdx) { // check whether GC content has changed
	    curGCIdx = cs.idx[base];
	    // update GC content dependent parameters,
	    // no new computation of SegProbs (empty interval [2,1]
	    updateToLocalGCEach(curGCIdx, 2, 1);
	}
	states[stateidx]->viterbiForwardAndSampling( viterbi, forward, state, base, 
						     doBacktracking, oli );
	if (!Constant::overlapmode && ((oli.base >= base && oli.state == state) || (oli.base > base + 10))) {
	  throw ProjectError("Viterbi got stuck at state " + itoa(state) + "and base " + itoa(base) +
			     ".\n(oli.state =" + itoa(oli.state) + ", oli.base="+ itoa(oli.base) +")");
	}
	State *temp = new State(oli.base+1, base, getStateType(stateidx));
	temp->setTruncFlag(base, oli.base, dnalen);
	viterbiPath->push(temp);
	if (Constant::overlapmode && oli.predEnd > -INT_MAX)
	  base = oli.predEnd;
	else 
	  base = oli.base;
	state = oli.state;
	oli.reset();
    }
    if (profileModel)
	profileModel->appendMatchesTo(viterbiPath->proteinMatches);
    return viterbiPath;
}

/*
 * NAMGene::doViterbiPiecewise
 * Iterate the viterbi algorithm on pieces of DNA small enough so that the DP matrices fit into memory.
 */
Transcript *NAMGene::doViterbiPiecewise(SequenceFeatureCollection& sfc, AnnoSequence *annoseq, Strand strand){
  list<AltGene> *geneList = new list<AltGene>;
  char *dna = annoseq->sequence;
  char *curdna;
  AnnoSequence *curAnnoSeq;
  char *seqname = annoseq->seqname;
  vector<Double> *origInitProbs;
  vector<Double> *origTermProbs;
  vector<Double> *synchProbs;
//  SequenceFeatureCollection *partSFC;
  static int geneid = 1; // made this static so gene numbering goes across sequences and is unique
  int transcriptid;
  int curdnalen;
  
  bool singlestrand = false; // singlestrand = no shadow states
  try {
    singlestrand = Properties::getBoolProperty("singlestrand");
  } catch (...) {}
  
  int maxstep = 1000000;
  int endPos, beginPos;
  int seqlen = strlen(dna);
  try {
    maxstep = Properties::getIntProperty( "maxDNAPieceSize" );
  } catch (...) {}
  if (maxstep < 1000) {
      cerr << "maxDNAPieceSize is too small: " << maxstep << endl;
      exit(1);
  }

  beginPos = 0;
  try {
    synchstate = Properties::getIntProperty( "/NAMGene/SynchState" );  
  } catch (...) {
    cerr << "No synch state specified!" << endl;
  }

  /*
   * If we can't do it in one piece - we must do multiple pieces.
   * At the cutting points only the synchstate is allowed. So
   * change the initial and terminal probabilities for the cutting points.
   */
  origInitProbs = new vector<Double>(statecount, 0.0);
  origTermProbs = new vector<Double>(statecount, 0.0);
  synchProbs = new vector<Double>(statecount, 0.0);
    
  for (int i=0; i<statecount; i++) {
    (*origInitProbs)[i] = initProbs[i];
    (*origTermProbs)[i] = termProbs[i];
    if (i == synchstate)
      (*synchProbs)[i] = 1.0;
    else
      (*synchProbs)[i] = 0.0;
  }
 
  /*
   * loop over the pieces
   */
  do {
     endPos = getNextCutEndPoint(dna, beginPos, maxstep, sfc);
     if (show_progress) 
	 cerr << "examining piece " << beginPos + annoseq->offset + 1 << ".." 
	      << endPos + annoseq->offset + 1<< " (" << (endPos-beginPos+1)
	      << " bp)" << endl;
#ifdef DEBUG
      cout << "# examining piece " << beginPos + annoseq->offset + 1<< ".." << endPos + annoseq->offset + 1 << " (" << (endPos-beginPos+1) << " bp)" << endl;
#endif
    curAnnoSeq = new AnnoSequence();
    curdnalen = endPos-beginPos+1;
    curAnnoSeq->length = curdnalen;
    curdna = newstrcpy(dna+beginPos, curdnalen);
    curAnnoSeq->sequence = curdna;
    curAnnoSeq->offset = annoseq->offset + beginPos;

    /*
     * Set the initial and terminal probabilities.
     */
    for (int i=0; i<statecount; i++) {
      if (beginPos == 0)
	initProbs[i] = (*origInitProbs)[i];
      else 
	initProbs[i] = (*synchProbs)[i];
      if (endPos == seqlen-1)
	termProbs[i] = (*origTermProbs)[i];
      else 
	termProbs[i] = (*synchProbs)[i];
    }

    list<AltGene> *pieceGenes = new list<AltGene>;

    if (!singlestrand) {
	SequenceFeatureCollection partSFC(sfc, beginPos, endPos);
	pieceGenes = getStepGenes(curAnnoSeq, partSFC, strand);
    } else {
	if (strand == plusstrand || strand == bothstrands){
	    SequenceFeatureCollection partSFC(sfc, beginPos, endPos);
	    pieceGenes = getStepGenes(curAnnoSeq, partSFC, plusstrand);
	}
	if (strand == minusstrand || strand == bothstrands) {
	    list<AltGene> *pieceForwardGenes, *pieceReverseGenes;
	    char *reverseDNA = reverseComplement(curdna);
	    delete [] curAnnoSeq->sequence;
	    curAnnoSeq->sequence = reverseDNA;
	    SequenceFeatureCollection partSFC(sfc, beginPos, endPos, true);
	    pieceForwardGenes = getStepGenes(curAnnoSeq, partSFC, plusstrand);
	    pieceReverseGenes = reverseGeneList(pieceForwardGenes, endPos-beginPos);
	    pieceGenes->splice(pieceGenes->end(), *pieceReverseGenes);
	}
    }
    pieceGenes->sort();
    // shift gene coordinates, set sequence name, gene and transcript names
    for (list<AltGene>::iterator agit = pieceGenes->begin(); agit != pieceGenes->end(); ++agit){
	agit->shiftCoordinates(beginPos + annoseq->offset);
	agit->seqname = seqname;
     
	if (Constant::uniqueGeneId){
	    char buf[38];
	    sprintf(buf, "%.30s.g%d", seqname, geneid);
	    agit->id = buf;
	} else {
	    agit->id = "g"+itoa(geneid);
	} // TEMPorarily uncommented
      
	agit->sortTranscripts();

	transcriptid = 1;
	for(list<Transcript*>::iterator it = agit->transcripts.begin();it != agit->transcripts.end(); ++it ) {
	    (*it)->seqname = seqname;
	    (*it)->id = "t" + itoa(transcriptid);
	    (*it)->geneid = agit->id;
	    transcriptid++;
	}
	geneid++;
    }
    // print the genes
    printGeneList(pieceGenes, annoseq, Constant::codSeqOutput, Constant::proteinOutput, sfc.collection->hasHintsFile);

    // append pieceGenes to geneList
    geneList->splice(geneList->end(), *pieceGenes);

    
	
    delete curAnnoSeq;
    beginPos = endPos + 1;
  } while (beginPos < seqlen);

  /*
   * Reset the original initial and terminial probabilities.
   */
  for (int i=0; i<statecount; i++) {
    initProbs[i] = (*origInitProbs)[i];
    termProbs[i] = (*origTermProbs)[i];
  }
  delete origInitProbs;
  delete origTermProbs;
  delete synchProbs;
  if (geneList->empty() && !Constant::MultSpeciesMode)
      cout << "# (none)" << endl;
  return getPtr(geneList);
}

/*
 * NAMGene::getStepGenes
 */
list<AltGene> *NAMGene::getStepGenes(AnnoSequence *annoseq, SequenceFeatureCollection& sfc, Strand strand, bool onlyViterbi){
    list<AltGene> *genesAllHints, *genesPartialHints, *genes;
    list<list<AltGene> *> *results;
    const char *dna = annoseq->sequence;
    char *curdna;
    int curdnalen;
    int n = strlen(dna);
    bool uniqueCDS; 
    try {
	uniqueCDS = Properties::getBoolProperty("uniqueCDS");
    } catch(...) {
	uniqueCDS = false;
    }

    sfc.setSeqLen(n);
    sfc.computeHintedSites(dna);
    sfc.makeGroups();
    sfc.prepareLocalMalus(dna);
    //sfc->printGroups(); //AUSGABE
    StateModel::setSFC(&sfc);
    /*
     * First make a prediction run on the whole sequence using all hints at the same time.
     */
    genesAllHints = findGenes(dna, strand, onlyViterbi);
#ifdef DEBUG
    // printGeneList(genesAllHints, annoseq, false, false, false);
#endif
    results = new list<list<AltGene> *>;
    results->push_back(genesAllHints);
    sfc.createPredictionScheme(genesAllHints);
    PredictionScheme *scheme  = sfc.predictionScheme;
#ifdef DEBUG
    //scheme->print(annoseq->offset);
#endif
    /*
     * Then make further runs with certain hintgroups deactivated.
     */
#ifdef DEBUG
    int numruns = scheme->predictionRuns.size();
    cout << "Make " << numruns << " additional prediction runs." << endl; 
#endif

    for (list<PredictionRun>::iterator rit = scheme->predictionRuns.begin(); rit != scheme->predictionRuns.end(); ++rit){
#ifdef DEBUG
	cout << "Make run ";
	rit->print();
#endif
	curdnalen = rit->end - rit->begin + 1;
	curdna = newstrcpy(dna + rit->begin, curdnalen);
	sfc.setActiveFlag(rit->omittedGroups, false);
	sfc.shift(-rit->begin);
	sfc.setSeqLen(curdnalen);
	sfc.computeHintedSites(curdna);
	sfc.prepareLocalMalus(curdna);
        // TODO: set initial and terminal probs to allow partial genes again.
	genesPartialHints = findGenes(curdna, strand, onlyViterbi); 
	for (list<AltGene>::iterator agit = genesPartialHints->begin(); agit != genesPartialHints->end(); ++agit)
	    agit->shiftCoordinates(rit->begin);
#ifdef DEBUG
	printGeneList(genesPartialHints, annoseq, false, false, false);
#endif
	results->push_back(genesPartialHints);
	sfc.setSeqLen(n);
	sfc.shift(rit->begin); // shift back
	sfc.setActiveFlag(rit->omittedGroups, true);
	delete [] curdna;
    }
    genes = sfc.joinGenesFromPredRuns(results, maxtracks, uniqueCDS);
    // TODO: delete results 
    postProcessGenes(genes, annoseq); // truncate masked UTRs
    return genes;
}


/*
 * NAMGene::findGenes
 * onlyViterbi: if true, only the Viterbi transcripts are output
 * numSample: number of sample iterations, if onlyViterbi is true, these are used
 * to compute the a posteriori probabilities of the Viterbi transcripts, exons, introns.
 *
 */

list<AltGene> *NAMGene::findGenes(const char *dna, Strand strand, bool onlyViterbi){
  list<AltGene> *agl;
  list<Transcript*> alltranscripts;
  list<Transcript*> *filteredTranscripts = new list<Transcript*>;
  list<Transcript*>::iterator geneit1, geneit2;
  if (sampleiterations < 1)
      sampleiterations = 1;
  list<AltGene>::iterator agit;
  Transcript *genes = NULL, *g;
  StatePath *viterbiPath;
  StatePath *condensedViterbiPath;
  
  // compute the viterbi and forward table, main work done here
  viterbiAndForward(dna, profileModel);
#ifdef DEBUG
  cerr << "After viterbi: average load is " << viterbi.load() << ",\n"
       << "                      used are " << viterbi.used_load() << " of " << statecount << " states.\n";
  // cerr << "      Forward: average load is " << forward.load() << ",\n"
  //      << "                      used are " << forward.used_load() << " of " << statecount << " states.\n";
  cerr << "globalThreshUsedCount=" << globalThreshUsedCount << endl;
  globalThreshUsedCount=0;
#endif
 

  /*
   * add the viterbi transcripts to the list of genes
   */
  viterbiPath = getViterbiPath(dna, "");
  //getPathEmiProb(viterbiPath, dna); // for testing
  condensedViterbiPath = StatePath::condenseStatePath(viterbiPath);
  //cout << "Viterbi path:" << endl;
  //viterbiPath->print(); // for testing
  //cout << "Condensed Viterbi path:" << endl;
  //condensedViterbiPath->print(); // for testing
  //printDPMatrix();
  genes = condensedViterbiPath->projectOntoGeneSequence("g");
  delete viterbiPath;

  // clear viterbi matrix (not needed anymore) to save heap space;
  // the length will stay unchanged, but the entries at each base
  // position are deleted (this way we allow viterbiForwardAndSampling
  // to access the values so it won't cause errors; since in sampling
  // mode the values are ignored, we can assume from now on that they
  // are all 0).
  int dnalen = viterbi.size();
  viterbi.assign(dnalen, statecount); 
#ifdef DEBUG
  cerr << "Now we deleted the viterbi matrix.\n";
#endif

  for (g = genes; g != NULL; g = g->next) {
      g->apostprob = 1.0;        // add the viterbi transcripts with a very light weight
      g->setStatePostProbs(1.0); // just to ensure that they are present when sampling
      g->setSampleCount(1);
      g->hasProbs = true;
      g->throwaway = false;
      g->viterbi = true;
      alltranscripts.push_back(g);
  }

  int progress, oldprogress=0;
  if (sampleiterations > 1) {
      if (show_progress) {
	  cerr << "sampling algorithm progress:\n[%]: ";
	  oldprogress = 0;
      }
    /*
     * Sample and add the sampled genes to the list of genes
     */
    StatePath *sampledPath, *condensedsampledPath;
    for (int i=0; i < sampleiterations-1; i++) {
#ifdef DEBUG
	cerr << "Sample iteration " << i << endl;
#endif
      if (show_progress) {
	  progress = 100*i/(sampleiterations-2);
	  while (progress > oldprogress){
	      oldprogress++;
	      if (oldprogress%10<9 && oldprogress%10>0)
		  cerr << ".";
	      else if (oldprogress%10==0)
		  cerr << oldprogress;
	  }
      }
      // sample the transcripts 
      sampledPath = getSampledPath(dna, "");
      condensedsampledPath = StatePath::condenseStatePath(sampledPath);
      // condensedsampledPath->print(); // for testing
      char gr[9];
      sprintf(gr, "s%d-", (i+1));
      genes = condensedsampledPath->projectOntoGeneSequence(gr);
      delete sampledPath;
      delete condensedsampledPath;
      
      //cout << "i=" << i << endl;
      // store all sampled transcripts
      for (g = genes; g != NULL; g = g->next) {
	  //  g->printGFF();
        // initialize, all apostprobs and counts must be set to 1 
	g->apostprob = 1.0;
	g->setStatePostProbs(1.0);
	g->setSampleCount(1);
	g->hasProbs = true;
	g->viterbi = false;
	if (!alternatives_from_sampling)
	  g->throwaway = true; // no alternatives requested, throw sampled gene away later
	alltranscripts.push_back(g);
      }
    } // for i<sampleiterations
    if (show_progress)
	cerr << endl;
        
    alltranscripts.sort(ptr_comparison<Transcript>());
    // now remove multiple copies and increase apostprob instead
    for (geneit1 = alltranscripts.begin(); geneit1 != alltranscripts.end();){
	geneit2 = geneit1;
	for (geneit2++; geneit2 != alltranscripts.end() && (*geneit2)->geneBegin() == (*geneit1)->geneBegin();)
	    if (**geneit1 == **geneit2){
		(*geneit1)->throwaway &= (*geneit2)->throwaway; // throw away only if all versions should be thrown away.
		(*geneit1)->viterbi |= (*geneit2)->viterbi;
		(*geneit1)->addSampleCount(1);
		(*geneit1)->apostprob += 1.0;
		(*geneit1)->addStatePostProbs(1.0);
		delete *geneit2; // free space for duplicated transcript
		geneit2 = alltranscripts.erase(geneit2); // remove deleted pointer from the list as well
	    } else
		geneit2++;
	geneit1++;
    }
   
    /*
     * compute the posterior probabilities of the states and transcripts
     */
    for (geneit1 = alltranscripts.begin(); geneit1 != alltranscripts.end(); geneit1++){
	geneit2 = geneit1;
	geneit2++;
	for (; geneit2 != alltranscripts.end() && (*geneit2)->geneBegin() <= (*geneit1)->geneEnd(); geneit2++)
	    (*geneit1)->updatePostProb(*geneit2);
    }
    for (geneit1 = alltranscripts.begin(); geneit1 != alltranscripts.end(); geneit1++){
	(*geneit1)->normPostProb(sampleiterations); // +1 wegen Viterbipfad
    }
  }
  /*
   * filter transcripts by probabilities, strand
   */ 
  filterGenePrediction(alltranscripts, *filteredTranscripts, dna, strand, noInFrameStop, minmeanexonintronprob, minexonintronprob);

  // determine transcripts with maximum expected accuracy criterion
  if (Constant::MultSpeciesMode){
      // store alltranscripts in member variable of NAMGene
      setAllTranscripts(filteredTranscripts);
      agl = new list<AltGene>;
  } else if (mea_prediction){
      agl = groupTranscriptsToGenes(getMEAtranscripts(*filteredTranscripts, dna));
  } else { //filter transcripts by maximum track number 
      agl = groupTranscriptsToGenes(*filteredTranscripts);
  }
  
  if (sampleiterations>1) {
    /*
     * compute apostprob of genes. The posterior 'probability' of a gene is the 
     * expected number of transcripts overlaping the region of the gene on 
     * the given strand. This may now be larger than 1 (changed March 15, Mario).
     */
    for (agit = agl->begin(); agit != agl->end(); ++agit) {
	for (list<Transcript*>::iterator tit = alltranscripts.begin(); tit != alltranscripts.end() && (*tit)->geneBegin() <= agit->maxcodend; ++tit){
	    if (agit->overlaps(*tit))
		agit->apostprob += (*tit)->apostprob;
	}
    }
  } else { // no sampling, set hasProbs to false
    for (agit = agl->begin(); agit != agl->end(); ++agit) {
      agit->hasProbs = false;
      for (list<Transcript*>::iterator git = agit->transcripts.begin(); git != agit->transcripts.end(); ++git){
	(*git)->hasProbs = false;
	(*git)->setStateHasScore(false);
      }
    }
  }

  delete condensedViterbiPath;
 
#ifdef DEBUG
  cerr << "At the end of findGenes: average load is " << viterbi.load() << ",\n"
       << "                                used are " << viterbi.used_load() << " of " << statecount << " states.\n";
  cerr << "                Forward: average load is " << forward.load() << ",\n"
       << "                                used are " << forward.used_load() << " of " << statecount << " states.\n";
#endif
  return agl;
  // TODO: the *filteredTranscripts has to be deleted eventually after cmpgenepred has called .getAllTranscripts()
}


/* 
 * NAMGene::getNextCutEndPoint 
 * Find a suitable point to cut dna into pieces for separate analysis. 
 * Cutting point should be at most maxstep from dna + beginPos, as large as possible  
 * but also likely to be in the intergenic region (synch state). 
 */ 
int NAMGene::getNextCutEndPoint(const char *dna, int beginPos, int maxstep, SequenceFeatureCollection& sfc){
  int restlen = strlen(dna+beginPos);
  char *curdna; 
  int cutendpoint=0;
  int examChunkSize=50000;
  int examIntervalStart, examIntervalEnd;
  SequenceFeatureCollection *partSFC;
  list<Feature> *groupGaps;
  list<Feature>::iterator lastirit;
  
  if (examChunkSize < 0.2*maxstep)
    examChunkSize = (int) (0.2*maxstep);
  if (examChunkSize > 150000)
    examChunkSize = 150000;

  if (restlen <= maxstep){
    return beginPos + restlen-1;
  } else {
      /*
       * try to find a good breakpoint
       * 1st try: window of size examChunkSize around last groupgap in range or at end of range if groupgap does not exist.
       *          If have two genes with intergenic region in groupgap: take center of largest such ir.    
       * 2nd try: Double the size of the window. Try second to last groupgap.
       *          If have two genes with intergenic region in groupgap: take center of largest such ir.
       *          Take internal predicted intergenic region with largest overlap to a groupgap if exists.
       *          Take the larger predicted marginal intergenic region in groupgap.
       *          Give up: Take right boundary of range.
       */
      
      groupGaps = new list<Feature>;
      if (sfc.groupGaps) {
	  for (list<Feature>::iterator it = sfc.groupGaps->begin(); it != sfc.groupGaps->end() && it->start <= beginPos+maxstep; ++it)
	      if (it->end > beginPos && it->start <= beginPos+maxstep)
		  groupGaps->push_back(*it);
      } else {
	  groupGaps->push_back(Feature(1, strlen(dna)-1, irpartF, bothstrands, -1, "groups"));
      }

      int centerLastGap;
      if (!groupGaps->empty()) {
	  lastirit = groupGaps->begin();
	  while (lastirit != groupGaps->end() && lastirit->start < beginPos + maxstep)
	      ++lastirit;
	  if (lastirit != groupGaps->begin())
	      --lastirit;
	  if (lastirit->end - lastirit->start < examChunkSize)
	      centerLastGap = (lastirit->end + lastirit->start) / 2;
	  else 
	      centerLastGap = lastirit->end - examChunkSize/2; // for huge gaps, take the right subwindow of size examChunkSize
      } else {
	  centerLastGap = beginPos + maxstep - 1;
      }
      if (examChunkSize > maxstep) {
	  examIntervalStart = beginPos;
	  examIntervalEnd = beginPos + maxstep - 1;
      } else {
	  examIntervalStart = centerLastGap - examChunkSize/2;
	  examIntervalEnd   = centerLastGap + examChunkSize/2;
      
	  if (examIntervalEnd >= beginPos + maxstep){ // too far right: shift left
	      examIntervalStart -= (examIntervalEnd - (beginPos + maxstep-1));
	      examIntervalEnd   = beginPos + maxstep - 1;
	  }
	  if (examIntervalStart < beginPos) {// too far left: shift right
	      examIntervalEnd += beginPos - examIntervalStart;
	      examIntervalStart = beginPos;
	  }
      }

      /*
       * 1st try: window of size examChunkSize around last groupgap in range or at end of range if groupgap does not exist.
       *          If have two genes with intergenic region in groupgap: take center of largest such ir.
       */
      StatePath *viterbiPath; 
      StatePath *condensedViterbiPath;

      curdna = newstrcpy(dna + examIntervalStart, examIntervalEnd-examIntervalStart+1);
      partSFC = new SequenceFeatureCollection(sfc, examIntervalStart, examIntervalEnd + 10000);
      partSFC->setSeqLen(examIntervalEnd-examIntervalStart+1);
      partSFC->computeHintedSites(curdna);
      partSFC->prepareLocalMalus(curdna);
      StateModel::setSFC(partSFC);
      viterbiAndForward(curdna);      // do not use protein profile here
      viterbiPath = getViterbiPath(curdna, "temp");
      delete partSFC; //achtung, das muß nach getViterbiPath stehen
      delete [] curdna;
      condensedViterbiPath = StatePath::condenseStatePath(viterbiPath);
      //cout << "getNextCutEndPoint viterbi path:" << endl;
      //condensedViterbiPath->print();
      cutendpoint = tryFindCutEndPoint(condensedViterbiPath, examIntervalStart, examIntervalEnd, groupGaps, true);
      delete viterbiPath;
      delete condensedViterbiPath;
      if (cutendpoint != -1) { // found cutendpoint in first try
	  //cout << "# found cutendpoint in first try : " << cutendpoint << endl;
      } else {
	  /* 
	   * 2nd try: Double the size of the window. Try second to last groupgap.
	   *          Take internal predicted intergenic region with largest overlap to a groupgap if exists.
	   *          Take the larger predicted marginal intergenic region in groupgap.
	   *          Take the largest intergenic region (don't care about groupgaps).
	   *          Give up: Take right boundary of range.
	   */
	  examChunkSize *=2;
	  if (examChunkSize > maxstep)
	      examChunkSize = maxstep;
	  if (!groupGaps->empty()) {
	      // take the last but one group gap if existent and if the last group gap is not larger than the previous exam interval
	      if (lastirit->start >= examIntervalStart && lastirit != groupGaps->begin())
		  lastirit--;

	      if (lastirit->end - lastirit->start < examChunkSize)
		  centerLastGap = (lastirit->end + lastirit->start) / 2;
	      else 
		  centerLastGap = lastirit->end - examChunkSize/2; // for huge gaps, take the right subwindow of size examChunkSize
	  } else {
	      centerLastGap = beginPos + maxstep - 1;
	  }
	  examIntervalStart = centerLastGap - examChunkSize/2;
	  examIntervalEnd   = centerLastGap + examChunkSize/2;
	  
	  if (examIntervalEnd >= beginPos + maxstep){ // too far right: shift left
	      examIntervalStart -= (examIntervalEnd - (beginPos + maxstep-1));
	      examIntervalEnd   = beginPos + maxstep - 1;
	  }
	  if (examIntervalStart < beginPos){ // too far left: shift right
	      examIntervalEnd += beginPos - examIntervalStart;
	      examIntervalStart = beginPos;
	  }

	  curdna = newstrcpy(dna + examIntervalStart, examIntervalEnd-examIntervalStart+1);
	  partSFC = new SequenceFeatureCollection(sfc, examIntervalStart, examIntervalEnd + 10000);
	  partSFC->setSeqLen(examIntervalEnd-examIntervalStart+1);
	  partSFC->computeHintedSites(curdna);
	  partSFC->prepareLocalMalus(curdna);
	  StateModel::setSFC(partSFC);
	  viterbiAndForward(curdna); // do not use protein profile here
	  viterbiPath = getViterbiPath(curdna, "temp");
	  delete partSFC; // note: this has to be called AFTER getViterbiPath
	  delete [] curdna;
	  condensedViterbiPath = StatePath::condenseStatePath(viterbiPath);
	  cutendpoint = tryFindCutEndPoint(condensedViterbiPath, examIntervalStart, examIntervalEnd, groupGaps, true);
	  if (cutendpoint == -1) {
	      cutendpoint = tryFindCutEndPoint(condensedViterbiPath, examIntervalStart, examIntervalEnd, groupGaps, false);
	      if (cutendpoint == -1) {
		  cutendpoint = tryFindCutEndPoint(condensedViterbiPath, examIntervalStart, examIntervalEnd, NULL, false);
		  if (cutendpoint == -1) {
		      cutendpoint = beginPos + maxstep - 1;
		  }
	      }
	  }
	  delete viterbiPath; 
	  delete condensedViterbiPath; 
      }
      
      delete groupGaps;
    
  } 
  if (cutendpoint <= beginPos + 0.05 * maxstep || cutendpoint <= beginPos + 5000) // move by at least 5% and by at least 5000bp
      cutendpoint = beginPos + maxstep - 1;
  return cutendpoint; 
} 

/*
 * NAMGene::tryFindCutEndPoint
 * 
 * returns the center of the largest intersection of a predicted intergenic region with a groupgap
 * only considers internal intergenic regions if onlyInternalIR=true, except if a marginal IR is at least half the size of the range 
 * If groupGaps is NULL, don't intersect with groupgaps. Simply find the largest IR.
 * return -1. if no such IR found.
 * irsize is the size of the likely intergenic region
 */ 

int NAMGene::tryFindCutEndPoint(StatePath *condensedExamPath, int examIntervalStart, int examIntervalEnd, list<Feature> *groupGaps, bool onlyInternalIR){
    // determine the internal intergenic region that has largest overlap with a groupGap (if any)
    int irbegin, irend;
    int maxirbegin = -1, maxirend = -1;
    State *st = condensedExamPath->first;
    bool groupGapsNULL = (groupGaps == NULL);
    if (groupGapsNULL){
	groupGaps = new list<Feature>;
	groupGaps->push_back(Feature(0, INT_MAX, irpartF, bothstrands, -1, "groups"));
    }
    if (st) {
	while (st) {
	    if(st->type == igenic){ // interior intergenic region
		irbegin = examIntervalStart + st->begin;
		irend = examIntervalStart + st->end;
		// intersect intergenic region with groupgaps
		int lgbegin=-1, lgend=-1;
		for (list<Feature>::iterator it = groupGaps->begin(); it != groupGaps->end(); ++it){
		    if (it->start<irbegin && it->end<=irend && it->end>=irbegin && it->end-irbegin> lgend-lgbegin) { 
			// predicted ir          -----------------
			// groupgap          ------------
			// intersection          --------
			lgbegin = irbegin;
			lgend = it->end;
		    } else if (it->start < irbegin && it->end>irend && irend-irbegin > lgend-lgbegin) {
			// predicted ir          -----------------
			// groupgap          -----------------------
			// intersection          -----------------
			lgbegin = irbegin;
			lgend = irend;
		    } else if (it->start > irbegin && it->end < irend && it->end-it->start > lgend-lgbegin) {
			// predicted ir          -----------------
			// groupgap                 ---------       
			// intersection             ---------     
			lgbegin = it->start;
			lgend = it->end;
		    } else if (it->start >= irbegin && it->start <= irend && it->end >= irend && irend-it->start > lgend-lgbegin) {
			// predicted ir          -----------------
			// groupgap                           ---------       
			// intersection                       ----
			lgbegin = it->start;
			lgend = irend;
		    } 
		}
		
		if (lgend-lgbegin > maxirend-maxirbegin &&
		    ((st != condensedExamPath->first && st->next != NULL) || // internal intergenic region
		     !onlyInternalIR || // allow marginal ir when onlyInternalIR is false
		     // or when the size is at least half the size of the examined range
		     lgend - lgbegin > (examIntervalEnd - examIntervalStart)/2)){
		    maxirbegin = lgbegin;
		    maxirend = lgend;
		}
	    }
	    st = st->next;
	} 
    }
    if (groupGapsNULL)
	delete groupGaps;

    if (maxirend-maxirbegin>0) {
	return (maxirend + maxirbegin)/2;
    } else {
	return -1;
    }
}

/*
 * NAMGene::getTrainViterbiPath
 * 
 * Does the whole job of finding the most likely path for training.
 * No splitting into chunks. No alternatives. Not using protein profile.
 */
StatePath* NAMGene::getTrainViterbiPath(const char *dna, SequenceFeatureCollection *sfc){
    bool crf = inCRFTraining;
    inCRFTraining = false; // this way Viterbi will be much more efficient because of use of static variables
    sfc->computeHintedSites(dna);
    StateModel::setSFC(sfc);
        viterbiAndForward(dna); // do not use protein profile here
    inCRFTraining = crf;
    return getViterbiPath(dna, "");
}


/*
 * NAMGene::getEmissionProbability()
 * 
 * returns the probability of the emission in the last run of the forward algorithm
 */ 

Double NAMGene::getEmissionProbability(){
    if (!needForwardTable)
	return 0.0;
    Double alpha = 0.0;
    const ViterbiColumnType& lastCol = forward[forward.size()-1];
    for (int i=0; i < statecount; i++) {
	alpha += lastCol[i] * termProbs[i];
    }
    return alpha;
}

/*
 * NAMGene::getPathEmiProb
 * compute the common probability of the path and the emission
 *
 * Only count the emissionprobabilities in the interval [countStart, coundEnd]
 * | state1 | | state3 | | | | |        state8                   | | | | state 12   | | |
 *                                        |                                           |  
 *                                      countStart                                countEnd
 *                                        <--      partProb   --->
 *                             <---        emi                --->
 */

Double NAMGene::getPathEmiProb(StatePath *path, const char *dna, SequenceFeatureCollection& sfc, int countStart, int countEnd){
    State *curstate, *laststate = NULL;
    Double pathemi = 1.0, emi;
    int stateidx;
    StateModel::setSFC(&sfc);
    if (!cs.dna || strcmp(dna,cs.dna))
      cs.computeStairs(dna); //(re)initialize GC content classes
    curGCIdx = -1; // initialize with invalid GC content class
    if (countEnd == -1) { // count all emission probabilities
	countStart = -100;
	countEnd = strlen( dna ) + 100;
    }
    for (curstate = path->first; curstate != NULL; curstate = curstate->next){
	stateidx = getStateIndex(curstate->type);
	int pos = curstate->end;
	if (pos > cs.n-1)
	    pos = cs.n-1; // some states, e.g. exons, may extend past the sequence end
	if (pos < 0)
	    pos = 0; // can happen with a start codon at the very beginning of the sequence
	if (cs.idx[pos] != curGCIdx) {
          curGCIdx = cs.idx[pos];
	  updateToLocalGCEach(curGCIdx, 2, 1); // update GC content dependent parameters
	}
	emi = states[stateidx]->emiProbUnderModel(curstate->begin, curstate->end);
	curstate->prob = emi;
	if (curstate->begin >= countStart && curstate->end <= countEnd){ // fragment lies completely in count region
	    pathemi *= emi;
	    //cout << "emi " << curstate->begin << " " <<  curstate->end << " " << emi << endl;
	} else if (!(curstate->end < countStart || curstate->begin > countEnd)){
	    // nontrivial overlap, only count a fraction corresponding to the overlap
	    int ovlp = ((curstate->end > countEnd)? countEnd : curstate->end) -
		((curstate->begin > countStart)? curstate->begin : countStart) + 1;
	    Double partProb = exp((emi.log() * ovlp / (curstate->end - curstate->begin + 1)));
	    pathemi *= partProb;
	} // otherwise ignore this emission probability
	
	if (laststate && curstate->begin >= countStart && curstate->end <= countEnd) {
	    // multiply also by the transition probability of the transition "laststate->curstate"
	    pathemi *= transitions[getStateIndex(laststate->type)][stateidx];
	} else { 
	    // first state: multiply also by the initial probability
	    /*
	     * HACK: if the first state is an exon with a negative start position use a transition from
	     * igenic region (=0) to this exon also.
	     */
	    if (initProbs[stateidx] > 0.0 && curstate->begin >= countStart && curstate->end <= countEnd)
		pathemi *= initProbs[stateidx];
	    //else // assume viterbi was correct, not correct, TODO 
	    //pathemi *= initProbs[0] * transitions[0][curstate->type]; //igenic region -> curstate
	}
	laststate = curstate;
    }
    // multiply with terminal probability
    if (laststate->begin >= countStart && laststate->end <= countEnd){
	pathemi *= termProbs[getStateIndex(laststate->type)];
    }
    path->pathemiProb = pathemi;
    return pathemi;
}

void NAMGene::readTransAndInitProbs( ){
  ifstream istrm;
  string speciesValue = Properties::getProperty(SPECIES_KEY);
  if (!inCRFTraining)
      cout << "# " << speciesValue << " version.";
  try {
      // look if there is a species specific transition file and take it
      string fname = Constant::fullSpeciesPath() + speciesValue + "_" + Properties::getProperty(TRANSFILE_KEY);
      istrm.open(fname.c_str());
      if (istrm.peek() == EOF) // file doesn't exist
	  throw ProjectError(fname + " doesn't exist");
      if (!inCRFTraining)
	  cout << " Using species specific transition matrix: " << fname;
  } catch (...){
      istrm.clear();
      // use default transition parameter file
      string fname = Constant::modelPath() + Properties::getProperty(TRANSFILE_KEY);
      istrm.open(fname.c_str());
      if (!inCRFTraining)
	  cout << " Using default transition matrix.";
  }
  cout << endl;

  if( istrm ){
    int i, j, k;
    int count;
    int numStartPossible=0, numEndPossible=0;
    transitions.assign(statecount, statecount);
    initProbs.assign(statecount, 0.0);
    termProbs.assign(statecount, 0.0);
    stateReachable.assign(statecount, true);

    /* First determine the number of states specified in the file.
     * Throw an exception if this value isn't equal to statecount from the 
     * properties.
     */
    istrm >> comment >> count;
    if( count != statecount ){
      istrm.close();
      throw NAMGeneError( "Incorrect state count in transition file!!!" );
    }

    /*
     * read in the initial probabilities
     * THIS IS A HACK. Put this into a configuration file.
     */
    istrm >> goto_line_after( "[Initial]" ) >> comment >> numStartPossible;
    for (k=0; k<numStartPossible; k++ ) {
      istrm >> comment >> i;
      istrm >> initProbs[i]; 
    }
    istrm >> goto_line_after( "[Terminal]" ) >> comment >> numEndPossible;
    for (k=0; k < numEndPossible; k++ ) {
      istrm >> comment >> i;
      istrm >> termProbs[i];
    }
	
    /*
     *  read in the transition probabilities
     */
    istrm >> goto_line_after( "[Transition]" );
    while( istrm ){
      istrm >> comment >> i >> j ;
      if (i >= statecount || j >= statecount){
	cerr << "transition:" << i << "->" << j << endl;
	throw NAMGeneError( "State number in transition file too large!");
      }
      istrm >> transitions[i][j];
    }

    istrm.close();
  } else {
    throw NAMGeneError( "Could't open the file with transition probabilities." );
  }
}

void NAMGene::readOvlpLenDist( ){
  if (!Constant::overlapmode)
    return;
  // this function is only relevant for prokaryotes (--genemodel=bacterium)
  ifstream istrm;
  string speciesValue = Properties::getProperty(SPECIES_KEY);
  try {
      // look if there is a species specific ovlp_len.pbl file and take it
      string fname = Constant::fullSpeciesPath() + speciesValue + "_" + OVLPLENFILE;
      istrm.open(fname.c_str());
      if (istrm.peek() == EOF) // file doesn't exist
	  throw ProjectError(fname + " doesn't exist");
      cout << "# Using species specific overlap length distribution: " << fname << endl;
  } catch (...){
      istrm.clear();
      // use default overlap length distribution file
      string fname = Constant::modelPath() + OVLPLENFILE;
      istrm.open(fname.c_str());
      cout << "# Using default overlap length distribution file." << endl;
  }
  
  if (Constant::maxOvlp >= Constant::min_coding_len){
      cerr << string("WARNING: maxOvlp(") +  itoa(Constant::maxOvlp) +
	  string(") is not smaller than /Constant/min_coding_len(")
	  + itoa(Constant::min_coding_len) + ").\n" +
	  "This may result in nested genes, either properly nested, or several genes with the same stop codon.\n" + 
	  "You may want to change this by decreasing maxOvlp or increasing /Constant/min_coding_len." << endl;
  }
  if (istrm) {
    Constant::head2tail_ovlp.assign(Constant::maxOvlp+1, 0.0);
    istrm >>  goto_line_after( "[HEAD2TAIL]" );
    istrm >> comment;
    for (int i=0; i <= Constant::maxOvlp; i++)
      istrm >> Constant::head2tail_ovlp[i];
    Constant::head2head_ovlp.assign(Constant::maxOvlp+1, 0.0);
    istrm >>  goto_line_after( "[HEAD2HEAD]" );
    istrm >> comment;
    for (int i=0; i <= Constant::maxOvlp; i++)
      istrm >> Constant::head2head_ovlp[i];
    Constant::tail2tail_ovlp.assign(Constant::maxOvlp+1, 0.0);
    istrm >>  goto_line_after( "[TAIL2TAIL]" );
    istrm >> comment;
    for (int i=0; i <= Constant::maxOvlp; i++)
      istrm >> Constant::tail2tail_ovlp[i];
    istrm.close();
  } else {
    throw NAMGeneError( "Could't open the file with the overlap length distribution." );
  }
}

#define EPSILON_CONSISTENCY 0.00000001
void NAMGene::checkProbsConsistency( ){
  Double sum;
  int i,j;
    
  /*
   * check whether initProbs really is a vector of probabilities.
   */
  sum=0;
  for (i=0; i< statecount; i++) {
    sum += initProbs[i];
  }
#ifdef DEBUG
  if (sum > Double(1.0) + EPSILON_CONSISTENCY || sum < Double(1.0) - EPSILON_CONSISTENCY)
    cerr << "Warning: Sum of all initial probabilities is 1 + delta, delta=" << setprecision(7) << sum-1.0 << ". Not equal to 1." << endl;
#endif
  
  /*
   * check whether termProbs really is a vector of probabilities.
   */
  sum=0;
  for (i=0; i< statecount; i++) {
    sum += termProbs[i];
  }
#ifdef DEBUG
  if (sum > Double(1.0) + EPSILON_CONSISTENCY || sum < Double(1.0) - EPSILON_CONSISTENCY)
    cerr << "Warning: Sum of all terminal probabilities is 1 + delta, delta=" << setprecision(7) << sum-1.0 << ". Not equal to 1." << endl;
#endif

  /*
   * check whether transitions really is a vector of probabilities.
   */
  for (i=0; i<statecount; i++) {  
    sum=0;
    for (j=0; j<statecount; j++) {
      sum += transitions[i][j];
      if (transitions[i][j]<0)
	  cerr << "Error: transition probability from state " << i << " to " << j << " is negative:" << transitions[i][j] << endl;
    }
#ifdef DEBUG
    if (sum > Double(1.0) + EPSILON_CONSISTENCY || sum < Double(1.0) - EPSILON_CONSISTENCY) {
	cerr << "Warning: Sum of transition probabilities out of state "<< i << " is 1 + delta, delta=" << setprecision(7) << sum-1.0
	     << ". Not equal to 1." << endl;
	cerr << sum << "\n";
	cerr << "Transitions are: " << transitions[i][0];
	for (j=1; j<statecount; j++) 
	    cerr << ", " << transitions[i][j];
	cerr << "\n";
    }
#endif
  }
  for (j=0; j<statecount; j++) {
    sum=0;
    for (i=0; i<statecount; i++) {
      sum += transitions[i][j];
    }
#ifdef DEBUG
    if (sum < EPSILON_CONSISTENCY) {
      cerr << "Warning: Sum of all transition probabilities into state "<< j << " is 0." << endl;
    }
#endif
  }
}

void NAMGene::computeReachableStates( ){
    bool expanded = true;
    bool allReachable;
    for (int i=0; i < statecount; i++)
	stateReachable[i] = (initProbs[i] > 0.0);

    while (expanded) {
	expanded = false;
	for (int i=0; i<statecount; i++) {
	    for (int j=0; j<statecount; j++) {
		if (stateReachable[i] && !stateReachable[j] && transitions[i][j] > 0.0) {
		    stateReachable[j] = expanded = true;
		}
	    }
	}
    }
    allReachable = true;
    for (int i=0; i < statecount; i++)
	allReachable &= stateReachable[i];
#ifdef DEBUG
    if (!allReachable) {
	cerr << "The following states are not reachable: ";
	for (int i=0; i < statecount; i++)
	    if (!stateReachable[i])
		cerr << i << ", ";
	cerr << endl << "If this is not on purpose you should check the transition matrix and the initial probabilities." << endl;
    }
#endif
}

void NAMGene::createStateModels( ){
    states.resize( statecount );
    StateModel::resetModelCounts();
    
    for( int i = 0; i < statecount; ++i ){
	/*
	 * Determine the names of the states and create the state objects
	 */ 
	const char* statename = Properties::getProperty("/NAMGene/state", i);
	states[i] = StateModel::newStateModelPtr(statename);
    }
}

void NAMGene::initAlgorithms(){
  StateModel::resetPars();
  for ( int i = 0; i < statecount; ++i ){
      states[i]->initAlgorithms(transitions, i);
  }
  /*
   * give the states the possibly implicitly changed transition matrix above 
   * in order to contruct the lists of predecessors
   *
   */
  for( int i = 0; i < statecount; i++ )
    states[i]->initPredecessors(transitions, i);
}

void NAMGene::updateToLocalGCEach(int idx, int from, int to){
    StateModel::updateToLocalGC(idx, from, to); 
    for ( int i = 0; i < statecount; ++i )
	states[i]->updateToLocalGCEach(transitions, i);
    /*
     * give the states the possibly implicitly changed transition matrix above 
     * in order to contruct the lists of predecessors
     *
     */
    for( int i = 0; i < statecount; i++ )
	states[i]->initPredecessors(transitions, i);

}

void NAMGene::prepareModels(const char* dna, int dnalen) {
    StateModel::prepareViterbi(dna, dnalen, stateMap);
#ifdef DEBUG
    PP::ExonScorer::transitions = &transitions;
#endif
    PP::SubstateModel::setStateStrands(stateMap);
}

/*
 * NAMGene::setPathAndProb
 * compute the emission probability of the path induced by the annotated genes
 * for scoreTx option to augustus
 */

void NAMGene::setPathAndProb(AnnoSequence *annoseq, FeatureCollection &extrinsicFeatures){
    while (annoseq){
	annoseq->anno->path = StatePath::getInducedStatePath(annoseq->anno->genes, annoseq->length);
	SequenceFeatureCollection &sfc = extrinsicFeatures.getSequenceFeatureCollection(annoseq->seqname);
	sfc.prepare(annoseq, false);
	sfc.computeHintedSites(annoseq->sequence);
	sfc.prepareLocalMalus(annoseq->sequence);
	prepareModels(annoseq->sequence, annoseq->length);
	annoseq->anno->emiProb = getPathEmiProb(&(*(annoseq->anno->path)), &(*(annoseq->sequence)), sfc);
	annoseq = annoseq->next;
    }
}

void NAMGene::printDPMatrix(){
    cout << "----- Viterbi matrix: -----" << endl;
    for (int t=0; t <viterbi.size(); t++){
	cout << t;
	for (int s = 0; s < statecount; s++)
	    cout << "\t" << stateTypeIdentifiers[stateMap[s]] << ": " << viterbi[t][s];
	cout << endl;
    }
}
 
