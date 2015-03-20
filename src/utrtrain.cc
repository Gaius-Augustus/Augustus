/**********************************************************************
 * file:    utrtrain.cc
 * licence: 
 *          
 * descr.:  training of untranslated region model
 * authors: Mario Stanke (mario@gobics.de)
 * 
 * date    |   author      |  changes 
 * --------|---------------|------------------------------------------ 
 * 01.02.06| Mario Stanke  | creation of the file
 * 27.03.06| Mario Stanke  | 5' UTR intron
 * 15.05.06| Mario Stanke  | 3' UTR training
 **********************************************************************/

#include "utrmodel.hh"

// project includes
#include "commontrain.hh"
#include "gene.hh"
#include "motif.hh"
#include "gene.hh"
#include "properties.hh"
#include "intronmodel.hh" // so content model can be reused here

// standard C/C++ includes
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>


/*
 * UtrModel::buildModel
 */
void UtrModel::buildModel( const AnnoSequence* annoseq, int parIndex){

    if (verbosity)
	cout << "** building model for untranslated regions. *UTR*" << endl;
    if (!annoseq->anno) 
	throw ProjectError("Tried to train Utrmodel with 0 sequences.");
    
    tssMotifTATA = new Motif(tss_end+tss_start, 0, 1);                  // winsize memory pseudocount neighbors
    tssMotif     = new Motif(tss_end+tss_start, 0, 1);                  // winsize memory pseudocount neighbors
    ttsMotif     = new Motif(Constant::d_polyasig_cleavage, tts_motif_memory, 1, 1);   // winsize memory pseudocount neighbors (memory 1 better than 0)
    tataMotif    = new Motif(tata_start+tata_end, 0, tata_pseudocount); // winsize memory pseudocount neighbors
    
    if (!hasLenDist) {
	lenCount5Single.assign(exonLenD+1,0);
	lenCount5Initial.assign(exonLenD+1,0);
	lenCount5Internal.assign(exonLenD+1,0);
	lenCount5Terminal.assign(exonLenD+1,0);
	numHuge5Terminal = numHuge5Initial = numHuge5Internal = numHuge5Single = 0;
	
	lenCount3Single.assign(exonLenD+1,0);
	lenCount3Initial.assign(exonLenD+1,0);
	lenCount3Internal.assign(exonLenD+1,0);
	lenCount3Terminal.assign(exonLenD+1,0);
	numHuge3Terminal = numHuge3Initial = numHuge3Internal = numHuge3Single = 0;
    }
    // Initialize global count variables
    initCountVars( );
    
    buildTSSModel(annoseq);
    buildTTSModel(annoseq);
    buildProbabilities(annoseq);
    

    if (!hasLenDist && verbosity > 0) {
	cout << "5' single:  " << num5Single << endl;
	cout << "5' initial: " << num5Initial << endl;
	cout << "5' internal:" << num5Internal << endl;
	cout << "5' terminal:" << num5Terminal << endl;
	cout << "5' intron:  " << num5Introns << endl;
	
	cout << "3' single:  " << num3Single << endl;
	cout << "3' initial: " << num3Initial << endl;
	cout << "3' internal:" << num3Internal << endl;
	cout << "3' terminal:" << num3Terminal << endl;
	cout << "3' intron:  " << num3Introns << endl;
    }
    if (verbosity>2){
	cout << "number of sequences used for training of:\ntssMotif: " << tssMotif->numSeqs << endl 
	     << "tssMotifTATA: " << tssMotifTATA->numSeqs << endl 
	     << "tataMotif: " << tataMotif->numSeqs << endl;
	cout << "distance from start of first tata motif relative to tss:" << endl;
	cout << "pos\tfreq" << endl;
	for (int i=0; i <= d_tss_tata_max - d_tss_tata_min; i++)
	    cout << -d_tss_tata_max + i << "\t" << distCountTata[i] << endl;
    }
    tssMotif->makeProbs();
    tssMotifTATA->makeProbs();
    ttsMotif->makeProbs();
    tataMotif->makeProbs();
    if (!hasLenDist) 
	computeLengthDistributions();
    hasLenDist = true;
    lastParIndex = parIndex;
}

/*
 * UtrModel::buildProbabilities
 */
void UtrModel::buildProbabilities(const AnnoSequence* annoseq){
    const Gene *curgene;
    const AnnoSequence *as = annoseq;
    while (as){
	sequence = as->sequence;
	curgene = dynamic_cast<const Gene*> (as->anno->genes); // assume here that sequences with multiple genes have already been split
	if (curgene){
	    gweight = curgene->weight;
	    if (curgene->utr5exons || curgene->utr3exons){
		try {
		    processStates(curgene);
		} catch (UtrModelError e) {
		    cerr << "Sequence " << curgene->seqname << ":" << endl;
		    cerr << e.getMessage() << endl;
		}
	    }
	
	    /*
	      if (curgene->utr3exons){
	  
	      if (curgene->utr3exons->next == NULL)
	      cout << "3utr no intron. length " << curgene->utr3exons->length() << endl;
	      else {
	      cout << "3utr intron. distance " << curgene->utr3exons->end - curgene->lastExon()->end << endl;
	      if (curgene->utr3exons->end - curgene->lastExon()->end > 50)
	      cout << "gene disobeys 50 nt rule " << curgene->geneid << endl;
	      }
	      }*/
	}
	as = as->next;
    }
    
    /*
     * determine emission probabilities
   */
    if (verbosity>2)
	cout << "number of bases used for training tss upwindow markov chain: " << tssup_gesbasen << endl;
    // 5' utr emission 

    vector<Double> patternprobs(utr5_emicount.size());
    makeProbsFromCounts(patternprobs, utr5_emicount, k, utr_patpseudo, false);
    utr5_emiprobs.probs.resize(utr5_emicount.size());
    utr5_emiprobs.order = k;
    computeEmiFromPat(patternprobs, utr5_emiprobs.probs, k);    

    // 5' utr initial emission (inital exon and single exon)
    patternprobs.resize(utr5init_emicount.size());
    makeProbsFromCounts(patternprobs, utr5init_emicount, k, utr_patpseudo, false);
    utr5init_emiprobs.probs.resize(utr5init_emicount.size());
    utr5init_emiprobs.order = k;
    computeEmiFromPat(patternprobs, utr5init_emiprobs.probs, k);    

    // 3' utr emission
    patternprobs.resize(utr3_emicount.size());
    makeProbsFromCounts(patternprobs, utr3_emicount, k, utr_patpseudo, false);
    utr3_emiprobs.probs.resize(utr3_emicount.size());
    utr3_emiprobs.order = k;
    computeEmiFromPat(patternprobs, utr3_emiprobs.probs, k);    
}

/*
 * ===[ UtrModel::registerPars  ]===================================
 */
void UtrModel::registerPars (Parameters* parameters){
  try {
    Constant::CRFtrainUTR = Properties::getBoolProperty("CRFtrainUTR");
  } catch (...) {} 
  if (!parameters || !Constant::CRFtrainUTR)
    return;

  for (int idx=0; idx < Constant::decomp_num_steps; idx++){
    parameters->addMMGroup(&GCutr5_emiprobs[idx]);
    parameters->addMMGroup(&GCutr5init_emiprobs[idx]);
    parameters->addMMGroup(&GCutr3_emiprobs[idx]);
  }
}


/*
 * UtrModel::buildTSSModel
 */
void UtrModel::buildTSSModel(const AnnoSequence* annoseq){
    const Gene *curgene;
    int dnalen;
    while (annoseq){
	curgene = dynamic_cast<Gene*> (annoseq->anno->genes); // assume here that sequences with multiple genes have already been split
	if (curgene && curgene->utr5exons){
	    try {
		dnalen = strlen(annoseq->sequence);
		int tsspos = curgene->utr5exons->begin;
		if ((tsspos - Constant::tss_upwindow_size) > 0 && (tsspos + tss_end - 1 < dnalen)){
		    //cout << "good TSS: " << curgene->seqname << " "; 
		    processTSS(annoseq->sequence + tsspos - Constant::tss_upwindow_size);
		}
	    } catch (UtrModelError e) {
		cerr << "Sequence " << curgene->seqname << ":" << endl;
		cerr << e.getMessage() << endl;
	    }
	}
	annoseq = annoseq->next;
    }
    /*
     * determine emission probabilities of tss_upwindow
     */

    vector<Double> patternprobs(tssup_emicount.size());
    makeProbsFromCounts(patternprobs, tssup_emicount, tssup_k, tssup_patpseudo, false);
    tssup_emiprobs.resize(tssup_emicount.size());
    computeEmiFromPat(patternprobs, tssup_emiprobs, tssup_k);
}

/*
 * UtrModel::buildTTSModel
 * 
 * In training the tts is fixed and the polyA signal found.
 * In prediction the polyA signal (aataaa motif) is fixed and the tts at a fixed relative position.
 */
void UtrModel::buildTTSModel(const AnnoSequence* annoseq){
  Gene *curgene;
  double threshold = 0.01; // pattern must be at least this frequent among polyasig_consensus or 1-base variants to be used
  aataaa_probs.assign(POWER4TOTHE(aataaa_boxlen), 0.0);
  aataaa_count.assign(POWER4TOTHE(aataaa_boxlen), 0);
  Seq2Int s2i(aataaa_boxlen);
  int bestpn = -1;
  int count_found=0, count_notfound=0;
  int *numposfound = new int[d_polya_cleavage_max-d_polya_cleavage_min+1]; // number of patterns found at this position
  for (int i= 0; i<d_polya_cleavage_max-d_polya_cleavage_min+1; i++)
    numposfound[i]=0;
  // this is actually trained later
  // for now, this table defines a priority list of all 1-base variants of aataaa or tgtaac
  // This is ONLY for the case in which more than one of those patterns occur in the window.
  if (polyasig_consensus == "aataaa"){
      aataaa_probs[s2i("aataaa")] = 0.651;
      aataaa_probs[s2i("attaaa")] = 0.170;
      aataaa_probs[s2i("tataaa")] = 0.035;
      aataaa_probs[s2i("agtaaa")] = 0.030;
      aataaa_probs[s2i("aatata")] = 0.021;
      aataaa_probs[s2i("cataaa")] = 0.016;
      aataaa_probs[s2i("gataaa")] = 0.016;
      aataaa_probs[s2i("aataca")] = 0.014;
      aataaa_probs[s2i("aagaaa")] = 0.013;
      aataaa_probs[s2i("aatgaa")] = 0.010;
      aataaa_probs[s2i("aataga")] = 0.009;
      aataaa_probs[s2i("actaaa")] = 0.008;
      aataaa_probs[s2i("aaaaaa")] = 0.001;
      aataaa_probs[s2i("aacaaa")] = 0.001;
      aataaa_probs[s2i("aatcaa")] = 0.001;
      aataaa_probs[s2i("aattaa")] = 0.001;
      aataaa_probs[s2i("aataac")] = 0.001;
      aataaa_probs[s2i("aataag")] = 0.001;
      aataaa_probs[s2i("aataat")] = 0.001;
  } else if (polyasig_consensus == "tgtaac") {
      aataaa_probs[s2i("tgtaac")] = 0.387;
      aataaa_probs[s2i("tgtaaa")] = 0.212;
      aataaa_probs[s2i("tgtaag")] = 0.206;
      aataaa_probs[s2i("tgtaat")] = 0.195;
  } else {
      aataaa_probs[s2i(polyasig_consensus.c_str())] = 0.01; // this needs to be extended appropriately
  }
  while (annoseq){
      curgene = dynamic_cast<Gene*> (annoseq->anno->genes); // assume here that sequences with multiple genes have already been split
      if (curgene && curgene->complete3utr && curgene->transend >= 0){
	  // search for best matching pattern
	  bool found = false;
	  int bestpos = 999;
	  // string upregion(annoseq->sequence+curgene->transend - d_polya_cleavage_min - aataaa_boxlen + 1, d_polya_cleavage_max - d_polya_cleavage_min + 1);
	  // cout << "tts upregion\t" << upregion << endl;
	  for (int pos = curgene->transend - d_polya_cleavage_min - aataaa_boxlen + 1; pos >= curgene->transend - d_polya_cleavage_max - aataaa_boxlen + 1; pos--) {
	      if (pos >= 0 && pos < annoseq->length - aataaa_boxlen) {
		  try {
		      int pn = s2i(annoseq->sequence+pos);
		      if ((found && (aataaa_probs[pn] > aataaa_probs[bestpn])) || 
			  (!found && aataaa_probs[pn] > 0.0)) {
			  bestpn = pn;
			  bestpos = pos - (curgene->transend - d_polya_cleavage_max - aataaa_boxlen + 1); // base after the polyA signal
			  found = true;
		      }
		  } catch (...) {}
	      }
	  }
	  if (found) {
	      numposfound[bestpos]++;
	      aataaa_count[bestpn]++;
	      ttsMotif->addSequence(annoseq->sequence+curgene->transend -d_polya_cleavage_max + 1 + bestpos, curgene->weight);
	      count_found++;
	  } else {
	      count_notfound++;
	  }
      }
      annoseq = annoseq->next;
  }
  
  int numpatterns=0;
  for (int i=0; i<POWER4TOTHE(aataaa_boxlen); i++)
    if (aataaa_count[i]>0)
      numpatterns++;
  int *goodpatterns = new int[numpatterns];
  int p=0;
  for (int i=0; i<POWER4TOTHE(aataaa_boxlen); i++)
    if (aataaa_count[i]>0)
      goodpatterns[p++]=i;
  // sort patterns by frequency
  for (int i=0; i<numpatterns-1; i++)
    for(int j=i+1;j<numpatterns; j++)
      if (aataaa_count[goodpatterns[i]] < aataaa_count[goodpatterns[j]]){
	// swap i and j
	int k = goodpatterns[j];
	goodpatterns[j] = goodpatterns[i];
	goodpatterns[i] = k;
      }

  if (verbosity > 1) {
    cout << "polyadenylation signal:" << endl;
    cout << "In " << count_found << " cases " << polyasig_consensus << " or a variant was found\nending in the window [" << -d_polya_cleavage_max << "," << -d_polya_cleavage_min <<
      "] relative to the tts.\nIn " << count_notfound << " cases no such pattern was found." << endl;
    if (verbosity == 1)
      cout << "Increase /UtrModel/verbosity to see pattern frequencies and their position." << endl;
  }
  if (verbosity > 2) {
    cout << "pattern frequencies" << endl;
    bool notused = false;
    
    cout << "used:" << endl;
    for (int i=0; i<numpatterns; i++) {
      if (aataaa_count[goodpatterns[i]] < threshold * count_found && !notused){
	cout << "not used (because appears less than "<< 100*threshold <<"%):" << endl;
	notused = true;
      }
      cout << s2i.inv(goodpatterns[i]) << " " << aataaa_count[goodpatterns[i]] << endl;
    }
    cout << "position frequency (w.r.t. base right after polyA signal: tts-/UtrModel/d_polyasig_cleavage+1" << endl;
    for (int i=0; i< d_polya_cleavage_max-d_polya_cleavage_min+1; i++){
      cout << i-d_polya_cleavage_max << "   \t" << numposfound[i] << endl;
    }
  }
  int new_count_found=0;
  for (int i=0; i<numpatterns; i++) {
    if (aataaa_count[goodpatterns[i]] < threshold * count_found)
      aataaa_count[goodpatterns[i]] = 0; // too rare is just as good as not at all
    new_count_found += aataaa_count[goodpatterns[i]];
  }
  for (int i=0; i<aataaa_probs.size(); i++) 
    aataaa_probs[i] = Double(aataaa_count[i])/new_count_found;
  delete [] numposfound;
  delete [] goodpatterns;
}

/*
 * UtrModel::processTSS
 * IR---|--------|**tata*a****|---------|-----|----|--5'UTR---...
 *      ^<----------- tss_upwindow_size ----->^
 *      |start                                |TSS 
 *
 */
void UtrModel::processTSS(const char* start){
  int winlen = Constant::tss_upwindow_size + tss_end;
  char *tss = newstrcpy(start, winlen);
  bool tataprom;
  int tatapos, tatawinbegin = Constant::tss_upwindow_size - d_tss_tata_max;
  tatapos = findTATA(tss + tatawinbegin, d_tss_tata_max - d_tss_tata_min);
  tataprom = (tatapos>-1);
  //cout << "tss " << tss;
 
  if (tataprom) {
      if (tatapos < 0 || tatapos > d_tss_tata_max - d_tss_tata_min)
	  throw  UtrModelError("processTSS: tata box out of range");
      distCountTata[tatapos]++;
  }
  if (tataprom){
    tssMotifTATA->addSequence(start + Constant::tss_upwindow_size - tss_start, 1);
    tataMotif->addSequence(start + tatawinbegin + tatapos - tata_start, 1);
    processTssupSequence(start, start + tatawinbegin + tatapos - tata_start - 1);
    processTssupSequence(start + tatawinbegin + tatapos + tata_end, start + Constant::tss_upwindow_size - tss_start - 1);
  } else {
    tssMotif->addSequence(start + Constant::tss_upwindow_size - tss_start, 1);
    processTssupSequence(start, start + Constant::tss_upwindow_size - tss_start -1);
  }
  //cout << endl;
}

/*
 * UtrModel::processStates
 */
void UtrModel::processStates( const Gene* gene ){
    State *exon, *lastexon;
    /* 
     *    5' UTR
     */
    if (gene->utr5exons) {
	exon = gene->utr5exons;
	if( !exon->next ){
	    try{
		process5SingleExon( exon, (gene->exons != NULL));
	    } catch (UtrModelError e) {
		if (verbosity)
		    cerr << "gene " << gene->geneid  << " transcr. " << gene->id << " in sequence " << gene->seqname << ": " << e.getMessage() << endl;
	    }
	} else { // 5' initial, internal, terminal exons and introns
	    try {
		process5InitialExon( exon );
	    } catch( UtrModelError e ){
		if (verbosity)
		    cerr << "gene " << gene->geneid << " transcr. " << gene->id << " in sequence " << gene->seqname << ": " << e.getMessage() << endl;
	    }
	    lastexon = exon;
	    exon = exon->next;
	    while( exon->next ){
		try {
		    process5InternalExon( exon );
		    process5Intron(lastexon->end + 1, exon->begin - 1);
		} catch( UtrModelError e ){
		    if (verbosity)
			cerr << "gene " << gene->geneid << " transcr. " << gene->id << " in sequence " << gene->seqname << ": " << e.getMessage() << endl;
		}
		lastexon = exon;
		exon = exon->next;
	    } 
	    try {
		process5TerminalExon( exon );
		process5Intron(lastexon->end + 1, exon->begin - 1);
	    } catch( UtrModelError e ){
		if (verbosity)
		    cerr << "gene " << gene->geneid << " transcr. " << gene->id << " in sequence " << gene->seqname << ": " << e.getMessage() << endl;
	    }
	}
    }
    
    /* 
     *    3' UTR
     */
    if (gene->utr3exons) {
	exon = gene->utr3exons;
	if( !exon->next ){
	    try{
		process3SingleExon( exon, (gene->exons != NULL) && gene->complete3utr);
	    } catch (UtrModelError e) {
		if (verbosity)
		    cerr << "gene " << gene->geneid  << " transcr. " << gene->id << " in sequence " << gene->seqname << ": " << e.getMessage() << endl;
	    }
	} else { // 3' initial, internal, terminal exons and introns
	    try {
		process3InitialExon( exon );
	    } catch( UtrModelError e ){
		if (verbosity)
		    cerr << "gene " << gene->geneid << " transcr. " << gene->id << " in sequence " << gene->seqname << ": " << e.getMessage() << endl;
	    }
	    lastexon = exon;
	    exon = exon->next;
	    while( exon->next ){
		try {
		    process3InternalExon( exon );
		    process3Intron(lastexon->end + 1, exon->begin - 1);
		} catch( UtrModelError e ){
		    if (verbosity)
			cerr << "gene " << gene->geneid << " transcr. " << gene->id << " in sequence " << gene->seqname << ": " << e.getMessage() << endl;
		}
		lastexon = exon;
		exon = exon->next;
	    } 
	    try {
		process3TerminalExon( exon , gene->complete3utr);
		process3Intron(lastexon->end + 1, exon->begin - 1);
	    } catch( UtrModelError e ){
		if (verbosity)
		    cerr << "gene " << gene->geneid << " transcr. " << gene->id << " in sequence " << gene->seqname << ": " << e.getMessage() << endl;
	    }
	}
    }
}

/*
 * UtrModel::process5SingleExon
 */
void UtrModel::process5SingleExon( const State* exon, bool withLen) {
  int length = exon->end-exon->begin+1;
  if (length < 1)
    throw UtrModelError("Single 5'UTR training exon too short.");
  process5InitSequence(sequence + exon->begin + tss_end + k, sequence + exon->end - Constant::trans_init_window);
  if (!hasLenDist && withLen) {
    if (length <= exonLenD) {
      lenCount5Single[length] += 1;
    } else {
      numHuge5Single += 1;
    }
    num5Single += 1;
  }
}

/*
 * UtrModel::process5InitialExon
 */
void UtrModel::process5InitialExon( const State* exon, bool withLen) {
  int length = exon->end-exon->begin+1;
  if (length < 1)
    throw UtrModelError("Initial 5'UTR training exon too short.");
  process5InitSequence(sequence + exon->begin + tss_end + k, sequence + exon->end - Constant::dss_start);
  if (!hasLenDist && withLen) {
    if (length <= exonLenD) {
      lenCount5Initial[length] += 1;
    } else {
      numHuge5Initial += 1;
    }
    num5Initial += 1;
  }
}

/*
 * UtrModel::process5InternalExon
 */
void UtrModel::process5InternalExon( const State* exon) {
  int length = exon->end-exon->begin+1;
  if (length < 1)
    throw UtrModelError("Internal 5'UTR training exon too short.");
  process5Sequence(sequence + exon->begin + Constant::ass_end + k, sequence + exon->end - Constant::dss_start);
  
  if (!hasLenDist) {
    if (length <= exonLenD) {
      lenCount5Internal[length] += 1;
    } else {
      numHuge5Internal += 1;
    }
    num5Internal += 1;
  }
}

/*
 * UtrModel::process5TerminalExon
 */
void UtrModel::process5TerminalExon( const State* exon) {
   int length = exon->end-exon->begin+1;
  if (length < 1)
      throw UtrModelError("Terminal 5'UTR training exon too short.");
  process5Sequence(sequence + exon->begin + Constant::ass_end + k, sequence + exon->end - Constant::trans_init_window);
  
  if (!hasLenDist) {
    if (length <= exonLenD) {
      lenCount5Terminal[length] += 1;
    } else {
      numHuge5Terminal += 1;
    }
    num5Terminal += 1;
  }
}

/*
 * UtrModel::processIntron
 */
void UtrModel::process5Intron( int begin, int end){
  
  /* Testausgabe
// utr splice site profiles are very much like the others
  static char *dssfenster = new char[18];
  strncpy( dssfenster, sequence + begin-6, 17 );
  dssfenster[18]='\0';
  cout << "utr dss-window "  << dssfenster << endl;
  static char *assfenster = new char[36];
  strncpy( assfenster, sequence + end-30, 35 );
  assfenster[35]='\0';
  cout << "utr ass-window "  << assfenster << endl;
  */
    if (!hasLenDist)
	num5Introns++;
}


/*
 * UtrModel::process3SingleExon
 */
void UtrModel::process3SingleExon( const State* exon, bool withLen) {
  int length = exon->end-exon->begin+1;
  if (length < 1)
    throw UtrModelError("Single 3'UTR training exon too short.");
  process3Sequence(sequence + exon->begin + k, sequence + exon->end);
  if (!hasLenDist && withLen) {
    if (length <= exonLenD && length > 2) {
	lenCount3Single[length] += 1;
    } else if (length > exonLenD) {
	numHuge3Single += 1;
    }
    num3Single += 1;
  }
}

/*
 * UtrModel::process3InitialExon
 */
void UtrModel::process3InitialExon( const State* exon, bool withLen) {
  int length = exon->end-exon->begin+1;
  if (length < 1)
    throw UtrModelError("Initial 3'UTR training exon too short.");
  process3Sequence(sequence + exon->begin + k, sequence + exon->end - Constant::dss_start);
  if (!hasLenDist && withLen) {
    if (length <= exonLenD) {
      lenCount3Initial[length] += 1;
    } else {
      numHuge3Initial += 1;
    }
    num3Initial += 1;
  }
}

/*
 * UtrModel::process3InternalExon
 */
void UtrModel::process3InternalExon( const State* exon) {
  int length = exon->end-exon->begin+1;
  if (length < 1)
    throw UtrModelError("Internal 3'UTR training exon too short.");
  process3Sequence(sequence + exon->begin + Constant::ass_end + k, sequence + exon->end - Constant::dss_start);
  
  if (!hasLenDist) {
    if (length <= exonLenD) {
      lenCount3Internal[length] += 1;
    } else {
      numHuge3Internal += 1;
    }
    num3Internal += 1;
  }
}

/*
 * UtrModel::process3TerminalExon
 */
void UtrModel::process3TerminalExon( const State* exon, bool withLen) {
   int length = exon->end-exon->begin+1;
  if (length < 1)
      throw UtrModelError("Terminal 3'UTR training exon too short.");
  process3Sequence(sequence + exon->begin + Constant::ass_end + k, sequence + exon->end);
  
  if (!hasLenDist && withLen) {
    if (length <= exonLenD) {
      lenCount3Terminal[length] += 1;
    } else {
      numHuge3Terminal += 1;
    }
    num3Terminal += 1;
  }
}

/*
 * UtrModel::processIntron
 */
void UtrModel::process3Intron( int begin, int end){
    if (!hasLenDist)
	num3Introns++;
}


/*
 * UtrModel::computeLengthDistributions
 */
void UtrModel::computeLengthDistributions( ){

    Smooth sm(minwindowcount, slope_of_bandwidth);
  
    /*
     * compute length distributions
     *
     */
    lenDist5Single.resize(max_exon_length+1);
    lenDist5Initial.resize(max_exon_length+1);
    lenDist5Internal.resize(max_exon_length+1);
    lenDist5Terminal.resize(max_exon_length+1);
 
    lenDist3Single.resize(max3singlelength+1);
    lenDist3Initial.resize(max_exon_length+1);
    lenDist3Internal.resize(max_exon_length+1);
    lenDist3Terminal.resize(max3termlength+1);
  
    // compute mean exon lengths
    if (verbosity) {
	long int sum;
	int num, i;

	cout << "5' single, initial, internal, terminal mean exon lengths :" << endl;
	sum = 0; num = 0;
	for (i=0; i<= exonLenD; i++){
	    sum += i*lenCount5Single[i]; num += lenCount5Single[i];
	}
	if (num>0)
	    cout << sum/num << "\t";
	else 
	    cout << "n.a." << "\t";
	sum = 0; num = 0;
	for (i=0; i<= exonLenD; i++){
	    sum += i*lenCount5Initial[i]; num += lenCount5Initial[i];
	}
	if (num>0)
	    cout << sum/num << "\t";
	else 
	    cout << "n.a." << "\t";
	sum = 0; num = 0;
	for (i=0; i<= exonLenD; i++){
	    sum += i*lenCount5Internal[i]; num += lenCount5Internal[i];
	}
	if (num>0)
	    cout << sum/num << "\t";
	else 
	    cout << "n.a." << "\t";
	sum = 0; num = 0;
	for (i=0; i<= exonLenD; i++){
	    sum += i*lenCount5Terminal[i]; num += lenCount5Terminal[i];
	}
	if (num>0)
	    cout << sum/num << endl;
	else 
	    cout << "n.a." << endl;

	cout << "3' single, initial, internal, terminal mean exon lengths :" << endl;
	sum = 0; num = 0;
	for (i=0; i<= exonLenD; i++){
	    sum += i*lenCount3Single[i]; num += lenCount3Single[i];
	}
	if (num>0)
	    cout << sum/num << "\t";
	else 
	    cout << "n.a." << "\t";
	sum = 0; num = 0;
	for (i=0; i<= exonLenD; i++){
	    sum += i*lenCount3Initial[i]; num += lenCount3Initial[i];
	}
	if (num>0)
	    cout << sum/num << "\t";
	else 
	    cout << "n.a." << "\t";
	sum = 0; num = 0;
	for (i=0; i<= exonLenD; i++){
	    sum += i*lenCount3Internal[i]; num += lenCount3Internal[i];
	}
	if (num>0)
	    cout << sum/num << "\t";
	else 
	    cout << "n.a." << "\t";
	sum = 0; num = 0;
	for (i=0; i<= exonLenD; i++) {
	    sum += i*lenCount3Terminal[i]; num += lenCount3Terminal[i];
	}
	if (num>0)
	    cout << sum/num << endl;
	else 
	    cout << "n.a." << endl;
    }

    if (numHuge5Single>0 && verbosity)
	cout << numHuge5Single << " 5' single exons were " << " longer than " << exonLenD << endl;
    if (numHuge5Initial>0 && verbosity)
	cout << numHuge5Initial << " 5' initial exons were " << " longer than " << exonLenD << endl;
    if (numHuge5Internal>0 && verbosity)
	cout << numHuge5Internal << " 5' internal exons were " << " longer than " << exonLenD << endl;
    if (numHuge5Terminal>0 && verbosity)
	cout << numHuge5Terminal << " 5' terminal exons were " << " longer than " << exonLenD << endl;

    if (numHuge3Single>0 && verbosity)
	cout << numHuge3Single << " 3' single exons were " << " longer than " << exonLenD << endl;
    if (numHuge3Initial>0 && verbosity)
	cout << numHuge3Initial << " 3' initial exons were " << " longer than " << exonLenD << endl;
    if (numHuge3Internal>0 && verbosity)
	cout << numHuge3Internal << " 3' internal exons were " << " longer than " << exonLenD << endl;
    if (numHuge3Terminal>0 && verbosity)
	cout << numHuge3Terminal << " 3' terminal exons were " << " longer than " << exonLenD << endl;
   

    sm.smoothCounts(lenCount5Single, lenDist5Single);
    sm.smoothCounts(lenCount5Initial, lenDist5Initial);
    sm.smoothCounts(lenCount5Internal, lenDist5Internal);
    sm.smoothCounts(lenCount5Terminal, lenDist5Terminal);
    scaleDblVector(lenDist5Single, Double(num5Single-numHuge5Single) / num5Single);
    scaleDblVector(lenDist5Initial, Double(num5Initial-numHuge5Initial) / num5Initial);
    scaleDblVector(lenDist5Internal, Double(num5Internal-numHuge5Internal) / num5Internal);
    scaleDblVector(lenDist5Terminal, Double(num5Terminal-numHuge5Terminal) / num5Terminal);

    sm.smoothCounts(lenCount3Single, lenDist3Single);
    sm.smoothCounts(lenCount3Initial, lenDist3Initial);
    sm.smoothCounts(lenCount3Internal, lenDist3Internal);
    sm.smoothCounts(lenCount3Terminal, lenDist3Terminal);
    scaleDblVector(lenDist3Single, Double(num3Single-numHuge3Single) / num3Single);
    scaleDblVector(lenDist3Initial, Double(num3Initial-numHuge3Initial) / num3Initial);
    scaleDblVector(lenDist3Internal, Double(num3Internal-numHuge3Internal) / num3Internal);
    scaleDblVector(lenDist3Terminal, Double(num3Terminal-numHuge3Terminal) / num3Terminal);

    fillTailsOfLengthDistributions();
}


/*
 * UtrModel::processSequence
 */
void UtrModel::process5InitSequence( const char* start, const char* end){
  Seq2Int s2i(k+1);
  for( ; start <= end; start++){
      try {
	  utr5init_emicount[s2i(start-k)] += gweight;
	  utr5init_gesbasen += gweight;
      } catch (InvalidNucleotideError e) {
      }
  }
}

/*
 * UtrModel::process5Sequence
 */
void UtrModel::process5Sequence( const char* start, const char* end){
  Seq2Int s2i(k+1);
  for( ; start <= end; start++){
      try {
	  utr5_emicount[s2i(start-k)] += gweight ;
	  utr5_gesbasen += gweight;
      } catch (InvalidNucleotideError e) {
      }  
  }
}

/*
 * UtrModel::process3Sequence
 */
void UtrModel::process3Sequence( const char* start, const char* end){
  Seq2Int s2i(k+1);
  for( ; start <= end; start++){
      try {
	  utr3_emicount[s2i(start-k)] += gweight ;
	  utr3_gesbasen += gweight;
      } catch (InvalidNucleotideError e) {
      }  
  }
}

/*
 * UtrModel::processTssupSequence
 */
void UtrModel::processTssupSequence( const char* start, const char* end){
  Seq2Int s2i(tssup_k+1);
  for( ; start <= end; start++){
    try {
      tssup_emicount[s2i(start-tssup_k)] += 1 ;
      tssup_gesbasen += 1;
    } catch (InvalidNucleotideError e) {
    }  
  }
}

/*
 * storeGCPars
 * store GC dependent model probabilities in the respective array variables
 */
void UtrModel::storeGCPars(int idx){
  try {
    if (!Constant::utr_option_on)
      return;
  } catch(ProjectError e) {
    return;
  }
  if (!IntronModel::GCemiprobs || utr5_emiprobs.probs.empty() || utr3_emiprobs.probs.empty() || utr5init_emiprobs.probs.empty() || !tssMotif || !ttsMotif || !tssMotifTATA || !tataMotif)
    return;
  
  /* Set emission probabilities to an actual mixture with the intron model emission probabilities
   * This is done at two different times:
   * The old code (utr{5,3}patternweight) does this after reading it from the *_utr_probs.pbl file.
   * This newer code does this right after training, so it can be used as a more neutral start for CRF training.
   */
  for (int i=0; i<POWER4TOTHE(k+1); i++) {
    utr5init_emiprobs.probs[i] = utr5init_emiprobs.probs[i] * utr5prepatternweight + IntronModel::GCemiprobs[idx].probs[i] * (1.0-utr5prepatternweight);
    utr5_emiprobs.probs[i]     = utr5_emiprobs.probs[i]     * utr5prepatternweight + IntronModel::GCemiprobs[idx].probs[i] * (1.0-utr5prepatternweight);
    utr3_emiprobs.probs[i]     = utr3_emiprobs.probs[i]     * utr3prepatternweight + IntronModel::GCemiprobs[idx].probs[i] * (1.0-utr3prepatternweight);
  }

  GCutr5init_emiprobs[idx] = utr5init_emiprobs;
  GCutr5_emiprobs[idx] = utr5_emiprobs;
  GCutr3_emiprobs[idx] = utr3_emiprobs;
  std::stringstream out1, out2, out3;
  out1 << "utr5 emiprob gc" << (idx+1);
  GCutr5_emiprobs[idx].setName(out1.str());
  out2 << "utr5init emiprob gc" << (idx+1);
  GCutr5init_emiprobs[idx].setName(out2.str());
  out3 << "utr3 emiprob gc" << (idx+1);
  GCutr3_emiprobs[idx].setName(out3.str());
  GCtssup_emiprobs[idx] = tssup_emiprobs;
  GCtssMotif[idx] = *tssMotif;
  GCttsMotif[idx] = *ttsMotif;
  GCtssMotifTATA[idx] = *tssMotifTATA;
  GCtataMotif[idx] = *tataMotif;
}

/*
 * UtrModel::printProbabilities
 */
void UtrModel::printProbabilities( int idx, BaseCount *bc, const char* suffix ){
    try{
      const char* fname = Properties::getProperty("/UtrModel/outfile"); 
      string filename = Constant::fullSpeciesPath() + fname;
      if (suffix)
	  filename += suffix;
      if (idx < 0) { // do nothing but deleting the outfile
	ofstream ofstrm(filename.c_str(), ios::trunc);
	ofstrm.close();
	return;
      }
      ofstream ofstrm(filename.c_str(), ios::app);
  
      if( ofstrm ){
	  if (verbosity)
	  cout << "Writing UTR model parameters [" << idx+1 << "] to file " << filename << "." << endl;
	  if (idx == 0 ) {
	      ofstrm << "# UTR model parameters\n# begin of content independent part" << endl;
	      ofstrm << "\n# Length distributions\n[UTRLENGTH]" << endl;
	      ofstrm << "# maximal individually stored length probability d=\n" << exonLenD << endl;
	      ofstrm << "# slope of smoothing bandwidth =\n" << slope_of_bandwidth << endl;
	      ofstrm << "# smoothing minwindowcount =\n" << minwindowcount << endl;
	      ofstrm << "# length 5' sing  5' init  5' int  5' term  3' sing  3' init  3' int  3' term " << endl;
	      ofstrm << "# total number of exons of above types" << endl;
	      ofstrm << setw(10) << num5Single << setw(10) << num5Initial << setw(10)  << num5Internal <<setw(10) << num5Terminal << setw(10)
		     << num3Single << setw(10) << num3Initial << setw(10)  << num3Internal <<setw(10) << num3Terminal  << endl;
	      ofstrm << "# number of exons exceeding length d=" << exonLenD << endl;
	      ofstrm << setw(10) << numHuge5Single << setw(10) << numHuge5Initial << setw(10) << numHuge5Internal << setw(10) << numHuge5Terminal << setw(10)
		     << numHuge3Single << setw(10) << numHuge3Initial << setw(10) << numHuge3Internal << setw(10) << numHuge3Terminal << endl;
	      
	      ofstrm << "# 1000 P(len=k), k=0,1,..., " << exonLenD << endl;
	      for( int i = 0; i <= exonLenD; i++ ){
		  ofstrm << i 
		   << "\t" << 1000 * lenDist5Single[i] << "\t" << 1000 * lenDist5Initial[i] 
			 << "\t" << 1000 * lenDist5Internal[i] << "\t" << 1000 * lenDist5Terminal[i] 
			 << "\t" << 1000 * lenDist3Single[i] << "\t" << 1000 * lenDist3Initial[i] 
			 << "\t" << 1000 * lenDist3Internal[i] << "\t" << 1000 * lenDist3Terminal[i] 
			 << endl;
	      }

	      ofstrm << "#\n# AATAAA/TGTAA box probabilities\n"
		     << "[AATAAA]\n";
	      ofstrm << "# Size of vector\n" << aataaa_probs.size() << endl;
	      ofstrm << "# Probabilities\n";
	      Seq2Int s2i(aataaa_boxlen);
	      for( int i = 0; i < aataaa_probs.size(); i++ ){
		  if( aataaa_probs[i] > 0.0 )
		      ofstrm << s2i.inv(i) << "\t" << aataaa_probs[i] << endl;
	      }

	      ofstrm << "# end of content independent part" << endl;
	  }
		
	  ofstrm << "[" << idx+1 << "]" << endl;
	  if (bc)
	      ofstrm << "# (a,c,g,t)= " << *bc << endl;

	  ofstrm << "#\n# Probabilities file for the UTR model\n#\n#"
		 << endl;

	  // save emission probabilities: single and initial UTR exons
	  ofstrm << "\n#\n# The emission probabilities of single and initial 5' UTR exons\n#" << endl;
	  ofstrm << "[EMISSION-5INITIAL]" << endl;
	  ofstrm << "# size of the emission vector\n" << GCutr5init_emiprobs[idx].probs.size() << endl;
	  ofstrm << "#k=\n" << k << endl;
	  ofstrm << "# patpseudo : pseudocount for sequence patterns\n" << utr_patpseudo << endl;
	  Seq2Int s2i(k+1);
	  for( int i = 0; i != GCutr5init_emiprobs[idx].probs.size(); ++i )
	      ofstrm << s2i.inv(i) << '\t' << GCutr5init_emiprobs[idx].probs[i] << endl;

	  /************ frequencies fo the (k+1) tuples
	  ofstrm << "\n# 5' initial and single exon patterns:" << endl;
	  for( int j = 0; j < utr5init_emicount.size(); j++ )
	      ofstrm << "#\t" << s2i.inv(j) << "\t" << utr5init_emicount[j] << endl;
	  //  *****************************************/

	  // save emission probabilities: internal and terminal UTR exons
	  ofstrm << "\n#\n# The emission probabilities of internal and terminal 5' UTR exons\n#" << endl;
	  ofstrm << "[EMISSION-5]" << endl;
	  ofstrm << "# size of the emission vector\n" << GCutr5_emiprobs[idx].probs.size() << endl;
	  ofstrm << "#k=\n" << k << endl;
	  ofstrm << "# patpseudo : pseudocount for sequence patterns\n" << utr_patpseudo << endl;
	  for( int i = 0; i != GCutr5_emiprobs[idx].probs.size(); ++i ){
	      ofstrm << s2i.inv(i) << '\t' << GCutr5_emiprobs[idx].probs[i] << endl;
	  }

	  /************* frequencies fo the (k+1) tuples
	  ofstrm << "\n# 5' internal and terminal exon patterns:" << endl;
	  for( int j = 0; j < utr5_emicount.size(); j++ )
	      ofstrm << "#\t" << s2i.inv(j) << "\t" << utr5_emicount[j] << endl;
	  //  *****************************************/

	  // save emission probabilities: 3' UTR exons
	  ofstrm << "\n#\n# The emission probabilities of 3' UTR exons\n#" << endl;
	  ofstrm << "[EMISSION-3]" << endl;
	  ofstrm << "# size of the emission vector\n" << GCutr3_emiprobs[idx].probs.size() << endl;
	  ofstrm << "#k=\n" << k << endl;
	  ofstrm << "# patpseudo : pseudocount for sequence patterns\n" << utr_patpseudo << endl;
	  for( int i = 0; i != GCutr3_emiprobs[idx].probs.size(); ++i ){
	      ofstrm << s2i.inv(i) << '\t' << GCutr3_emiprobs[idx].probs[i] << endl;
	  }

	  /************* frequencies fo the (k+1) tuples
	  ofstrm << "\n# 3' exon patterns:" << endl;
	  for( int j = 0; j < utr3_emicount.size(); j++ )
	      ofstrm << "#\t" << s2i.inv(j) << "\t" << utr3_emicount[j] << endl;
	      //  *****************************************/

	  // save emission probabilities: tss upwindow
 
	  ofstrm << "\n#\n# The emission probabilities of the tss upwindow\n#" << endl;
	  ofstrm << "[EMISSION-TSSUPWIN]" << endl;
	  ofstrm << "# size of the emission vector\n" << GCtssup_emiprobs[idx].size() << endl;
	  ofstrm << "#tssup_k=\n" << tssup_k << endl;
	  ofstrm << "# patpseudo : pseudocount for sequence patterns\n" << tssup_patpseudo << endl;
	  Seq2Int s2iup(tssup_k+1);
	  for( int i = 0; i != GCtssup_emiprobs[idx].size(); ++i ){
	      ofstrm << s2iup.inv(i) << '\t' << GCtssup_emiprobs[idx][i] << endl;
	  }

	  ofstrm << "\n# motif around the TSS of TATA-less promoters\n[TSSMOTIF]" << endl;
	  GCtssMotif[idx].write(ofstrm);
	  ofstrm << "\n# motif around the TSS of TATA promoters\n[TSSMOTIFTATA]" << endl;
	  GCtssMotifTATA[idx].write(ofstrm);
	  ofstrm << "\n# tata box motif \n[TATAMOTIF]" << endl;
	  GCtataMotif[idx].write(ofstrm);
	  ofstrm << "\n# motif after polyA signal\n[TTSMOTIF]" << endl;
	  GCttsMotif[idx].write(ofstrm);
	  ofstrm.close();
      }
    } catch( ProjectError ) { 
	if (verbosity)
	    cout << "Utr model parameters not saved." << endl;
    }
}
