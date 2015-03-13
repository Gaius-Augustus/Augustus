/**********************************************************************
 * file:    introntrain.cc
 * licence: 
 *          
 * descr.:  
 * authors: Mario Stanke, mario@gobics.de
 *
 * date    |   author      |  changes 
 * --------|---------------|------------------------------------------
 * 02.01.06| Mario Stanke  | creating the file
 * 11.05.06| Mario Stanke  | fixed error in case no intron is longer than d
 **********************************************************************/

#include "intronmodel.hh"

// project includes
#include "gene.hh"
#include "properties.hh"
#include "merkmal.hh"
#include "projectio.hh"
#include "commontrain.hh"

// standard C/C++ includes
#include <fstream>
#include <climits>
#include <sstream>

/*
 * IntronModel::buildModel
 */
void IntronModel::buildModel( const AnnoSequence* annoseq, int parIndex){
    if (verbosity)
	cout << "** building model for introns. *INTRON*" << endl;
  
   if (!annoseq->anno)
	throw ProjectError("Tried to train Intronmodel with 0 sequences.");

    assMotif = new Motif(ass_upwindow_size, ass_motif_memory, 2, ass_motif_radius);            // second parameter: k=3
 
    /*
     * Initialize global count variables
     */ 
    initCountVars( );
    
    storeIntronLengths(annoseq);
    if (verbosity)
	printLengthQuantiles(); // just for info   

    buildLenDist(annoseq);
    buildProbabilities(annoseq);
    lastParIndex = parIndex;

    if (verbosity) {
	cout << "number of bases used for training: " << gesbasen << endl;
	cout << "end *INTRON*" << endl;
    }
}

/*
 * ===[ IntronModel::registerPars  ]===================================
 */
void IntronModel::registerPars (Parameters* parameters){
  int ssalpha = 10;
  try {
    ssalpha = Properties::getIntProperty("/IntronModel/ssalpha");
  } catch (...){}
  
  if (!parameters)
	return;
  if (Constant::CRFtrainSS){
    if (Constant::ass_maxbinsize > 0)
      parameters->addMMGroup(&assBinProbs, ssalpha);
    if (Constant::dss_maxbinsize > 0)
      parameters->addMMGroup(&dssBinProbs, ssalpha);
  }
  if (Constant::CRFtrainIntron)
    for (int idx=0; idx < Constant::decomp_num_steps; idx++)
      parameters->addMMGroup(&GCemiprobs[idx]);
}

void IntronModel::storeIntronLengths ( const AnnoSequence* annoseq){
    introns  = 0;

    /*
     * Determine the maximal length of an intron, allocate a vector
     * and store the intron lengths. 
     */
    int maxintronlength = 0, minintronlength = INT_MAX, len;
    State *tintron;
    for (AnnoSeqGeneIterator geneit = AnnoSeqGeneIterator(annoseq); geneit.hasMoreElements(); ++geneit) {
	tintron = geneit.gene->introns;
	while (tintron) {
	    len = tintron->length();
	    //cout << "ilen= " << len << endl;
	    if (len > maxintronlength) {
		maxintronlength = len;
	    }
	    if (len < minintronlength) {
		minintronlength = len;
	    }
	    tintron = tintron->next;
	}
    }
    if (verbosity)
	cout << "minimal / maximal length of a training-intron: " 
	     << minintronlength << " / " << maxintronlength << endl;
    intlencount.assign(maxintronlength + 1, 0);
    for (AnnoSeqGeneIterator geneit = AnnoSeqGeneIterator(annoseq); geneit.hasMoreElements(); ++geneit) {
	tintron = geneit.gene->introns;
	while (tintron) { 
	    intlencount[tintron->length()] += 1; //geneit.annoseq->weight;
	    introns += 1; //geneit.annoseq->weight;
	    tintron = tintron->next;
	}
    }
}


void IntronModel::initCountVars( ){
  /*
   * Initialize global count variables
   */
  if (!hasSpliceSites){ // in CRF training the counts are reused for multiple gc content classes
    asscount.assign(POWER4TOTHE(Constant::ass_size()), 0 );
    dsscount.assign(POWER4TOTHE(Constant::dss_size()), 0 );
    c_ass = 0;
    c_dss = 0;
  }
  
  emicount.assign( POWER4TOTHE( k+1 ), 0 );
  gesbasen = 0;  
}


/*
 * IntronModel::printLengthQuantiles
 */
void IntronModel::printLengthQuantiles(){
    /*
     * take the shortest q% and the longest q% of the introns as initial
     * set for training 'short' and 'long' introns, respectively.
     */
  double q = .333;
    int lensum, len, size = intlencount.size(), d_short, d_long;
    for (len=0, lensum=0; (lensum < q * introns) && (len <= size); len++){
	lensum += intlencount[len];
    }
    d_short = len;
    for (len=size-1, lensum=0; (lensum < q * introns) && (len >= 0); len--){
	lensum += intlencount[len];
    }
    d_long = len;
    if (verbosity)
	cout << 100*q << "% and " << 100*(1-q) << "% quantile of the intron len. dist.: " 
	     << d_short << ", " << d_long << endl;
}


/*
 * IntronModel::processSequence
 */
void IntronModel::processSequence( const char* start, const char* end){
    Seq2Int s2i(k+1);
    for( ; start <= end; start++){
	try {
	    emicount[s2i(start-k)] += gweight ;
	    emicount[s2i.rc(start-k) ] += gweight ;
	    gesbasen += 2*gweight;
	} catch (InvalidNucleotideError e) {
	}  
    }
}



void IntronModel::readSpliceSites(){
    string filename = Constant::fullSpeciesPath();
    bool withMotif = false; //default, best for human
    try {
	filename += Properties::getProperty("/IntronModel/splicefile");
    } catch (...){
	if (verbosity) 
	    cout << "no splice file read" << endl;
	return;
    }
    try {
	withMotif = Properties::getBoolProperty("/IntronModel/sf_with_motif");
    } catch (...){}
    ifstream istrm(filename.c_str());
    char buff[400];
    int dssoffset=39;
    int assoffset=42;
    char *ass, *dss;
    if (istrm) {
	istrm >> goto_line_after( "[DSS-OFFSET]" ); 
	istrm >> comment >> dssoffset;
	istrm >> goto_line_after( "[ASS-OFFSET]" ); 
	istrm >> comment >> assoffset;
	do {
	    istrm.getline( buff, 399 );
	    if (strncmp (buff, "dss ", 4) == 0) {
		dss = buff + 4;
		processDSS(dss, dssoffset);
	    } else if (strncmp (buff, "ass ", 4) == 0) {
		ass = buff + 4;
		processASS(ass, assoffset, withMotif);
	    }
	  } while( istrm ); 
	istrm.clear( ios::failbit );
    } else {
	cerr << "Intronmodel: Couldn't open file " << filename << endl;
    }
}



/*
 * IntronModel::processASS
 * pos is the position of the first nucleotide in the exon
 * xxxxxxxxxxxxxxxxxx pp ag | ppp
 *
 *   intron                 | exon  
 */
void IntronModel::processASS(const char* dna, int pos, Boolean withMotif){
    /*
     * examine the region upstream of the splice site
     *
     * at most minIntronSize - Constant::ass_start-ASS_MIDDLE-DSS_MIDDLE-Constant::dss_end =
     * minIntronSize - 10
     * the shortest intron in the drosophila training set was 44bp
     */
    
    if (withMotif && pos - ASS_MIDDLE - Constant::ass_start - ass_upwindow_size >= assMotif->k) 	
	assMotif->addSequence(dna + pos - ASS_MIDDLE - Constant::ass_start - ass_upwindow_size, gweight);
    
    if (hasSpliceSites)
	return;
    static char *astr = new char[Constant::ass_size()+1]; // waive the delete
    Seq2Int s2i(Constant::ass_size());
    if(!onASS(dna + pos - ASS_MIDDLE))
        throw IntronModelError("ASS error! Expected 'ag' but found '" + 
			       string(dna, pos - ASS_MIDDLE,2) + 
			       "' at position " + itoa(pos - ASS_MIDDLE) + ".");
    

    strncpy( astr, dna + pos - ASS_MIDDLE - Constant::ass_start, Constant::ass_start );
    strncpy( astr + Constant::ass_start, dna + pos, Constant::ass_end );
    astr[Constant::ass_size()] = '\0';
    try {
	asscount[s2i(astr)]++;
	c_ass += 1;
    } catch (InvalidNucleotideError e){
	//cerr << "Acceptor splice site " << astr << " contains an invalid character." << endl;
    }
   
    /*/ output for seqlogo of ass
    static char *assfenster = new char[51];
    strncpy( assfenster, dna + pos-1-44, 50 );
    assfenster[50]='\0';
    cout << "cds ass "  << assfenster << endl;*/
}

/*
 * IntronModel::processDSS
 * dss is the last position of the exon
 *           pppp | gt ppppp
 *           exon |    intron
 */
void IntronModel::processDSS( const char* dna, int pos){ 
    if (hasSpliceSites || pos < Constant::dss_start + 1){
	return;
    }
    
    static char *dstr = new char[Constant::dss_size()+1]; // forget the delete
    Seq2Int s2i(Constant::dss_size());
    if (!onGenDSS(dna+pos+1))
	throw IntronModelError("DSS error! Expected 'gt'" +
			       string(Constant::dss_gc_allowed ?  " or 'gc'" : "") +
			       ", but found '" + string(dna, pos+1, 2) + "' at position " +
			       itoa(pos+1) + ".");
    else {
	strncpy( dstr, dna + pos - Constant::dss_start + 1, Constant::dss_start );
	strncpy( dstr + Constant::dss_start, dna + pos + 3, Constant::dss_end );
	dstr[Constant::dss_size()] = '\0';
	try{
	    dsscount[s2i(dstr)]++;
	    c_dss += 1;
	} catch (InvalidNucleotideError e) {
	    //cerr << "Donor splice site " << dstr << " contains an invalid character." << endl;
	}
	/*/output for seqlogo of ass 
	static char *dssfenster = new char[11];
	strncpy( dssfenster, dna + pos+1-4, 10);
	dssfenster[10]='\0';
	cout << "cds dss "  << dssfenster << endl;*/
    }
}



/*
 * IntronModel::buildProbabilities
 */
void IntronModel::buildProbabilities(const AnnoSequence* annoseq){
    
    Seq2Int s2i(k+1);
    State* in;
    State* exon, *lastexon;
    const Gene *curgene;
    seqProb(-1,0);
    int numErrorSplicesites = 0;
    /*
     * process the genes
     */

    const AnnoSequence *as = annoseq;
    while (as){
      sequence = as->sequence;
      curgene = dynamic_cast<const Gene*> (as->anno->genes); // assume here that sequences with multiple genes have already been split
      if (curgene) {
	  gweight = curgene->weight;
	  if (curgene->exons){
	      // upstream UTR introns and igenic region
	      processSequence( sequence + k, sequence + curgene->geneBegin()-1-20); // could be more accurate here 
	      if (curgene->utr5exons) {
		  // UTR introns
		  lastexon = curgene->utr5exons;
		  exon = lastexon->next;
		  while (exon) {
		      processSequence (sequence + lastexon->end + 1 + k + Constant::dss_end + DSS_MIDDLE, 
				       sequence + exon->begin - 1 - Constant::ass_start - ASS_MIDDLE 
				       - Constant::ass_upwindow_size);
		      lastexon = exon;
		      exon = exon->next;
		  }
	      }
	  
	      // the intronic sequences (separately, so splice site exceptions are ignored here
	      /*in = curgene->introns;
		while( in ){
		processSequence(sequence + in->begin + Constant::dss_end + DSS_MIDDLE + k,
		sequence + in->end - ASS_MIDDLE - Constant::ass_start -
		Constant::ass_upwindow_size);
		in = in->next;
		}*/
	      // splice sites
	      in = curgene->introns;
	      while( in ) {
		  try {
		      processDSS(sequence, in->begin-1);
		      processASS(sequence, in->end + 1, true);
		  } catch (IntronModelError e) {
		      numErrorSplicesites++;
		      if (verbosity) {
			  if (numErrorSplicesites <= 20) { 
			      cerr << "Sequence " << curgene->seqname << ":" << endl;
			      cerr << e.getMessage() << endl;
			  } 
			  if (numErrorSplicesites == 20) 
			      cerr << "further detailed output of splice site errors supressed." << endl;
		      }
		  }
		  in = in->next;
	      }
	      // downstream UTR introns and igenic region
	      if (curgene->utr3exons) {
		  // UTR introns
		  lastexon = curgene->utr3exons;
		  exon = lastexon->next;
		  while (exon) {
		      processSequence (sequence + lastexon->end + 1 + k + Constant::dss_end + DSS_MIDDLE, 
				       sequence + exon->begin - 1 - Constant::ass_start - ASS_MIDDLE 
				       - Constant::ass_upwindow_size);
		      lastexon = exon;
		      exon = exon->next;
		  }
	      }
	      processSequence( sequence + curgene->geneEnd() +  1 + k,
			       sequence + strlen(sequence) - 1 );
	  }
      }
      as = as->next;
    }
    
    /*
     * determine emission probabilities
     */
    {
    vector<Double> patternprobs(emicount.size());
    makeProbsFromCounts(patternprobs, emicount, k, patpseudo, false);
    emiprobs.probs.resize(emicount.size());
    emiprobs.order = k;
    computeEmiFromPat(patternprobs, emiprobs.probs, k);    
    }

    /*
     * read in the other splice sites
     */
    gweight = 1;
    readSpliceSites();
    
    if (numErrorSplicesites > 5 && verbosity)
	cerr << "number of non " << (Constant::dss_gc_allowed ? "GT/AG or GC/AG" : "GT/AG") 
	     << " introns: " << numErrorSplicesites << endl;
    
    /*
     * Compute the frequencies of the acceptor splice sites
     * Use pseudocounts.
     */
    // Important: If this computation is changed then readProbabilities needs to be changed also
    
    assprobs.resize( asscount.size() );
    Double together = Double(c_ass) +  asscount.size() * asspseudo;
    for( int i = 0; i < asscount.size(); i++ )
	assprobs[i] = (Double(asscount[i]) + asspseudo) / together;
    
    // make ass CRF features
    if (Constant::CRFtrainSS && Constant::ass_maxbinsize > 0){
      assBinProbs.reset();
      for (int i=0; i < assprobs.size(); i++ )
        assBinProbs.addProb(assprobs[i]);
      assBinProbs.trainBins(Constant::ass_maxbinsize);
      cout << "ass bin boundaries" << endl;
      assBinProbs.printBoundaries();
    }

    /*
     * Compute frequencies of the donor splice sites using pseudocounts.
     */
    makeDSSProbs();

    /*
     * Initialize the ASS Motif
     */
    
    assMotif->makeProbs(); // this is GC-content dependent
    hasSpliceSites = true;
}

void IntronModel::buildLenDist(const AnnoSequence* annoseq){
    /*
     *  Determine the maximal length of an intron, allocate a vector for the lengths
     */
    int maxintronlength = 0, len;
    State * tintron;
    for (AnnoSeqGeneIterator geneit = AnnoSeqGeneIterator(annoseq); geneit.hasMoreElements(); ++geneit) {
	const Transcript* curgene = geneit.gene;
	tintron = curgene->introns;
	while (tintron) {
	    len = tintron->length();
	    if (len < 0)
		throw IntronModelError("Intron has negative length.");
	    if (len > maxintronlength) {
		maxintronlength = len;
	    }
	    tintron = tintron->next;
	}
	curgene = curgene->next;
    }
    if (maxintronlength < d){
	if (verbosity)
	    cout << "WARNING: /Intronmodel/d larger than maximal intron length in training set (" << maxintronlength << ")." << endl
		 << "If you expect this to be because of an unrepresentative training set or an insufficient sample size, " << endl
		 << "then decrease the parameter /Intronmodel/d. I will go ahead and assume there are no introns longer than " << d << endl;
    }

    
    // store all intron lengths in intlencount
    intlencount.assign(maxintronlength + 1, 0);   
    for (AnnoSeqGeneIterator geneit = AnnoSeqGeneIterator(annoseq); geneit.hasMoreElements(); ++geneit) {
	const Transcript* curgene = geneit.gene;
	tintron = curgene->introns;
	while (tintron) {
	    intlencount[tintron->length()] += 1; //curgene->weight;
	    tintron = tintron->next;
	}
	curgene = curgene->next;
    }
    

    Smooth sm(minwindowcount, slope_of_bandwidth);
    sm.smoothCounts(intlencount, lenDist, d+1);
    scaleDblVector(lenDist, 1.0);
    //d = sm.geoCutOff(intlencount, lenDist);

    if (verbosity)
	cout << "d = " << d << endl;
    int sumOfLong = 0;
    introns_d = 0;
    for (int i = 0; i <= d; i++) {
	introns_d += intlencount[i];
    }
    for (int i = d+1; i < intlencount.size(); i++) {
	sumOfLong += intlencount[i] * i;
    }

    if (sumOfLong > 0)
	mal = sumOfLong / (introns - introns_d)-d;
    else 
	mal = 10; // some small number
    probShortIntron = Double(1.0)/(Double(1.0) + mal * lenDist[d]);
    //probShortIntron = (double) introns_d/introns;
    //mal = Double(introns-introns_d)/(lenDist[d]*Double(introns_d));

    if (verbosity && introns > 0) {
	cout << 100*introns_d/introns << "% of the introns were of length at most " << d << endl;
	cout << "In the model " << 100 * probShortIntron << "% of the introns are of length at most " << d << endl;
	cout << "P(len = d  ) = " << lenDist[d] * probShortIntron << endl;
	cout << "P(len = d+1) = " 
	     << (Double(1.0)-probShortIntron) / mal << endl;    
	cout << "model:mean length additional to " << d << " : " << mal << endl;
	cout << "real :mean length additional to " << d << " : " << sumOfLong / (introns - introns_d)-d << endl;
    }
/*
    Smooth sm(minwindowcount, slope_of_bandwidth);
    int sumOfLong = 0;
    d = sm.geoCutOff(intlencount, lenDist);
    
    if (verbosity)
	cout << "computed d = " << d << endl;
    introns_d = 0;
    for (int i = 0; i <= d; i++) {
	introns_d += intlencount[i];
    }
    for (int i = d+1; i < intlencount.size(); i++) {
	sumOfLong += intlencount[i] * i;
    }
    cout << "#introns_d/#introns = " << introns_d << " / " << introns 
	 << " = " << (double) introns_d/introns << endl;
    
    // Information reg. calculation of geom. Intron transition probability
     
    probShortIntron = (double) introns_d/introns;
    
    if (introns > introns_d) {
	//mal = Double(introns-introns_d)/(lenDist[d]*Double((double)introns_d));
	mal = sumOfLong / (introns - introns_d)-d;
	if (verbosity){
	    cout << "weighed mean additional intron length of the " << introns - introns_d 
		 << " introns longer than d = " << d <<" :" 
		 << mal << endl;
	    cout << "P(len=d+1)=" << Double(introns-introns_d)/introns/mal  
		 << endl;
	    cout << setprecision(2) << 100 * probShortIntron 
		 <<"% out of the weighed " << introns << " introns were <= " << d << endl;
	}
    }
    else if (verbosity)
	    cout << "All weighed " << introns << " introns were of length <= " << d << endl;

    if (verbosity)
	cout << "P(len=d) = " << lenDist[d]*Double((double)introns_d/introns) << endl;*/
}



/*
 * IntronModel::makeDSSProbs
 */

void IntronModel::makeDSSProbs ( ){

    // give all neighbor-patterns a little of the weight
    dssprobs.assign( dsscount.size(), 0.0);
    for (int i=0; i < dsscount.size(); i++ ){
	dssprobs[i] += dsscount[i] + dsspseudo;
	for (int j=0; j < Constant::dss_size(); j++) {
	    int nj = (i % POWER4TOTHE(j+1))/POWER4TOTHE(j);
	    for (int nn=0; nn<4; nn++){
		int pn = i +(nn-nj)*POWER4TOTHE(j);
		dssprobs[pn] += dssneighborfactor * (dsscount[i] + dsspseudo );
	    }
	}
    } 

    Double normsum = 0.0;
    for (int i=0; i < dssprobs.size(); i++ )
	normsum += dssprobs[i];
    
    if (verbosity>1) 
      cout << "normsum = " << normsum << ",  c_dss=" << c_dss << endl;
    for (int i=0; i < dssprobs.size(); i++ )
	dssprobs[i] /= normsum;

    // make dss CRF features
    if (Constant::CRFtrainSS && Constant::dss_maxbinsize > 0){
      dssBinProbs.reset();
      for (int i=0; i < dssprobs.size(); i++ )
	dssBinProbs.addProb(dssprobs[i]);
      dssBinProbs.trainBins(Constant::dss_maxbinsize);
      cout << "dss bin boundaries" << endl;
      dssBinProbs.printBoundaries();
    }
}

/* 
 * storeGCPars
 * store GC dependent model probabilities in the respective array variables
 */
void IntronModel::storeGCPars(int idx){
  GCemiprobs[idx] = emiprobs;
  std::stringstream out;
  out << "intron emiprob gc" << (idx+1);
  GCemiprobs[idx].setName(out.str());
  if (assMotif)
    GCassMotif[idx] = *assMotif;
  GCprobShortIntron[idx] = probShortIntron;
  GCmal[idx] = mal;
}

/*
 * IntronModel::printProbabilities
 */
void IntronModel::printProbabilities( int idx, BaseCount *bc, const char* suffix ){
   try{     
      const char* fname = Properties::getProperty("/IntronModel/outfile"); 
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
	  cout << "Writing intron model parameters [" << idx+1 << "] to file " << filename << "." << endl;
	if (idx == 0 ) { 
	  Seq2Int s2i_a(Constant::ass_size());
	  Seq2Int s2i_d(Constant::dss_size());

	  ofstrm << "#intron model parameters\n# begin of content independent part" << endl;
	  ofstrm << "#\n# ASS probabilities\n"
		 << "#only nonpseudocount values are shown\n[ASS]\n";
	  Double mincount;
	  ofstrm << "# Size of vector\n" << assprobs.size() << endl;
	  ofstrm << "# c_ass (ASS count)\n" << c_ass << endl;
	  ofstrm << "# asspseudocount (added to all possible patterns, no matter if they occur)\n" 
		 << asspseudo << endl;
	  ofstrm << "# Probabilities * 1000\n";
	  // don't save the values resulting from the pseudocounts
	  // mincount is between the pseudocount value and the smallest 'real' value
	  mincount = Double(asspseudo+ 0.5)/(Double(c_ass) + assprobs.size() * asspseudo);
	  for( int i = 0; i < assprobs.size(); i++ ){
	    if( assprobs[i] > mincount )
	      ofstrm << s2i_a.inv(i) << "\t" << 1000 * assprobs[i] << endl;
	  }
	  if (assBinProbs.nbins > 0){
	      ofstrm << "\n#\n# ASS binned probabilities\n[ASSBIN]\n";
	      assBinProbs.write(ofstrm);
	  }
		
	  ofstrm << "\n\n#\n# DSS probabilities\n#only nonpseudocount values are shown\n[DSS]\n";
	  ofstrm << "# Size of vector\n" << dssprobs.size() << endl;
	  ofstrm << "# c_dss (DSS count)\n" << c_dss << endl;
	  ofstrm << "# dsspseudocount (added to all possible patterns, no matter if they occur)\n" 
		 << dsspseudo << endl;
	  ofstrm << "# Probabilities * 1000\n";
	  for( int i = 0; i < dssprobs.size(); i++ ){
	    ofstrm << s2i_d.inv(i) << "\t" << 1000 * dssprobs[i] << endl;
	  }
	  if (dssBinProbs.nbins > 0){
	      ofstrm << "\n#\n# DSS binned probabilities\n[DSSBIN]\n";
	      dssBinProbs.write(ofstrm);
	  }
		    
	  // save the length distribution
	  //--------------------------------------------
	  ofstrm << "\n#\n# Length probabilities\n#" << endl;
	  ofstrm << "[LENGTH]" << endl;
	  ofstrm << "# The 'd' variable\n" << d << endl;
	  ofstrm << "# The length probabilities from '0' to 'd' (*1000) "
		 << endl;
	  /*for( Vector<Double>::iterator it = intlenprobs.begin();
	    it != intlenprobs.end();
	    ++it )
	    ofstrm  << *it << endl;
	  */
	  for (int i=0; i< lenDist.size(); i++) {
	    ofstrm << 1000 * lenDist[i] << endl;
	  }

	  ofstrm << "# end of content independent part" << endl;		    
	}
		
	ofstrm << "[" << idx+1 << "]" << endl;
	if (bc)
	  ofstrm << "# (a,c,g,t)= " << *bc << endl;

	ofstrm << "#\n# Probabilities file for the intron model\n#\n#"
	       << endl;

	// Transition probabilities into an intron
	// Exceptionally managed by intron model.
	// ---------------------------------------
	ofstrm << "# Transition probabilities\n#" << endl;
	ofstrm << "[TRANSITION]" << endl;
	ofstrm << "# the probability of an intron of length at most d\n" 
	       << GCprobShortIntron[idx] << endl;
	ofstrm << "# mean additional length of introns with length > d\n"
	       << GCmal[idx] << endl;
		
	// save the emission probabilities
	//-----------------------------------------------
	int precision = LLDouble::getOutputPrecision();
	LLDouble::setOutputPrecision(4);
	ofstrm << "\n#\n# The emission probabilities of introns\n#" << endl;
	ofstrm << "[EMISSION]" << endl;
	ofstrm << "# size of the emission vector\n" << GCemiprobs[idx].probs.size() << endl;
	ofstrm << "#k=\n" << k << endl;
	ofstrm << "# patpseudo : pseudocount for sequence patterns\n" << patpseudo << endl;
	Seq2Int s2i(k+1);
	for( int i = 0; i != GCemiprobs[idx].probs.size(); ++i )
	  ofstrm << s2i.inv(i) << '\t' << GCemiprobs[idx].probs[i] << endl;
	
	LLDouble::setOutputPrecision(precision);

	// ************* frequencies fo the (k+1) tuples
	ofstrm << "\n# patterns:" << endl;
	for( int j = 0; j < emicount.size(); j++ )
	  ofstrm << "#\t" << s2i.inv(j) << "\t" << emicount[j] << endl;
	  //  *****************************************

	ofstrm << "\n# motif upstream of acceptor splice site\n[ASSMOTIF]" << endl;
	GCassMotif[idx].write(ofstrm);
				
	ofstrm.close();
      }
    } catch( ProjectError ) { 
      if (verbosity)
	cout << "Intron model parameters not saved." << endl;
    }
}
