/**********************************************************************
 * file:    merkmal.cc
 * licence: Artistic Licence, see file LICENCE.TXT or
 *          http://www.opensource.org/licenses/artistic-license-2.0.php
 * descr.:  Features (German: Merkmale) for the Conditional Random Field
 * authors: Mario Stanke, mario@gobics.de
 *
 * date    |   author      |  changes
 * --------|---------------|------------------------------------------
 * 27.07.08| Mario Stanke  | creation of the classes
 **********************************************************************/

#include "merkmal.hh"

// project includes
#include "geneticcode.hh"
#include "projectio.hh"
#include "namgene.hh"
#include "evaluation.hh"
#include "extrinsicinfo.hh"
#include "gene.hh"
#include "statemodel.hh"
#include "utrmodel.hh"

// standard C/C++ includes
#include <fstream>


vector<AnnoSequence*> CRF::evalAnnoSeqs;
FeatureCollection *CRF::evalExtrinsicFeatures = NULL;
int CRF::statecount = 0;
StateModel **CRF::states = NULL;

void MMGroup::registerPars(Parameters* parameters, double alphavalue){
    this->parameters = parameters;

    for (int idx=0; idx < getNumPars(); idx++)
	parameters->addMerkmal(getFactor(idx), alphavalue);
}

void MMGroup::addCount(int index, int count) {
    if (parameters)
	parameters->addCount(offset+index, count);
}

Double* PatMMGroup::getFactor(int index){
    return &probs[index];
}

string PatMMGroup::verbalDescription(int index){
    ostringstream descr;
    descr << parname;
    descr << " ";
    descr << Seq2Int(order+1).INV(index);
    return descr.str();
}

Double PatMMGroup::getMinProb(float qthresh){
  return quantile(probs, qthresh);
}

void FramedPatMMGroup::getFramePat(int index, int &frame, int &pattern){
    frame = 0;
    while (frame < 3 && index >= probs[frame].size()){
	index -= probs[frame].size();
	frame++;
    }
    if (frame > 2)
	throw ProjectError("FramedPatMMGroup: Index out of range.");
    pattern = index;
}

Double* FramedPatMMGroup::getFactor(int index){
    int frame, pattern;
    getFramePat(index, frame, pattern);
    return &probs[frame][pattern];
}

string FramedPatMMGroup::verbalDescription(int index){
    int frame, pattern;
    getFramePat(index, frame, pattern);
    ostringstream descr;
    descr << parname;
//    descr << " idx ";
//    descr << index;
    descr << " frame ";
    descr << frame;
    descr << " ";
    descr << Seq2Int(order+1).INV(pattern);
    return descr.str();
}


/*
 * BinnedMMGroup::trainBins
 *
 * Find bin boundaries and starting average values by computing the quantiles of the origprobs
 */

void BinnedMMGroup::trainBins(int maxNumBins){
    int i, j, n;
    int minBinSize;
    Double sum;

    avprobs.clear();
    bb.clear();
    if (maxNumBins < 1){
	nbins=0; // use no binning
	return;
    }
    if (origprobs.size() < 1){
	throw ProjectError("Tried to train binned feature group without training data.");
    }
    std::sort(origprobs.begin(), origprobs.end());

    minBinSize = origprobs.size() / maxNumBins;
    if (minBinSize < 1)
	minBinSize = 1;
    nbins = 0;
    i = j = 0;
    
    while (i < origprobs.size()){
	sum = 0.0;
	n = 0; // number of items in current bin
	while(i < origprobs.size() && (n < minBinSize || origprobs[i] == origprobs[i-1])){
	    sum += origprobs[i];
	    n++;
	    i++;
	}
	avprobs.push_back(sum / n); // (arithmetic average for now)
	if (i < origprobs.size()) {
	    bb.push_back((origprobs[i] + origprobs[i-1]) / 2); // warning origprobs[i] and origprobs[i-1] could be almost equal
	}
	nbins++;
    }
    monotonify();
}

void BinnedMMGroup::monotonify(){
    if (monotonic == 0)
	return;
    // if an adjacent pair of probs violates the monotonicity criteria, set both to their average
    for (int j=0; j < nbins-1; j++){
	if (monotonic * (avprobs[j+1] - avprobs[j]) < 0){
	    cout << "monotonify: set " << verbalDescription(j) << "= " << avprobs[j] << " and "
		 << verbalDescription(j+1) << "= " << avprobs[j+1] << " to their average." << endl;
	    avprobs[j] = avprobs[j+1] = (avprobs[j] + avprobs[j+1]) / 2;
	}
    }
    // it will be closer to monotinicity now, but it is not monotone in general, TODO if time
}

void BinnedMMGroup::printBoundaries(){
    for (int j=0; j < nbins-1; j++)
	cout << "bb[" << j << "]= " << bb[j] << endl;
}

int BinnedMMGroup::getIndex(Double p){
    // make a binary search to find the bin, in which p belongs
    int a=0, b=nbins-1, m;
    while (a < b){
	m = (a+b)/2;
	if (p < bb[m])
	    b = m;
	else
	    a = m+1;
    }
    if (a > 0 && a < nbins -1 && !(bb[a-1] <= p && bb[a] > p))
	cerr << "getIndex error " << p  << " " << bb[a-1] << " " << bb[a] << endl;
    return a;
}

void BinnedMMGroup::write(ostream &out){
    out << "# number of bins:\n" << nbins << endl;
    out << "# monotonic: ";
    if (monotonic>0)
	out << "increasing";
    else if (monotonic<0)
	out << "decreasing";
    else 
	out << "no";
    out << endl;
    out << "# lower bin boundary <= average probability" << endl;
    out << "\t" << avprobs[0] << endl;
    for (int i=1; i < nbins; i++)
	out << bb[i-1] << "\t" << avprobs[i] << endl;
}

void BinnedMMGroup::read(ifstream &in){
    if (!in)
	throw ProjectError("Error when reading in BinnedMMGroup.");
    in >> comment >> nbins;
    if (nbins < 1)
	throw ProjectError("BinnedMMGroup:: less than one bin.");
    avprobs.resize(nbins);
    bb.resize(nbins-1);
    in >> comment >> avprobs[0];
    for (int i=1; i < nbins; i++)
	in >> bb[i-1] >> avprobs[i];
}

Double BinnedMMGroup::getMinProb(float qthresh){
  if (nbins==0)
    return quantile(origprobs, qthresh);
  else
    return quantile(avprobs, qthresh);
}

Double* BinnedMMGroup::getFactor(int index){
    return &avprobs[index];
}

string BinnedMMGroup::verbalDescription(int index){
    ostringstream descr;
    descr << parname;
    descr << " bin[";
    if (index > 0)
	descr << bb[index-1];
    else
	descr << "-infty";
    descr << ", ";
    if (index < nbins-1)
	descr << bb[index];
    else
	descr << "+infty";
    descr << ")";
    return descr.str();
}

int Parameters::addMMGroup(MMGroup* mmgroup, double alphavalue){
    int offset = weights.size();
    // check consistency of offset 
    if (mmgroups && mmgroups->size() > 0){
	MMGroup* lastMMGroup = (*mmgroups)[mmgroups->size()-1];
	if (offset !=  lastMMGroup->getOffset() + lastMMGroup->getNumPars())
	   throw ProjectError("Parameters::addMMGroup: inconsistency detected");
    }
    mmgroup->registerPars(this, alphavalue);
    if (!mmgroups)
	mmgroups = new vector<MMGroup*>;
    mmgroups->push_back(mmgroup);
    mmgroup->setOffset(offset);
    return offset;
}

void Parameters::addMerkmal(Double* parptr, double alphavalue){
    if (!parptrs)
	parptrs = new vector<Double*>;
    parptrs->push_back(parptr);
    weights.push_back(parptr->log());
    alpha.push_back(alphavalue);
    counts.push_back(0);
}

void Parameters::addWeights(vector<double> &h){
    if (h.size() != weights.size())
	throw ProjectError ("Parameters::addWeights: vector sizes not matching.");
    for (int i=0; i < h.size(); i++)
	weights[i] += h[i];
}

void Parameters::setWeights(vector<double> &v){
    if (v.size() != counts.size())
	throw ProjectError ("Parameters::addWeights: vector sizes not matching.");
    weights = v;
}

void Parameters::smoothFeatures(){
    for (int i=0; i<mmgroups->size(); i++)
	(*mmgroups)[i]->smooth();
}

void Parameters::updatePars(){
    if (!parptrs)
	return;
    if (parptrs->size() != weights.size())
	throw ProjectError("Unequal parameter vectors in Parameters::updatePars.");
    for (int idx=0; idx < weights.size(); idx++)
	(*(*parptrs)[idx]) = LLDouble::exp(weights[idx]);
}

void Parameters::updateWeights(){
    if (!parptrs)
        return;
    if (parptrs->size() != weights.size())
        throw ProjectError("Unequal parameter vectors in Parameters::updatePars.");
    for (int idx=0; idx < weights.size(); idx++)
	weights[idx] = ((*parptrs)[idx])->log();
}

string Parameters::verbalDescription(int index){
    if (!mmgroups)
	return "";
    for (int i=0; i<mmgroups->size(); i++){
	int offset = (*mmgroups)[i]->getOffset();
	if ( offset <= index && offset + (*mmgroups)[i]->getNumPars() > index)
	    return (*mmgroups)[i]->verbalDescription(index-offset);
    }
    return "error";
}



/*
 * CRF::train
 */
void CRF::train(Parameters* startPars, vector<AnnoSequence*> trainGenes, FeatureCollection &extrinsicFeatures){
    /*
     * check whether there are hints for the evaluation set
     * TODO
     */
    evalExtrinsicFeatures = new FeatureCollection();


    /*
     * decide on the training method
     */
    string method;
    try {
	method = Properties::getProperty("CRF_TRAIN");
    } catch  (...) {
        method = "OLM";
    }
    // train parameters
    if (method == "IIS") {
	cout << "use IIS algorithm (Improved Iterative Scaling) for training" << endl;
	improvedIterativeScaling(startPars, trainGenes, extrinsicFeatures);
    } else {
	cout << "use OLM algorithm (Online Large Margin) for training" << endl;
	onlineLargeMarginTraining(startPars, trainGenes, extrinsicFeatures);
    }
}

/*
 * onlineLargeMarginTraining (OLM)
 * 
 */
void CRF::onlineLargeMarginTraining(Parameters* parameters, vector<AnnoSequence*> trainGenes, FeatureCollection &extrinsicFeatures){
    inCRFTraining = true;
    NAMGene namgene;
    Double viterbiPathEmiProb, correctPathEmiProb;
    double correctScore, viterbiScore;
    double loss; // value of loss function
    AnnoSequence* as;
    vector<int> a; // count vector of correct features minus viterbi features
    vector<double> h(parameters->size(), 0.0); // weight difference vector
    vector<double> v(parameters->size(), 0.0); // accumulated weights for averaging at the end
    int N; // number of rounds
    int m; // training set size
    int numiterations=0, u; // at most Nm
    bool termination=false;
    double target;
    StatePath *viterbiPath, *correctPath, *condensedViterbiPath=NULL;//, *condensedCorrectPath=NULL;
    int chunksize = 100; // consider this many training genes to be one training example (will be increased a little below to fit)
    int maxTrainSeqs = -1; // use only this many training genes for CRF iterations per round (-1: use all)
    vector<double> bs;
    vector< vector<int> > viterbiCounts, correctCounts; // holds all feature counts of a chunk
    vector<int> chunkViterbiCounts(parameters->size(), 0), chunkCorrectCounts(parameters->size(), 0);
    vector<SequenceFeatureCollection*> sfcs; // sequence feature collection for each training sequence
    UtrModel::setTtsSpacing(1); // so the true parses are possible in the model

// temp for optimizing parameters
    double scorediffweight = 0.05;
    double lossweight = 2.0;
    try {
        scorediffweight = Properties::getdoubleProperty("scorediffweight");
    } catch (...){}
    try {
        lossweight = Properties::getdoubleProperty("lossweight");
    } catch (...){}
    cout << "scorediffweight= " << scorediffweight << "\t lossweight= " << lossweight << endl << flush;

    
    /*
     * preparations
     */
    try {
	N = Properties::getIntProperty("CRF_N");
	if (N < 1){
	    cerr << "CRF_N < 1 not allowed. Please use --CRF=0 to turn off CRF training. Will use CRF_N = 1..." << endl;
	    N = 1;
	}
    } catch  (...) {
	N = 3;
    }

    /*
     * precompute and save correct paths
     * removed since this takes 72 bytes for every state (basically 70Mb per mega base)
     */
    //vector<StatePath*> correctPaths;
    m = trainGenes.size();
    chunksize = (chunksize > m) ? m : (int)((float)m/(int)(m/chunksize) + .99999); // increase chunksize so the chunks are all equal size and fit in training set size evenly
    cout << "training set size: " << m << " chunk size: " << chunksize << endl;
    if (evalAnnoSeqs.size() <= 0)
      cout << "have no evaluation set" << endl;
    vector<int> numCorrect;
    numCorrect.resize(m, 0);
    bs.resize(chunksize, 0.0);
    viterbiCounts.resize(chunksize);
    correctCounts.resize(chunksize);

    for (int round = 1; round <= N && !termination; round++){
	/*
	 * evaluate parameters
	 */
	if (evalAnnoSeqs.size() > 0){ // have evaluation set
            UtrModel::setTtsSpacing(10); 
	    cout << "Evaluating on " << evalAnnoSeqs.size() << " sequences " << endl;
	    Evaluation* eval = predictAndEvaluate(evalAnnoSeqs, *evalExtrinsicFeatures);
	    target = (3*eval->nukSens + 2*eval->nukSpec + 4*eval->exonSens + 3*eval->exonSpec + 2*eval->geneSens + 1*eval->geneSpec)/15;
	    cout << "Evaluation before round " << round << " target = " << target << endl;
	    eval->finishEvaluation();
	    eval->print();
	    // TODO Kriterien zum Abbruch des Trainings, z.B. round > 2 und wurde nicht mehr besser
	    UtrModel::setTtsSpacing(1);
	    termination = false;
	    delete eval;
	}
	if (termination){
	    cout << "Terminated training algorithm" << endl;
	    break;
	}
	cout << "*** CRF training round " << round << " of " << N <<  " ***" << endl;
	for (int i=0; i < m && (maxTrainSeqs < 0 || i < maxTrainSeqs); i += chunksize) {
	    // reset chunk counters
	    bs.assign(chunksize, 0.0);
	    // loop over chunk	   
	    for (u = 0; u < chunksize && i + u < m; u++) {
		as = trainGenes[i+u];
		cout << "using training sequence " << as->seqname << " (length=" << as->length << ")." << endl;
		if (round == 1) {
		    SequenceFeatureCollection &sfc = extrinsicFeatures.getSequenceFeatureCollection(as->seqname);
		    if (sfcs.size() != i+u) {
			throw ProjectError("inconsistency with index in sfcs");
		    }
		    sfc.prepare(as, true);
		    sfcs.push_back(&sfc);
		}
	    
		/*
		 * compute Viterbi path
		 */
		viterbiPath = namgene.getTrainViterbiPath(as->sequence, sfcs[i+u]);
		condensedViterbiPath = StatePath::condenseStatePath(viterbiPath);
		//condensedViterbiPath->print();
		Gene* viterbiGenes = (Gene*) condensedViterbiPath->projectOntoGeneSequence("v");
		if (viterbiGenes && round < 4)
		    printGeneList(viterbiGenes, as, false, false);
		/*
		 * correct path
		 */
		correctPath = StatePath::getInducedStatePath(as->anno->genes, as->length, false);
		//condensedCorrectPath = StatePath::condenseStatePath(correctPath);
		//condensedCorrectPath->print();
		if (round < 4)
		    as->anno->printGFF();

		/*
		 * compare Viterbi with correct path
		 */
		if (*correctPath == *viterbiPath){
		    cout << "Gene(s) correctly predicted in " << as->seqname << ". I don't do nothin." << endl;
		    numCorrect[u+i]++;
		} else {
		    cout << "Viterbi wrong in " << as->seqname << "." << endl;
		    // testing for now: Only count features falling in range of annotated gene(s)
		    int countFrom = as->anno->genes->geneBegin()-10000;
		    int countTo = as->anno->genes->geneEnd()+10000; //TODO: deal with multi-gene sequences
		    cout << "counting features only in region " << countFrom << ".." << countTo << endl;
		    StateModel::setCountRegion(countFrom, countTo);
		    parameters->resetCounts();
		    viterbiPathEmiProb = namgene.getPathEmiProb(viterbiPath, as->sequence, *sfcs[i+u], countFrom, countTo);
		    viterbiCounts[u] = parameters->getCounts();
		    viterbiScore = viterbiPathEmiProb.log();
		    cout << "Viterbi pathEmiProb=" << viterbiPathEmiProb << " score = " << setprecision(6) << viterbiScore << endl;
		    if (round == 1){
			cout << "counts of viterbi path" << endl;
			parameters->print<int>(viterbiCounts[u], 5, true);
		    }
		    
		    parameters->resetCounts();
		    correctPathEmiProb = namgene.getPathEmiProb(correctPath, as->sequence, *sfcs[i+u], countFrom, countTo);
		    correctCounts[u] = parameters->getCounts();
		    correctScore = correctPathEmiProb.log();
		    cout << "correct pathEmiProb=" << correctPathEmiProb << " score = " << setprecision(6) << correctScore << endl;
		    if (round == 1){
			cout << "counts of correct path" << endl;
			parameters->print<int>(correctCounts[u], 5, true);
		    }
		    parameters->resetCounts();
		    
		    if (!(correctPathEmiProb > 0.0)){
			cout << "Correct path IMPOSSIBLE in model. Can't do anything about this by only changing weights." << endl;
		    } else {
			/*
			 * compute loss function that measures how far Viterbi is off
			 */
			loss = lossFct(as->anno->genes, viterbiGenes); // TODO
			if (viterbiScore < correctScore){
			    cout << "viterbi score < correct score" << endl;
			}
			bs[u] = scorediffweight * ((viterbiScore > correctScore)? (viterbiScore - correctScore) : 0) + lossweight * loss;
			cout << "b = scorediffweight * ((viterbiScore > correctScore)? (viterbiScore - correctScore) : 0) + lossweight * loss = " 
			     << scorediffweight << " * ( " << viterbiScore << " - " << correctScore << ") + " << lossweight 
			     << " * " << loss << " = " << bs[u] << endl;
		    }
		} // *correctPath != *viterbiPath
		cout << "correct: " << numCorrect[u+i] << " of " << N << endl;
		// clean up
		delete correctPath;
		delete viterbiPath;
		if (condensedViterbiPath)
		    delete condensedViterbiPath;
	    } // end of loop over chunk elements
	    
	    vector<double> bsc = capOutliers(bs);
	    h.assign(h.size(), 0.0);
	    cout << "updating parameters after chunk processing" << endl;
	    for (int uu = 0; uu < bsc.size(); uu++){
		cout << "chunk item no " << uu << "\tcapped b = " << bsc[uu];
		if (bsc[uu] > 0.0){
		    int numfeat = 0;
		    a = correctCounts[uu];
		    for (int j=0; j < a.size(); j++)
			a[j] -= viterbiCounts[uu][j];
		    double normOfA = 0.0;
		    for (int j=0; j < a.size(); j++){
			normOfA += parameters->getAlpha(j) * a[j]*a[j];
			if (a[j] != 0)
			    numfeat++;
		    }
		    if (normOfA > 0.0) {
			double scalefactor = bsc[uu] / normOfA;
			for (int j=0; j < h.size(); j++){
			  // This is a hack against the following problem caused by the fact that currently not all features are CRF-trained
			  // Consider a difference only in a small initial exon. Then the only ss/exon feature difference is one 
			  // false dss feature, which gets all the change (in one case something like a factor of 65)
			  double change = scalefactor * parameters->getAlpha(j) * a[j];
			  if (change > 0.5) // = factor of 2
			    change = 0.5;
			  if (change < -0.5)
                            change = -0.5;
			  h[j] += change;
			}
		    }
		    cout << "\tnumfeat=" << numfeat << "\tnormOfA=" << normOfA;
		}
		cout << endl;
	    }
	    /*
	     * update weights ...
	     */
	    vector<double> startWeights = parameters->getWeights();     //TEMP
	    parameters->addWeights(h);
	    parameters->updatePars();
	    parameters->smoothFeatures();
	    parameters->updateWeights();
	    vector<double> endWeights = parameters->getWeights();       //TEMP
	    CRF::compareWeights(parameters, startWeights, endWeights, 30);       //TEMP
	    if (false){
		cout << "updated Parameters:" << endl;
		parameters->print(150, false);
	    }
	    
	    /*
	     * for testing: compare the scores of viterbi and the correct path under the new weights, b should be 0
	     *
	     viterbiPathEmiProb = namgene.getPathEmiProb(viterbiPath, as->sequence, sfc[i]);
	     viterbiScore = viterbiPathEmiProb.log();
	     correctPathEmiProb = namgene.getPathEmiProb(correctPath, as->sequence, sfc[i]);
	     correctScore = correctPathEmiProb.log();
	     b = viterbiScore - correctScore + loss;
	     cout << "with new pars: b = viterbiScore - correctScore + loss" << " = " << viterbiScore
	     << " - " << correctScore << " + " << loss << " = " << b << endl;
	     if (abs(b) > 0.001)
	     throw ProjectError("updating not correct");
	    */
	    
	    // add up the weight vector
	    for (int j=0; j < v.size(); j++)
		v[j] += parameters->getWeight(j);
	    numiterations++;
	
	    // ******** test evaluation
	    if (evalAnnoSeqs.size() > 0 && (numiterations <=200  || 400 < chunksize || ((int) numiterations/chunksize) % 4 == 0)){// 200 was 500, 400 was 100
		cout << "Evaluating on " << evalAnnoSeqs.size() << " sequences " << endl;
		UtrModel::setTtsSpacing(10);
		Evaluation* eval = predictAndEvaluate(evalAnnoSeqs, *evalExtrinsicFeatures);
		target = (3*eval->nukSens + 2*eval->nukSpec + 4*eval->exonSens + 3*eval->exonSpec + 2*eval->geneSens + 1*eval->geneSpec)/15;
		cout << "Evaluation after iteration " << numiterations << " and training sequence " << i + u << ". target = " << target << endl;
		eval->finishEvaluation();
		eval->print();
		UtrModel::setTtsSpacing(1);
		delete eval;
	    }
	    
	    //parameters->setWeights(startWeights); // undo changes, temp for testing
	    //parameters->updatePars();

	} // end of loop over training examples
	if (states && round < N){ // save intermediate parameters to file
	    char* suffix = new char[4];
	    sprintf(suffix, "%d", round);
	    cout << "saving intermediate parameters " << suffix << endl;
	    ContentDecomposition cd;
	    BaseCount bc;
	    for( int i = 0; i < statecount; i++ ) // first delete all files with this suffix
	      states[i]->printProbabilities( -1, &bc, suffix );
	    for (int idx=0; idx < cd.n; idx++) {
	      bc = cd.getBaseCount(idx);
	      for( int i = 0; i < statecount; i++ )
		states[i]->printProbabilities(idx, &bc, suffix);
	    }
	}
    } // end of loop over rounds
    cout << "N*m = " << N << " * " << m << " = " << N*m << ". Number of iterations: " << numiterations << endl;
    //for (int i=0; i < v.size(); i++)
    //	v[i] /= numiterations;
    cout << "using last parameters after training" << endl;
    //    parameters->setWeights(v);
    //    parameters->updatePars();
    // clean up
    // for (int i=0; i < correctPaths.size(); i++)
    // delete correctPaths[i];
    inCRFTraining = false;
}

/*
 * improvedIterativeScaling
 * Improved Iterative Scaling Training Algorithm (IIS)
 * Maximum Likelihood without regularization
 */
void CRF::improvedIterativeScaling(Parameters* parameters, vector<AnnoSequence*> trainGenes, FeatureCollection &extrinsicFeatures){
    inCRFTraining = true;
    NAMGene namgene;
    namgene.setNeedForwardTable(true); // for sampling
    int m=0; // training set size, only feasible
    int sample=100; // number of sample iterations per sequence
    double up = 0.0; // upper limit on the feature sum
    double score, totalScore = 0.0, target;
    AnnoSequence* as;
    int N; // number of rounds
    vector<int> counts;
    StatePath *sampledPath, *correctPath, *condensedSampledPath=NULL;//, *condensedCorrectPath=NULL;
    Double pathEmiProb;
    bool termination = false;
    vector<int> Nf(parameters->size(), 0); // actual count of features in training data
    vector<int> Sf(parameters->size(), 0); // sampled count of features in training data
    vector<double> Ef(parameters->size(), 0.0); // expected count of features in training data
    vector<bool> feasible; // for each training gene: whether it is possible in the model
    vector<SequenceFeatureCollection*> sfcs; // sequence feature collection for each training sequence
    vector<double> h; // weight difference vector
    /*
     * preparations
     */
    try {
        N = Properties::getIntProperty("CRF_N");
    } catch  (...) {
        N = 1;
    }
    if (evalAnnoSeqs.size() <= 0)
      cout << "Have no evaluation set." << endl;
    /*
     * iterate over optimization rounds
     */
    for (int round = 1; round <= N && !termination; round++){
	/*
         * evaluate parameters
         */
        if (evalAnnoSeqs.size() > 0){ // have evaluation set
            Evaluation* eval = predictAndEvaluate(evalAnnoSeqs, *evalExtrinsicFeatures);
            target = (3*eval->nukSens + 2*eval->nukSpec + 4*eval->exonSens + 3*eval->exonSpec + 2*eval->geneSens + 1*eval->geneSpec)/15;
            cout << "Evaluation before round " << round << " target = " << target << endl;
            eval->print();
            termination = false;
            delete eval;
        }
        if (termination){
            cout << "Terminated training algorithm" << endl;
            break;
        }

	/*
	 * loop over the training examples
	 */
	cout << "*** CRF(IIS) training round " << round << " of " << N <<  " ***" << endl;
	for (int k=0; k < parameters->size(); k++) // reset sampled counts to 0
	    Sf[k] = 0;
	totalScore = 0.0;
	for (int i=0; i < trainGenes.size(); i++){
	    as = trainGenes[i];
	    if (round > 1 && feasible.size() >i && !feasible[i]) {
		cout << "skipping infeasible seq " << as->seqname << endl;
		continue;
	    }
	    //cout << "computing expected feature counts for training example " << as->seqname << " (length=" << as->length << ")." << endl;
	    if (round == 1) {
		SequenceFeatureCollection &sfc = extrinsicFeatures.getSequenceFeatureCollection(as->seqname);
		if (sfcs.size() != i) {
		    throw ProjectError("inconsistency with index in sfcs");
		}
		sfc.prepare(as, true);
		sfcs.push_back(&sfc);
	    }
	    StatePath *viterbiPath = namgene.getTrainViterbiPath(as->sequence, sfcs[i]);
	    delete viterbiPath; // we actually don't need this in this algorithm, but namgene needs the forward table..
	    
	    /*
             * in the first round compute the feature sum for all training examples
             * this is done simultaneously because viterbi/forward/ automatically initializes for getPathEmiProb
             */
            if (round == 1 || feasible[i]){
                correctPath = StatePath::getInducedStatePath(as->anno->genes, as->length, false);
                //condensedCorrectPath = StatePath::condenseStatePath(correctPath);
                //cout << "correctPath:" << endl; condensedCorrectPath->print();
                parameters->resetCounts();
		//cout << "getting path emi prob" << endl << flush;
                pathEmiProb = namgene.getPathEmiProb(correctPath, as->sequence, *sfcs[i]);
		delete correctPath;
		//delete condensedCorrectPath;
                //cout << "correct pathEmiProb=" << pathEmiProb << endl;
		if (pathEmiProb > 0.0){
		    if (round == 1) {
			feasible.push_back(true);
			m++;
			if (as->length > up)
			    up = as->length;
			vector<int> counts = parameters->getCounts();
			for (int k=0; k < parameters->size(); k++)
			    Nf[k] += counts[k];
		    }
		    score = pathEmiProb.log();
		    totalScore += score;
		    cout << "training seq no " << m << " " << as->seqname << ", score= " << score << endl;
                } else {
		    if (round == 1) {
			cout << "training gene structure on " << as->seqname << " not feasible. Skipping this example from now on." << endl;
			feasible.push_back(false);
		    }
		}
            }
	    if (feasible[i]){
		for (int j=0; j < sample; j++){
		    sampledPath = namgene.getSampledPath(as->sequence, as->seqname);
		    condensedSampledPath = StatePath::condenseStatePath(sampledPath);
		    //cout << j+1 << "th sampled path:" << endl; condensedSampledPath->print();
		    parameters->resetCounts();
		    namgene.getPathEmiProb(sampledPath, as->sequence, *sfcs[i]);
		    vector<int> counts = parameters->getCounts();
		    for (int k=0; k < parameters->size(); k++)
			Sf[k] += counts[k];
		    delete sampledPath;
		    delete condensedSampledPath;
		}
	    }
	}
	if (round == 1){
	    cout << "training set has size " << m << endl;
	    cout << "total feature counts of all features: " << endl;
	    parameters->print<int>(Nf, true, true);
	    cout << "up = uppper limit for feature count of all training sequences and any gene structure = " << up << endl;
	}
	cout << "total Score of all training genes: " << setprecision(6) << totalScore << endl;
	/*
	 * compare expected counts with actual counts
	 */
	// cout << "sampled count of features:" << endl; parameters->print<int>(Sf, true, true);
	cout << "expected count of features:" << endl;
	for (int k=0; k < parameters->size(); k++)
	    Ef[k] = (double) Sf[k] / sample;
        parameters->print<double>(Ef, true, true);
	cout << "comparison actual -> expected feature counts:" << endl;
	CRF::compareWeights<int,double>(parameters, Nf, Ef);
	/*
	 * update weights
	 */
	if (round <= 50)
	    up = 300;
	else 
	    up = 500;
	h.resize(parameters->size());
	for (int k=0; k < h.size(); k++){
	    // TODO: regularization
	    if (Ef[k] == 0.0 && Nf[k] > 0){
		h[k] = 0.1; //arbitrary
	    } else if (Nf[k] == 0) {
		h[k] = -0.1;
	    } else {
		h[k] = log(Nf[k]/Ef[k])/up;
	    }
	    cout << "h[" << k << "] = " << h[k];
	    if (k%10 == 0 || k == h.size()-1)
		cout << endl;
	    else
		cout << "  ";
	}
	vector<double> startWeights = parameters->getWeights(); //TEMP
	parameters->addWeights(h);
	parameters->smoothFeatures();
	parameters->updatePars();
	vector<double> endWeights = parameters->getWeights();   //TEMP
	CRF::compareWeights(parameters, startWeights, endWeights, 50);   //TEMP
	cout << "updated Parameters:" << endl;
	parameters->print(50, false);
	// ******** test evaluation
	if (evalAnnoSeqs.size() > 0 && round == N){
	    cout << "Evaluating on " << evalAnnoSeqs.size() << " sequences " << endl;
	    Evaluation* eval = predictAndEvaluate(evalAnnoSeqs, *evalExtrinsicFeatures);
	    target = (3*eval->nukSens + 2*eval->nukSpec + 4*eval->exonSens + 3*eval->exonSpec + 2*eval->geneSens + 1*eval->geneSpec)/15;
	    cout << "Evaluation after last round " << round << " target = " << target << endl;
	    eval->print();
	    delete eval;
	}
    }
    inCRFTraining = false;
}


double CRF::lossFct(Transcript* correct, Transcript* predicted){
    Evaluation eval;
    eval.addToEvaluation(predicted, correct, bothstrands);
    cout << "nukFP=" << eval.nukFP << " nukFPinside=" << eval.nukFPinside << " nukFN=" << eval.nukFN << " nucUFPinside=" << eval.nucUFPinside << " nucUFN=" << eval.nucUFN << endl;
    return 0.001 * eval.nukFPinside + 0.001 * eval.nukFN // coding bases
      + 0.0004 * eval.nucUFPinside + 0.0004 * eval.nucUFN; // noncoding bases
}



/*
 * capOutliers
 * replace bs[u] by min{bs[u], cap}, where cap is the qth Quantile
 * this weighs down the effect of extreme outliers which could be misannotated or chanceless to be predicted correctly 
 */
vector<double> CRF::capOutliers(vector<double> bs){
    double q = 0.9, cap;
    vector<double> t = bs, bsc(bs.size(), 0.0);
    try {
      q = Properties::getdoubleProperty("capthresh");
    } catch (...){
      q = 0.9;
    }
    if (bs.size() > 0){
	std::sort(t.begin(), t.end());
	cap = t[(int) (q * (bs.size()-1))];
	cout << "capOutliers: " << 100*q << "% quantile= " << cap << endl;
	//cout << "t\tbs\tbsc" << endl;
	for (int i=0; i < bs.size(); i++){
   	    bsc[i] = (bs[i] <= cap)? bs[i] : 0; // this was set to cap instead of 0 before
	    if (bsc[i] < bs[i]){
		cout << "capped outlier " << bs[i] << " --> " << bsc[i] << endl;
	    }
	    //cout << t[i] << "\t" << bs[i] << "\t" << bsc[i] << endl;
	}
    }
    return bsc;
}
