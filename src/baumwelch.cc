/**********************************************************************
 * file:    baumwelch.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  Implementation of a special case of Baum-Welch parameter estimation
 * authors: Mario Stanke, mario@gobics.de
 *
 * date    |   author      |  changes
 * --------|---------------|------------------------------------------
 * 15.08.02|               |
 **********************************************************************/

#include "baumwelch.hh"

// project includes
#include "projectio.hh"
#include "statemodel.hh"

// standard C/C++ includes
#include <iomanip>

/*
 * BaumWelch::initialTypeParameters()
 * estimate the content models from the input data together with the type
 * information.
 */

void BaumWelch::initialTypeParameters(Double patpseudocount){
    
    int i, t;  
    Seq2Int s2i(k+1);
    
    /*
     * Initialize the data structures
     */
    if (!models) {
	models = new vector<ContentModel>(modelTypes);
	for (t=0; t < modelTypes; t++) {
	    (*models)[t].init(k, numFrames);
	}
    }
    if (!modelTypeProbs) {
	modelTypeProbs = new vector<Double>(modelTypes, 0.0);
    }
    Matrix< vector<Integer> > patCounts;
    patCounts.assign(modelTypes, numFrames, vector<Integer>(POWER4TOTHE(k+1)));
    vector<Integer> modelTypeCount(modelTypes);
    
    TrainingData *td = inputSeqs;
    while (td) {
	// determine type i
	if (td->type < 0) 
	    i=0;
	else 
	    i = td->type;
	if (i <0 || i>modelTypes )
	    cerr << "i out of range: " << i << endl;

	for (int a=0; a < td->seqLen - k; a++) {
	    try{
		patCounts[i][modm(td->frame+a+k, numFrames)][s2i(td->seq + a)] += 
		    td->weight;
	    } catch (InvalidNucleotideError e) {}
	}
	/*
	 * count the expected number of starts
	 */
	modelTypeCount[i] += td->weight;
	td = td->next;
    }

    for (i=0; i < modelTypes; i++) {
	for (int f=0; f < numFrames; f++) {
	    vector<Double> probs(POWER4TOTHE(k+1));
	    cout << "contentModel["<<i<<"], frame = " << f << endl;
	    StateModel::determineShortPatterns(patCounts[i][f], k, 200);
	    StateModel::makeProbsFromCounts(probs, patCounts[i][f], 
					    k, patpseudocount, true);
	    (*models)[i].setPatProb(f, probs);
	}
    }

    Double sum = 0.0;
    for (i=0; i < modelTypes; i++) {
	sum += modelTypeCount[i];
    }
    for (i=0; i < modelTypes; i++) {
	(*modelTypeProbs)[i] = modelTypeCount[i];
	if (sum > 0.0) 
	    (*modelTypeProbs)[i] /= sum;
    }
    cout << "modelTypeProbs[0]= "  <<(*modelTypeProbs)[0]  << ", modelTypeProbs[1] =" 
	 <<(*modelTypeProbs)[1] << endl;
}


/*
 * BaumWelch::initialRandomParameters()
 * generate random pattern distributions and start distribution
 */

void BaumWelch::initialRandomParameters(){
    int i, t;  
    /*
     * Initialize the data structures
     */
    if (!models) {
	models = new vector<ContentModel>(modelTypes);
	for (t=0; t < modelTypes; t++) {
	    (*models)[t].init(k, numFrames);
	}
    }
    if (!modelTypeProbs) {
	modelTypeProbs = new vector<Double>(modelTypes);
    }
     
     /*
      * Set the probabilities of the types
      */
    
     for (i=0; i < modelTypes; i++) {
	 (*modelTypeProbs)[i]= (double) 1 / modelTypes; // first they are equally likely
     }

     /*
      * Set the probabilities of the patterns
      */
     for (t=0; t < modelTypes; t++) {
	 (*models)[t].setRandomPatProbs();
     }
}

/*
 * void BaumWelch::initialExonTrainingData(const Gene* geneList)
 * makes a list of TrainingData out of the coding regions
 * 
 */

void BaumWelch::initialExonTrainingData(AnnoSequence *seqList){
    AnnoSequence *curSeq = seqList;
    int numData=0;

    numData = 0;
    while( curSeq ){
	Gene *curgene = NULL;
	if (curSeq->anno) 
	    curgene = curSeq->anno->genes;
	while(curgene) {
	    TrainingData *td = new TrainingData();
	    td->weight = curgene->weight; 
	    td->seqLen = curgene->clength - STOPCODON_LEN; // exclude the stop codon
	    td->seq = new char[td->seqLen+1];
	    td->frame = 0;                                 // begin with the start codon
	    td->type = 1;
	    td->name = newstrcpy(curgene->seqname);
	    State *exon = curgene->exons;
	    int curpos = 0;	     
	    while (exon) {
		//if (exon->contentType != 1) 
		//   td->type = exon->contentType;
		int cplength = exon->length();
		if (!exon->next)                // last exon
		    cplength -= STOPCODON_LEN;  // don't copy stop codon
		//assert (curpos + cplength <= td->seqLen);
		strncpy(td->seq + curpos, curSeq->sequence + exon->begin, cplength);
		curpos += cplength;
		exon = exon->next;
	    }
	    td->seq[td->seqLen]='\0';
	    addTrainingData(td);
	    curgene = curgene->next;
	}
	curSeq = curSeq->next;
    }
}


/*
 * void BaumWelch::initialIntronTrainingData(const Gene* geneList)
 * makes a list of TrainingData out of the introns
 * 
 */

void BaumWelch::initialIntronTrainingData(const Gene* geneList){
    const Gene *curgene = geneList;
    int numData;

    /* determine the number of genes containing introns
    while (curgene) {
	if (curgene->exons && curgene->exons->next)
	    numData++;
	curgene = curgene->next;
	}*/

    numData = 0;
    curgene = geneList;
    while( curgene ){
	if (curgene->exons && curgene->exons->next) {  // gene has introns
	    
	    TrainingData *td = new TrainingData();
	    // determine the length of all introns together
	    State *exon = curgene->exons;
	    td->seqLen = 0;
	    while (exon) {
		if (exon->next) {
		    td->seqLen += exon->next->begin - exon->end - 1;
		}
		exon = exon->next;
	    }
	    td->weight = curgene->weight; 
	    td->seq = new char[td->seqLen+1];
	    td->frame = 0;
	    exon = curgene->exons;
	    int curpos = 0;	     
	    while (exon) {
		if (exon->next) {
		    int cplength = exon->next->begin - exon->end - 1;
		    //assert (curpos + cplength <= td->seqLen);
		    // uncommented because member sequence doesn't exist anymore in class gene
		    //strncpy(td->seq + curpos, curgene->sequence + exon->end + 1, cplength);
		    curpos += cplength;
		}
		exon = exon->next;
	    }
	    td->seq[td->seqLen]='\0';
	    addTrainingData(td);
	}
	curgene = curgene->next;
    }
}

void  BaumWelch::classifyTrainingData(){
    int bestModel;
    Double bestProb, curProb;
    vector<int> classcount(modelTypes, 0);

    TrainingData *td = inputSeqs;
    while (td) {
	bestProb = 0.0;
	bestModel = -1;
	for (int i=0; i<modelTypes; i++) {
	    curProb = (*models)[i].seqProbUnderModel(td->seq, td->seqLen, td->frame);
	    if (curProb > bestProb) {
		bestProb = curProb;
		bestModel = i;
	    }
	}
	td->type = bestModel;
	classcount[bestModel]++;
	td = td->next;
    }
    cerr << "classification: " << classcount << endl;
}


/*
 * BaumWelch::reestimate()
 * make one step in the Baum-Welch parameter reestimation
 * return the likelihood of the observation under the old model
 * P(s_1) * P(s_2) * ... * P(s_n)
 */

Double BaumWelch::reestimate(Double patpseudocount){
    int i;
    Seq2Int s2i(k+1);
    Double odds = 1.0;    // probability of the observation
    vector<Double> p_s_given_m(modelTypes);
    vector<Double> p_m_given_s(modelTypes);
    Double p_s;           // probability of a certain sequence
    vector<Matrix<Double> > newpatprobs(modelTypes, Matrix<Double>(numFrames, POWER4TOTHE(k+1)));;
    vector<Double> newmodelTypeProbs;

    TrainingData *td = inputSeqs;
    while (td) {
	//cout << "---------- j=" << j << endl;
	// determine the probability of the j-th training data under each model i
	for (i=0; i<modelTypes; i++) {
	    p_s_given_m[i] = (*models)[i].seqProbUnderModel(td->seq, td->seqLen, td->frame);
	    //cout << "p_s_given_m[" << i << "]=" << p_s_given_m[i] << ", ";
	}
	//cout << endl;
	// determine the overall probability of the j-th training data
	p_s = 0.0;
	for (i=0; i<modelTypes; i++) {
	    p_s += (*modelTypeProbs)[i] * p_s_given_m[i];
	}
	//cout << "p_s = " << p_s << endl;
	if (!(p_s>0.0)) {
	    cerr << "Couldn't reestimate: Emission was impossible under model" << endl;
	    return 0.0;
	} 
	// determine the probability of model i under the j-th training data
	// using Bayes formula
	for (i=0; i<modelTypes; i++) {
	    p_m_given_s[i] = (*modelTypeProbs)[i] * p_s_given_m[i] / p_s;
	    //cout << p_m_given_s[i] << "  ";
	}
	//cout << endl;
	

	/*
	 * count the expected number of patterns using the p_m_given_s[i]
	 */	
	for (i=0; i<modelTypes; i++) {
	    for (int a=0; a < td->seqLen - k; a++) {
		try{
		    newpatprobs[i][modm(td->frame+a+k, numFrames)][s2i(td->seq + a)] += 
			td->weight * p_m_given_s[i];
		} catch (InvalidNucleotideError e) {}
	    }
	}

	/*
	 * count the expected number of starts using the p_m_given_s[i]
	 */
	for (i=0; i<modelTypes; i++) {
	    newmodelTypeProbs[i] += td->weight * p_m_given_s[i];
	}
	odds *= p_s;
	td = td->next;
    }

    /*
     * Now, add the pseudocounts one to each pattern in each frame in each model.
     */
    if (patpseudocount > 0.0) 
	for (i=0; i<modelTypes; i++) 
	    for (int pn=0; pn < newpatprobs[i].getRowSize(); pn++) 
		for (int f=0; f < numFrames; f++) 
		    newpatprobs[i][f][pn] += patpseudocount;

    
    // normalize
    for (i=0; i < modelTypes; i++) {
	for (int f=0; f < numFrames; f++) {
	    Double sum = 0.0;
	    for (int pn=0; pn < newpatprobs[i].getRowSize(); pn++) {
		sum += newpatprobs[i][f][pn];
	    }
	    for (int pn=0; pn < newpatprobs[i].getRowSize(); pn++) {
		newpatprobs[i][f][pn] = newpatprobs[i][f][pn]/sum;
		if (newpatprobs[i][f][pn]< 1e-9 && newpatprobs[i][f][pn] > 0) {
		    cerr << "i= " << i << " f=" << f << " p = " << s2i.inv(pn) 
			 << ", v= " << 
			newpatprobs[i][f][pn] << endl;
		}
	    }
	}
	(*models)[i].setPatProb(newpatprobs[i]);
    }
    Double sum = 0.0;
    for (i=0; i < modelTypes; i++) {
	sum += newmodelTypeProbs[i];
    }
    for (i=0; i < modelTypes; i++) {
	(*modelTypeProbs)[i] = newmodelTypeProbs[i] / sum;
    }

    // free memory
    // TODO
    // cout << "odds = " << odds << endl;
    return odds;
}

ContentModel BaumWelch::getContentModel(int type){
    if (type <0 || type >= modelTypes) {
	throw ProjectError("BaumWelch::getContentModel: unknown type");
    }
    return (*models)[type];
}
 
/*
 * BaumWelch::printThisContentModel(ostream &out)
 */ 
   
void BaumWelch::printAllContentModels(ostream &out) {
    printContentModels(out, models, modelTypeProbs, k);
}

/*
 * BaumWelch::printInputSeqs(ostream &out)
 */ 
void BaumWelch::printInputSeqs(ostream &out){
    out << setw(12) << "basecounts" << ", "
	<< setw(7) << "seqLen" << ", "
	<< setw(2) << "fr" << ", "
	<< setw(2) << "ty" << ", "
	<< setw(2) << "wg" << ", "
	<< "sequence" << endl << endl;
    BaseCount *bc;

    char *seq;
    TrainingData *td = inputSeqs;
    while (td) {
	seq = new char[td->seqLen+1];
	strncpy(seq, td->seq, td->seqLen);
	seq[td->seqLen]='\0';
	bc = new BaseCount(td->seq, td->seqLen);
	out << *bc << ", "
	    << setw(7) << td->seqLen << ", "
	    << setw(2) << td->frame << ", "
	    << setw(2) << td->type << ", "
	    << setw(2) << td->weight << ", "
	    << seq << endl;
	td = td->next;
	delete bc;
	delete seq;
    }
}

/*
 * readAllContentModels
 */

void BaumWelch::readAllContentModels(istream &in){

    in >> comment >> numFrames;
    in >> comment >> modelTypes;
    
    //cerr << "numFrames= " << numFrames << endl;
    //cerr << "modelTypes= " << modelTypes << endl;
    //cerr << "k=" << k << endl;

    // check whether the dimensions are correct
    
    if (modelTypeProbs && (*modelTypeProbs).size() != modelTypes){
	cerr << "BaumWelch::readAllContentModels: wrong dimension of modelTypeProbs\n";
	delete modelTypeProbs;
	//assert(modelTypeProbs==NULL);
    }
    if (modelTypeProbs == NULL) {
	modelTypeProbs = new vector<Double>(modelTypes);	
    }
    
    if (models && (*models).size() != modelTypes){
	cerr << "BaumWelch::readAllContentModels: wrong dimension of models\n";
	delete models;
	//assert(models==NULL);
    }
    if (models == NULL) {
    	models = new vector<ContentModel>(modelTypes);	
	for (int t=0; t < modelTypes; t++) {
	    (*models)[t].init(k, numFrames);
	}
    }
    
    // read modelTypeProbs
    in >> comment;
    for (int t=0; t < modelTypes; t++) {
	in >> (*modelTypeProbs)[t];
    }
    
    // read content models
    in >> comment;   
    Double dbl;
    Seq2Int s2i_dummy(k+1);
    for (int i=0; i < (*models)[0].getNumPatterns(); i++) {
	in >> comment;
	if (i != s2i_dummy.read(in))
	    cerr << "BaumWelch::readAllContentModels: wrong order of patterns\n";
	for (int t=0; t < modelTypes; t++) {
	    for (int f=0; f < numFrames; f++) {
		in >> dbl;
		(*models)[t].getPatProb(f, i) = dbl / 1000;
		
	    }
	}
    } 
}


/*
 * printContentModel(...)
 */

void printContentModels(ostream& out, vector<ContentModel> *models,
		       vector<Double> *modelTypeProbs, int k){

    Seq2Int s2i(k+1);
    
    if (modelTypeProbs->size()==0 || models->size()==0) {
	out << "Could not write content model probabilities. No data." << endl;
	return;
    }
    
    int numTypes = modelTypeProbs->size();
    int numFrames = (*models)[0].getNumFrames();

    out << "# number of reading frames\n" << numFrames << endl;
    out << "# number of types\n" << modelTypeProbs->size() << endl;
    out << "# probabilities of the types" << endl;
    cout << setw(k+1) << " ";
    for (int t=0; t < numTypes; t++) {
	out << setw(12*numFrames+2) << (*modelTypeProbs)[t];
    }
    out << endl;
    out << "# 1000 * probabilities of the patterns" << endl;
    for (int i=0; i < (*models)[0].getNumPatterns(); i++) {
	out << s2i.inv(i) << "  ";
	for (int t=0; t < numTypes; t++) {
	    for (int f=0; f < numFrames; f++) {
		out << setw(9) << 1000 * (*models)[t].getPatProb(f, i);
	    }
	    if (t < numTypes) {
		out << "  ";
	    }
	}
	out << endl;
    }
    out << "\n# mean (a,c,g,t)-contents of the models:" << endl;
    out << "#" << setw(k) << " ";
    for (int t=0; t < numTypes; t++) 
	out << "    " << (*models)[t].getMeanContent();
    out << endl;
}
