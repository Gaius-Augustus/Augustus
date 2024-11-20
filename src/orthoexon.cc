/*
 * orthoexon.cc
 *
 * License: Artistic License, see file LICENSE.TXT or 
 *          https://opensource.org/licenses/artistic-license-1.0
 * 
 * Description: maintains orthologous exons for comparative gene prediction
 */

#include "orthoexon.hh"
#include "graph.hh"

#include <fstream>
#include <iostream>
#include <cmath>

const char* phyleticPatternIdentifiers[6]={"0", "1", "-", "_", "g", "l"};

OrthoExon::OrthoExon(int_fast64_t k, size_t n)
    : key(k), omega(-1.0), Eomega(-1.0), VarOmega(-1.0), leftBoundaryExtOmega(-1.0),
      rightBoundaryExtOmega(-1.0), leftBoundaryIntOmega(-1.0), rightBoundaryIntOmega(-1.0),
      probClamsa(-1.0), intervalCount(0), intervalCountClamsa(0), subst(-1), cons(-1.0),
      leftCons(-1.0), rightCons(-1.0), diversity(-1.0), leftBoundaryEbony(-1.0), rightBoundaryEbony(-1.0) {
    orthoex.resize(n);
    orthonode.resize(n);
    weights.resize(n,0);
    labels.resize(n,3); // initialize all species as "unaligned" (label 3)
}

double OrthoExon::getLeftBoundaryEbony() const {
    return leftBoundaryEbony; // Return the value of leftBoundaryEbony
}
double OrthoExon::getRightBoundaryEbony() const {
    return rightBoundaryEbony; // Return the value of rightBoundaryEbony
}

StateType OrthoExon::getStateType() const{
    for (int s = 0; s < orthoex.size(); s++)
	if (orthoex[s])
	    return orthoex[s]->getStateType();
    return TYPE_UNKNOWN;
}

int OrthoExon::numExons() const{
    int k=0;
    for (size_t s = 0; s < orthoex.size(); s++)
	if (orthoex[s])
	    k++;
    return k;
}

bool OrthoExon::exonExists(int pos) const{
    if(pos >= orthoex.size() || pos < 0)
	throw ProjectError("Internal error in OrthoExon::exonExists: out of bound access");
    if(orthoex[pos])
	return true;
    return false;	    
}

/*
 * setPresent() and setAbsent()
 * initialization of 'labels'
 */

void OrthoExon::setPresent(bit_vector v){
    for(int i=0; i<v.size(); i++){
        if(v[i] == 1)
            labels[i] = 0;
    }
}

void OrthoExon::setAbsent(bit_vector v){
    for(int i=0; i<v.size(); i++){
        if(v[i] == 1)
            labels[i] = 2;
    }
}

/*
 * phyletic patterns
 *-----------------------
 * label of Species i in {0,1,2,3,4,5} equiv. {"0","1","-","_","g","l"}
 * 
 * 0 - EC present in i, but not predicted as exon
 * 1 - EC present and predicted as exon in i
 * 2 - EC absent in i, but alignment present
 * 3 - alignment not present in i
 * 4 - EC present in i and "gained" through dual decomposition
 * 5 - EC present in i and "lost" through dual decomposition
 */
string OrthoExon::getPhyleticPattern() const {
    string phyletic_pattern = "";
    for(size_t pos = 0; pos < labels.size(); pos++){
	phyletic_pattern+=phyleticPatternIdentifiers[labels[pos]];
    }
    return phyletic_pattern;
}

void OrthoExon::setPhyleticPattern(map<int, list<int> > &pp_init, map<int, list<int> > &pp_opt){
            
    int numSpecies = labels.size();
    vector<int> labels_opt(numSpecies,0);
    vector<int> labels_init(numSpecies,0);

    map<int, list<int> >::iterator it = pp_opt.find(ID);
    if(it != pp_opt.end()){
	for(list<int>::iterator lit = it->second.begin(); lit != it->second.end(); lit++)
	    labels_opt[*lit]=1;
    }
    it = pp_init.find(ID);
    if(it != pp_init.end()){
	for(list<int>::iterator lit = it->second.begin(); lit != it->second.end(); lit++)
	    labels_init[*lit]=1;
    }
    for(size_t pos = 0; pos < labels.size(); pos++){
	if(exonExists(pos)){
	    if(labels_init[pos] == 0 && labels_opt[pos] == 1)
		labels[pos]=4; // exon gain through dual decomposition
	    else if(labels_init[pos] == 1 && labels_opt[pos] == 0)
		labels[pos]=5; // exon lost through dual decomposition
	    else{ // labels_init[pos] == labels_opt[pos]
		labels[pos]=labels_opt[pos]; 
	    }
	}
    }
}

void OrthoExon::setTree(PhyloTree* t) {
    tree = t;
    setDiversity(tree->getDiversity());
}
/*
vector<int> OrthoExon::getRFC(vector<int> offsets) const{
    vector<int> rfc;
    for (size_t s = 0; s < orthoex.size(); s++){
	if (orthoex[s] == NULL)
	    rfc.push_back(-1);
	else
	    rfc.push_back((offsets[s] + orthoex[s]->getFirstCodingBase()) % 3);
    }
    return rfc;
}
*/
void OrthoExon::setOmega(vector<double>* llo, CodonEvo* codonevo , bool oeStart){
    if(oeStart){
      loglikOmegaStarts.push_back(*llo);
      //cout<<"set Omega at oeStart: "<<getAliStart()<<":"<<getAliEnd()<<":"<<getStateType()<<"\t(omega, omega squared, count) = "<<"("<<omega<<", "<<omegaSquared<<", "<<omegaCount<<")"<<endl;
    }else{
      
      vector<double> loglikOmegas;
      if(!loglikOmegaStarts.empty()){
	loglikOmegas = loglikOmegaStarts.front();
	loglikOmegaStarts.pop_front();
      }else{
	if(!llo->empty()){
	  loglikOmegas = *llo;
	}
      }
      
        //calculate posterior mean of omega
        int k = llo->size();

	if(*llo == loglikOmegas)
	  k = 0;

	double currOmega;
	double currVarOmega;

	if(k == 0){ // no likelihood was calculated
	  currOmega = -1;
	  currVarOmega = -1;
	  omega = -1;
	  storeOmega(currOmega, currVarOmega);
	  //cerr<<"ortho exon "<<this->ID<<" has no omega"<<endl;
	  return;
	}
	//cout<<"number of omegas: "<<k<<endl;
	if(loglikOmegas.size() != k)
	    loglikOmegas.resize(k,0.0);
	//cout<<"log likelihood of omega"<<endl;
	for (int u=0; u < k; u++){
	  // sum of likelihoods in intervall OE start and OE end
	  loglikOmegas[u] = (*llo)[u] - loglikOmegas[u];
	  //	  cout<<"omega "<<codonevo->getOmega(u)<<"\t"<<loglikOmegas[u]<<endl;
	}
	double maxloglik = *max_element(loglikOmegas.begin(), loglikOmegas.end());
	double sum = 0;
	vector<double> postprobs(k, 0.0);

	for (int u=0; u < k; u++){
	  sum += postprobs[u] = exp(loglikOmegas[u] - maxloglik) * codonevo->getPrior(u);
	}
	
	//cout << "posterior distribution and prior of omega" << endl;
	for (int u=0; u < k; u++){
	  //	  cout << codonevo->getOmega(u) << "\t" << postprobs[u] <<"\t"<<codonevo->getPrior(u)<< endl;
	  postprobs[u] /= sum;
	}
	currOmega = 0;

        //cout<<"---------------------------------------------------------------------------"<<endl;
	//cout<<"wi\t\tloglikOmegas\tmaxloglik\tpostprobs/sum\tprior\tsum\texp(loglik - maxloglik)"<<endl;
 
        for (int u=0; u < k; u++){
	  currOmega += postprobs[u] * codonevo->getOmega(u);
	  //cout<<codonevo->getOmega(u)<<"\t\t"<<loglikOmegas[u]<<"\t\t"<<maxloglik<<"\t\t"<<postprobs[u]<<"\t\t"<<codonevo->getPrior(u)<<"\t\t"<<sum<<"\t\t"<<exp(loglikOmegas[u] - maxloglik)<<endl;                                                     
        }


	//cout<<"curOmega: "<<currOmega<<endl;
	//cout<<"set Omega at oeEnd: "<<getAliStart()<<":"<<getAliEnd()<<":"<<getStateType()<<"\t(omega, omega squared, count) = "<<"("<<omega<<", "<<omegaSquared<<", "<<omegaCount<<")"<<" Eomega: "<<currOmega<<endl;
	double omega_maxML = 0;
	currVarOmega = 0.0;
	for (int u=0; u < k; u++){
	    currVarOmega += postprobs[u] * pow(codonevo->getOmega(u) - currOmega, 2);
	    if(loglikOmegas[u] == maxloglik){
	      omega_maxML = codonevo->getOmega(u);
	    }
	}
	//cout<<"curVarOmega: "<<currVarOmega<<endl;
	if(omega_maxML > 0)
	  omega = omega_maxML;

	storeOmega(currOmega, currVarOmega);
	loglikOmegas.clear();
    }
}

void OrthoExon::storeOmega(double currOmega, double currVarOmega){
  
  switch(intervalCount){
  case 0: leftBoundaryExtOmega = currOmega;
    break;
  case 1: leftBoundaryIntOmega = currOmega;
    break;
  case 2: Eomega = currOmega; VarOmega = currVarOmega;
    break;
  case 3:rightBoundaryIntOmega = currOmega;
    break;
  case 4: rightBoundaryExtOmega = currOmega;
    break;
  default: throw ProjectError("Error in setOmega(): too many intervals were calculated.");
    break;
  }
  intervalCount++;
}


// clamsa related code
void OrthoExon::setClamsa(const cumValues &cv, CodonEvoDiscr* codonevodiscr , bool oeStart){
    const vector<double>* llo = &cv.logliks;
    
    if (oeStart){
    	loglikClamsaStarts.push_back(*llo);
    } else {
        // at right boundary finish aggregating computations
      	vector<double> loglikClamsa;

      	if (!loglikClamsaStarts.empty()){
            loglikClamsa = loglikClamsaStarts.front();
            loglikClamsaStarts.pop_front();
      	} else {
            if (!llo->empty()){
                loglikClamsa = cv.logliks;
            }
    	}
		
        // calculate mean of likelihood over all codons
        int k = llo->size();

        if (*llo == loglikClamsa)
            k = 0;

        double clamsaProb;

        if (k == 0){ // no likelihood was calculated
            clamsaProb = -1;
            // omega = -1; no equivalent for clamsa 
            storeClamsa(clamsaProb);
            //cerr<<"ortho exon "<<this->ID<<" has no omega"<<endl;
            return;
        }
        cout << "number of rate matrices: " << k << endl;
        if (loglikClamsa.size() != k)
            loglikClamsa.resize(k, 0.0);
        cout << "log likelihoods for each rate matrix:" << endl;
        for (int u=0; u < k; u++){
            // sum of likelihoods in intervall OE start and OE end
            loglikClamsa[u] = (*llo)[u] - loglikClamsa[u];
            cout << u << "\t" << loglikClamsa[u] << endl;
        }

        cout << "OrthoExon::setClamsa, numSites=" << cv.numSites << endl;
        // assign the mean of loglikClamsa computed for each model in clamsa
        for (int u=0; u < k; u++){
            loglikClamsa[u] /= cv.numSites;
            cout << "loglikClamsa[" << u << "]=\t" << loglikClamsa[u] << "\t";
        }
        cout << endl;
        cout << "cumValues alignment:" << endl;
        for (int i=0; i < cv.rows.size(); i++)
            cout << cv.rows[i] << endl;
        cout << "len ID=" << ID << "\t" << cv.rows[0].size() << endl;
        clamsaProb = codonevodiscr->getProb(loglikClamsa);
        storeClamsa(clamsaProb);
        // cout << "---- stored clamsa for OE" << ID << "\t" << getClamsaProb() << endl;
        loglikClamsa.clear();
    }
}


// clamsa related code
void OrthoExon::setClamsa2(const cumValues &cv, CodonEvoDiscr* codonevodiscr){
    const vector<double> &llo = cv.logliks;		
    // calculate mean of likelihood over all codons
    int k = llo.size();
    double clamsaProb;

    if (k == 0){ // no likelihood was calculated
        clamsaProb = -1;
    } else {
        vector<double> loglikClamsa(k, 0.0);
        // cout << "OrthoExon::setClamsa2, numSites=" << cv.numSites << endl;
        for (int u=0; u < k; u++){
            loglikClamsa[u] = llo[u] / cv.numSites;       
            cout << u << "\t" << loglikClamsa[u] << endl;
        }
        clamsaProb = codonevodiscr->getProb(loglikClamsa);
    }
    storeClamsa(clamsaProb);
    // cout << "---- stored clamsa for OE" << ID << "\t" << getClamsaProb() << endl;
}


void OrthoExon::storeClamsa(double currClamsa){
    probClamsa = currClamsa;
}


void OrthoExon::setSubst(int subs, bool oeStart){
  if(oeStart){
    if(intervalCount == 1)
	subst = subs;
  }else{
    if(subs != -1 && intervalCount == 3){
      if(subst >= 0)
	subst = subs - subst;
      else
	subst = subs - subst - 1;
    }
    else{
      if(subst >= 0 && intervalCount == 3)
	throw ProjectError("Error in setSubs(): numSubs was defined at OE start but is not at OE end!");
    }
  }
}


double logit(double p) {
    if (p == -1.0) { // Special case: return 0 for p == -1.0
        return 0;
    }
    if (p == 1.0) {
        return std::log(0.999 / (1 - 0.999));
    } 
    if (p == 0.0) {
        return std::log(0.001 / (1 - 0.001));
    }
    if (p < 0.0 || p > 1.0) {
        throw std::domain_error("Probability p must be in the range [0, 1].");
    }
    return std::log(p / (1 - p));
}

double OrthoExon::getLogRegScore() const{
    // pre-definitions for the boundary feature
    double b_l;
    double b_r;
    if (getLeftIntOmega() > 0 && getLeftExtOmega() > 0 && getRightIntOmega() > 0 && getRightExtOmega() > 0){
        b_l = log(getLeftIntOmega()) - log(getLeftExtOmega());
        b_r = log(getRightIntOmega()) - log(getRightExtOmega());
    } else {
        b_l = 0;
        b_r = 0;
    }

    // GM added the following to include up to 21 features
    double addScore = 0.0;
#if TESTING
    addScore = 
        Constant::ex_sc[31]*hasOmega()*getLeftExtOmega()	
        //'leftBoundaryExtOmega', 	14
        + Constant::ex_sc[32]*hasOmega()*getRightExtOmega()	
        //'rightBoundaryExtOmega'		15
        + Constant::ex_sc[33]*hasOmega()*getLeftIntOmega()	
        // 'leftBoundaryIntOmega'		16
        + Constant::ex_sc[34]*hasOmega()*getRightIntOmega()	
        // 'rightBoundaryIntOmega'		17
        + Constant::ex_sc[35]*hasConservation()*getLeftConsScore()
        //	'leftBoundaryCons'			18
        + Constant::ex_sc[36]*hasConservation()*getRightConsScore();
    //	'rightBoundaryCons'			19
#endif
    
    double clamsaScore = getClamsaScore();
    double score;
    std::cout << "addScore: " << addScore << std::endl;

    std::cout << "clamsaScore: " << clamsaScore << std::endl;

    std::cout << "Eomega: " << Eomega << std::endl;

    std::cout << "VarOmega: " << VarOmega << std::endl;

    std::cout << "cons: " << cons << std::endl;

    std::cout << "containment: " << containment << std::endl;

    std::cout << "diversity: " << diversity << std::endl;

    std::cout << "numExons: " << numExons() << std::endl;

    std::cout << "orthoex.size(): " << orthoex.size() << std::endl;

    std::cout << "b_l: " << b_l << std::endl;

    std::cout << "b_r: " << b_r << std::endl;

    std::cout << "hasOmega: " << hasOmega() << std::endl;

    std::cout << "hasVarOmega: " << hasVarOmega() << std::endl;

    std::cout << "hasConservation: " << hasConservation() << std::endl;

    std::cout << "hasContainment: " << hasContainment() << std::endl;

    std::cout << "hasDiversity: " << hasDiversity() << std::endl;
    std::cout << "leftBoundaryEbony: " << leftBoundaryEbony << std::endl;
    std::cout << "rightBoundaryEbony: " << rightBoundaryEbony << std::endl;
    std::cout << "; oeScore: " << std::endl;

    int hasLeftEbony = (leftBoundaryEbony >= 0);                                                                                                                                       
    int hasRightEbony = (rightBoundaryEbony >= 0); 
    double leftEbonyScore = logit(leftBoundaryEbony);
    double rightEbonyScore = logit(rightBoundaryEbony);
    std::cout << "leftEbonyScore: " << leftEbonyScore << std::endl;
    std::cout << "rightEbonyScore: " << rightEbonyScore << std::endl;
    
    score = addScore + Constant::ex_sc[6]  * Eomega * hasOmega()
             + Constant::ex_sc[7]  * VarOmega * hasVarOmega()
             + Constant::ex_sc[8]  * cons * hasConservation()
             + Constant::ex_sc[9]  * containment * hasContainment()
             + Constant::ex_sc[10] * diversity * hasDiversity()
             + Constant::ex_sc[11] * numExons()
             + Constant::ex_sc[15] * numExons() / orthoex.size()
             + Constant::ex_sc[13] * cons * diversity * hasConservation() * hasDiversity() 
             + Constant::ex_sc[14] * Eomega * hasOmega() * diversity * hasDiversity()
             - Constant::ex_sc[1]  * hasOmega()
             + Constant::ex_sc[16] * ( b_l * exp(Constant::lambda*b_l) + b_r * exp(Constant::lambda*b_r) ) / ( exp(Constant::lambda*b_l) + exp(Constant::lambda*b_r) )
             - Constant::ex_sc[2] // for being a HECT
             + Constant::ex_sc[17] * clamsaScore;

#ifdef EBONY
    if (Constant::ebonyScores) { 
    score = score
             + Constant::ex_sc[19] * score  // EBONY-specific score theta-one multiplying by oldscore                                                                                                   
             + Constant::ex_sc[20] * hasRightEbony // EBONY-specific score hasRightEbony                                                                                                                
             + Constant::ex_sc[21] * rightEbonyScore// EBONY-specific score rightEbonyScore                                                                                          
             + Constant::ex_sc[22] * hasLeftEbony// EBONY-specific score hasLeftEbony                                                                                                                   
             + Constant::ex_sc[23] * leftEbonyScore; // EBONY-specific score lefttEbonyScore  
    }
#endif
    return score;         
}

ostream& operator<<(ostream& ostrm, const OrthoExon &oe){
    ostrm << stateTypeIdentifiers[oe.getStateType()];
    for (int s = 0; s < oe.orthoex.size(); s++){
	ExonCandidate *ec = oe.orthoex[s];
	if (ec)
	    ostrm << "\t" << ec->begin+1 << "\t" << ec->end - ec->begin + 1;
	else
	    ostrm << "\t" << 0 << "\t" << 0 << "\t";
    }
    return ostrm;
}

istream& operator>>(istream& istrm, OrthoExon& ex_tuple){

    string exontype;
    long int begin, length;

    istrm >> exontype;
    for (int i = 0; i < ex_tuple.orthoex.size(); i++){
	istrm >> begin >> length;
	if (begin != 0 && length != 0){
	    ExonCandidate *exoncand = new ExonCandidate(toExonType(exontype.c_str()), begin-1, begin+length-2);
	    ex_tuple.orthoex[i] = exoncand;
	}
    }
    return istrm;
}
