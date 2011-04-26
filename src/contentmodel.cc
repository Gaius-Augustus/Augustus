/**********************************************************************
 * file:    contentmodel.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  Model typical sequences by pattern frequencies
 * authors: Mario Stanke, mario@gobics.de
 *
 * date    |   author      |  changes
 * --------|---------------|------------------------------------------
 * 19.08.02| Mario Stanke  |  creation of the file
 **********************************************************************/

#include "contentmodel.hh"

/*
 * sort an array increasingly p from index a to index b (included).
 */
void QuickSort(Double *p, int a, int b);

/*
 * Double seqProbUnderModel(char *s, int len, int frame);
 * return the probability of sequence s of lenth len under this model
 * assume seq starts with reading frame f
 */
Double ContentModel::seqProbUnderModel(const char *seq, int len, int frame) const {
    if (len <= k) {
	return 1.0;
    }
    int f;
    Double p, pShortPattern;
    Seq2Int s2i(k+1), s2is(k);

    // begin pattern
    f = modm(frame + k, getNumFrames()); 
    try {
	p = patprob[f][s2i(seq + k - k)];
    } catch (InvalidNucleotideError) {
	p = POWER4TOTHE(k+1);
    }
    for (int i = k+1; i < len; i++) {
	int fr = modm(frame + i, getNumFrames());
	try {
	    pShortPattern = 
		patprob[fr][s2is(seq + i-1 - (k-1)) * 4]
		+ patprob[fr][s2is(seq + i-1 - (k-1)) * 4 + 1]
		+ patprob[fr][s2is(seq + i-1 - (k-1)) * 4 + 2]
		+ patprob[fr][s2is(seq + i-1 - (k-1)) * 4 + 3];
	} catch (InvalidNucleotideError e) {
	    pShortPattern = POWER4TOTHE(k+1);
	}
	if (!(pShortPattern>0.0)) {
	    return 0.0;
	}
	try {
	    p *= patprob[fr][s2i(seq + i - k)] / pShortPattern;
	} catch (InvalidNucleotideError e) {
	    p = .25;
	}
    }
    return p;
}
/*
 * char *ContentModel::generateRandSeq(int n)
 * generate a sequence of length n>k under this probabilistic content model
 * generate start pattern according to pattern probs
 * and then nucleotide by nucleotide according to conditional probabilities
 * start with frame 'startFrame'
 */

string ContentModel::generateRandSeq(int n, int startFrame) const {
    if (n <= k) {
	return NULL;
    }
    int i; 
    Seq2Int s2i(k+1);

    // the k+1 first nucleotides at once
    int patternindex = randomIndex(patprob[modm(startFrame+k, getNumFrames())], patprob.getRowSize());
    string erg = s2i.inv(patternindex);
    // the rest one by one
    Double emiprobs[4];
    char nucs[4]={'a', 'c', 'g' ,'t'};
    int nucIdx;
    for (i=k+1; i<n; i++) {
	for (int j=0; j<4; j++) {
	    erg += nucs[j];
	    emiprobs[j] = patprob[modm(startFrame + i, getNumFrames())][s2i(erg.c_str()+i-k)]; 
	}
	nucIdx = randomIndex(emiprobs, 4);
	erg[i] = nucs[nucIdx];
    }
    // TODO: tolower(erg);
    return erg;
}

/*
 * int ContentModel::randomIndex(Double *p, int size)
 * normalize the vector p to a probability vector and 
 * return a random index from this distribution on {0,1, ..., size}
 */

int ContentModel::randomIndex(const Double* p, int count){
    Double sum=0.0;
    for (int i=0; i < count; i++) {
	sum += p[i];
    }
    if (!(sum>0.0)) {
	throw ProjectError("ContentModel::randomIndex: probability vector was 0");
    }
    Double r = sum * ((double) rand()/RAND_MAX);
    sum = 0.0;
    int i = 0;
    while (!(sum > r) && i < count){
	sum += p[i++];
    }
    return i-1;
}

/*
 * Double getMeanContent()
 * returns the mean  content of (a,c,g,t) of this model
 */


BaseCount ContentModel::getMeanContent(){
    Seq2Int s2i(k+1);
    BaseCount bc(0,0,0,0);
    bc.ra = bc.rc = bc.rg = bc.rt = 0.0;
    for (int f=0; f<getNumFrames(); f++) {
	for(int pn=0; pn < getNumPatterns(); pn++) {
	    string s = s2i.inv(pn);
	    for (int c=0; c < k+1; c++) {
		switch (s[c]) {
		    case 'A': case 'a':
			bc.ra += patprob[f][pn].doubleValue();
			break;	
		    case 'C': case 'c':
			bc.rc += patprob[f][pn].doubleValue();
			break;
		    case 'G': case 'g':
			bc.rg += patprob[f][pn].doubleValue();
			break;
		    case 'T': case 't':
			bc.rt += patprob[f][pn].doubleValue();
			break;
		    default:
			throw ProjectError("BaumWelch::getMeanContent:"
					       "unknown nucleotide");
		}
	    }
	}
    }
    bc.ra /= ((double) getNumFrames()*(k+1));
    bc.rc /= ((double) getNumFrames()*(k+1));
    bc.rg /= ((double) getNumFrames()*(k+1));
    bc.rt /= ((double) getNumFrames()*(k+1));
      
    return bc;
}

void ContentModel::setRandomPatProbs(){
    int f,i;
    for (f=0; f < getNumFrames(); f++) {
	Double sum = 0.0;
	for (i=0; i < getNumPatterns(); i++) {
	    patprob[f][i] = 30 + (double) 10 * rand() / RAND_MAX;
	    sum += patprob[f][i];
	}
	for (i=0; i < getNumPatterns(); i++) {
	    patprob[f][i] /= sum;
	}
    }
}

void ContentModel::setPatProb(const Matrix<Double>& newpatprob) { 
    if(newpatprob.getColSize() != getNumFrames() ||
       newpatprob.getRowSize() != getNumPatterns()) {
	throw ProjectError("ContentModel::setPatProb: wrong dimensions");
    }
    patprob = newpatprob;
}

void ContentModel::setPatProb(int frame, const vector<Double>& newRow) {
    if(frame<0 || frame>=getNumFrames() ||
       newRow.size() != getNumPatterns()) 
	throw ProjectError("ContentModel::setPatProb: wrong dimensions");
    for (int pn=0; pn < getNumPatterns(); pn++)
	patprob[frame][pn] = newRow[pn];
}


/*
 * ContentModel::print()
 */

void ContentModel::print(){
    Seq2Int s2i(k+1);
    
    cout << "# number of reading frames\n" << getNumFrames() << endl;
    cout << "# k\n" << k << endl;
    
    for (int i=0; i < getNumPatterns(); i++) {
	cout << s2i.inv(i) << "  ";
	for (int f=0; f < getNumFrames(); f++) {
	    cout << setw(12) << patprob[f][i];
	} 
	cout << endl;
    }
    cout << endl;

    cout << "\n# mean (a,c,g,t)-content:" << endl;
    BaseCount bc = getMeanContent();
    cout << "#   " << bc << endl;
}

/*
 * Double ContentModel::DiscrimDifference(ContentModel *a, ContentModel *b, int n);
 * Vaguely: The DiscrimDifference is the Probability of determining the right model
 * when given a sequence of length n under one of the two models.
 * We assume that the decision is based on a Neyman-Pearson-Test:
 * For some threshold c we s
 * say 'a' iff P_b(seq)/P_a(seq) > c
 * It is max_c min{ P_a( P_b(seq)/P_a(seq))<=c, P_b( P_b(seq)/P_a(seq)>c)}
 */

double ContentModel::DiscrimDifference(const ContentModel *a, const ContentModel *b, 
				       int n, int fa, int fb){
    /*
     *
     */
    
    const int numRepeats = 1000; // number of repeats of the simulation
    Double *sa, *sb;             // scores of samples of model a and b, respectively
    Double pa, pb;
    /*
     * score = P_b(seq) / P_a(seq)
     */
    sa = new Double[numRepeats]; 
    sb = new Double[numRepeats];
    
    for (int i=0; i < numRepeats; i++) {
	// emission of model a
	string emis = a->generateRandSeq(n, fa);
	pa = a->seqProbUnderModel(emis.c_str(), n, fa);
	pb = b->seqProbUnderModel(emis.c_str(), n, fb);
	if (pa>0.0) {
	    sa[i] = pb/pa;
	} else {
	    sa[i] = numeric_limits<double>::infinity();
	}
	// emission of model b
	emis = b->generateRandSeq(n, fb);
	pa = a->seqProbUnderModel(emis.c_str(), n, fa);
	pb = b->seqProbUnderModel(emis.c_str(), n, fb);
	if (pa>0.0) {
	    sb[i] = pb/pa;
	} else {
	    sb[i] = Double::getMaxDouble();
	}
    }
 
    /*
     * Now determine a cutoff value c for the score such that
     * |#{i | sa[i]>c} - #{i | sb[i]<=c}| <= 1 
     * estimate the DiscrimDifference by 
     * min{ #{i | sa[i]<=c}, #{i | sb[i]>c}} / numRepeats
     */
    
    QuickSort(sa, 0, numRepeats-1);
    QuickSort(sb, 0, numRepeats-1);
  
    
/*
    cerr <<" sa sorted:" << endl;
    for (int i=0; i<numRepeats; i++) {
	cerr << sa[i] << ", ";
    }
    cerr << endl;

    
    cerr <<"sb sorted:" << endl;
    for (int i=0; i<numRepeats; i++) {
	cerr << sb[i] << ", ";
    }
    cerr << endl;
*/
    int k = numRepeats-1;
    while (k>0 && (sa[k] > sb[numRepeats-1 - k])) {
	k--;
    }
    return (double) (k+1) / numRepeats;
}


void QuickSort(Double *p, int a, int b) {
    if (a >= b)
	return;
    Double pivot = p[a], temp;
    bool partitioned = false;
    int i=a-1, j=b+1;
    while(!partitioned) {
	do {j--;} while(p[j] > pivot);
	do {i++;} while(p[i] < pivot);
	if (i<j) {      // exchange p[i] and p[j]
	    temp = p[i];
	    p[i] = p[j];
	    p[j] = temp;
	} else 
	    partitioned = true;
    }
    QuickSort(p, a, j);
    QuickSort(p, j+1, b);
}

