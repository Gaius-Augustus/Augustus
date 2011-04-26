/*****************************************************************************\
 * Filename : ContentModel.hh
 * Author   : Mario Stanke
 * Project  : HMM
 * Version  : 0.2
 *
 *
 * Description: Model typical sequences by pattern frequencies
 *
 *
 * Date                |   Author              |  Changes
 *---------------------|-----------------------|---------------------------------
 * Die Jun 25 2002     | Mario Stanke          | creation of the file
\******************************************************************************/

#ifndef __CONTENTMODEL_HH
#define __CONTENTMODEL_HH

// project includes
#include "types.hh"
#include "motif.hh"

// standard C/C++ includes
#include <string>
#include <iostream>
#include <iomanip>

using namespace std;

class ContentModel {
public:
    /*
     * default contructor
     */
    ContentModel(){}
    
    /*
     * contructor
     */
    ContentModel(int k, int f = 1){
	init(k,f);
    };

     void init(int k, int f){
	this->k = k;
	patprob.assign(f, POWER4TOTHE(k+1));
    }
 
    void setRandomPatProbs();
    string generateRandSeq(int n, int startFrame = 0) const;
    static int randomIndex(const Double* p, int count);
    Double seqProbUnderModel(const char *seq, int len, int frame) const;
    BaseCount getMeanContent();
    Double getPatProb(int f, int pn) const{ return patprob[f][pn]; }
//     const Double* getPatProb(int f) const {return (*patprob)[f];}
    vector<Double> copyPatProb(int f) const {return patprob.getRow(f);}
    void setPatProb(const Matrix<Double>& newpatprob);
    void setPatProb(int frame, const vector<Double>& newRow);
//     void setPatProbClone(const Matrix<Double>& newpatprob);
    int getNumFrames() const { return patprob.getColSize(); }
    int getNumPatterns() const { return patprob.getRowSize(); }
    void print();
    static double DiscrimDifference(const ContentModel *a, const ContentModel *b, 
				    int n, int fa=0, int fb=0);

private:
    int k;
    Matrix<Double> patprob;
};



#endif    //__CONTENTMODEL_HH
