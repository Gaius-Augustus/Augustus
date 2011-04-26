/*****************************************************************************\
 * Filename : baumwelch.hh
 * Author   : Mario Stanke
 * Project  : HMM
 * Version  : 0.1
 *
 *
 * Description: Implementation of a special case of Baum-Welch parameter estimation
 *
 *
 * Date       |   Author              |  Changes
 *------------|-----------------------|------------------------------------------
 * 15.08.2002 | Stanke                | creation of the class
\******************************************************************************/

#ifndef _BAUMWELCH_HH
#define _BAUMWELCH_HH

// project includes
#include "types.hh"
#include "gene.hh"
#include "contentmodel.hh"

// standard C/C++ includes
#include <iostream>


/*
 * class TrainingData
 * a sequence together with some information
 */
class TrainingData {
public:
    TrainingData(){
	type = -1;            // type usually not known
	next = NULL;
	name = NULL;
    }

    char *seq;
    int seqLen;
    int frame;   // reading Frame of the first nucleotide
    int type;
    int weight;
    char *name;
    TrainingData *next;
};


class BaumWelch {
public:

    BaumWelch(int k, int numFrames, int modelTypes = 2){
	this->k = k;
	this->numFrames = numFrames;
	this->modelTypes = modelTypes;
	models = NULL;
	modelTypeProbs = NULL;
	inputSeqs = NULL;

    }
    void initialTypeParameters(Double patpseudocount = 0.0);
    void initialRandomParameters();
    void initialExonTrainingData(AnnoSequence *seqList);
    void initialIntronTrainingData(const Gene* geneList); 
    void classifyTrainingData();  
    void setK(int k){
	if (k != this->k) {
	    delete models;
	    models = NULL;
	    delete modelTypeProbs;
	    modelTypeProbs = NULL;
	    this->k = k;
	}
    };
    void setTrainingData(TrainingData *inputSeqs){
	if (this->inputSeqs) {
	    // TODO
	}
	this->inputSeqs = inputSeqs;
    }
    Double reestimate(Double patpseudocount = 0.0);
    void printAllContentModels(ostream &out);
    void readAllContentModels(istream &in);
    void printInputSeqs(ostream &out);
    TrainingData* getInputSeqs(){return inputSeqs;};
    ContentModel getContentModel(int type);
    int getK() { return k;}
    int getModelTypes() { return modelTypes;}
    int getNumFrames() { return numFrames;}
    Double getModelTypeProb(int i) { return (*modelTypeProbs)[i];}
    void addTrainingData(TrainingData *td){
	td->next = inputSeqs;
	inputSeqs = td;
	numData++;
    }
private:
    //Double seqProbUnderModel(Matrix<Double> &patprob, char *s, int len, int frame); 

private:
    int k;                  // order of the HMM
    int modelTypes;         // number of different models
    int numFrames;          // number of different reading frames = periodicity
    int numData;            // number of input sequences

    /* For each model type we have a matrix which holds in row f (f=0,1,2, frame) 
     * and row p (p = 1..4^(k+1), number of patrern) the probability of that pattern
     * ending in reading frame f in the model of the respective type.
     */
    vector<ContentModel> *models;
    vector<Double> *modelTypeProbs;       // a priori probability of that model
    TrainingData *inputSeqs;
};


void printContentModels(ostream& out, vector<ContentModel> *models,
			vector<Double> *modelTypeProbs, int k);

#endif
