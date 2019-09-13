/*
 * codonMSA.cc
 *
 * License: Artistic License, see file LICENSE.TXT or 
 *          https://opensource.org/licenses/artistic-license-1.0
 */

#include "codonMSA.hh"
#include "exonmodel.hh"
#include <iostream>
#include <unistd.h>

using namespace std;

CodonMSA::CodonMSA(string codonAliFilename, double branchlength){

  aliLen = 0;
  ctree=NULL;

  if(Properties::hasProperty(TREE_KEY)){
    string treeFilename =  Properties::getProperty(TREE_KEY);
    ctree = new PhyloTree(treeFilename);
    readAlignment(codonAliFilename);
  }else{
    readAlignment(codonAliFilename);
    ctree = new PhyloTree(speciesNames,branchlength);
  }

  vector<string> species;
  ctree->getSpeciesNames(species);
  //for(int i = 0; i<species.size();i++)
  //cout << "species: " << species[i] << " speciesNames: " << speciesNames[i] << endl;

  if(species != speciesNames)
    throw ProjectError("inconsistent species name vectors");

  /*  double ctree_scaling_factor = 1; // scaling factor to scale branch lengths in codon tree to one codon substitution per time unit
  try {
    ctree_scaling_factor = Properties::getdoubleProperty("/CompPred/scale_codontree");
  } catch (...) {
    ctree_scaling_factor = 1;
  }
  if(ctree_scaling_factor <= 0.0){
    cerr << "No negative scaling factor allowed. /CompPred/scale_codontree must be a positive real number. Will use =1." << endl;   
    ctree_scaling_factor=1;
  }
  
  ctree->scaleTree(ctree_scaling_factor); // scale branch lengths to codon substitutions                                             
  */  

vector<double> ct_branchset;
  ctree->getBranchLengths(ct_branchset);
  int k; // number of omegas
  try {
    k = Properties::getIntProperty("/CompPred/num_omega");
  }catch(...){
    k = 20;
  }
  // TODO: codonusage

  //BaseCount::init();
    
  //PP::initConstants();
  //  NAMGene namgene; // creates and initializes the states                                                                            
  StateModel::readAllParameters(); // read in the parameter files: species_{igenic,exon,intron,utr}_probs.pbl                       

  double *pi = ExonModel::getCodonUsage();
  /**********
  double pi[64];
  cout << "pi: ";
  for (int i=0; i<64; i++){
    pi[i] = (double)1/64;
    cout << pi[i] << " ";
  }
  cout << endl;
  ***********/
  codonevo.setKappa(4.0);
  codonevo.setPi(pi);
  codonevo.setBranchLengths(ct_branchset, 25);
  codonevo.setOmegas(k);
  codonevo.setPrior(0.5);
  //cout << "Omegas, for which substitution matrices are stored:" << endl;
  //codonevo.printOmegas();
  if(Constant::useAArates){
    codonevo.setAAPostProbs();
  }
  /*cout << "Omegas, for which substitution matrices are stored:" << endl;                                                          
    codonevo.printOmegas();*/
  codonevo.computeLogPmatrices();

  // gsl_matrix *P = codonevo.getSubMatrixLogP(0.3, 0.25);                                                                          
  // printCodonMatrix(P);                                                                                                           
  //GeneMSA::setCodonEvo(&codonevo);
}

void CodonMSA::readAlignment(string filename){
  
  string rowseq;
  string speciesName;
  int c = 0;

  if(ctree != NULL){
    ctree->getSpeciesNames(speciesNames);
    aliRows.resize(speciesNames.size(),"");
  }
  
  ifstream Alignmentfile;
  Alignmentfile.open(filename.c_str(), ifstream::in);
  if (!Alignmentfile) {
    string errmsg = "Could not open the codon alignment file " + filename + ".";
    throw PropertiesError(errmsg);
  }
  while (!Alignmentfile.eof()) {
    try {
      Alignmentfile >> speciesName >> rowseq;
      if(!Alignmentfile.eof()){
	speciesName.erase(0,1);
	if(ctree == NULL){
	  aliRows.push_back(rowseq);
	  speciesNames.push_back(speciesName);
	}else{
	  int pos = find(speciesNames.begin(), speciesNames.end(), speciesName) - speciesNames.begin();
	  if(pos >= speciesNames.size()){
	    throw ProjectError( string("Species ") + speciesName + " not found in phylogenetic tree.");
	  }
	  if(aliRows[pos] == "")
	    aliRows[pos] = rowseq;
	  else
	    cerr << "Warning: Multiple sequences for species " << speciesName << " found! Paralogs will be ignored. Run ESPOCA with a gene tree that includes all paralogs or use the default star tree if paralogs shall be considered." << endl;
	  if(c == 0)
	    refSpeciesIdx = pos;
	  c++;
	}
      }
      // cout << "species name: " << speciesName << "\t rowseq: " << rowseq << endl;
    } catch (std::ios_base::failure &e) {
      throw ProjectError(string("Could not open file ") + filename + ". Make sure this is not a directory.\n");
    }
    if(aliLen){
      if(aliLen != rowseq.length()){
	throw ProjectError("codon alignment rows have different size.");
      }
    }else{
      aliLen = rowseq.size();
    }
  }
}


void CodonMSA::printOmegaStats(){

  codonevo.graphOmegaOnCodonAli(aliRows, ctree, refSpeciesIdx);
}


