/**********************************************************************
 * file:    codonMSA.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  datastructure for the codon alignment
 * author:  Lizzy Gerischer
 *
 * date    |   author           |  changes
 * --------|--------------------|------------------------------------------
 * 14.07.16| Lizzy Gerischer    | creation of the file
 **********************************************************************/

#include "codonMSA.hh"
#include "exonmodel.hh"
#include "randseqaccess.hh"
#include "geneMSA.hh"
#include "liftover.hh"
#include <iostream>
#include <unistd.h>
#include <cstring>

using namespace std;


CodonMSA::CodonMSA(string codonAliFilename, double branchlength){

  ctree = NULL;
  seqname = "chr1";

  if(Properties::hasProperty(TREE_KEY)){
    string treeFilename =  Properties::getProperty(TREE_KEY);
    ctree = new PhyloTree(treeFilename);
    ctree->getSpeciesNames(speciesNames);
    readAlignment(codonAliFilename);
  }else{
    readAlignment(codonAliFilename);
    ctree = new PhyloTree(speciesNames,branchlength);
    ctree->getSpeciesNames(speciesNames);
  }

  //for(int i = 0; i<species.size();i++)
  //cout << "species: " << species[i] << " speciesNames: " << speciesNames[i] << endl;

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


  MemSeqAccess msa(speciesSeqs, seqLengths, speciesNames);
  RandSeqAccess *rsa = &msa;
  rsa->printStats();
  GeneMSA geneRange(rsa, alignment);
  geneRange.setTree(ctree);
  geneRange.printStats();

  vector<map<int_fast64_t, ExonCandidate*> > exoncands(speciesNames.size());
  vector<map<int_fast64_t, ExonCandidate*> > addECs(speciesNames.size());

  vector<AnnoSequence*> seqranges(speciesNames.size());
  
  cout << "****************" << endl;     
  for(map<string, char*>::iterator sit=speciesSeqs.begin(); sit!=speciesSeqs.end(); sit++){
    string name = sit->first;
    name.erase(0,5);
    map<string, int>::iterator slit = seqLengths.find(name);
    if(slit == seqLengths.end())
      cout << "seqlength of species " << sit->first << " not found!!!" << endl;
    cout << "species: " << sit->first << " length: " << slit->second << " sequence: ";
    char *c = sit->second;
    for(int i=0; i<slit->second; i++)
      cout << c[i];
    cout << endl;
  }
  cout << "****************" << endl;

  for(int s = 0; s < speciesNames.size(); s++){
    map<string, char*>::iterator sit = speciesSeqs.find(seqname+"."+speciesNames[s]);
    if(sit == speciesSeqs.end())
      cout << "sequence of species " << speciesNames[s] << " not found!!!" << endl;
    map<string, int>::iterator slit = seqLengths.find(speciesNames[s]);
    if(slit == seqLengths.end())
      cout << "seqlength of species " << sit->first << " not found!!!" << endl;

    cout << "speciesname: " << sit->first << " length: " << slit->second; 
    cout << " sequence: ";
    char *e = sit->second;
    for(int i=0; i<slit->second; i++)
      cout << e[i];
    cout << endl;
    if(sit != speciesSeqs.end()){
      cout << " found in speciesSeqs" << endl;
      geneRange.createExonCands(s, sit->second, exoncands[s], addECs[s]);
      seqranges[s] = new AnnoSequence();
      seqranges[s]->sequence = sit->second;
      map<string, int>::iterator slit = seqLengths.find(speciesNames[s]);
      if(slit != seqLengths.end())
	seqranges[s]->length = slit->second;
    }else{
      cout << " NOT found in speciesSeqs" << endl;
    }
  }

  vector<int> offsets = geneRange.getOffsets();
  LiftOver lo(geneRange.getAlignment(), offsets);
  map<int_fast64_t, list<pair<int,ExonCandidate*> > > alignedECs; // hash of aligned ECs 
  lo.projectToAli(addECs,alignedECs);
  addECs.clear();

  lo.projectToGenome(alignedECs, seqranges, exoncands, true);

  geneRange.setExonCands(exoncands);
  exoncands.clear(); 

  //geneRange.printExonCands();

  ExonEvo evo(2);
  vector<double> branchset;
  ctree->getBranchLengths(branchset);
  evo.setBranchLengths(branchset);
  evo.computeLogPmatrices();
  
  list<OrthoExon> hects;
  geneRange.createOrthoExons(hects, alignedECs, &evo);

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
  cout << "setKappa" << endl;
  codonevo.setKappa(4.0);
  cout << "setPi" << endl;
  codonevo.setPi(pi);
  cout << "setBranchLengths" << endl;
  codonevo.setBranchLengths(ct_branchset, 25);
  cout << "setOmega" << endl;
  codonevo.setOmegas(k);
  cout << "setPrior" << endl;
  codonevo.setPrior(0.5);
  //cout << "Omegas, for which substitution matrices are stored:" << endl;
  //codonevo.printOmegas();
  if(Constant::useAArates){
    codonevo.setAAPostProbs();
  }
  /*cout << "Omegas, for which substitution matrices are stored:" << endl;                                                          
    codonevo.printOmegas();
  */
  cout << "computeLogPmatrices" << endl;;
  codonevo.computeLogPmatrices();
 
  // gsl_matrix *P = codonevo.getSubMatrixLogP(0.3, 0.25);                                                                          
  // printCodonMatrix(P);                                                                                                           
  
  GeneMSA::setCodonEvo(&codonevo);

  // get features
  geneRange.computeOmegasEff(hects, seqranges, ctree, NULL);

  geneRange.printOrthoExons(hects);
  
}



void CodonMSA::readAlignment(string alignFilename){ 

  string rowseq, buffer;
  string completeName;
  map<string, size_t> notExistingSpecies;
  string speciesName;
  string seqID;
   int numSpecies;

  ifstream Alignmentfile;

  if(ctree==NULL){
    // read file to find out number of species
    numSpecies = 0;
    Alignmentfile.open(alignFilename.c_str(), ifstream::in);
    if (!Alignmentfile) {
      string errmsg = "Could not open the alignment file " + alignFilename + ".";
      throw PropertiesError(errmsg);
    }
    while(!Alignmentfile.eof()){
      Alignmentfile >> buffer;
      if(buffer[0] == '>')
	numSpecies++;
    }
    Alignmentfile.close();
  }else{
    numSpecies = ctree->numSpecies();
  }
  Alignmentfile.open(alignFilename.c_str(), ifstream::in);
  if (!Alignmentfile) {
    string errmsg = "Could not open the alignment file " + alignFilename + ".";
    throw PropertiesError(errmsg);
  }

  alignment =  new Alignment(numSpecies);
  int numSpeciesFound = 0;

  while (!Alignmentfile.eof()) {
    try {
      Alignmentfile >> buffer;
    } catch (std::ios_base::failure e) {
      throw ProjectError(string("Could not open file ") + alignFilename + ". Make sure this is not a directory.\n");
    }
      // loop over the lines of the alignment block, don't search any further if enough species were found                                                                                                                                                                   
      while (!Alignmentfile.eof()) {
	
	if(buffer[0] == '>'){

	  // store last sequence
	  
	  if(rowseq != ""){
	    addAliRow(speciesName, rowseq, numSpeciesFound, &notExistingSpecies);
	  }
	  speciesName = buffer.substr(1);
	  rowseq = "";
	} else {
	  rowseq += buffer;
	}
	Alignmentfile >> buffer;
      }
      // store last entry
      addAliRow(speciesName, rowseq, numSpeciesFound, &notExistingSpecies);

      if (numSpeciesFound <= 1) // do not consider alignments with less than 2 rows
	delete alignment;
  }
  // clean up
  Alignmentfile.close();
  if (!notExistingSpecies.empty()){
    cerr << "Warning: Species ";
    for (map<string,size_t>::iterator it = notExistingSpecies.begin(); it != notExistingSpecies.end(); ++it)
      cerr << it->first << " ";
    cerr << ((notExistingSpecies.size() > 1)? "are": "is") << " not included in the target list of species. These alignment lines are ingored." << endl;
  }

  cout << *alignment << endl;

}

void CodonMSA::addAliRow(string speciesName, string rowseq, int &numSpeciesFound, map<string, size_t> *notExistingSpecies){

  AlignmentRow *row; 

  if (!alignment->aliLen)
    alignment->aliLen = rowseq.length();
  else if (alignment->aliLen != rowseq.length()) {
    throw ProjectError("Error in multiFasta. Alignment row does not agree in length.");
  }

  row = new AlignmentRow (seqname, 1, plusstrand, rowseq);
  
  int index = ctree->getLeaf(speciesName)->getIdx();
  cout << "index of species " << speciesName << " is " << index << endl;
  
  if (index >= 0) { // species name in the white list                                                                                                                                                                                                                    
    if (!(alignment->rows[index])){ // first row for this species in this block                                                                                                                                                                                          
      alignment->rows[index] = row; // place at the right position                                                                                                                                                                                                       
      int seqlen = alignment->rows[index]->getSeqLen();
      char* seq = new char[seqlen+1];
      rowseq = removeGaps(rowseq);
      if(seqlen != rowseq.size())
	throw ProjectError("Error in CodonMSA::addAliRow: sequence length conflict.");
      copy(rowseq.begin(), rowseq.end(), seq);
      seq[seqlen] = '\0';
            
      speciesSeqs.insert(pair<string, char*>(seqname+"."+speciesName, seq));
      seqLengths.insert(pair<string, int>(speciesName, seqlen));
    
      cout << "length of sequence of species " << speciesName << " is " << alignment->rows[index]->getSeqLen() << endl;
  
      // store chrLen and check whether consistent with previous chrLen                                                                                                                                                                                                  
      numSpeciesFound++;
    } else {
      // multiple rows of same species in same block                                                                                                                                                                                                                     
      throw ProjectError("Error in multiFasta. Multiple entries of same species exist.");
    }
  } else {
    cout << "Species " << speciesName << " not in tree" << endl;
    notExistingSpecies->insert(pair<string, size_t>(speciesName, 1));
    delete row;
  }
}


string CodonMSA::removeGaps(string row){

  row.erase(std::remove(row.begin(), row.end(), '-'), row.end());
  
  for(int i = 0; i < row.size(); i++){
    row[i] = tolower(row[i]);
  }
  return row;
}


void CodonMSA::printOmegaStats(){

  codonevo.graphOmegaOnCodonAli(aliRows, ctree, refSpeciesIdx);
}


