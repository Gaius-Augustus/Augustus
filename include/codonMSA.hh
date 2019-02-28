/**********************************************************************
 * file:    codonMSA.hh
 * licence: Artistic Licence, see file LICENCE.TXT or
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  datastructure for the codon alignment
 * author:  Lizzy Gerischer
 *
 * date    |   author           |  changes
 * --------|--------------------|------------------------------------------
 * 14.07.16| Lizzy Gerischer    | creation of the file
 **********************************************************************/


#ifndef _CODONMSA
#define _CODONMSA

#include <fstream>

#include "properties.hh"
#include "phylotree.hh"
#include "contTimeMC.hh"
#include "alignment.hh"
#include "orthoexon.hh"
#include "geneMSA.hh"

class CodonMSA{
public:
  CodonMSA(string filename, double branchlength);
  ~CodonMSA(){
    delete ctree;
    for(auto it = speciesSeqs.begin(); it != speciesSeqs.end(); it++)
      delete[] it->second;
    delete geneRange;
  }
  
  void readAlignment(string filename);
  void addAliRow(string speciesName, string rowseq, int &numSpeciesFound, map<string, size_t> *notExistingSpecies);
  string removeGaps(string row);
  void printOmegaStats();

  vector<string> aliRows;
  map<string, char*> speciesSeqs; // string is species abbreviation and char* is a pointer to the actual sequence of DNA+Gap characters 
  map<string, int> seqLengths; // contains the sequence length for each species
  vector<string> speciesNames; // name abbreviations of all species
  Alignment *alignment;
  int refSpeciesIdx; // reference Species index as it appears in speciesNames
  CodonEvo codonevo;
  PhyloTree* ctree;
  string seqname;
  list<OrthoExon> hects;
  GeneMSA* geneRange;
};

#endif
