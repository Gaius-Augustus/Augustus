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

class CodonMSA{
public:
  CodonMSA(string filename, double branchlength);
  ~CodonMSA(){
    delete ctree;
    for(auto it = speciesSeqs.begin(); it != speciesSeqs.end(); it++)
      delete[] it->second;
  }
  
  void readAlignment(string filename);
  void addAliRow(string speciesName, string rowseq, int &numSpeciesFound, map<string, size_t> *notExistingSpecies);
  string removeGaps(string row);
  void printOmegaStats();

  vector<string> aliRows;
  map<string, char*> speciesSeqs;
  map<string, int> seqLengths;
  vector<string> speciesNames;
  Alignment *alignment;
  int refSpeciesIdx;
  CodonEvo codonevo;
  PhyloTree* ctree;
  string seqname;
};

#endif
