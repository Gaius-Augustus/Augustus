/*
 * codonMSA.hh
 *
 * License: Artistic License, see file LICENSE.TXT or 
 *          https://opensource.org/licenses/artistic-license-1.0
 */

#ifndef _CODONMSA
#define _CODONMSA

#include <fstream>

#include "properties.hh"
#include "phylotree.hh"
#include "contTimeMC.hh"
#include "codonevo.hh"

/**
 * @brief datastructure for the codon alignment
 * 
 * @author Lizzy Gerischer
 */
class CodonMSA{
public:
  CodonMSA(string filename, double branchlength);
  ~CodonMSA(){
    delete ctree;
  }
  
  void readAlignment(string filename);
  void printOmegaStats();

  vector<string> aliRows;
  vector<string> speciesNames;
  size_t aliLen;
  int refSpeciesIdx;
  CodonEvo codonevo;
  PhyloTree* ctree;
};

#endif
