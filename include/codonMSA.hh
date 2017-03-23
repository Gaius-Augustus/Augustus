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

class CodonMSA{
public:
  CodonMSA(string filename);
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
