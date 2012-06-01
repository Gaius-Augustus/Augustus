/**********************************************************************
 * file:    compgenepred.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  comparative gene prediction on multiple species
 * authors: Mario Stanke
 *
 * date    |   author      |  changes
 * --------|---------------|------------------------------------------
 * 07.03.12| Mario Stanke  | creation of the file
 **********************************************************************/

#include "compgenepred.hh"
#include "orthograph.hh"
#include "orthoexon.hh"

CompGenePred::CompGenePred(){
  if (!Constant::speciesfilenames.empty()) {
    rsa = new MemSeqAccess();
  } else {
    rsa = new DbSeqAccess();
  }
}

void CompGenePred::start(){

  OrthoGraph::tree = new PhyloTree(Constant::treefile);  //has to be initialized before OrthoGraph

  OrthoGraph::numSpecies = OrthoGraph::tree->species.size();

  OrthoGraph::initOutputFiles();

  //GenomicMSA msa;
  //GeneMSA * geneMSA = msa.getNextGene() ;
  //while(geneMSA){
  //geneMSA = msa.getNextGene()

  // TODO: OrthoGraph orthograph(rsa, orthoSeqRanges, GeneMSA *geneMSA);
  OrthoGraph orthograph(rsa);

  // iterative optimization
  orthograph.optimize();

  // transfer longest paths to genes + filter + ouput
  orthograph.outputGenes(minusstrand);

  // }

  OrthoGraph::closeOutputFiles();

  // free memory space of tree
  delete OrthoGraph::tree;

}
