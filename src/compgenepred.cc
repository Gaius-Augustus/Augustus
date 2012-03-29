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
#include "randseqaccess.hh"
#include "namgene.hh"
#include "statemodel.hh"
#include "extrinsicinfo.hh"

CompGenePred::CompGenePred(){
  if (!Constant::speciesfilenames.empty()) {
    rsa = new MemSeqAccess();
  } else {
    rsa = new DbSeqAccess();
  }
}

void CompGenePred::start(){
  // read in alignment, determine orthologous sequence fragments
  
  // determine object that holds a sequence range for each species
  // loop over species

  // first species (testing only)
  NAMGene namgene; // creates and initializes the states
  SequenceFeatureCollection sfc(NULL); // hints, empty for now, will later read in hints for sequence ranges from database
  StateModel::readAllParameters(); // read in the parameter files: species_{igenic,exon,intron,utr}_probs.pbl

  string speciesname = "dog";
  string chrName = "chr31";
  int start = 39798072;
  int end = 39798218;
  Strand strand = minusstrand;

  for (int s=0;s<2;s++) {
    AnnoSequence *seqRange = rsa->getSeq(speciesname, chrName, start+s*1000, end+s*1000, strand);
    if (seqRange){
      namgene.doViterbiPiecewise(sfc, seqRange, bothstrands); // builds graph for each species
      // G = namgene.buildgraph
      // G.addExons() orthologous 
      // add G for species to OrthoGraph
    } else {
      cerr << "random sequence access failed on " << speciesname << ", " << chrName << ", " << start << ", " <<  end << ", " <<  strand << endl;
    }
  } 
}
