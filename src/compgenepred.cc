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
#include "orthoexon.hh"
#include "extrinsicinfo.hh"
#include "mea.hh"

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

  map< vector<string>, list<OrthoExon> > all_orthoex = readOrthoExons(Constant::orthoexons);  //read in orthologous exons from file
  //writeOrthoExons(all_orthoex);
  
  OrthoGraph orthograph;

  vector<string> speciesname = OrthoExon::species;
  vector<string> chrName;
  chrName.push_back("chr21");
  chrName.push_back("chr17");
  chrName.push_back("chr31");
  int start[] = {44836600, 31983200, 39789600} ;
  int end[] = {44846200, 31992000, 39798300};
  Strand strand = minusstrand;
 
  for (int s = 0; s < 3; s++) {
    AnnoSequence *seqRange = rsa->getSeq(speciesname[s], chrName[s], start[s], end[s], strand);
    if (seqRange){
      //cout << seqRange->seqname << "\t" << seqRange->length << "\t" << seqRange->offset << endl;
      //cout << seqRange->sequence << endl;
      /*
       * build list of additional exoncandidates, which are inserted in the graph
       */
      list<Status> additionalExons;
      for(list<OrthoExon>::iterator it = all_orthoex[chrName].begin(); it !=  all_orthoex[chrName].end(); it++){
	if(it->orthoex[s] != NULL && it->orthoex[s]->begin >= start[s] &&  it->orthoex[s]->end <= end[s] ){
	  Status state = Status(CDS, it->orthoex[s]->begin-seqRange->offset, it->orthoex[s]->end-seqRange->offset, 0.0, it->orthoex[s]);
	  additionalExons.push_back(state);
	}
      }
      /*for(list<Status>::iterator it = additionalExons.begin(); it != additionalExons.end(); it++){
	cout << it->begin <<"\t"<< it->end <<"\t"<<((State*) it->item)->type << endl;
	}*/

      //namgene.doViterbiPiecewise(sfc, seqRange, bothstrands); // builds graph for each species

      list<Gene> *alltranscripts = NULL; //von namgene
      if(alltranscripts){
	if(!alltranscripts->empty()){
	  /*
	   * build datastructure for graph representation
	   * @stlist : list of all sampled states
	   */
	  list<Status> stlist;
	  buildDatastructure(alltranscripts, false, stlist);
	  //build graph
	  AugustusGraph *singleGraph = new AugustusGraph(&stlist, seqRange->length);
	  singleGraph->buildGraph(additionalExons);
	  //add graph for species to OrthoGraph
	  orthograph.addSingleGraph(speciesname[s], singleGraph);
	}
      }
    } else {
      cerr << "random sequence access failed on " << speciesname[s] << ", " << chrName[s] << ", " << start[s] << ", " <<  end[s] << ", " <<  strand << endl;
    }
  } 
}
