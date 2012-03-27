/**********************************************************************
 * file:    mea7.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  interface for MEA prediction in a graph  with seven neutral lines
 * authors: Stefanie König
 *
 * date    |   author      |  changes
 * --------|---------------|------------------------------------------
 * 16.18.12|Stefanie König | creation of the file
 **********************************************************************/


#include "mea.hh"
#include "mea7.hh"
#include "meaPath.hh"
#include "orthograph.hh"
#include "orthoexon.hh"
#include "phylotree.hh"
#include <iostream>


void getMEAtranscripts7(list<Gene> *meaGenes, list<Gene> *alltranscripts, int strlength){

  if(!alltranscripts->empty()){
   
    bool utr;
    try {
      utr = Properties::getBoolProperty("UTR");
    } catch (...) {
      utr = false;
    }

    if(utr){
      cerr<<"UTR prediction not possible in combination with option 'seven_neutralLines'"<<endl;
      utr = false;
    }

    list<Status> stlist;

    //builds datastructure needed for the graph representation
    buildDatastructure(alltranscripts, utr, stlist);

    //read in orthologous exons
    list<OrthoExon> all_orthoex = readOrthoExons(Constant::orthoexons);
    
    //in addition to the sampled exons, more exoncandidates can be added
    list<Status*> additionalExons;

    for(list<OrthoExon>::iterator it = all_orthoex.begin(); it != all_orthoex.end(); it++){
      additionalExons.push_back(it->orthoex[0]);
    }

    //build Graph
    OrthoGraph myGraph =  OrthoGraph(&stlist, strlength, additionalExons);

    //find shortest path
    MEApath path(myGraph.graphs[0]);
    path.findMEApath7();
    
    //convert back to augustus datastructure for genes
    getMeaGenelist7(path.getPath(), meaGenes);
  }   
}

void getMeaGenelist7(list<Node*> meaPath, list<Gene> *meaGenes){

  Node* current = meaPath.back();
  Node* head = meaPath.front();
  Node* predcurrent;

  Gene *currentGene = new Gene();

  while(current != head){

    while(current->item == NULL){
      if (current == head){
	goto end;
      }
      current=current->pred;  //exon1
    }
    State *ex = new State(*((State*)current->item));
    addExonToGene(currentGene, ex);
    if(current->pred->n_type  == IR){ //end of gene
      setGeneProperties(currentGene);
      meaGenes->push_front(*currentGene);
      delete currentGene;
      currentGene = new Gene();
      current=current->pred;
    }
    else{
      predcurrent=current->pred;
      while(predcurrent->item == NULL){
	if(predcurrent == head){
	  setGeneProperties(currentGene);
	  meaGenes->push_front(*currentGene);
	  delete currentGene;
	  goto end;
	}
       	predcurrent=predcurrent->pred; //exon2
      }
      addIntronToGene(currentGene, predcurrent, current); //add intron exon2->exon1
      current=predcurrent;
    } 
  }
 end:;
}
