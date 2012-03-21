#include "mea.hh"
#include "mea7.hh"
#include "meaPath.hh"
#include <iostream>


/*********************************************************************
 * 
 * date     |   author      |  description
 * ---------|---------------|-----------------------------------------
 * 16.18.12 |Stefanie König |  interface for MEA prediction in a graph
 *          |               |  with seven neutral lines
 *          |               |  can't handle overlapping genes (no
 *          |               |  back edges)
 *          |               |  UTR prediction is not possible
*********************************************************************/



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

    //in addition to the sampled exons, more exoncandidates can be added
    list<Status> additionalExons;

    //build Graph
    AugustusGraph myGraph(&stlist, strlength);
    myGraph.buildGraph(&additionalExons);

    //find shortest path
    MEApath path(&myGraph);
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
