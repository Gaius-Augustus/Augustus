#include <iostream>

#include "graph.hh"
#include "meaPath.hh"

using namespace std;

void MEApath::findMEApath(){

  getTopologicalOrdering();
  graph->addBackEdges();
  relax();   
  
  //backtracking
  backtracking();
  
  //for graphviz dot to draw path
  for(list<Node*>::iterator node=graph->nodelist.begin(); node!=graph->nodelist.end(); node++){
    bool nodeInPath = false;
    bool predInPath = false;
    for(list<Node*>::iterator pathNode=meaPath.begin(); pathNode!=meaPath.end(); pathNode++){
      if((*node)->pred == *pathNode)
	predInPath = true;
      if(*node == *pathNode)
	nodeInPath = true;
      if(predInPath && nodeInPath)
	goto nextNode;
    }
    (*node)->pred = NULL;
    nextNode:;
    }
  //  graph->printGraph("/fs1.data/lizzy/species/tribolium/param_training/train_log/tmp_opt_triboliumMEA/MEA_graph.dot");
}

/*
 * deapth first search
 */

void MEApath::dfs(Node *n){

  processed[graph->getKey(n)] = n;
  for(list<Edge>::iterator edge=n->edges.begin(); edge!=n->edges.end(); edge++){
    if(processed[graph->getKey(edge->to)]==0)
      dfs(edge->to);
  }
  topSort.push_back(n);
}

void MEApath::relax(){

  for(list<Node*>::iterator it=graph->nodelist.begin(); it!=graph->nodelist.end(); it++)
    (*it)->score = - numeric_limits<double>::max(); // set to minimum score

  graph->head->score = 0;  
  bool continueRelax = true;

  while(continueRelax){
    bool nothingChanged = true;
    for(int i = topSort.size()-1; i >= 0; i--){
      for(list<Edge>::iterator edge=topSort[i]->edges.begin(); edge!=topSort[i]->edges.end(); edge++){	 
	if(topSort[i]->score + edge->score > edge->to->score){
	  // update exon distance
	  edge->to->score = topSort[i]->score + edge->score;
	  edge->to->pred = topSort[i];
	  nothingChanged = false;
	}
      }    
    }
    /*
     * @MultSpeciesMode: in this mode, no backedges are added, the graph is truly acyclic. 
     * It is, therefore, sufficient to relax each node exactly one time.
     */
    if(Constant::MultSpeciesMode)
      cerr<<"WARNING! relax(): MultSpeciesMode is on!!!"<<endl;
    if(nothingChanged)
      continueRelax = false;
  }

  for(list<Node*>::iterator node = graph->nodelist.begin(); node!=graph->nodelist.end(); node++){   
    for(list<Edge>::iterator edge = (*node)->edges.begin(); edge != (*node)->edges.end(); edge++)
      if(edge->to->score < (*node)->score + edge->score)
	cerr<<"MEA (relax): wrong distance "<<edge->to->score<<" at: "<<edge->to->begin<<":"<<edge->to->end<<endl;
  }
}

void MEApath::getTopologicalOrdering(){

  for(list<Node*>::iterator node=graph->nodelist.begin(); node!=graph->nodelist.end(); node++){   
    if(processed[graph->getKey(*node)]==0)
      dfs(*node);
  }
}

void MEApath::backtracking(){

  Node *pos = graph->tail;
  meaPath.push_front(pos);
  while(pos->pred != NULL){
    meaPath.push_front(pos->pred);
    pos = pos->pred;
  }
}
