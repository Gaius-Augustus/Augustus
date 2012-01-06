#include <iostream>

#include "graph.hh"
#include "meaPath.hh"

using namespace std;

void MEApath::findMEApath(){

  getTopologicalOrdering();
  graph->addBackEdges();
 
  //graph->printGraph("flybase_graph.dot");


  relax();   
  
  // graph->printGraphToShell();
  //backtracking
  Node *pos = topSort[0];
  meaPath.push_front(pos);  
  while(pos->pred != NULL){
    meaPath.push_front(pos->pred); 
    pos = pos->pred;
  }
  
  /*
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
  graph->printGraph("2R.1000001-2000000last_section.dot"); */
}

void MEApath::getTopologicalOrdering(){

  for(list<Node*>::iterator node=graph->nodelist.begin(); node!=graph->nodelist.end(); node++){   
    if(processed[graph->getKey(*node)]==0)
      dfs(*node);
  }
}

void MEApath::dfs(Node *n){

  processed[graph->getKey(n)] = n;
  for(list<Edge>::iterator edge=n->edgeoffsets.begin(); edge!=n->edgeoffsets.end(); edge++){
    if(processed[graph->getKey(edge->to)]==0 )
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
      for(list<Edge>::iterator edge=topSort[i]->edgeoffsets.begin(); edge!=topSort[i]->edgeoffsets.end(); edge++){	 
	if(topSort[i]->score + edge->score > edge->to->score){
	  edge->to->score = topSort[i]->score + edge->score;
	  edge->to->pred = topSort[i];
	  nothingChanged = false;
	}
      }    
    }
    if(nothingChanged) 
      continueRelax = false;
  }
  for(list<Node*>::iterator node = graph->nodelist.begin(); node!=graph->nodelist.end(); node++)
    for(list<Edge>::iterator edge = (*node)->edgeoffsets.begin(); edge != (*node)->edgeoffsets.end(); edge++)
      if(edge->to->score < (*node)->score + edge->score)
	cerr<<"MEA (relax): wrong distance "<<edge->to->score<<" at: "<<edge->to->begin<<":"<<edge->to->end<<endl;
 
}
