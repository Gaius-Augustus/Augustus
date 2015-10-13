#include <iostream>
#include <limits>
#include <list>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>

#include "graph.hh"
#include "mea.hh"

using namespace std;

Graph::~Graph(){

  for(list<Node*>::iterator it=nodelist.begin(); it!=nodelist.end(); it++)
    delete *it;
}


void Graph::buildGraph(){

  vector<Node*> neutralLine; //represents the area of intergenic regions

  // get size of neutral line
  getSizeNeutralLine();

  for(int i=0; i<max-min+1; i++)
    neutralLine.push_back(NULL);
  
  head = new Node(-1,-1);
  nodelist.push_back(head);
  tail = new Node(max+1,max+1);
  nodelist.push_back(tail);

  calculateBaseScores();
 
  // add all existing exons and introns to the graph

  for(list<Status>::iterator it=statelist->begin(); it!=statelist->end(); it++){
    if(it->name == CDS || it->name == utr3 || it->name == utr5){     
      if(it->next != NULL){
	if(it->next->name == CDS || it->next->name == utr3 || it->next->name == utr5)
	  addPair(&(*it), it->next, neutralLine);
	else{
	  if(it->next->next != NULL)
	    addPair(&(*it), it->next->next, neutralLine);
	  else{
	    addExon(&(*it), neutralLine);
	    addEdgeToTail(&(*it));
	  }
	}
      }
      else{
	addExon(&(*it), neutralLine);
	addEdgeToTail(&(*it));
      }
    }
  } 
  createNeutralLine(neutralLine);
 
  // add edges of incomplete transcripts to the neutral line

  addEdgeFromHead(&statelist->front());
 
  for(list<Status>::iterator it=statelist->begin(); it!=statelist->end(); it++){
    if(it->next==NULL){      
      if(&(*it) != &statelist->back()){
	it++;
	addEdgeFromHead(&(*it));
	--it;
      }
    }
  }
  nodelist.sort(compareNodes);
  
  bool noCompatibleEdges;
  try {
    noCompatibleEdges = Properties::getBoolProperty("/MeaPrediction/no_compatible_edges");
  } catch (...) {
    noCompatibleEdges = false;
  }
  if(!noCompatibleEdges)
    addCompatibleEdges();

  // add node weight to edge weight
  addWeightToEdge();
  
  //sort edgeoffset lists

  for(list<Node*>::iterator node=nodelist.begin(); node!=nodelist.end(); node++)
    (*node)->edges.sort(compareEdges);
}


bool Graph::edgeExists(Node *e1, Node *e2){
  if(e1 == NULL || e2 == NULL)
    return false;
  for(list<Edge>::iterator it=e1->edges.begin(); it!=e1->edges.end(); it++){
    if(it->to == e2)
      return true;
  }
  return false;
}

Node* Graph::addExon(Status *exon, vector<Node*> &neutralLine){

  if(!alreadyProcessed(exon)){
    Node *ex = new Node(exon->begin, exon->end, setScore(exon), exon->item, sampled);
    nodelist.push_back(ex);
    addToHash(ex);
    if(exonAtGeneStart(exon)){
      // include edge from neutral line to exon

      Node *neut = new Node(ex->begin, ex->begin, 0.0, NULL, IR);
      Edge intron(ex, false);
      if(!alreadyProcessed(neut)){	
	neutralLine[ex->begin-min] = neut;
	neut->edges.push_back(intron);
	nodelist.push_back(neut);
	addToHash(neut);
      }
      else{
	getNode(neut)->edges.push_back(intron);
	delete neut;
      }
    }
    if(exonAtGeneEnd(exon)){
	// include edge from exon to neutral line      

      Node *neut = new Node(ex->end, ex->end, 0.0, NULL, IR);
      if(!alreadyProcessed(neut)){
	neutralLine[ex->end-min] = neut;
	Edge intron(neut, false);
	ex->edges.push_back(intron);
	nodelist.push_back(neut);
	addToHash(neut);
      }
      else{
	Edge intron(getNode(neut), false);
	ex->edges.push_back(intron);      
	delete neut;
      }
    }
    return ex;
  }
  return getNode(exon);
}

void Graph::addPair(Status *exon1, Status *exon2, vector<Node*> &neutralLine){

  Node *e1 = addExon(exon1, neutralLine);
  Node *e2 = addExon(exon2, neutralLine);
  if(exon1->next == exon2 && !edgeExists(e1,e2)){
    Edge in(e2, false);
    e1->edges.push_back(in);  
  }
  else if(exon1->next != exon2 && !edgeExists(e1,e2) && !mergedStopcodon(e1,e2)){
    Edge in(e2, false, setScore(exon1->next), exon1->next->item);
    e1->edges.push_back(in);  
  }
}

void Graph::createNeutralLine(vector<Node*> &neutralLine, double weight, bool onlyComplete){

    Node *pos = head;
    int i = 0;
    int n = neutralLine.size();

    for(int j=0; j<n; j++){
	if(neutralLine[j] != NULL){
	    if( pos!=head || (!onlyComplete || neutralLine[j]->n_type == IR) ){
		double w = weight *(j-i);
		Edge neut(neutralLine[j], true, w);
		pos->edges.push_back(neut);
	    }
	    pos = neutralLine[j];
	    i = j;
	}
    }
    if( !onlyComplete || pos->n_type == IR ){
	double w = weight *(n-i);
	Edge lastEdge(tail, true, w);
	pos->edges.push_back(lastEdge);
    }
}

/*
 * adds introns between exons under biologically reasonable conditions
 */

void Graph::addCompatibleEdges(){

  statelist->sort(compareStatus);
  map<string,Node*> processedStartNodes;
  map<string,Node*> processedEndNodes;
  
  for(list<Node*>::iterator e1=nodelist.begin(); e1!=nodelist.end(); e1++){
    if((*e1)->item != NULL && processedStartNodes[getKey(*e1)] == 0){
      processedStartNodes[getKey(*e1)] = *e1;
      int count = 0;
      processedEndNodes.clear();
      for(list<Node*>::iterator e2=e1; e2!=nodelist.end(); e2++){
	if((*e2)->item != NULL && processedEndNodes[getKey(*e2)] == 0){
	  processedEndNodes[getKey(*e2)] = *e2;
	 
	  if(compatible(*e1,*e2) && !edgeExists(*e1,*e2)){
	    insertIntron(*e1,*e2);
	    count++;
	  }
	  if(count>10 || (*e2)->begin - (*e1)->end > 5000) //TODO: max_intron_length restriction!
	    break;	 
	}	
      }
    }
  }
}

/*
 * adds back edges to the neutral line only if it does not cause non-neutral loops
 */

void Graph::addBackEdges(){
  // get neutral line nodes out of all nodes in the order from head to tail
  map<string,Node*> inQueue;
  list<Node*> neutralNodes;
  Node *pos = head;
  while(pos != tail){
	neutralNodes.push_back(pos);
    for(list<Edge>::iterator edge=pos->edges.begin(); edge!=pos->edges.end(); edge++){
      if(edge->neutral){
	pos = edge->to;
	break;
      }
      else if(&(*edge) == &pos->edges.back())
	cerr<<"ERROR: (mea) neutral line has gap!"<<endl;
    }
  }
  neutralNodes.push_back(tail);

  // for every neutral line node
  for(list<Node*>::iterator fromNeut=neutralNodes.begin(); fromNeut!=neutralNodes.end(); fromNeut++){  
    if(nonneutralIncomingEdge(*fromNeut)){
      bool nonNeutralLoop = false;
      // for every neutral line node (backwards) in front of the actual neutral line node (core node)
      for(list<Node*>::iterator toNeut=fromNeut; toNeut!=neutralNodes.begin(); toNeut--){   
	if(nonNeutralLoop)
	  break;
	nonNeutralLoop = true;
	queue<Node*> q;
	inQueue.clear();
	int nrNonNeutralEdges = 0;
    // figure out if there are edges to nonneutral nodes
	for(list<Edge>::iterator edge=(*toNeut)->edges.begin(); edge!=(*toNeut)->edges.end(); edge++)
	  if(!edge->neutral)
	    nrNonNeutralEdges++;

	// if there are edges to nonneutral nodes do the following to check if there is a nontrivial way to the core node
	if(nrNonNeutralEdges > 0 && *fromNeut != *toNeut){
	  q.push(*toNeut);
	  while(!q.empty()){
	    Node *pos = q.front();
	    q.pop();
	    for(list<Edge>::iterator edge=pos->edges.begin(); edge!=pos->edges.end(); edge++){  
	      if(inQueue[getKey(edge->to)] == 0){
		q.push(edge->to);
		inQueue[getKey(edge->to)] = edge->to;
	      }
          // if the programm finds a neutral line node which comes not later then the core node on this nontrivial way, no additional backadges will be added
	      if(pos->item != NULL && edge->to->item == NULL && edge->to->begin <= (*fromNeut)->begin)
		goto nextEdge;
	      // if the programm is already farer then the core node without finding a nontrivial way to a neutral line node, a backedge can be added
	      if(minInQueue(&q) > (*fromNeut)->begin){
		insertIntron(*fromNeut,*toNeut);
		nonNeutralLoop = false;
		goto nextEdge;
	      }
	    }	    
	  }
	}
      nextEdge:;
      }
    }
  }
}

void Graph::addBackEdgesComp(){
  // get neutral line nodes out of all nodes in the order from head to tail
  bool strongOverlap = false;
  list<Node*> neutralNodes;
  Node *pos = head;
  while(pos != tail){
	neutralNodes.push_back(pos);
    for(list<Edge>::iterator edge=pos->edges.begin(); edge!=pos->edges.end(); edge++){
      if((edge->to)->n_type == IR || (edge->to) == tail){
	  pos = edge->to;
	  break;
      }
      else if(&(*edge) == &pos->edges.back())
	cerr<<"ERROR: (mea) neutral line has gap!"<<endl;
    }
  }
  neutralNodes.push_back(tail);
  Node *lastSource = tail;
  tail->nextNontrivialNeutNode = tail;
  // for every neutral line node backwards, starting at tail, add backedges and, if necessary, auxiliary Nodes
  for (list<Node*>::reverse_iterator actSource=neutralNodes.rbegin(); actSource!=neutralNodes.rend(); actSource++){
    if (*actSource == head){break;}             // maybe change the for-loop, to insert this condition in it (without this condition, the program adds unneeded auxiliary nodes on the head)
    (*actSource)->nextNontrivialNeutNode = lastSource->nextNontrivialNeutNode;
    queue<Node*> exonQueue;
    size_t count = 0;
    list<Node*> actStopExons;
    // for every "nonneutral" (start-)Edge from source go every possible new way (dont go to a node, if you already found it from another source)
    for (list<Edge>::iterator startEdge=(*actSource)->edges.begin(); startEdge!=(*actSource)->edges.end(); startEdge++){
      if ((startEdge->to)->prevNontrivialNeutNode == NULL && (startEdge->to)->n_type != IR){
        (startEdge->to)->prevNontrivialNeutNode = (*actSource);
        (startEdge->to)->index = count;
        exonQueue.push(startEdge->to);
        Node *actMinStopExon = NULL;            // save the actual minimum stop for this path; queue is NOT ordered with end positions, so it´s possible, that you found the minimum at the end
        // while there is still a new node (a node that was not already found) 
        while (!exonQueue.empty()){
          Node *pos=exonQueue.front();
          exonQueue.pop();
          // for every edge from startEdge go every possible new way (dont go to a node, if you already found it from another source); add a back edge if you found a way back to IR
          for (list<Edge>::iterator edge=pos->edges.begin(); edge!=pos->edges.end(); edge++){
            if (pos->n_type != IR && (edge->to)->n_type == IR){


              actStopExons.push_back(pos);
              if(!edgeExists(pos,lastSource)){
		// if strong overlap is wanted (one transcript can overlap a whole transcript)
		if (strongOverlap){
                  if ((edge->to)->begin > lastSource->begin){             // see YY     
		    edge->to = lastSource;                                // see YY
		  }                                                       // YY: if backedges should get a penalty, this lines has to be deleted
                }else{	// only weak overlap
		  Node *suitableTarget = lastSource;
		    while (suitableTarget->nextNontrivialNeutNode != NULL && suitableTarget->nextNontrivialNeutNode->end + 1 <= edge->to->end && !(suitableTarget == tail)){
		      Node* tempNode = NULL;
		      for (list<Edge>::iterator targetEdge=suitableTarget->edges.begin(); targetEdge!=suitableTarget->edges.end(); targetEdge++){
		        if ((targetEdge->to)->n_type == IR){
			  tempNode = targetEdge->to;
		        }
		        if ((targetEdge->to) == tail){
			  tempNode = edge->to;
		        }
		      }
		      if (tempNode != NULL){
			suitableTarget = tempNode;
		      }
		    }
		    if ((edge->to)->begin > suitableTarget->begin){
		      edge->to = suitableTarget;
		    }
		}
                // Edge edgeNew(lastSource,false, edge->score);         // see WW
                // pos->edges.push_back(edgeNew);                       // WW: these lines are needed, if backedges should get a penalty 
              }
              if (actMinStopExon){
                if (actMinStopExon->end > pos->end){actMinStopExon = pos;}
              }else{actMinStopExon = pos;}
            }


            if ((edge->to)->prevNontrivialNeutNode == NULL || ((edge->to)->prevNontrivialNeutNode == (*actSource) && (edge->to)->index != count)){
              (edge->to)->prevNontrivialNeutNode = (*actSource);
              (edge->to)->index = count;
              exonQueue.push(edge->to);
            }
          }
        }
	if (actMinStopExon){
          (startEdge->to)->nextNontrivialNeutNode = actMinStopExon;
	  if ((*actSource)->nextNontrivialNeutNode == NULL){
	    (*actSource)->nextNontrivialNeutNode = actMinStopExon;
	  }else{
	    if ((*actSource)->nextNontrivialNeutNode->end > actMinStopExon->end){
	      (*actSource)->nextNontrivialNeutNode = actMinStopExon;
	    }
	  }
	}
      }
      count++;
    }

    // if more then one exon starts from this source, add auxiliary nodes for every different transcript end (only if strongOverlap is true)
    if ((*actSource)->edges.size()>2 && strongOverlap){
      actStopExons.sort(compareNodeEnds);
      list<Node*> auxiliaryNodeList;
      // (again) for every "nonneutral" (start-)Edge from source add auxiliary nodes for every transcript with same start and different minStop
      for (list<Edge>::iterator startEdge=(*actSource)->edges.begin(); startEdge!=(*actSource)->edges.end(); startEdge++){
        if ((startEdge->to)->n_type == IR){continue;}                       // this condition can be extended, if auxNodes gets a new type (!=IR)
        int auxEnd;
        if ((startEdge->to)->nextNontrivialNeutNode){
          auxEnd = ((startEdge->to)->nextNontrivialNeutNode)->end;
        }else{
          auxEnd = tail->end;
        }
        Node* auxiliaryNode = NULL;
        // if there is already a auxNode with the same End, use it instead of creating a new one
        for (list<Node*>::iterator auxNode=auxiliaryNodeList.begin(); auxNode!=auxiliaryNodeList.end(); auxNode++){
          if ((*auxNode)->end == auxEnd){auxiliaryNode = (*auxNode);}
        }
        // create a new auxNode if necessary
        if(auxiliaryNode == NULL){
          auxiliaryNode = new Node((*actSource)->begin,auxEnd);
          auxiliaryNode->n_type=IR;                             // maybe another node_type later
          auxiliaryNodeList.push_back(auxiliaryNode);
        }
        // add edge from auxiliaryNode to startEdge->to
        auxiliaryNode->prevNontrivialNeutNode = (*actSource);
        Edge edgeNew(startEdge->to,false);
        auxiliaryNode->edges.push_back(edgeNew);
        // delete the startEdge (this edge is not forbidden but also not necessary and this step reduces the size of the graph)
        startEdge = (*actSource)->edges.erase(startEdge);
        startEdge--;
      }
      int lastEnd = (*actSource)->begin + 1;
      auxiliaryNodeList.sort(compareNodeEnds);
      Node* lastAuxNode = (*actSource);
      // for every auxiliary Node, add edges from every actual (found from the actual source) stop exon if this auxiliary node leads to the nearest next stop exon
      for (list<Node*>::iterator auxNode=auxiliaryNodeList.begin(); auxNode!=auxiliaryNodeList.end(); auxNode++){
        for (list<Node*>::iterator stopExon=actStopExons.begin(); stopExon!=actStopExons.end(); stopExon++){
          if((*stopExon)->end >= lastEnd && (*stopExon)->end < (*auxNode)->end){
            //double edgeScore = 0;                             // see XX
            for (list<Edge>::iterator edge=(*stopExon)->edges.begin(); edge!=(*stopExon)->edges.end(); edge++){
              if ((edge->to)->n_type == IR){
                if ((edge->to)->begin > (*auxNode)->begin){
                  //edgeScore = edge->score;                    // see XX
                  edge->to = (*auxNode);                        // if backedges should get a penalty, this line has to be deleted
                }
              }
            }
            // Edge edgeNew((*auxNode),false, edgeScore);       // see XX
            // (*stopExon)->edges.push_back(edgeNew);           // XX: these lines are needed, if backedges should get a penalty 
          }
        }
        // connect the auxiliary nodes in the sorted order, beginning from the neutral line node at this position
        Edge edgeNew((*auxNode),false);
        lastAuxNode->edges.push_back(edgeNew);

        lastAuxNode = (*auxNode);
        nodelist.push_back(*auxNode);
        lastEnd = (*auxNode)->end;
      }

      // add Edge from the last auxiliary node to the neutral line
      Edge edgeNew(lastSource,false);
      lastAuxNode->edges.push_back(edgeNew);

      for (list<Node*>::iterator auxNode=auxiliaryNodeList.begin(); auxNode!=auxiliaryNodeList.end(); auxNode++){
        (*auxNode)->end = (*auxNode)->begin;
      }
    }
    lastSource = (*actSource);
  }
}

void Graph::tarjan(){
  for (list<Node*>::iterator node=nodelist.begin(); node!=nodelist.end(); node++){
    (*node)->label = 0;            // if label is 1, node is in stack
    (*node)->index = 0;
  }
  stack<Node*> S;
  size_t index = 1;

  for (list<Node*>::iterator node=nodelist.begin(); node!=nodelist.end(); node++){
    if ((*node)->index == 0){
      tarjanAlg(*node, S, index);
    }
  }
}

void Graph::tarjanAlg(Node* from, stack<Node*> &S, size_t &index){
  from->index = index;
  from->lowlink = index;
  index++;
  S.push(from);
  from->label = 1;
  for (list<Edge>::iterator edge=from->edges.begin(); edge!=from->edges.end(); edge++){
    if ((edge->to)->index == 0){
      tarjanAlg(edge->to, S, index);
      from->lowlink = std::min(from->lowlink, (edge->to)->lowlink);
    }else if ((edge->to)->label == 1){
      from->lowlink = std::min(from->lowlink, (edge->to)->index);
    }
  }
  if (from->lowlink == from->index){
    Node* actNode = NULL;
    list<Node*> component;
    do{
      actNode = S.top();
      S.pop();
      actNode->label = 0;
      component.push_back(actNode);
    }while (from != actNode);
    if (component.size() > 1){cerr << "There is a loop in the graph of the cgp. It is possible, that the result is not the optimum!" << endl;}
  }
}

int Graph::minInQueue(queue<Node*> *q){

  int min = q->front()->begin;
  int n = q->size();
  for(int i=0; i<n; i++){
    Node *node = q->front();
    if(node->begin < min)
      min = node->begin;
    q->pop();
    q->push(node);
  }
  return min;
}

bool Graph::nonneutralIncomingEdge(Node *exon){

  for(list<Node*>::iterator ex=nodelist.begin(); ex!=nodelist.end(); ex++){
    if((*ex)->begin > exon->begin)
      return false;
    for(list<Edge>::iterator edge=(*ex)->edges.begin(); edge!=(*ex)->edges.end(); edge++){
      if(!edge->neutral && edge->to == exon)
	return true;
    }
  }
  return false;
}

void Graph::printGraphToShell(){
  cout<<"****************GRAPH******************"<<endl<<endl;;
for(list<Node*>::iterator node = nodelist.begin(); node != nodelist.end(); node++){
  cout<<"-------------------------------------------------------"<<endl<<setw(10)<<"Node";
  if((*node)->item != NULL) 
    cout<<((State*)(*node)->item)->type;
  cout<<":"<<setw(10)<<(*node)->begin<<setw(10)<<(*node)->end<<setw(10)<<(*node)->score;
    if((*node)->item == NULL)
      cout<<setw(10)<<"neutral";
    for(list<Edge>::iterator edge = (*node)->edges.begin(); edge != (*node)->edges.end(); edge++){
      cout<<endl<<setw(10);
      if(edge->item != NULL)
	cout<<((State*)edge->item)->type;
      cout<<"----->"<<setw(10)<<edge->to->begin<<setw(10)<<edge->to->end<<setw(10)<<edge->to->score;
      if(edge->neutral)
	cout<<setw(10)<<"neutral";
      if(edge->to->pred == *node)
	cout<<setw(10)<<"pred";
      cout<<endl<<setw(10)<<edge->score<<endl;
    }
  }
}

void Graph::getSizeNeutralLine(){

  max = 0; 
  min = numeric_limits<int>::max();
  for(list<Status>::iterator it=statelist->begin(); it!=statelist->end(); it++){
    if(it->name >= CDS && it->name <intron){
      if(it->end > max)
	max = it->end;
      
      if(it->begin < min)
	min = it->begin;    
    }
  }
}

void Graph::addWeightToEdge(){

  //adds node weight to all outgoing edges

  for(list<Node*>::iterator node = nodelist.begin(); node != nodelist.end(); node++){
    for(list<Edge>::iterator edge = (*node)->edges.begin(); edge != (*node)->edges.end(); edge++){
      edge->score += (*node)->score;
    }
  }
}

/*
 * program specific function implementation
 */


bool AugustusGraph::exonAtGeneStart(Status *st){
  
  StateType Type = ((State*)st->item)->type;

  return( (utr==false && ( Type == singleG || Type == initial0 || Type == initial1 || Type == initial2 || Type == rsingleG || Type == rterminal0 || Type == rterminal1 || Type == rterminal2 )) || (utr==true && (Type == utr5single || Type == utr5init || Type == rutr3single || Type == rutr3term)) );      
}

bool AugustusGraph::exonAtGeneEnd(Status *st){

  StateType Type = ((State*)st->item)->type;

  return((utr==false && (Type == singleG || Type == terminal || Type == rsingleG || Type == rinitial)) || (utr==true && (Type == utr3single || Type == utr3term || Type == rutr5single || Type == rutr5init)) );
}

bool AugustusGraph::exonAtCodingStart(Node *st){
 
  StateType Type = ((State*)st->item)->type;

  return(Type >= singleG && Type < internal0) || Type == rsingleG || (Type >= rterminal0 && Type <=rterminal2);
}


bool AugustusGraph::exonAtCodingEnd(Node *st){

  StateType Type = ((State*)st->item)->type;

  return(Type == singleG || Type == terminal || Type == rsingleG || Type == rinitial);
}

// generates an identification key
string AugustusGraph::getKey(Node *n){

  if(n->item == NULL)
    return (itoa(n->begin) +  ":neutral");
  else
    return (itoa(n->begin) + ":" + itoa(n->end) + ":" + itoa( (int)((State*)n->item)->type )); 
}

string AugustusGraph::getKey(Status *st){

  return (itoa(st->begin) + ":" + itoa(st->end) + ":" + itoa( (int)((State*)st->item)->type )); 
}

string AugustusGraph::getKey(State *st){

  return (itoa(st->begin) + ":" + itoa(st->end) + ":" + itoa( (int)(st->type)) ); 
}

string AugustusGraph::getKey(ExonCandidate* exoncand){
  
  return (itoa(exoncand->begin) + ":" + itoa(exoncand->end) + ":" + itoa( ((ExonCandidate*)exoncand)->getStateType() )); 

}
double AugustusGraph::getIntronScore(Status *predExon, Status *nextExon){

  if(nextExon != NULL){
    for(list<Status>::iterator st=statelist->begin(); st!=statelist->end(); st++){
      if(&(*st) == nextExon && &(*st) != &statelist->front()){
	if((--st)->next != NULL)
	  return setScore(&(*st));
	else
	  return 0;
      }
    }
  }
  else if(predExon != NULL){
    for(list<Status>::iterator st=statelist->begin(); st!=statelist->end(); st++){
      if(&(*st) == predExon && st->next != NULL){
	st++;
	return setScore(&(*st));
      }
    }
  }
  return 0;
}

//for incomplete features at the beginning of a sequence 
void AugustusGraph::addEdgeFromHead(Status *exon){

  if(exon != NULL){
    if(exon->name >= CDS && exon->name <= utr5){

      if(!edgeExists(head,getNode(exon)) && !exonAtGeneStart(exon)){
     
	Edge toExon(getNode(exon),false,getIntronScore(NULL,exon));
	head->edges.push_back(toExon);	
      }
    }
    else if(exon->next != NULL){
      addEdgeFromHead(exon->next);
    }
  }
}

//for incomplete features at the end of a sequence
void AugustusGraph::addEdgeToTail(Status *exon){
  
  Node *ex = getNode(exon);
 
  if(!edgeExists(ex,tail)){
    if(!exonAtGeneEnd(exon))
      {
	Edge toTail(tail,false,getIntronScore(exon,NULL));
	ex->edges.push_back(toTail);
      }
  }
}

void Graph::insertIntron(Node *exon1, Node *exon2){

  // insert neutral edge 
  if(exon1->item==NULL && exon2->item==NULL){    
    Edge intr(exon2);
    exon1->edges.push_back(intr);
  }
  //insert compatible edge
  else if(!mergedStopcodon(exon1, exon2)){
    Edge intr(exon2, false);
    exon1->edges.push_back(intr);
  }
}

bool AugustusGraph::compatible(Node *exon1, Node *exon2){

  if(exon1->item != NULL && exon2->item != NULL){
    StateType type1 = ((State*)exon1->item)->type; 
    StateType type2 = ((State*)exon2->item)->type;

    return(isCodingExon(type1) && isCodingExon(type2) && exon1->end < exon2->begin && sameStrand(type1,type2) && sameReadingFrame(exon1,exon2) && !exonAtCodingEnd(exon1) && !exonAtCodingStart(exon2))

      ||

      (exon1->end == exon2->begin-1 && sameStrand(type1,type2) && (((type1 == utr5single || type1 == utr5term) && type2 >= singleG && type2 < internal0) || ((type1 == rutr3init || type1 == rutr3single) && ((type2 > rinternal2 && type2 <= rterminal2) || type2 ==rsingleG)) || ((type1 == singleG || type1 == terminal) && (type2 == utr3single || type2 == utr3init)) || ((type1 == rsingleG || type1 == rinitial) && (type2 == rutr5single || type2 == rutr5term)) ))

      ||

      (((is3UTRExon(type1) && is3UTRExon(type2) &&
 (((type1 == utr3init || type1 == utr3internal) && (type2 == utr3internal || type2 == utr3term)) ||
 ((type1 == rutr3term || type1 == rutr3internal) && (type2 == rutr3internal || type2 == rutr3init)))) 
||
(is5UTRExon(type1) && is5UTRExon(type2) &&
 (((type1 == utr5init || type1 == utr5internal) && (type2 == utr5internal || type2 == utr5term)) ||
 ((type1 == rutr5term || type1 == rutr5internal) && (type2 == rutr5internal || type2 == rutr5init))))) 
&& sameStrand(type1,type2) && exon1->end < exon2->begin);

  }
  else 
    return false;
}

bool AugustusGraph::sameStrand(StateType typeA, StateType typeB){
  return(((typeA < rsingleG && typeA >= singleG) && (typeB < rsingleG && typeB >= singleG)) || ((typeA >= rsingleG && typeA < intron_type) && (typeB >= rsingleG && typeB < intron_type)));
}

bool AugustusGraph::sameReadingFrame(Node *e1, Node *e2){

  /*
   * find framenumber of startnucleotide by using length and frame of e2
   * test if frame of e1 is equal to that framenumber 
   */
  
  StateType typeA = ((State*)e1->item)->type;
StateType typeB = ((State*)e2->item)->type;

  if(isOnFStrand(typeA) && isOnFStrand(typeB))
    return(((State*)e1->item)->frame() == mod3(((State*)e2->item)->frame() - (((State*)e2->item)->end-((State*)e2->item)->begin+1)%3));
  else
    return(((State*)e1->item)->frame() == mod3(((State*)e2->item)->frame() + (((State*)e2->item)->end-((State*)e2->item)->begin+1)%3));

}

void AugustusGraph::calculateBaseScores(){

  map<string,Status*> exonProcessed;
  for(list<Status>::iterator st=statelist->begin(); st!=statelist->end(); st++){
    if(exonProcessed[getKey(&(*st))] == 0){
      exonProcessed[getKey(&(*st))] = &(*st);
      for(int pos=st->begin; pos<=st->end; pos++)
	if(getBasetype(&(*st), pos)>=0){
	  baseScore[getBasetype(&(*st), pos)*seqlength + pos] += st->score;
	}
    }
  } 
}

/*
 * calculate scores for exons and introns.
 * the scoring function is peacewise linear with 5 partially fixed points
 * in getPoints the 2 Points used for the linear scoring function corresponding to the aposteriori probability of that state are identified.
 */

double AugustusGraph::setScore(Status *st){
  
  if(Constant::logreg){

    if(st->name >= CDS && st->name < intron){
      return (   Constant::lg_es[0]
	       + Constant::lg_es[1] * st->score
	       + Constant::lg_es[2] * getAvgBaseProb(st)
	       + Constant::lg_es[3] * log(st->end - st->begin + 1) );
    }else{
      return (   Constant::in_sc[0]
	       + Constant::in_sc[1] * st->score
	       + Constant::in_sc[2] * getAvgBaseProb(st)
	       + Constant::in_sc[3] * log(st->end - st->begin + 1) );
    }
  }else{

    double a1 = 0;
    double a2 = 0;
    double b1 = 0;
    double b2 = 0;
  
    if(st->name >= CDS && st->name < intron){
      double s_se = 0;
      getPoints(st,st->score,&a1,&a2,&b1,&b2);
      
      s_se = (st->score - a1) * (b2-a2)/(b1-a1) + a2;

      double s_be = 0;
      double p_b = 0;
    
      for(int pos = st->begin; pos<=st->end; pos++){
	if(getBasetype(st, pos)>=0){
	  p_b = baseScore[getBasetype(st, pos)*seqlength + pos];
	  getPoints(st,p_b,&a1,&a2,&b1,&b2);	  
	}
	if(Constant::MultSpeciesMode){
	  s_be += p_b - r_be;
	}
	else{
	  s_be += (p_b - a1) * (b2-a2)/(b1-a1) + a2;
	}
      }
      s_be /= st->end - st->begin + 1;
    
      // cout<<"exon\tlength: "<<st->end-st->begin+1<<"\tapostprob: "<<st->score<<"\tstate score: "<<s_se<<"\tscore: "<<alpha_e * s_se + s_be<<"\tpoints: ("<<a1<<","<<a2<<")\t("<<b1<<","<<b2<<")"<<endl;
      return alpha_e * s_se + s_be;
    }
    else{
      double s_bi = 0;
      double p_b = 0;
      for(int pos = st->begin; pos<=st->end; pos++){
	if(getBasetype(st, pos)>=0){
	  p_b = baseScore[getBasetype(st, pos)*seqlength + pos];
	  getPoints(st,p_b,&a1,&a2,&b1,&b2);
	}
	if(Constant::MultSpeciesMode){
	  s_bi += p_b - r_bi;
	}
	else{
	  s_bi += (p_b - a1) * (b2-a2)/(b1-a1) + a2;
	}
      }
      s_bi /= st->end - st->begin + 1;
      double s_si = 0;
      getPoints(st,st->score,&a1,&a2,&b1,&b2);
      s_si = (st->score - a1) * (b2-a2)/(b1-a1) + a2;
       
      // cout<<"intron\tlength: "<<st->end-st->begin+1<<"\tapostprob: "<<st->score<<"\tstate score: "<<s_si<<"\tscore: "<<alpha_i * s_si + s_bi<<"\tpoints: ("<<a1<<","<<a2<<")\t("<<b1<<","<<b2<<")"<<endl;
      return alpha_i * s_si + s_bi;
    }
  }
}

int AugustusGraph::getBasetype(Status *st, int pos){

  StateType Type = ((State*)st->item)->type;

  if(st->name==CDS){
    int frameAtPos = mod3(((State*)st->item)->frame() - (st->end-pos+1)%3);
    if(Type >= singleG && Type <= terminal){

      switch(frameAtPos){
      case 0:
	return 0; //CDS+0
      case 1:
	return 2; //CDS+1
      case 2:
	return 4; //CDS+2
      default:
	return -1;
      }
    }
    else{
      switch(frameAtPos){
      case 0:
	return 1; //CDS-0
      case 1:
	return 3; //CDS-1
      case 2:
	return 5; //CDS-2
      default:
	return -1;
      }
    }
  }
  // TODO: add 4 more cases: ncExon+ ncExon- ncIntron+ ncIntron-
  else if(st->name==utr3 || st->name==utr5){
    if(Type >= utr5single && Type <= utr3term)
      return 6; //UTR+
    else if(Type >= rutr5single && Type <= rutr3term)
      return 7; //UTR-
    else
      return -1;
  }
  
  else if(st->name >= intron){
    if((Type >= singleG && Type < rsingleG) || Type == intron_type )
      return 8; //intron+
    else
      return 9; //intron-
  }
  else
    return -1;
}

void AugustusGraph::getPoints(Status *st, double p, double *a1, double *a2, double *b1, double *b2){

  if(st->name >= CDS && st->name < intron){
    if(p < i1_e){
      *a1 = 0; *a2 = x0_e;
      *b1 = i1_e; *b2 = j1_e;
    } else if(p >= i1_e && p < y0_e){
      *a1 = i1_e; *a2 = j1_e;
      *b1 = y0_e; *b2 = 0;
    } else if(p >= y0_e && p < i2_e){
      *a1 = y0_e; *a2 = 0;
      *b1 = i2_e; *b2 = j2_e;
    } else {
      *a1 = i2_e; *a2 = j2_e;
      *b1 = 1; *b2 = x1_e;
    }
  } else {
    if(p < i1_i){
      *a1 = 0; *a2 = x0_i;
      *b1 = i1_i; *b2 = j1_i;
    } else if(p >= i1_i && p < y0_i){
      *a1 = i1_i; *a2 = j1_i;
      *b1 = y0_i; *b2 = 0;
    } else if(p >= y0_i && p < i2_i){
      *a1 = y0_i; *a2 = 0;
      *b1 = i2_i; *b2 = j2_i;
    } else {
      *a1 = i2_i; *a2 = j2_i;
      *b1 = 1; *b2 = x1_i;
    }
  }
}

/*
 * creates an input file for graphviz
 */

void AugustusGraph::printGraph2(string filename){

  ofstream file;
  try{
    file.open((filename).c_str());
  } catch (...) {
    cerr<<"AugustusGraph::printGraph() can't open file"<<endl;
    return;
  }
  file<<"digraph MEAgraph {\n";
  file<<"rankdir=LR;\n";

  file<<"\tnode[shape=box];\n";

  head->begin = 0;
  head->end = 0;

  for(list<Node*>::iterator pos = nodelist.begin(); pos != nodelist.end(); pos++){
    for(list<Edge>::iterator it = (*pos)->edges.begin(); it != (*pos)->edges.end(); it++){
      string name1 = "exon";
      string name2 = "exon";
      StateType type1;
      StateType type2;

      if((*pos)->item != NULL)
	type1 = ((State*)(*pos)->item)->type;
     else
       type1 = TYPE_UNKNOWN;
     if(it->to->item != NULL)
       type2 = ((State*)it->to->item)->type;
     else
       type2 = TYPE_UNKNOWN;

     if((*pos)->item != NULL)
	name1 = name1 + itoa(type1) + "_";
     name1 = name1 + itoa((*pos)->begin) + "_" + itoa((*pos)->end);
      if(it->to->item != NULL)
	name2 = name2 + itoa(type2) + "_";
      name2 = name2 + itoa(it->to->begin) + "_" + itoa(it->to->end);
     
      file<<name1<<"[";
      if((*pos)->begin == it->to->begin)
	file<<"label="<<(*pos)->begin<<",style=filled,fillcolor=yellow,";
      if((*pos)->end == it->to->end)
	file<<"style=filled,fillcolor=turquoise,";
      if((type1 >= singleG && type1 <=terminal) || (type1 >= rsingleG && type1 <= rterminal2))
	file<<"shape=ellipse,";
      file<<"];\n";

      file<<name2<<"[";
      if((*pos)->begin == it->to->begin)
	file<<"style=filled,fillcolor=yellow,";
      if((*pos)->end == it->to->end)
	file<<"label="<<it->to->end<<",style=filled,fillcolor=turquoise,";
      if((type2 >= singleG && type2 <=terminal) || (type2 >= rsingleG && type2 <= rterminal2))
	file<<"shape=ellipse,";
      file<<"];\n";
     
      
      if((*pos)==head)
	file<<name1<<"[style=filled,label=head];\n";
      if(it->to==tail)
	file<<name2<<"[style=filled,label=tail];\n";

      if((*pos)->begin == it->to->begin)
	file<<"{ rank=same; "<<name1<<";"<<name2<<";}\n";
      if((*pos)->end == it->to->end)
	file<<"{ rank=same; "<<name1<<";"<<name2<<";}\n";
      
     
      file<<name1<<"->"<<name2<<"[";
      if((*pos)->begin == (*pos)->end && it->to->begin == it->to->end){
	if((*pos)->begin <= it->to->begin){
	  if(it->to->pred == (*pos))
	    file<<"weight=100,color=blue,";
	  else
	    file<<"weight=100,color=red,";
	}
	else{
	  if(it->to->pred == (*pos))
	    file<<"weight=-100,color=blue,";
	  else
	    file<<"weight=-100,color=green,";
	}
      }
      else if(it->to->pred == (*pos))
	file<<"color=blue,";
      file<<"label=\""<<it->score<<"\"];\n";  
    }
  }
  file<<"}\n";
  file.close();

  head->begin = -1;
  head->end = -1;
}

void AugustusGraph::printGraph(string filename){

  map<string,Node*> inQueue;
  queue<Node*> q;
  head->begin = 0;
  head->end = 0;
  q.push(head);
  ofstream file;
  try{
    file.open((filename).c_str());
  } catch (...) {
    cerr<<"AugustusGraph::printGraph() can't open file"<<endl;
    return;
  }
  file<<"digraph MEAgraph {\n";
  file<<"rankdir=LR;\n";
 
  file<<"\tnode[shape=box];\n";
  while(!q.empty()){
    Node *pos = q.front();
    q.pop(); 
   
    for(list<Edge>::iterator it=pos->edges.begin(); it!=pos->edges.end(); it++){
     if(inQueue[getKey(it->to)] == 0){
	q.push(it->to);
	inQueue[getKey(it->to)] = it->to;
      }
     string name1 = "exon";
     string name2 = "exon";
     StateType type1;
     StateType type2;

     if(pos->item != NULL)
        type1 = ((State*)pos->item)->type;
     else
       type1 = TYPE_UNKNOWN;
     if(it->to->item != NULL)
       type2 = ((State*)it->to->item)->type;
     else
       type2 = TYPE_UNKNOWN;

      if(pos->item != NULL)
	name1 = name1 + itoa(type1) + "_";
      name1 = name1 + itoa(pos->begin) + "_" + itoa(pos->end);
      if(it->to->item != NULL)
	name2 = name2 + itoa(type2) + "_";
      name2 = name2 + itoa(it->to->begin) + "_" + itoa(it->to->end);
      
      file<<name1<<"[";
      if(pos->begin == it->to->begin)
	file<<"label="<<pos->begin<<",style=filled,fillcolor=yellow,";
      if(pos->end == it->to->end)
	file<<"style=filled,fillcolor=turquoise,";
      if((type1 >= singleG && type1 <=terminal) || (type1 >= rsingleG && type1 <= rterminal2))
	file<<"shape=ellipse,";
      file<<"];\n";

      file<<name2<<"[";
      if(pos->begin == it->to->begin)
	file<<"style=filled,fillcolor=yellow,";
      if(pos->end == it->to->end)
	file<<"label="<<it->to->end<<",style=filled,fillcolor=turquoise,";
      if((type2 >= singleG && type2 <=terminal) || (type2 >= rsingleG && type2 <= rterminal2))
	file<<"shape=ellipse,";
      file<<"];\n";
     
      
      if(pos==head)
	file<<name1<<"[style=filled,label=head];\n";
      if(it->to==tail)
	file<<name2<<"[style=filled,label=tail];\n";

      if(pos->begin == it->to->begin)
	file<<"{ rank=same; "<<name1<<";"<<name2<<";}\n";
      if(pos->end == it->to->end)
	file<<"{ rank=same; "<<name1<<";"<<name2<<";}\n";
      
     
      file<<name1<<"->"<<name2<<"[";
      if(pos->begin == pos->end && it->to->begin == it->to->end){
	if(pos->begin <= it->to->begin){
	   if(it->to->pred == pos)
	     file<<"weight=100,color=blue,";
	   else
	     file<<"weight=100,color=red,";
	}
	else{
	  if(it->to->pred == pos)
	     file<<"weight=-100,color=blue,";
	  else
	    file<<"weight=-100,color=green,";
	}
      }
      else if(it->to->pred == pos)
	file<<"color=blue,";
      file<<"label=\""<<it->score<<"\"];\n";  
      
    }
  }
  file<<"}\n";
  file.close();

  head->begin = -1;
  head->end = -1; 
}


bool compareNodes(Node *first, Node *second){
  return(first->begin <= second->begin);
}

bool compareEdges(Edge first, Edge second){
  return(first.to->begin <= second.to->begin);
}

bool compareNodeEnds(const Node *first, const Node *second)
{
  return (first->end <= second->end);
}

//casts node to StateType
StateType Node::castToStateType(){

  if(n_type == sampled || n_type == utrExon){
    return ((State*)item)->type;
  }
  else if(n_type == unsampled_exon){
    return ((ExonCandidate*)item)->getStateType();
  } 
  return TYPE_UNKNOWN;

}

string stateNameIdentifiers[NUM_STATENAMES]={
    "CDS","3'UTR","5'UTR","intron","intron","intron"
};

string nodeTypeIdentifiers[NUM_NODETYPES]=
  {"IR", "plus0", "plus1", "plus2", "minus0", "minus1", "minus2", "T_plus1", "TA_plus2", "TG_plus2", "T_minus1", "C_minus1", "YY_minus0",
   "ncintr", "rncintr",
   "utr5intr", "TLstart", "TLstop", "utr3intr", "rutr5intr", "rTLstart", "rTLstop", "rutr3intr", "utr",
   "sampled", "unsampled_exon"};

//print function for nodes and edges
ostream& operator<<(ostream& ostrm, Node *node){

    if(node->n_type >= IR && node->n_type < utrExon){
	ostrm << node->begin << "\t" << node->end << "\t" << nodeTypeIdentifiers[node->n_type]<< "\t";
    }
    else if(node->n_type >= utrExon){
	ostrm << node->begin << "\t" << node->end << "\t" << stateTypeIdentifiers[node->castToStateType()]<< "\t";
    }
    else if( node->begin == -1)
	ostrm << "head";
    else{
	ostrm << "tail";
    }
    return ostrm;
}

ostream& operator<<(ostream& ostrm, const Edge &edge){
  
  ostrm << (edge.to) << "\t" << edge.score;
  return ostrm;
}


bool AugustusGraph::mergedStopcodon(Node* exon1, Node* exon2){
    return mergedStopcodon(exon1->castToStateType(), exon2->castToStateType(), exon1->end, exon2->begin);
}

bool AugustusGraph::mergedStopcodon(Status* exon1, Status* exon2){
    if(exon1 && exon2)
	return mergedStopcodon(((State*)exon1->item)->type, ((State*)exon2->item)->type, exon1->end, exon2->begin);
    return false;
}

bool AugustusGraph::mergedStopcodon(StateType type1, StateType type2, int end1, int begin2){

    char joinedCodon[4] = "";

    if(isCodingExon(type1) && isCodingExon(type2)){
	if(type1 == initial1 || type1 == internal1 || type1 == rterminal1 || type1 == rinternal1 ){
	    strncat(joinedCodon, sequence + end1, 1);
	    strncat(joinedCodon, sequence + begin2, 2);
	}
	else if(type1 == initial2 || type1 == internal2 ||  type1 == rterminal0 || type1 == rinternal0 ){
	    strncat(joinedCodon, sequence + end1 - 1, 2);
	    strncat(joinedCodon, sequence + begin2, 1);
	}
	if(joinedCodon[0] != '\0'){
	    if (isOnFStrand(type1) && GeneticCode::isStopcodon(joinedCodon) ){
		return true;
	    }
	    if( !isOnFStrand(type1) && GeneticCode::isRCStopcodon(joinedCodon) ){
		return true;
	    }
	}
    }
    return false;
}

void Node::addWeight(float weight){

    for(list<Edge>::iterator edge = edges.begin(); edge != edges.end(); edge++){
	edge->score += weight;
    }
}

Edge* Node::getEdge(Node* succ){
    Edge* e = NULL;
    for(list<Edge>::iterator it=edges.begin(); it!=edges.end(); it++){
	if(it->to == succ){
	    e = &(*it);     
	    break;
	}
    }
    return e;
}

State* Node::getIntron(Node* succ){

    Edge* edge = getEdge(succ);
    State *intron = NULL;
    if(edge && edge->item != NULL){
	intron = new State(*((State*)edge->item));
    }
    return intron;
}

bool isTlstartOrstop(Status *predExon, Status *succExon){
    if(predExon && succExon){
	if( (predExon->isCDS() && succExon->isUTR()) || (succExon->isCDS() && predExon->isUTR()) )
	    return true;
    }
    return false;
}

float AugustusGraph::getAvgBaseProb(Status *st){
    if(st->isExon() || st->isIntron()){
	float prob = 0.0;
	for(int pos = st->begin; pos<=st->end; pos++){
	    if(getBasetype(st, pos)>=0){
		prob += baseScore[getBasetype(st, pos)*seqlength + pos];
	    }
	}
	return prob /= st->getLen();
    }
    return 0.0;
}
