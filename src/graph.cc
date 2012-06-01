#include <iostream>
#include <limits>
#include <list>
#include <fstream>
#include <sstream>
#include<iomanip>

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
    Node *ex = new Node(exon->begin, exon->end, setScore(exon), exon->item);
    nodelist.push_back(ex);
    addToHash(ex);

    if(exonAtGeneStart(exon)){
      // include edge from neutral line to exon

      Node *neut = new Node(ex->begin, ex->begin);
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

      Node *neut = new Node(ex->end, ex->end);
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
  else if(exon1->next != exon2 && !edgeExists(e1,e2)){
    Edge in(e2, false, setScore(exon1->next), exon1->next->item);
    e1->edges.push_back(in);  
  }
}

void Graph::createNeutralLine(vector<Node*> &neutralLine){

  Node *pos = head;
  int n = neutralLine.size();
  for(int i=0; i<n; i++){
    if(neutralLine[i] != NULL){
      Edge neut(neutralLine[i]);
      pos->edges.push_back(neut);
      pos = neutralLine[i];
    }
  }
  Edge lastEdge(tail);
  pos->edges.push_back(lastEdge);
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

  for(list<Node*>::iterator fromNeut=neutralNodes.begin(); fromNeut!=neutralNodes.end(); fromNeut++){  
    if(nonneutralIncomingEdge(*fromNeut)){
      bool nonNeutralLoop = false;
      for(list<Node*>::iterator toNeut=fromNeut; toNeut!=neutralNodes.begin(); toNeut--){   
	if(nonNeutralLoop)
	  break;
	nonNeutralLoop = true;
	queue<Node*> q;
	inQueue.clear();
	int nrNonNeutralEdges = 0;
	for(list<Edge>::iterator edge=(*toNeut)->edges.begin(); edge!=(*toNeut)->edges.end(); edge++)
	  if(!edge->neutral)
	    nrNonNeutralEdges++;
	
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
	      if(pos->item != NULL && edge->to->item == NULL && edge->to->begin <= (*fromNeut)->begin)
		goto nextEdge;
	      
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
  cout<<"-------------------------------------------------------"<<endl<<setw(10)<<"Node:"<<setw(10)<<(*node)->begin<<setw(10)<<(*node)->end<<setw(10)<<(*node)->score;
    if((*node)->item == NULL)
      cout<<setw(10)<<"neutral";
    for(list<Edge>::iterator edge = (*node)->edges.begin(); edge != (*node)->edges.end(); edge++){
      cout<<endl<<setw(10)<<"----->"<<setw(10)<<edge->to->begin<<setw(10)<<edge->to->end<<setw(10)<<edge->to->score;
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
  else{
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

double AugustusGraph::setScore(Status *st){

  if(st->name >= CDS && st->name < intron){
    double s_be = 0;
    for(int pos = st->begin; pos<=st->end; pos++)
      if(getBasetype(st, pos)>=0)
	s_be += baseScore[getBasetype(st, pos)*seqlength + pos] - r_be;
    return alpha_se * (st->score - r_se) + alpha_be * s_be;
  }
  else{
    double s_bi = 0;
    for(int pos = st->begin; pos<=st->end; pos++)
      if(getBasetype(st, pos)>=0)
	s_bi += baseScore[getBasetype(st, pos)*seqlength + pos] - r_bi;
    return alpha_si * (st->score - r_si) + alpha_bi * s_bi;
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


/*
 * creates an input file for graphviz
 */

void AugustusGraph::printGraph(string filename){

  map<string,Node*> inQueue;
  queue<Node*> q;
  head->begin = 0;
  head->end = 0;
  q.push(head);
  ofstream file;
  try{
    file.open(("/home/Lizzy/augustus_output/" + filename).c_str());
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

/*********************************************************************
 * 
 * date     |   author      |  added
 * ---------|---------------|-----------------------------------------
 * 16.01.12 |Stefanie König |  functions to construct graph
 *                          |  with seven neutral lines
*********************************************************************/

/*
 * builds Graph with seven neutral lines representing
 * the intergenetic region
 * the plus strand (3 neutral lines for each coding frame)
 * and the minus strand (3 neutral lines for each coding frame)
 */

string nodeTypeIdentifiers[NUM_NODETYPES]=
  {"IR", "plus0", "plus1", "plus2", "minus0", "minus1", "minus2", "sampled", "unsampled_exon"};

void SpeciesGraph::buildGraph(){

  vector< vector<Node*> > neutralLines; //represents the seven neutral lines
  
  getSizeNeutralLine();

  /*
   * extend NeutralLine in case that additional Exons are beyond the scope
   */
  if(!additionalExons.empty()){
    for(list<ExonCandidate*>::iterator it=additionalExons.begin(); it!=additionalExons.end(); it++){
      if((*it)->end > max)
	max = (*it)->end;	  
      if((*it)->begin < min)
	min = (*it)->begin;	
    }
  }

  vector<Node*> intergenetic;   
  neutralLines.push_back(intergenetic);
  vector<Node*> plus0;
  neutralLines.push_back(plus0);
  vector<Node*> plus1;
  neutralLines.push_back(plus1);
  vector<Node*> plus2;
  neutralLines.push_back(plus2);
  vector<Node*> minus0;
  neutralLines.push_back(minus0);
  vector<Node*> minus1;
  neutralLines.push_back(minus1);
  vector<Node*> minus2;
  neutralLines.push_back(minus2);

  for(int i=0; i<max-min+1; i++){
    for(int j=0; j<neutralLines.size(); j++){
      neutralLines.at(j).push_back(NULL);
    }
  }

  head = new Node(-1,-1); // initialize nodelist of the graph with head and tail
  nodelist.push_back(head);
  tail = new Node(max+1,max+1);
  nodelist.push_back(tail);

  calculateBaseScores();
  
  // add all sampled states to Graph (Augusutus information)

  bool gene_end = true;

  for(list<Status>::iterator it=statelist->begin(); it!=statelist->end(); it++){

    if(gene_end == true && it->name == intron){ //start -> intron  (transcript incomplete at the start)
      Node* first = addExon(it->next, neutralLines);
      addIntron(head, first, &(*it));
    }
    if(it->name == CDS){ // CDS -> CDS   or   CDS -> end
      gene_end = false;
      Node* e1 = addExon(&(*it), neutralLines);
      if(it->next != NULL){

	if(it->next->name == intron && it->next->next == NULL){ //CDS -> intron -> end (transcript incomplete at the end)
	  addIntron(e1, tail, it->next);
	}

	if(it->next != NULL && it->next->name == intron && it->next->next != NULL){ //CDS -> intron -> CDS
          Node* e2 = addExon(it->next->next, neutralLines);
	  addIntron(e1, e2, it->next);
	}
      }
    }
    if(it->next == NULL) // end of transcript
      gene_end = true;
  }

  //add additional Exoncandidates
  if(!additionalExons.empty()){
    for(list<ExonCandidate*>::iterator it = additionalExons.begin(); it!=additionalExons.end(); it++){
      addExon(*it, neutralLines);
    }
  }
  //create neutral lines by connecting neutral nodes in the vector (all entries in the Vector, which are not NULL)
  //edges directed from smaller positions to larger
  for(int j=0; j<neutralLines.size(); j++){
    createNeutralLine(neutralLines.at(j));
  }

  //add Node weight to edge weight and initialize node distances for DAG longest path algorithm
  addWeightToEdge();
  
  //find longest path in the graph
  topSort();

  //relax all nodes in topological order
  relax(head, tail);

  setPathLabels(head, tail, 1);

  printGraph(speciesname);

}

Node* SpeciesGraph::addExon(Status *exon, vector< vector<Node*> > &neutralLines){

  if(!alreadyProcessed(exon)){
    Node *ex = new Node(exon->begin, exon->end, setScore(exon), exon->item, sampled);
    nodelist.push_back(ex);
    addToHash(ex);
    addNeutralNodes(ex, neutralLines);    
    return ex;
  }
  return getNode(exon);
}

Node* SpeciesGraph::addExon(ExonCandidate *exon, vector< vector<Node*> > &neutralLines){

  if(!alreadyProcessed(exon)){
    Node *ex = new Node(exon->begin, exon->end, exon->score, exon, unsampled_exon);
    nodelist.push_back(ex);
    addToHash(ex);
    addNeutralNodes(ex, neutralLines);
  }
  return getNode(exon);
}

void SpeciesGraph::addNeutralNodes(Node *node,vector< vector<Node*> > &neutralLines){

  Node *neut_to = new Node(node->begin, node->begin, 0.0, NULL, (NodeType)fromNeutralLine(node) );
  Edge intron(node, false);
  if(!alreadyProcessed(neut_to)){
    neutralLines.at(fromNeutralLine(node)).at(node->begin-min) = neut_to;
    neut_to->edges.push_back(intron);
    nodelist.push_back(neut_to);
    addToHash(neut_to);
  }
  else{
    getNode(neut_to)->edges.push_back(intron);
    delete neut_to;
  }
  
  Node *neut_from = new Node(node->end, node->end, 0.0, NULL, (NodeType)toNeutralLine(node));
  if(!alreadyProcessed(neut_from)){
    neutralLines.at(toNeutralLine(node)).at(node->end-min) = neut_from;
    Edge intron(neut_from, false);
    node->edges.push_back(intron);
    nodelist.push_back(neut_from);
    addToHash(neut_from);
  }
  else{
    Edge intron(getNode(neut_from), false);
    node->edges.push_back(intron);      
    delete neut_from;
  }
}

void SpeciesGraph::addIntron(Node* exon1, Node* exon2, Status *intr){

#ifdef DEBUG
  //cout << "intron: "<< intr->begin <<"-"<<intr->end<<"\ttype: " << ((State*)intr->item)->type <<endl;
#endif

  if( !edgeExists(exon1,exon2) ){
    Edge in(exon2, false, setScore(intr), intr->item);
    exon1->edges.push_back(in);  
  }
}

/* returns
*  0 if edge: neutral Line IR --> exon start
*  1 if edge: neutral Line 0+ --> exon start
*  2 if edge: neutral Line 1+ --> exon start
*  3 if edge: neutral Line 2+ --> exon start
*  4 if edge: neutral Line 0- --> exon start
*  5 if edge: neutral Line 1- --> exon start
*  6 if edge: neutral Line 2- --> exon start
*/
int SpeciesGraph::fromNeutralLine(Node *node){

  
  StateType type = node->castToStateType();
  int frame = mod3(stateReadingFrames[type]);

  if( type == singleG || type == initial0 || type == initial1 || type == initial2 || type == rsingleG || type == rterminal0 || type == rterminal1 || type == rterminal2 ){
    return 0;
  }
  if(type == internal0 || type == internal1 || type == internal2 || type == terminal){
    return (mod3(frame - mod3(node->end - node->begin + 1)) + 1);
  
  }
  if(type == rinternal0 || type == rinternal1 || type == rinternal2 || type == rinitial){
    return (mod3(frame + mod3(node->end - node->begin + 1)) + 4);
  }
  return -1;
}

/* returns
*  0 if edge: exon end --> neutral Line IR 
*  1 if edge: exon end --> neutral Line 0+
*  2 if edge: exon end --> neutral Line 1+
*  3 if edge: exon end --> neutral Line 2+
*  4 if edge: exon end --> neutral Line 0-
*  5 if edge: exon end --> neutral Line 1-
*  6 if edge: exon end --> neutral Line 2-
*/
int SpeciesGraph::toNeutralLine(Node *node){

  StateType type = node->castToStateType();
  int frame = mod3(stateReadingFrames[type]);

  if( type == singleG || type == terminal || type == rsingleG || type == rinitial ){
    return 0;
  }
  if(type == internal0 || type == internal1 || type == internal2 || type == initial0 || type == initial1 || type == initial2){
    return (frame + 1);
  }
  if(type == rinternal0 || type == rinternal1 || type == rinternal2 || type == rterminal0 || type == rterminal1 || type == rterminal2){
    return (frame + 4);
  }
  return -1;
}

StateType Node::castToStateType(){

  if(n_type == sampled){
    return ((State*)item)->type;
  }
  else if(n_type == unsampled_exon){
    return ((ExonCandidate*)item)->getStateType();
  } 
  return TYPE_UNKNOWN;

}

void SpeciesGraph::printGraph(string filename){

  //creates inputfile for graphviz
  map<string,Node*> inQueue;
  queue<Node*> q;
  q.push(head);
  ofstream file;
  file.open(("/home/stefanie/augustus_output/" + filename).c_str());
  
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
       type1 = pos->castToStateType();
     else
       type1 = TYPE_UNKNOWN;
     if(it->to->item != NULL)
       type2 = it->to->castToStateType();
     else
       type2 = TYPE_UNKNOWN;

      if(pos->item != NULL)
	name1 = name1 + itoa(type1) + "_";
      name1 = name1 + itoa(pos->begin+1) + "_" + itoa(pos->end+1);
      if(it->to->item != NULL)
	name2 = name2 + itoa(type2) + "_";
      name2 = name2 + itoa(it->to->begin+1) + "_" + itoa(it->to->end+1);
      
      file<<name1<<"[";
      if(pos==head)
	file<<"style=filled,label=head";
      else if(pos->item == NULL)
	file<<"shape=point";
      else if(pos->item != NULL && it->to->n_type == IR)
	file<<"style=filled,fillcolor=turquoise,";
      else if(pos->n_type == unsampled_exon)
	file<<"style=filled,fillcolor=lightgrey,";
      file<<"];\n";

      file<<name2<<"[";
      if(it->to==tail)
	file<<"style=filled,label=tail";
      else if(it->to->item == NULL)
	file<<"shape=point";
      else if(pos->n_type == IR && it->to != NULL)
	file<<"style=filled,fillcolor=yellow,";
      else if(it->to->n_type == unsampled_exon)
	file<<"style=filled,fillcolor=lightgrey,";
      file<<"];\n";
      
      if( !(pos == head && it->to ==tail) ){
	file<<name1<<"->"<<name2<<"[";
	if(pos->begin == pos->end && it->to->begin == it->to->end){
	  if(pos == head || it->to == tail){
	    if(it->to->pred == pos && it->to->label == 1 && pos->label == 1)
	      file<<"color=blue,";
	    else
	      file<<"color=red,";
	  }
	  else{
	    if(it->to->pred == pos && it->to->label == 1 && pos->label == 1)
	      file<<"weight=100,color=blue,";
	    else
	      file<<"weight=100,color=red,";
	  }
	}
	else if(it->to->pred == pos && it->to->label == 1 && pos->label == 1)
	  file<<"color=blue,";

	if(pos == head && it->to->n_type >= IR && it->to->n_type <= minus2)
	  file<<"label=" << nodeTypeIdentifiers[it->to->n_type] << "];\n";
	else if(it->to == tail && pos->n_type >= IR && pos->n_type <= minus2)
	  file<<"label=" << nodeTypeIdentifiers[pos->n_type] << "];\n";
	else if( pos->n_type >= IR && pos->n_type <= minus2  &&  it->to->n_type >= IR && it->to->n_type <= minus2){
	  file<<"label=" << nodeTypeIdentifiers[pos->n_type] << "];\n";
	}
	else
	  file<<"label="<<it->score<<"];\n";  
      }
    }
  }
  file<<"}\n";
  file.close();
}

void SpeciesGraph::topSort(){

  map<string,Node*> processed;

  for(list<Node*>::iterator node=nodelist.begin(); node!=nodelist.end(); node++){   
    if(processed[getKey(*node)]==0){
      dfs(*node, processed);
    }
  }
}

void SpeciesGraph::dfs(Node *node, map<string,Node*> &processed){

  static Node* temp = NULL;

  processed[getKey(node)] = node;
  for(list<Edge>::iterator edge=node->edges.begin(); edge!=node->edges.end(); edge++){
    if(processed[getKey(edge->to)] == 0)
      dfs(edge->to, processed);
  }
  if(temp){
    node->topSort_next = temp;
  }
  temp = node;
}

void SpeciesGraph::relax(Node *begin, Node *end){

  Node *next = begin;
  
  //initialize
  while(next != end){
    next = next->topSort_next;
    next->score =  - numeric_limits<double>::max();
  }
  begin->score = 0;
  next = begin;
  while(next != end){
    for(list<Edge>::iterator edge = next->edges.begin(); edge != next->edges.end(); edge++){	 
      if(next->score + edge->score > edge->to->score){
	if(isReachable(begin, next) && isReachable(edge->to, end)){
	  // update exon distance
	  edge->to->score = next->score + edge->score;
	  edge->to->pred = next;
	  // cout << *edge << endl;
	}
      }
    }
    next = next->topSort_next;
  }
}

void SpeciesGraph::setPathLabels(Node *begin, Node *end, bool label){

  while(end != begin){
    end->label = label;
    end = end->pred;
  }
  end->label = label;
}

void MoveObject::initLocalHeadandTail(size_t step_size){

 Node *left = NULL;
 Node *right = NULL;

 if(nodes.empty() && edges.empty()){
    throw ProjectError("MoveObject::setLocalHead(): no nodes and no edges specified in MoveObject");
  }
  else if(edges.empty()){
    left = nodes.front().node;
    right = nodes.back().node;
  }
  else if(nodes.empty()){
    left =  nodes.front().node;
    right = edges.back().edge->to;
  }
  else{
    if(compareNodes( edges.back().edge->to, nodes.back().node)){
      right = nodes.back().node;
    }
    else{
      right = edges.front().edge->to;
    }
    left = nodes.front().node;
  }
 setLocalTail(right);
 setLocalHead(left);

 if(local_head == local_tail){
   goLeftOnPath(step_size);
   goRightOnPath(step_size);
 }

    
}

void  MoveObject::setLocalHead(Node *node){
 
  while(node->label != 1){
    node = node->pred;
  }
  local_head = node;
}

void MoveObject::goLeftOnPath(size_t step_size){

  Node *node = local_head;
  
  while(step_size > 0 && node != graph->head){
    node = node->pred;
    step_size--;
  }
  local_head = node;
}

void MoveObject::goRightOnPath(size_t step_size){

  Node *node = local_tail;

  while(step_size > 0 && node != graph->tail){
    for(list<Edge>::iterator edge = node->edges.begin(); edge != node->edges.end(); edge++){
      if(edge->to->label == 1 &&edge->to->pred == node){
	node = edge->to;
	break;
      }
    }      
    step_size--;
  }
  local_tail = node;
}


void MoveObject::setLocalTail(Node *node){

  list<Node*> queue;
  queue.push_back(node);
  map<string, Node*> processed;
  processed[graph->getKey(node)] = node;
  while(!queue.empty()){
    Node *next = queue.front();
    queue.pop_front();
    if (next->label == 1){ 
      local_tail = next;
      break;
    }
    else{
      for(list<Edge>::iterator edge = next->edges.begin(); edge != next->edges.end(); edge++){      
	if (processed[graph->getKey(edge->to)] == 0){
	  processed[graph->getKey(edge->to)] = edge->to;
	  queue.push_back(edge->to);
	}
      }
    }
  }
}

bool SpeciesGraph::isReachable(Node *node, Node *target){

  if(target == tail || node == head){
    return true;
  }
  else{

    list<Node*> queue;
    queue.push_back(node);
    map<string, Node*> processed;
    processed[getKey(node)] = node;
    while(!queue.empty()){
      Node *next = queue.front();
      queue.pop_front();
      if (next == target){
	return true;
      }
      else{
	for(list<Edge>::iterator edge = next->edges.begin(); edge != next->edges.end(); edge++){      
	  if (processed[getKey(edge->to)] == 0 && edge->to->begin <= target->begin){
	    processed[getKey(edge->to)] = edge->to;
	    queue.push_back(edge->to);
	  }
	}
      }
    }
  }
  return false;
}

string SpeciesGraph::getKey(Node *n){

  if(n->item == NULL)
    return (itoa(n->begin) +  ":" + itoa((n->n_type)));
  else
    return (itoa(n->begin) + ":" + itoa(n->end) + ":" + itoa( (int)(n->castToStateType()) )); 
}

double SpeciesGraph::setScore(Status *st){

  if(st->name >= CDS && st->name < intron){
    double s_be = 0;
    for(int pos = st->begin; pos<=st->end; pos++)
      if(getBasetype(st, pos)>=0)
	s_be += baseScore[getBasetype(st, pos)*seqlength + pos] - r_be;
    updateMaxWeight(alpha_se * (st->score - r_se) + alpha_be * s_be);
    return alpha_se * (st->score - r_se) + alpha_be * s_be;
  }
  else{
    double s_bi = 0;
    for(int pos = st->begin; pos<=st->end; pos++)
      if(getBasetype(st, pos)>=0)
	s_bi += baseScore[getBasetype(st, pos)*seqlength + pos] - r_bi;
    updateMaxWeight(alpha_si * (st->score - r_si) + alpha_bi * s_bi);
    return alpha_si * (st->score - r_si) + alpha_bi * s_bi;
  }
}

double SpeciesGraph::localChange(MoveObject *move){

  double local_score = - getScorePath(move->getHead(), move->getTail());  //old score of local path

  cout << "\nSpecies: " << speciesname <<"\n\nbefore local change" << endl;
  printPath(move->getHead(), move->getTail());

  move->addWeights();

  setPathLabels(move->getHead(), move->getTail(), 0); // reset path labels locally  
  relax(move->getHead(), move->getTail());            // relax edges
  setPathLabels(move->getHead(), move->getTail(), 1); // local backtracking 

  move->undoAddWeights();

  local_score += getScorePath(move->getHead(), move->getTail());  // new score of local path

  cout << "\nafter local change" << endl;
  printPath(move->getHead(), move->getTail());

  printGraph(speciesname + ".opt");
  return local_score;
}

void SpeciesGraph::undoChange(MoveObject *move){

  setPathLabels(move->getHead(), move->getTail(), 0); // reset path labels locally  
  relax(move->getHead(), move->getTail());            // relax edges
  setPathLabels(move->getHead(), move->getTail(), 1); // local backtracking 

}

double SpeciesGraph::getScorePath(Node *begin, Node *end){

  double score = 0;

  while(begin != end){
 
    for(list<Edge>::iterator edge = begin->edges.begin(); edge != begin->edges.end(); edge++){
      if(edge->to->label == 1 && edge->to->pred == begin){
	score += edge->score;
	begin = edge->to;
	break;
      }
    }   
  }
  return score;
}

void SpeciesGraph::printPath(Node *begin, Node *end){

   double score = 0;

   cout << begin->begin << "\t" << begin->end << "\t" << score << endl;
 
  while(begin != end){
    for(list<Edge>::iterator edge = begin->edges.begin(); edge != begin->edges.end(); edge++){
      if(edge->to->label == 1 && edge->to->pred == begin){
	score += edge->score;
	cout << edge->to->begin << "\t" << edge->to->end << "\t" << score  << endl;
	begin = edge->to;
	break;
      }
    }
  }   
}


void MoveObject::addWeights(){
  for (list<MoveNode>::iterator it = nodes.begin(); it != nodes.end(); it++){
    for (list<Edge>::iterator iter =  it->node->edges.begin(); iter != it->node->edges.end(); iter++){
      iter->score += it->weight;
    }
  }
  for (list<MoveEdge>::iterator iter = edges.begin(); iter != edges.end(); iter++){
    iter->edge->score += iter->weight;
  }
}
void MoveObject::undoAddWeights(){
  for (list<MoveNode>::iterator it = nodes.begin(); it != nodes.end(); it++){
    for (list<Edge>::iterator iter =  it->node->edges.begin(); iter != it->node->edges.end(); iter++){
      iter->score -= it->weight;
    }
  }
  for (list<MoveEdge>::iterator iter = edges.begin(); iter != edges.end(); iter++){
    iter->edge->score -= iter->weight;
  }
}


ostream& operator<<(ostream& ostrm, Node *node){

  if(node->n_type >= IR){
    ostrm << node->begin << "\t" << node->end << "\t" << nodeTypeIdentifiers[node->n_type]<< "\t" << node->score;
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


