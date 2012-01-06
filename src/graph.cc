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

  // get size of neutral line

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
	  addPair(&(*it), it->next);
	else{
	  if(it->next->next != NULL)
	    addPair(&(*it), it->next->next);
	  else{
	    addExon(&(*it));
	    addEdgeToTail(&(*it));
	  }
	}
      }
      else{
	addExon(&(*it));
	addEdgeToTail(&(*it));
      }
    }
  } 
  createNeutralLine();
 
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

  // add node weight to edge weight

  for(list<Node*>::iterator node = nodelist.begin(); node != nodelist.end(); node++){
    for(list<Edge>::iterator edge = (*node)->edgeoffsets.begin(); edge != (*node)->edgeoffsets.end(); edge++){
      edge->score += (*node)->score;
    }
  }
  // printGraph("graph_visualization.dot");

  nodelist.sort(compareNodes);
 
  addCompatibleEdges();

  //sort edgeoffset lists

  for(list<Node*>::iterator node=nodelist.begin(); node!=nodelist.end(); node++)
    (*node)->edgeoffsets.sort(compareEdges);

  //printGraph("graph_visualization_extended.dot");
 
  
}

bool Graph::edgeExists(Node *e1, Node *e2){
  if(e1 == NULL || e2 == NULL)
    return false;
  for(list<Edge>::iterator it=e1->edgeoffsets.begin(); it!=e1->edgeoffsets.end(); it++){
    if(it->to == e2)
      return true;
  }
  return false;
}

Node* Graph::getNode(Node *n){
  return existingNodes[getKey(n)];
}

Node* Graph::getNode(Status *st){
  return existingNodes[getKey(st)];
}

Node* Graph::addExon(Status *exon){

  if(!alreadyProcessed(exon)){
    Node *ex = new Node(exon->begin, exon->end, setScore(exon), exon->item);
    nodelist.push_back(ex);
    addToHash(ex);

    if(exonAtGeneStart(exon)){
      
      //cout<<"exon "<<exon->begin<<":"<<exon->end<<" at gene start"<<endl;

      Node *neut = new Node(ex->begin, ex->begin);
      Edge intron(ex, false);
      if(!alreadyProcessed(neut)){	
	neutralLine[ex->begin-min] = neut;
	neut->edgeoffsets.push_back(intron);
	nodelist.push_back(neut);
	addToHash(neut);
      }
      else{
	getNode(neut)->edgeoffsets.push_back(intron);
	delete neut;
      }
    }
    if(exonAtGeneEnd(exon)){
      
      //cout<<"exon "<<exon->begin<<":"<<exon->end<<" at gene end"<<endl;

      Node *neut = new Node(ex->end, ex->end);
      if(!alreadyProcessed(neut)){
	neutralLine[ex->end-min] = neut;
	Edge intron(neut, false);
	ex->edgeoffsets.push_back(intron);
	nodelist.push_back(neut);
	addToHash(neut);
      }
      else{
	Edge intron(getNode(neut), false);
	ex->edgeoffsets.push_back(intron);      
	delete neut;
      }
    }
    return ex;
  }
  return getNode(exon);
}

void Graph::addPair(Status *exon1, Status *exon2){

  Node *e1 = addExon(exon1);
  Node *e2 = addExon(exon2);
  if(exon1->next == exon2 && !edgeExists(e1,e2)){
    Edge in(e2, false);
    e1->edgeoffsets.push_back(in);  
  }
  else if(exon1->next != exon2 && !edgeExists(e1,e2)){
    Edge in(e2, false, setScore(exon1->next), exon1->next->item);
    e1->edgeoffsets.push_back(in);  
  }
}


void Graph::createNeutralLine(){

  Node *pos = head;
  int n = neutralLine.size();
  for(int i=0; i<n; i++){
    if(neutralLine[i] != NULL){
      Edge neut(neutralLine[i]);
      pos->edgeoffsets.push_back(neut);
      pos = neutralLine[i];
    }
  }
  Edge lastEdge(tail);
  pos->edgeoffsets.push_back(lastEdge);
}

void Graph::insertIntron(Node *exon1, Node *exon2){
  if(exon1->item==NULL && exon2->item==NULL){
    Edge intr(exon2);
    exon1->edgeoffsets.push_back(intr);
  }
  else{
    Edge intr(exon2,false);
    exon1->edgeoffsets.push_back(intr);
  }
}

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


void Graph::addBackEdges(){

  map<string,Node*> inQueue;
  list<Node*> neutralNodes;
  Node *pos = head;
  while(pos != tail){
	neutralNodes.push_back(pos);
    for(list<Edge>::iterator edge=pos->edgeoffsets.begin(); edge!=pos->edgeoffsets.end(); edge++){
      if(edge->neutral){
	pos = edge->to;
	break;
      }
      else if(&(*edge) == &pos->edgeoffsets.back())
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
	for(list<Edge>::iterator edge=(*toNeut)->edgeoffsets.begin(); edge!=(*toNeut)->edgeoffsets.end(); edge++)
	  if(!edge->neutral)
	    nrNonNeutralEdges++;
	
	if(nrNonNeutralEdges > 0 && *fromNeut != *toNeut){
	  q.push(*toNeut);
	  while(!q.empty()){
	    Node *pos = q.front();
	    q.pop();
	    for(list<Edge>::iterator edge=pos->edgeoffsets.begin(); edge!=pos->edgeoffsets.end(); edge++){  
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
    for(list<Edge>::iterator edge=(*ex)->edgeoffsets.begin(); edge!=(*ex)->edgeoffsets.end(); edge++){
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
    for(list<Edge>::iterator edge = (*node)->edgeoffsets.begin(); edge != (*node)->edgeoffsets.end(); edge++){
      cout<<endl<<setw(10)<<"----->"<<setw(10)<<edge->to->begin<<setw(10)<<edge->to->end<<setw(10)<<edge->to->score;
      if(edge->neutral)
	cout<<setw(10)<<"neutral";
      if(edge->to->pred == *node)
	cout<<setw(10)<<"pred";
      cout<<endl<<setw(10)<<edge->score<<endl;
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


string AugustusGraph::getKey(Node *n){

  if(n->item == NULL)
    return (itoa(n->begin) + ":neutral");
  else
    return (itoa(n->begin) + ":" + itoa(n->end) + ":" + itoa( (int)((State*)n->item)->type )); 
}

string AugustusGraph::getKey(Status *st){

  return (itoa(st->begin) + ":" + itoa(st->end) + ":" + itoa( (int)((State*)st->item)->type )); 
}

double AugustusGraph::getIntronScore(Status *predExon, Status *nextExon){

  if(nextExon != NULL){
    for(list<Status>::iterator st=statelist->begin(); st!=statelist->end(); st++){
      if(&(*st) == nextExon && &(*st) != &statelist->front()){
	if((--st)->next != NULL)
	  return st->score;
	else
	  return 0;
      }
    }
  }
  else if(predExon != NULL){
    for(list<Status>::iterator st=statelist->begin(); st!=statelist->end(); st++){
      if(&(*st) == predExon && st->next != NULL){
	st++;
	return st->score;
      }
    }
  }
  return 0;
}

void AugustusGraph::addEdgeFromHead(Status *exon){

  if(exon != NULL){
    if(exon->name >= CDS && exon->name <= utr5){

      if(!edgeExists(head,getNode(exon)) && !exonAtGeneStart(exon)){
     
	Edge toExon(getNode(exon),false,getIntronScore(NULL,exon));
	head->edgeoffsets.push_back(toExon);	
      }
    }
    else if(exon->next != NULL){
      addEdgeFromHead(exon->next);
    }
  }
}


void AugustusGraph::addEdgeToTail(Status *exon){
  
  Node *ex = getNode(exon);
 
  if(!edgeExists(ex,tail)){
    if(!exonAtGeneEnd(exon))
      {
	Edge toTail(tail,false,getIntronScore(exon,NULL));
	ex->edgeoffsets.push_back(toTail);
      }
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

      (((is3UTRExon(type1) && is3UTRExon(type2) && ((type1 != utr3term && type2 != utr3init) || (type1 != rutr3init && type2 != rutr3term))) || (is5UTRExon(type1) && is5UTRExon(type2) && ((type1 != utr5term && type2 != utr5init) || (type1 != rutr5init && type2 != rutr5term)))) && sameStrand(type1,type2) && exon1->end < exon2->begin && sameReadingFrame(exon1,exon2));

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
      if(getBasetype(&(*st), pos)>=0)
	s_be += baseScore[getBasetype(&(*st), pos)*seqlength + pos] - r_be;
    return alpha_se * (st->score - r_se) + alpha_be * s_be;
  }
  else{
    double s_bi = 0;
    for(int pos = st->begin; pos<=st->end; pos++)
      if(getBasetype(&(*st), pos)>=0)
	s_bi += baseScore[getBasetype(&(*st), pos)*seqlength + pos] - r_bi;
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


void AugustusGraph::printGraph(string filename){

  cout<<"print graph to file:"<<endl;
  //creates inputfile for graphviz
  map<string,Node*> inQueue;
  queue<Node*> q;
  head->begin = 0;
  head->end = 0;
  q.push(head);
  ofstream file;
  file.open(("/home/lizzy/augustus_output/" + filename).c_str());
  
  file<<"digraph MEAgraph {\n";
  file<<"rankdir=LR;\n";
 
  file<<"\tnode[shape=box];\n";
  while(!q.empty()){
    Node *pos = q.front();
    q.pop(); 
   
    for(list<Edge>::iterator it=pos->edgeoffsets.begin(); it!=pos->edgeoffsets.end(); it++){
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
}



bool compareNodes(Node *first, Node *second){
  return(first->begin <= second->begin);
}

bool compareEdges(Edge first, Edge second){
  return(first.to->begin <= second.to->begin);
}
