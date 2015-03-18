#include "mea.hh"
#include "meaPath.hh"
#include "evaluation.hh"

#include <iostream>
#include <iomanip>

using namespace std;

/*
 * accuracy criterion that evaluates on gene level (complex but therefore not exact)
 */
void getMEAtranscripts(list<Gene> *MEAtranscripts, Gene **sampledGeneStructures, int n, const char* dna){

  bool utr;
  try {
    utr = Properties::getBoolProperty("UTR");
  } catch (...) {
    utr = false;
  }
  double w_gene;
  double w_exon;
  double w_base;
  double w_utr;
 try {
    w_gene = Properties::getdoubleProperty("/MeaPrediction/weight_gene");
  } catch (...) {
    w_gene = 0.0;
  }
 try {
    w_exon = Properties::getdoubleProperty("/MeaPrediction/weight_exon");
  } catch (...) {
    w_exon = 0.0;
  }
 try {
    w_base = Properties::getdoubleProperty("/MeaPrediction/weight_base");
  } catch (...) {
    w_base = 0.0;
  }
 try {
    w_utr = Properties::getdoubleProperty("/MeaPrediction/weight_utr");
  } catch (...) {
    w_utr = 0.0;
  }

  double maxAcc=-1;//infty
  
  Gene *bestG = NULL;

  for (int j=0; j<n; j++){          
    double acc = 0.0;
    for (int i=0; i<n; i++){ 
       Evaluation eval;
       eval.addToEvaluation(sampledGeneStructures[i], sampledGeneStructures[j], bothstrands);
       acc += w_gene * (eval.geneSens + eval.geneSpec) + w_exon * (eval.exonSens + eval.exonSpec) + w_base * (eval.nukSens + eval.nukSpec);
      if(utr)
      acc += w_utr * (eval.UTRexonSens + eval.UTRexonSpec) + w_base * (eval.nucUSens + eval.nucUSpec);
    }    
    if (acc>maxAcc){
      maxAcc = acc;
      bestG = sampledGeneStructures[j];
    }
  }

  while(bestG){
    MEAtranscripts->push_back(*bestG);
    bestG = (Gene*) bestG->next;
  }
}

/*
 * MEA using graph representation
 */

list<Transcript*> &getMEAtranscripts(list<Transcript*> &alltranscripts, const char* dna){
    list<Transcript*> *meaGenes = new list<Transcript*>;

    if (!alltranscripts.empty()){

	bool utr;
	try {
	    utr = Properties::getBoolProperty("UTR");
	} catch (...) {
	    utr = false;
	}

	list<Status> stlist;
	/*
	 * builds datastructure needed for the graph representation
	 */
	buildStatusList(alltranscripts, utr, stlist);

	//build Graph
	AugustusGraph myGraph(&stlist, dna);
	myGraph.buildGraph();

	//myGraph.printGraph("test_graph.dot");

	//find shortest path
	MEApath path(&myGraph);
	path.findMEApath();
        
	getMeaGenelist(path.getPath(), meaGenes);
    } 
    return *meaGenes;
}

void buildStatusList(list<Transcript*> &alltranscripts, bool utr, list<Status> &stlist){
    list<Status> stateList;

    for (list<Transcript*>::iterator it = alltranscripts.begin();it != alltranscripts.end(); it++){
	addToList ((*it)->exons, CDS, &stateList);
	addToList ((*it)->introns, intron, &stateList);
	Gene *g = dynamic_cast<Gene*> (*it);
	if (utr && g){ // the transcript is coding and we have predicted UTRs
	    addToList(g->utr5exons, utr5, &stateList);
	    addToList(g->utr3exons, utr3, &stateList);
	    addToList(g->utr5introns, utr5Intron, &stateList);
	    addToList(g->utr3introns, utr3Intron, &stateList);
	}
	// orders list after genes and startpositions of states
	stateList.sort(compareStatus);
	list<Status>::iterator st = stateList.begin();
	while (st != stateList.end()){
	    st->next = &(*(++st));
	}
	stateList.back().next = NULL;
   
	stlist.splice(stlist.end(), stateList);
    }
} 

void printStatelist(list<Status> *stateList){

  cout<<"-------------State list START-----------------"<<endl;
  list<Status>::iterator da;
  for(da=stateList->begin();da!=stateList->end();da++){
    if(da->name==CDS) 
      cout<<setw(10)<<"CDS";
    if(da->name==intron)
      cout<<setw(10)<<"intron";
    if(da->name==utr3Intron)
      cout<<setw(10)<<"utr3intron";
    if(da->name==utr5Intron)
      cout<<setw(10)<<"utr5intron";
    if(da->name==utr3)
      cout<<setw(10)<<"utr3";
    if(da->name==utr5)
      cout<<setw(10)<<"utr5";
    cout<<"-"<<stateTypeIdentifiers[((State*)da->item)->type]<<"\t"<<da->begin<<"\t"<<da->end<<"\t"<<da->score;
    
    if(da->next==NULL)
      cout<<"\tgene end";
    cout<<endl;
    }
  cout<<"-----------State list END----------------"<<endl;
}

void addToList(State *state, Statename name, list<Status> *slist){
    while (state){
	if (state->end >= state->begin){ // UTR exon can have length 0 when start codon comes right after splice site
	    Status someState(name, state->begin, state->end, (double)state->apostprob, state);
	    slist->push_back(someState);
	}
	state = state->next;
    }
}

bool compareStatus(Status first, Status second){
  return (first.begin < second.begin || (first.begin == second.begin && first.end < second.end));
}


/*
 * transfer nodelist of the graph representation to gene list for the AUGUSTUS output
 */
void getMeaGenelist(list<Node*> meaPath, list<Transcript*> *meaGenes){
    Transcript *currentGene = new Gene();
    for (list<Node*>::reverse_iterator node = meaPath.rbegin(); node != meaPath.rend(); node++){
	if ((*node)->item != NULL){
	    State *ex = new State(*((State*)(*node)->item));
	    addExonToGene(currentGene, ex);
	    if ((*node)->pred == NULL)
		cerr<<"ERROR in getMeaGenelist(): node in meaPath has no predecessor"<<endl;
	    if ((*node)->pred->item != NULL){
		if ((*node)->pred->end != (*node)->begin-1){
		    addIntronToGene(currentGene, (*node)->pred, *node);
		}
	    } else {
		setGeneProperties(currentGene);
		meaGenes->push_front(currentGene);
		currentGene = new Gene();
	    }
	}
    } 
}

void addExonToGene(Transcript *tx, State *exon){

    exon->next = NULL;

    if(isCodingExon(exon->type) || isNcExon(exon->type)){
	if(tx->exons == NULL)
	    tx->exons = exon;     
	else{
	    exon->next = tx->exons;
	    tx->exons = exon;
	}
    }
    else{
	Gene *gene = dynamic_cast<Gene *> (tx);
	if (gene){
	    if(is5UTRExon(exon->type)){
		if(gene->utr5exons == NULL)
		    gene->utr5exons = exon;
		else{
		    exon->next = gene->utr5exons;
		    gene->utr5exons = exon;
		}
	    }
	    else if(is3UTRExon(exon->type)){
		if(gene->utr3exons == NULL)
		    gene->utr3exons = exon;
		else{
		    exon->next = gene->utr3exons;
		    gene->utr3exons = exon;
		}
	    }
	}
    }
}

void addIntronToGene(Transcript* tx, Node* predExon, Node* succExon){
    Edge* intron = NULL;
    for(list<Edge>::iterator edge = predExon->edges.begin(); edge != predExon->edges.end(); edge++){
	if (edge->to == succExon){
	    intron = &(*edge);     
	    break;
	}
    }
    State* intr;
    if (intron != NULL && intron->item != NULL){
	intr = new State(*((State*)intron->item));
    } else {
	intr = new State(predExon->end+1, succExon->begin-1, getIntronStateType((State*)predExon->item,(State*)succExon->item));
    }
    addIntronToGene(tx, intr);
}

void addIntronToGene(Transcript* tx, State *intr){
    
    intr->next = NULL;
    if(isCodingIntron(intr->type) || isNcIntron(intr->type) || intr->type == intron_type || intr->type == rintron_type){
	if(tx->introns == NULL)
	    tx->introns = intr;
	else{
	    intr->next = tx->introns;
	    tx->introns = intr; 
	}
    }
    else{
	Gene *gene = dynamic_cast<Gene *> (tx);
	if (gene){	
	    if(is5UTRIntron(intr->type)){
		if(gene->utr5introns == NULL)
		    gene->utr5introns = intr;
		else{
		    intr->next = gene->utr5introns;
		    gene->utr5introns = intr;
		}
	    }
	    else if(is3UTRIntron(intr->type)){
		if(gene->utr3introns == NULL)
		    gene->utr3introns = intr;
		else{
		    intr->next = gene->utr3introns;
		    gene->utr3introns = intr;
		}
	    }
	}
    }
}

StateType getIntronStateType(State *exon1, State *exon2){

  if(exon1->type >= utr5single && exon1->type <= utr5term && exon2->type >= utr5single && exon2->type <= utr5term)
    return utr5intron;
  if(exon1->type >= rutr5single && exon1->type <= rutr5term && exon2->type >= rutr5single && exon2->type <= rutr5term)
    return rutr5intron;
  if(exon1->type >= utr3single && exon1->type <= utr3term && exon2->type >= utr3single && exon2->type <= utr3term)
    return utr3intron;
  if(exon1->type >= rutr3single && exon1->type <= rutr3term && exon2->type >= rutr3single && exon2->type <= rutr3term)
    return rutr3intron;
  if((exon1->type >= singleG && exon1->type <= terminal) ||(exon2->type >= singleG && exon2->type <= terminal))
    return intron_type;
  if((exon1->type >= rsingleG && exon1->type <= rterminal2) ||(exon2->type >= rsingleG && exon2->type <= rterminal2))
    return rintron_type;

  return TYPE_UNKNOWN;
}

void setGeneProperties(Transcript *tx){

  bool utr;
  
  try {
    utr = Properties::getBoolProperty("UTR");
  } catch (...) {
      utr = false;
  }

  tx->source = "AUGUSTUS";

  if(tx->exons)
      tx->strand = isOnFStrand(tx->exons->type)? plusstrand : minusstrand;

  int transStart, transEnd, codlength = 0;
  int codStart = 0, codEnd = 0;
  State *currState, *rcurrState;

  if(tx->exons != NULL){
    codStart = tx->exons->begin;
    codEnd = tx->exons->end;
    currState = tx->exons;
    while(currState){
      codlength += currState->length();
      if(currState->begin < codStart)
	codStart = currState->begin;
      if(currState->end > codEnd)
	codEnd = currState->end;
      currState = currState->next;
    }    
  }
  Gene *gene = dynamic_cast<Gene *> (tx);
  if(!gene){
      tx->transstart = codStart;
      tx->transend = codEnd;
      // check if transcript is truncated
      if(tx->exons->type != ncsingle &&  tx->exons->type != rncsingle){ // multi-exon-gene
	  if (tx->exons->type != ncinit && tx->exons->type != rncterm) // left truncated
	      tx->complete = false;
	  State *lastExon = tx->exons;
	  while(lastExon->next)
	      lastExon = lastExon->next;
	  if (lastExon->type != ncterm && lastExon->type != rncinit) // right truncated
	      tx->complete=false;
      }
  }
  else{
      gene->codingstart = codStart;
      gene->codingend = codEnd;
      gene->length = codEnd-codStart+1;
      if(!utr){
	  gene->transstart = codStart;
	  gene->transend = codEnd;
      }
      else{	
	  if(gene->strand == plusstrand){
	      currState = gene->utr5exons;
	      rcurrState = gene->utr3exons;
	      if(currState != NULL)
		  transStart = currState->begin;
	      else{
		  transStart = codStart;
	      }
	      if(rcurrState != NULL)
		  transEnd = rcurrState->end;
	      else{
		  transEnd = codEnd;
	      }
	  }
	  else{      
	      currState = gene->utr3exons;     
	      rcurrState = gene->utr5exons;
	      if(currState != NULL)
		  transStart = currState->begin;
	      else{
		  transStart = codStart;
	      }
	      if(rcurrState != NULL)
		  transEnd = rcurrState->end;
	      else{
		  transEnd = codEnd;
	      }
	  }
	  while(currState){
	      if(currState->begin < transStart)
		  transStart = currState->begin;
	      currState = currState->next;
	  }
	  while(rcurrState){
	      if(rcurrState->end > transEnd)
		  transEnd = rcurrState->end;
	      rcurrState = rcurrState->next;
	  }
	  gene->transstart = transStart;
	  gene->transend = transEnd;
      } 
      if(utr){
	  if(gene->strand == plusstrand){
	      if(gene->utr5exons != NULL && (gene->utr5exons->type == utr5internal || gene->utr5exons->type == utr5term))
		  gene->complete5utr = false;
	      if(gene->utr3exons != NULL){
		  currState = gene->utr3exons;
		  while(currState->next)
		      currState = currState->next;
		  if(currState->type == utr3init || currState->type == utr3internal)
		      gene->complete3utr = false;
	      }
	      else
		  gene->complete3utr = false;
	  }
	  else{
	      if(gene->utr3exons != NULL && (gene->utr3exons->type == rutr3internal || gene->utr3exons->type == rutr3init))
		  gene->complete3utr = false;
	      if(gene->utr5exons != NULL){
		  rcurrState = gene->utr5exons;
		  while(rcurrState->next)
		      rcurrState = rcurrState->next;
		  if(rcurrState->type == rutr5internal || rcurrState->type == rutr5term)
		      gene->complete5utr = false;
	      }
	      else
		  gene->complete5utr = false;
	  }
      }
      if(gene->exons != NULL){
	  
	  // determine coding length of aGene
	  gene->clength = 0;
	  for (State* ee = gene->exons; ee != NULL; ee = ee->next)
	      gene->clength += ee->length();
	  
	  State *lastExon = gene->exons;
	  while(lastExon->next)
	      lastExon = lastExon->next;
	  if(gene->strand == plusstrand)
	      gene->frame = mod3(gene->exons->frame() - gene->exons->length());
	  else{
	      gene->frame = mod3(gene->exons->frame() + gene->exons->length());
	      gene->frame = mod3(gene->frame - gene->clength + 1);
	  }
	  if(gene->exons->truncated == TRUNC_LEFT || !isFirstExon(gene->exons->type) || lastExon->truncated == TRUNC_RIGHT || !isLastExon(lastExon->type)){
	      gene->complete = false;
	  } 
      }  
      
#ifdef DEBUG
      cerr<<"################################################\n";
      cerr<<"# (MEA) gene properties\n";
      cerr<<"################################################\n";
      
      if(gene && gene->complete)
	  cerr<<"gene complete"<<endl;
      cerr<<"coding start : "<<gene->codingstart<<endl;
      cerr<<"coding end   : "<<gene->codingend<<endl;
      cerr<<"coding length: "<<gene->clength<<endl;
      cerr<<"first exon   : "<<gene->exons->begin<<":"<<gene->exons->end<<" Frame "<<gene->exons->frame()<<endl;
      cerr<<"------------------------------------------------"<<endl;
      cerr<<"trans start  : "<<gene->transstart<<endl;
      cerr<<"trans end    : "<<gene->transend<<endl;
      
      cerr<<"################################################\n";
#endif  
  }
}
