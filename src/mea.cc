#include "mea.hh"
#include "meaPath.hh"
#include "evaluation.hh"

#include <iostream>
#include <iomanip>

using namespace std;
void getMEAtranscripts(list<Gene> *MEAtranscripts, Gene **sampledGeneStructures, int n, int strlength){

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
	
      // TODO: parameter training to optimize accuracy criterion

      acc += w_gene * (eval.geneSens + eval.geneSpec) + w_exon * (eval.exonSens + eval.exonSpec) + w_base * (eval.nukSens + eval.nukSpec);
      if(utr)
      acc += w_utr * (eval.UTRexonSens + eval.UTRexonSpec) + w_base * (eval.nucUSens + eval.nucUSpec);
    }    
    if (acc>maxAcc){
      maxAcc = acc;
      bestG = sampledGeneStructures[j];
    }
  }
  //cout<<"maximum Accuracy: "<<maxAcc<<"\t best Genestructure: "<<bestG->id<<endl;

  while(bestG){
    MEAtranscripts->push_back(*bestG);
    bestG = bestG->next;
  }
}

void getMEAtranscripts(list<Gene> *meaGenes, list<Gene> *alltranscripts, int strlength){

  cerr<<"----entering getMEAtranscripts()----"<<endl;

  if(!alltranscripts->empty()){
   
    bool utr;
    try {
      utr = Properties::getBoolProperty("UTR");
    } catch (...) {
      utr = false;
    }

    /*
     * evaluates accuracy for each gene
     */
    /*
    list<AltGene> *genes = groupTranscriptsToGenes(alltranscripts);
    for(list<AltGene>::iterator altg=genes->begin(); altg!=genes->end(); altg++){

      for(list<Gene*>::iterator gPred=altg->transcripts.begin(); gPred!=altg->transcripts.end(); gPred++){
	for(list<Gene*>::iterator gAnno=altg->transcripts.begin(); gAnno!=altg->transcripts.end(); gAnno++){
	  Evaluation eval;
	  addToEvaluation(*gPred, *gAnno, bothstrands);
	  acc += w_gene * (eval.geneSens + eval.geneSpec) + w_exon * (eval.exonSens + eval.exonSpec) + w_base * (eval.nukSens + eval.nukSpec);
      if(utr)
      acc += w_utr * (eval.UTRexonSens + eval.UTRexonSpec) + w_base * (eval.nucUSens + eval.nucUSpec);
	}
	acc *= (*gPred)->apostprob;
      }
    }
    */
   
    /*
     * builds datastructure needed for the graph representation
     */
    list<Status> stateList, stlist;

    cerr<<"generating state list"<<endl;
    for(list<Gene>::iterator it=alltranscripts->begin();it!=alltranscripts->end();it++){
      addToList(it->exons,CDS,&stateList);
      addToList(it->introns,intron,&stateList);
      if(utr){
	addToList(it->utr5exons,utr5,&stateList);
	addToList(it->utr3exons,utr3,&stateList);
	addToList(it->utr5introns,utr5Intron,&stateList);
	addToList(it->utr3introns,utr3Intron,&stateList);
      }
      // orders list after genes and startpositions of states
      stateList.sort(compareStatus);
      list<Status>::iterator st = stateList.begin();
      while(st != stateList.end()){
	st->next = &(*(++st));
      }
      stateList.back().next = NULL;
   
      stlist.splice(stlist.end(),stateList);
   
    }
    // printStatelist(&stlist);

    cerr<<"initializing graph object"<<endl;
    //build Graph
    AugustusGraph myGraph(&stlist, strlength);

    cerr<<"finding shortest path"<<endl;
    //find shortest path
    MEApath path(&myGraph);
    
    getMeaGenelist(path.getPath(), meaGenes);
  }   
} 

void printStatelist(list<Status> *stateList){

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
    cout<<"-"<<((State*)da->item)->type<<"\t"<<da->begin<<"\t"<<da->end<<"\t"<<da->score;
    
    if(da->next==NULL)
      cout<<"\tgene end";
    cout<<endl;
    }
}

void addToList(State *state, Statename name, list<Status> *slist){
 
  while(state && state->end >= state->begin){
    Status someState(name, state->begin, state->end, (double)state->apostprob, state);
    slist->push_back(someState);
    state = state->next;
  }
}

bool compareStatus(Status first, Status second){
  return (first.begin < second.begin || (first.begin == second.begin && first.end < second.end));
}


/*
 * transfer nodelist of the graph representation to gene list for the AUGUSTUS output
 */

void getMeaGenelist(list<Node*> meaPath, list<Gene> *meaGenes){
    
  Gene *currentGene = new Gene();

  for(list<Node*>::reverse_iterator node=meaPath.rbegin(); node!=meaPath.rend(); node++){
    if((*node)->item != NULL){
      State *ex = new State(*((State*)(*node)->item));
      addExonToGene(currentGene, ex);
      if((*node)->pred->item != NULL){
	if((*node)->pred->end != (*node)->begin-1){
	  addIntronToGene(currentGene, (*node)->pred, *node);
	}
      }
      else{
	setGeneProperties(currentGene);
	meaGenes->push_front(*currentGene);
	delete currentGene;
	currentGene = new Gene();
      }
    }
  } 
  delete currentGene;  
}

void addExonToGene(Gene *gene, State *exon){

  exon->next = NULL;

  if(isCodingExon(exon->type)){
    if(gene->exons == NULL)
      gene->exons = exon;     
    else{
      exon->next = gene->exons;
      gene->exons = exon;
    }
  }
  else if(is5UTRExon(exon->type)){
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

void addIntronToGene(Gene* gene, Node* predExon, Node* succExon){
 
  Edge* intron = NULL;
  for(list<Edge>::iterator edge=predExon->edgeoffsets.begin(); edge!=predExon->edgeoffsets.end(); edge++){
    if(edge->to == succExon){
      intron = &(*edge);     
      break;
    }
  }
  State* intr;
  if(intron != NULL && intron->item != NULL){
    intr = new State(*((State*)intron->item));
  }
  else{
    intr = new State(predExon->end+1, succExon->begin-1, getIntronStateType((State*)predExon->item,(State*)succExon->item));
  }
  intr->next = NULL;
  if(isCodingIntron(intr->type) || intr->type == intron_type || intr->type == rintron_type){
    if(gene->introns == NULL)
      gene->introns = intr;
    else{
      intr->next = gene->introns;
      gene->introns = intr;
    }
  }
  else if(is5UTRIntron(intr->type)){
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

void setGeneProperties(Gene *gene){

  bool utr;

  try {
    utr = Properties::getBoolProperty("UTR");
  } catch (...) {
    utr = false;
  }

  gene->source = "AUGUSTUS";

  if(gene->exons->type >= rsingleG && gene->exons->type <= rutr3term)
    gene->strand = minusstrand;

  int transStart, transEnd, codlength = 0;
  int codStart = 0, codEnd = 0;
  State *currState, *rcurrState;

  if(gene->exons != NULL){
    codStart = gene->exons->begin;
    codEnd = gene->exons->end;
    currState = gene->exons;
    while(currState){
      codlength += currState->length();
      if(currState->begin < codStart)
	codStart = currState->begin;
      if(currState->end > codEnd)
	codEnd = currState->end;
      currState = currState->next;
    }    
  }
  if(!utr){
    transStart = codStart;
    transEnd = codEnd;
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
  } 
  gene->transstart = transStart;
  gene->transend = transEnd;
  gene->codingstart = codStart;
  gene->codingend = codEnd;
  gene->length = codEnd-codStart+1;
  gene->clength = codlength;

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
    State *lastExon = gene->exons;
    while(lastExon->next)
      lastExon = lastExon->next;
    int codingFrame;
    if(gene->strand == plusstrand)
      codingFrame = mod3(gene->exons->frame() - (gene->exons->end - gene->exons->begin + 1) % 3);
    else
      codingFrame = mod3(gene->exons->frame() + (gene->exons->end - gene->exons->begin + 1) % 3);    
    if(gene->exons->truncated == TRUNC_LEFT || !isFirstExon(gene->exons->type) || lastExon->truncated == TRUNC_RIGHT || !isLastExon(lastExon->type)){
	gene->complete = false;
	gene->frame = codingFrame;
    } 
  }  
 
#ifdef DEBUG
  cerr<<"################################################\n";
  cerr<<"# (MEA) gene properties\n";
  cerr<<"################################################\n";

  if(gene->complete)
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
