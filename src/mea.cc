#include "mea.hh"
#include "evaluation.hh"

#include <iostream>

using namespace std;
list<Gene>* getMEAtranscripts(Gene **sampledGeneStructures, int n, int strlength){

  cout<<"# entering MEA with sampled Gene Structures!---------------"<<endl; 

  double maxAcc=-1;//infty
  
  Gene *bestG = NULL;
 
  for (int j=0; j<n; j++){
   
    double acc = 0.0;

    for (int i=0; i<n; i++){
      Evaluation eval;
      
      eval.addToEvaluation(sampledGeneStructures[i], sampledGeneStructures[j]);
      
      // TODO: parameter for sensity to optimize accuracy criterion

	acc += (eval.exonSens + eval.exonSpec) / 2;
      
    }
    acc /= n;
    // cout<<"Accuracy of "<<sampledGeneStructures[j]->id<<" : "<<acc<<endl;
    if (acc>maxAcc){
      maxAcc = acc;
      bestG = sampledGeneStructures[j];
    }
  }
  cout<<"maximum Accuracy: "<<maxAcc<<"\t best Genestructure: "<<bestG->id<<endl;

  for(Gene *g=bestG; g!=NULL; g=g->next){
    g->printGFF();
    cout<<"---------------------"<<endl;
  }
 
  cout<<"# leaving MEA!-------------"<<endl;
 

 return geneToList(bestG);
}



list<Gene>* getMEAtranscripts(list<Gene> *alltranscripts, int strlength){

  cout<<"!!! entering MEA with alltranscripts"<<endl<<"---------------------------------------"<<endl;

  bool utr;
  try {
    utr = Properties::getBoolProperty("UTR");
  } catch (...) {
    utr = false;
  }

  list<Gene> *MEAtranscripts = alltranscripts;
  
/*
  list<Status> stateList;

  list<Gene>::iterator it;

  int noExons, noIntrons, noUTR;
  noExons = noIntrons = noUTR = 0;

  for(it=alltranscripts->begin();it!=alltranscripts->end();it++){
   
    addToList(it->exons,"exons",noExons,stateList);
    addToList(it->introns,"introns",noIntrons,stateList);
    if(utr){
      addToList(it->utr5exons,"UTR",noUTR,stateList);
      addToList(it->utr3exons,"UTR",noUTR,stateList);
      addToList(it->utr5introns,"UTR",noUTR,stateList);
      addToList(it->utr3introns,"UTR",noUTR,stateList);
    }
  }
  
  list<Status>::iterator da;
  for(da=stateList.begin();da!=stateList.end();da++){
    cout<<da->statename<<"\t"<<da->begin<<"\t"<<da->end<<"\t"<<da->score;
    cout<<endl;
  }
  */
  //build Graph
  //find shortest path


  return MEAtranscripts;
  
} 



Gene* listToGene(list<Gene> *genelist){

  Gene *g = &(*(genelist->begin()));

  for( list<Gene>::iterator it = genelist->begin(); it != genelist->end(); it++){
    g->next = &(*it);
    g = g->next;
  }
  return &(*(genelist->begin()));
}

list<Gene>* geneToList(Gene * genes){

  list<Gene> *genelist = new list<Gene>;
 
  while(genes){
    genelist->push_back(*genes);
    genes = genes->next;
  }
  return genelist;  
}

/*void addToList(State *st, string name, int &noStates, list<Status> &slist){

  for(State *state=st; state!=NULL; state=state->next){

    Status someState(name, state->begin, state->end, (double)state->apostprob, state);
    slist.push_back(someState);
    noStates++;
  }
  }*/
