#include "mea.hh"
#include "evaluation.hh"

#include<iostream>

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

  list<Gene> *MEAtranscripts = alltranscripts;

  list<Gene>::iterator it;

  for(it=MEAtranscripts->begin();it!=MEAtranscripts->end();it++){

    cout<<"****************************"<<endl;
    for(State *ex=it->exons; ex!=NULL; ex=ex->next){
      cout<<"EXONS:  number of samples: "<<ex->sampleCount<<"\t apost.prob: "<<ex->apostprob<<"\t interval: "<<ex->begin<<" : "<<ex->end<<endl;
    }
    for(State *in=it->introns; in!=NULL; in=in->next){
    cout<<"INTRONS:  number of samples: "<<in->sampleCount<<"\t apost.prob: "<<in->apostprob<<"\t intervall: "<<in->begin<<" : "<<in->end<<endl;
    }
    cout<<"GENE: "<<it->id<<"\t aposteriori prob: "<<it->apostprob<<"\t Interval: "<<it->transstart<<"-"<<it->codingstart<<"-"<<it->codingend<<"-"<<it->transend;
        
    if(it->complete){
      cout<<"\t coding region complete";
    }
    if(it->viterbi){
      cout<<"\t viterbi gene";
    }
    if(it->strand == plusstrand){
      cout<<"\t+++";
    }
    cout<<endl;
  }
  
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
