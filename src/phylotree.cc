/**********************************************************************
 * file:    phylotree.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  representation of the phylogenetic relationship of the species
 *          in comparative gene prediction
 * authors: Stefanie König
 *
 * date    |   author      |  changes
 * --------|---------------|------------------------------------------
 * 09.03.12|Stefanie König | creation of the file
 **********************************************************************/


#include "../include/phylotree.hh"
#include "parser/parser.h"
#include <queue>
#include <cmath>
#include <iostream>
#include <fstream>


void Treenode::printNode() const {
  cout<<"Species: "<<this->species<<" Distance: "<<this->dist_to_parent;
  if (this->parent == NULL){
    cout<<" root ";
  }
  if (this->children.empty()){
    cout<<" leaf "<<endl;
  }
  else {
    cout<<" Children :";
    for(list<Treenode*>::const_iterator node = this->children.begin(); node != this->children.end(); node++){
      cout<<(*node)->species<<endl;
    }
  }
}

double Treenode::calculateAlphaScore(bool label){
  
  double alpha_score = 1.0;
  for(list<Treenode*>::iterator it = this->children.begin(); it != this->children.end(); it++){

    alpha_score*=((P(label, 0, (*it)->dist_to_parent)*(*it)->alpha.at(0)) + (P(label, 1, (*it)->dist_to_parent)*(*it)->alpha.at(1)));
  }
  return alpha_score;

}

double P(bool label1, bool label2, double dist){  //the substitution model: Neyman-two-state model

  const double mu = 2; //TODO: command-line parameter

  if(label1 == label2){
    return (0.5)*(1+exp(-mu*dist));
  }
  else{
    return (0.5)*(1-exp(-mu*dist));
  }
}

void PhyloTree::printTree() const {
  for(list<Treenode*>::const_iterator node = this->treenodes.begin(); node != this->treenodes.end(); node++){
      (*node)->printNode();
  }
}

PhyloTree::PhyloTree(string filename){

  filebuf fb;
  fb.open(filename.c_str(),ios::in);
  if (fb.is_open()){
    istream istrm(&fb);
    Parser parser(&treenodes, istrm);  //define an object of the Parser class
#ifndef DEBUG
    parser.setDebug(false);
#endif
    parser.parse(); //start parsing
    fb.close();
  }
  else
    throw ProjectError("PhyloTree::PhyloTree: Could not open this file!");
}

PhyloTree::~PhyloTree(){
  for(list<Treenode*>::iterator it = treenodes.begin(); it != treenodes.end(); it++){
    delete *it;
  }
}

void PhyloTree::printWithGraphviz(string filename) const {

  //creates inputfile für graphviz
  queue<Treenode*> tree;
  tree.push(this->treenodes.back());
  int i=0;
  int j=-1;
  ofstream file;
  file.open(filename.c_str());
  
  file<<"digraph Tree {\n";
  file<<"rankdir=TB;\n";
  file<<"\tnode[shape=circle];\n";
  file<<"edge [arrowhead=none];\n";

  file<<i<<"[shape=point];\n";

  while(!tree.empty()){
    Treenode *next=tree.front();
    tree.pop();
    j++;
    for(list<Treenode*>::iterator node = next->children.begin(); node != next->children.end(); node++){
      tree.push((*node));
      i++;
      if( !(*node)->children.empty()){
	file<<i<<"[shape=point];\n";
      }
      else{
	file<<i<<"[label="<<(*node)->species<<",shape=ellipse];\n";
      }
      file<<j<<"->"<<i<<"[label="<<(*node)->dist_to_parent<<"];\n";
    }
  }
  file<<"}\n";
  file.close();
}

double PhyloTree::pruningAlgor(const OrthoExon &orthoex,const OrthoGraph &orthograph){

  double tree_score;
  string key = orthoex.getKey(orthograph);
  if(cache::inHash(key)){
    cache::incrementCounter(key);
    tree_score = cache::getScore(key);
  }
  else{
    for(list<Treenode*>::iterator it = this->treenodes.begin(); it != this->treenodes.end(); it++){
      if((*it)->isLeaf()){
	/*
	 * initialization
	*/
	if(key.at(OrthoExon::getVectorPositionSpecies((*it)->species)) == '2'){     // no Exon present -> sum over all possible labels (0,1)
	  cout<<(*it)->species<<"\t"<<(OrthoExon::getVectorPositionSpecies((*it)->species))<<"\t"<<"NULL"<<endl;
	  (*it)->alpha.at(0) = 1; 
	  (*it)->alpha.at(1) = 1;
	}
	else{
	  if(key.at(OrthoExon::getVectorPositionSpecies((*it)->species)) == '1'){   // Exon present and also part of the path
	    cout<<(*it)->species<<"\t"<<(OrthoExon::getVectorPositionSpecies((*it)->species))<<"\t"<<"1"<<endl;
	    (*it)->alpha.at(0) = 0; 
	    (*it)->alpha.at(1) = 1;
	  }
	  if(key.at(OrthoExon::getVectorPositionSpecies((*it)->species)) == '0'){   // Exon present, but not part of the path
	    cout<<(*it)->species<<"\t"<<(OrthoExon::getVectorPositionSpecies((*it)->species))<<"\t"<<"0"<<endl;
	    (*it)->alpha.at(0) = 1; 
	    (*it)->alpha.at(1) = 0;
	  }
	}
      }
      /*
	computation of the alpha values for the interior nodes
       */
      else{
	(*it)->alpha.at(0) = (*it)->calculateAlphaScore(0);
	(*it)->alpha.at(1) = (*it)->calculateAlphaScore(1);
      }
    }
    /*
      computation of the overall tree score
     */
    tree_score  = ( 0.5*(this->treenodes.back()->alpha.at(0) + this->treenodes.back()->alpha.at(1)) ); //equilibrium frequencies: (0.5, 0.5)
    cache::addToHash(key, tree_score);

    // #ifdef DEBUG
  cout<<"#####################################################################\n";
  cout<<"# tableau Prunning Alogrithm\n";
  cout<<"#####################################################################\n";
  cout<<"node\t\t"<<"label 0\t\t"<<"label 1\n";
  for(list<Treenode*>::iterator it = this->treenodes.begin(); it != this->treenodes.end(); it++){
     cout<<(*it)->species<<"\t\t"<<(*it)->alpha.at(0)<<"\t\t"<<(*it)->alpha.at(1)<<"\n";
    }
    cout<<"#####################################################################\n";
    //#endif
  }
  return tree_score;
}
