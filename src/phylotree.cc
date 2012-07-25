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


#include "phylotree.hh"
#include "properties.hh"
#include <queue>
#include <cmath>
#include <iostream>
#include <fstream>

#ifdef CPP0X
#include "parser/parser.h"
#endif



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

double Treenode::calculateAlphaScore(bool label, ExonEvoModel &evo){
  
    double alpha_score = 1.0;
    for(list<Treenode*>::iterator it = this->children.begin(); it != this->children.end(); it++){

	alpha_score*=((evo.P(label, 0, (*it)->dist_to_parent)*(*it)->alpha.at(0)) + (evo.P(label, 1, (*it)->dist_to_parent)*(*it)->alpha.at(1)));
    }
    return alpha_score;

}

void PhyloTree::printTree() const {
    for(list<Treenode*>::const_iterator node = this->treenodes.begin(); node != this->treenodes.end(); node++){
	(*node)->printNode();
    }
}

size_t PhyloTree::getVectorPositionSpecies(string name) {
    for(size_t pos=0; pos < species.size(); pos++){
	if (species.at(pos) == name){
	    return pos;
	}
    }
    return species.size();
}

PhyloTree::PhyloTree(string filename){

#ifdef CPP0X

    filebuf fb;
    fb.open(filename.c_str(),ios::in);
    if (fb.is_open()){
	istream istrm(&fb);
	Parser parser(&treenodes, &species, istrm);  //define an object of the Parser class
#ifndef DEBUG
	parser.setDebug(false);
#endif
	/*
	 * start parsing
	 * @error_message: if 0 parsing was successful, if 1 input syntax is wrong
	 */
	int error_message = parser.parse();
	fb.close();
	if(error_message == 1){
	    throw ProjectError("the parsing of " + filename + " has been unsuccessful. Please check, whether the syntax of your input file is correct" );
	}
    }
    else
	throw ProjectError("PhyloTree::PhyloTree: Could not open this file!");

#else
    throw ProjectError("Comparative gene prediction not possible with this compiled version. Please recompile with flag CPP0X.");
#endif
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

double PhyloTree::pruningAlgor(string labelpattern){
 
    for(list<Treenode*>::iterator it = this->treenodes.begin(); it != this->treenodes.end(); it++){
	if((*it)->isLeaf()){
	    /*
	     * initialization
	     */
	    if(labelpattern.at(getVectorPositionSpecies((*it)->species)) == '2'){     // no exon present -> sum over all possible labels (0,1)
		(*it)->alpha.at(0) = 1; 
		(*it)->alpha.at(1) = 1;
	    }
	    else{
		if(labelpattern.at(getVectorPositionSpecies((*it)->species)) == '1'){   // exon present and also part of the path;
		    (*it)->alpha.at(0) = 0; 
		    (*it)->alpha.at(1) = 1;
		}
		if(labelpattern.at(getVectorPositionSpecies((*it)->species)) == '0'){   // exon present, but not part of the path
		    (*it)->alpha.at(0) = 1; 
		    (*it)->alpha.at(1) = 0;
		}
	    }
	}
	/*
	  computation of the alpha values for the interior nodes
	*/
	else{
	    (*it)->alpha.at(0) = (*it)->calculateAlphaScore(0, this->evo);
	    (*it)->alpha.at(1) = (*it)->calculateAlphaScore(1, this->evo);
	}
    }
    /*
      computation of the overall tree score: log Probability multiplied with some constant factor
    */
    double tree_score  = log( ( this->evo.getEquilibriumFreq(0) * this->treenodes.back()->alpha.at(0) ) + ( this->evo.getEquilibriumFreq(1) * this->treenodes.back()->alpha.at(1)) ) * this->evo.getPhyloFactor();
   

    /*#ifdef DEBUG
    cout<<"#####################################################################\n";
    cout<<"# tableau pruning alogrithm\n";
    cout<<"#####################################################################\n";
    cout<<"node\t\t"<<"label 0\t\t"<<"label 1\n";
    for(list<Treenode*>::iterator it = this->treenodes.begin(); it != this->treenodes.end(); it++){
	cout<<(*it)->species<<"\t\t"<<(*it)->alpha.at(0)<<"\t\t"<<(*it)->alpha.at(1)<<"\n";
    }
    cout<<"#####################################################################\n";
    #endif*/
    return tree_score;
  
}

ExonEvoModel::ExonEvoModel(){
  
    try {
	mu  = Properties::getdoubleProperty("/CompPred/exon_loss");
    } catch (...) {
	mu  = 2.0;
    }
    try {
	lambda = Properties::getdoubleProperty("/CompPred/exon_gain");
    } catch (...) {
	lambda = 2.0;
    }
    if(mu <= 0.0 || lambda <= 0.0){
	throw ProjectError("the rates for exon loss/gain have to be positive");
    }
    try {
	phylo_factor  = Properties::getdoubleProperty("/CompPred/phylo_factor");
    } catch (...) {
	phylo_factor = 1.0;
    }
    if(phylo_factor <= 0.0){
	throw ProjectError("phylogenetic factor needs to be real positive number");
    }
}

double ExonEvoModel::P(bool label1, bool label2, double dist) const {  // the substitution model

    if(label1 == false && label2 == true){ // P(0 -> 1)
	return (lambda / (lambda + mu)) * (1 - exp(-(mu + lambda) * dist));
    }
    else if(label1 == true && label2 == false){ // P(1 -> 0)
	return (mu / (lambda + mu)) * (1 - exp(-(mu + lambda) * dist));
    }
    else if(label1 == false && label2 == false){ // P(0 -> 0)
	return (1 - (lambda / (lambda + mu)) * (1 - exp(-(mu + lambda) * dist)));
    }
    else{ // P(1 -> 1)
	return  (1 - (mu / (lambda + mu)) * (1 - exp(-(mu + lambda) * dist)));
    }
    
}

double ExonEvoModel::getEquilibriumFreq(bool label) const {

    if(label == 1){
	return (lambda / (lambda + mu));
    }
    else{
	return (mu / (lambda + mu));
    }  
}
