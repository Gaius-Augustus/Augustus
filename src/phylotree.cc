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
#include "contTimeMC.hh"
#include "orthoexon.hh"
#include <queue>
#include <cmath>
#include <iostream>
#include <fstream>

#ifdef COMPGENEPRED
#include "parser/parser.h"
#endif

void Treenode::printNode() const {
    cout<<"Species: "<<this->species<<" Distance: "<<this->distance;
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

void Treenode::addChild(Treenode *child){
    child->parent=this;
    this->children.push_back(child);
}

void Treenode::makeRoot(){
    this->parent=NULL;
    this->addDistance(0.0);
}

void Treenode::resizeTable(int size, double val){
    table.clear();
    table.resize(size,val);
}

void PhyloTree::printTree() const {
    for(list<Treenode*>::const_iterator node = this->treenodes.begin(); node != this->treenodes.end(); node++){
	(*node)->printNode();
    }
}

int PhyloTree::findIndex(string name) const {
    for(size_t pos=0; pos < species.size(); pos++){
	if (species.at(pos) == name){
	    return pos;
	}
    }
    return species.size();
}

PhyloTree::PhyloTree(string filename){

#ifdef COMPGENEPRED
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
	 * if error_message=0 parsing was successful, if error_message=1 input syntax is wrong
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
    throw ProjectError("Comparative gene prediction not possible with this compiled version. Please recompile with flag COMPGENEPRED set in common.mk.");
#endif
}

PhyloTree::PhyloTree(const vector<string> &speciesnames, double b){
    this->species=speciesnames;
    Treenode *root = new Treenode;
    for(int i=0;i<speciesnames.size(); i++){
	Treenode *leaf = new Treenode(speciesnames[i],b);
	root->addChild(leaf);
	this->treenodes.push_back(leaf);
    }
    this->treenodes.push_back(root);
}

//copy constructor
PhyloTree::PhyloTree(const PhyloTree &other){
    this->species=other.species;
    
    if(!other.treenodes.empty()){
	Treenode *root=other.treenodes.back();
	if(!root->isRoot()){
	    throw ProjectError("last node in tree ist not the root");
	}
	Treenode *root_copy = new Treenode();
	this->treenodes.push_front(root_copy);
	list<Treenode*> stack;
	list<Treenode*> stack_copy;
	stack.push_front(root);
	stack_copy.push_front(root_copy);
	while(!stack.empty()){
	    Treenode *p = stack.front();
	    stack.pop_front();
	    Treenode *p_copy = stack_copy.front();
	    stack_copy.pop_front();
	    if(!p->isLeaf()){
	      for(list<Treenode*>::reverse_iterator it = p->children.rbegin(); it != p->children.rend(); it++){
		    Treenode *copy = new Treenode((*it)->getSpecies(), (*it)->getDist());
		    p_copy->addChild(copy);
		    this->treenodes.push_front(copy);
		    stack_copy.push_front(copy);
		    stack.push_front(*it);
		}
	    }
	}
    }
}

PhyloTree::~PhyloTree(){
    for(list<Treenode*>::iterator it = treenodes.begin(); it != treenodes.end(); it++){
	delete *it;
    }
    treenodes.clear();
    species.clear();
}

void PhyloTree::printWithGraphviz(string filename) const {

    //creates inputfile für graphviz
    queue<Treenode*> tree;
    if(!treenodes.empty()){
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
		    file<<i<<"[label="<<(*node)->getSpecies()<<",shape=ellipse];\n";
		}
		file<<j<<"->"<<i<<"[label="<<(*node)->getDist()<<"];\n";
	    }
	}
	file<<"}\n";
	file.close();
    }
}
double PhyloTree::pruningAlgor(string labelpattern, Evo *evo, int u){
    vector<int> tuple;
    for (int i=0; i<labelpattern.length() ; i++){
    	tuple.push_back(atoi((labelpattern.substr(i,1)).c_str()));
    }
    return pruningAlgor(tuple, evo, u);
}


double PhyloTree::pruningAlgor(vector<int> &tuple, Evo *evo, int u){

    int states = evo->getNumStates();
    
    for(list<Treenode*>::iterator node = treenodes.begin(); node != treenodes.end(); node++){
	if((*node)->isLeaf()){
	    // initialization
	    int c =tuple.at(findIndex((*node)->getSpecies()));
	    if(c >= states || c < 0){    
		(*node)->resizeTable(states,1);  // in the case of unknown characters, we sum over all possibilities
	    }
	    else{
		(*node)->resizeTable(states,0);
		(*node)->setTable(c,1); 
	    }
	}
	else{
	    //recursion for the interior nodes
	    (*node)->resizeTable(states);	
	    for(int i=0; i<states; i++){
		double score = 1.0;
		for(list<Treenode*>::iterator it = (*node)->children.begin(); it != (*node)->children.end(); it++){
		    double sum=0;
		    gsl_matrix *P = evo->getSubMatrixP(u,(*it)->getDist());
		    for(int j=0; j<states; j++){
			sum+=gsl_matrix_get(P,i,j) * (*it)->getTable(j);
		    }
		    score*=sum; 
		}
		(*node)->setTable(i,score);	
	    }
	}
    }
    //printRecursionTable();

    //in the root, we take the weighted average over all states
    double tree_score=0;
    for(int i=0; i<states; i++){
	tree_score += (evo->getPi(i) * treenodes.back()->getTable(i));
    }
    return log(tree_score);
  
}

void PhyloTree::printRecursionTable() const{

    cout<<"+-------------------------------------------------------------------+\n";
    cout<<"|                       recursion table                             |\n";
    cout<<"+-------------------------------------------------------------------+\n";
    cout<<left;
    cout<<setw(15)<<"node";
    for(int i=0; i<(*treenodes.begin())->table.size(); i++){
	cout<<"state "<<setw(9)<<i;
    }
    cout<<endl;
    for(list<Treenode*>::iterator it = treenodes.begin(); it != treenodes.end(); it++){
	cout<<setw(15)<<(*it)->getSpecies();
	for(int i=0; i<(*it)->table.size(); i++){
	    cout<<setw(15)<<(*it)->getTable(i);
	}
	cout<<endl;
    }
    cout<<"+-------------------------------------------------------------------+\n";

}

void PhyloTree::getBranchLengths(vector<double> &branchset) const {

  for(list<Treenode*>::iterator node = treenodes.begin(); node != treenodes.end(); node++){
	if(!((*node)->isRoot()))
	    branchset.push_back((*node)->getDist());
    }
}

Treenode *PhyloTree::getLeaf(string species) const {
  for(list<Treenode*>::iterator node = treenodes.begin(); node != treenodes.end(); node++){
	if( (*node)->getSpecies() == species)
	    return (*node);
    }
    throw ProjectError("species " + species + " not in tree") ;
}

void PhyloTree::drop(string species){

    Treenode *node=getLeaf(species);
    drop(node);
}

void PhyloTree::drop(Treenode *node, ExonEvo *evo){
    
    if(node->isLeaf()){
	Treenode *p=node->getParent();
	if(p){
	    p->removeChild(node);
	    if(p->children.size()==1){
		collapse(p,evo);
	    }
	}
	treenodes.remove(node);
	delete node;
    }
}

void PhyloTree::collapse(Treenode *node, ExonEvo *evo){

    Treenode *child = node->children.front();

    if(!(node->isRoot())){
	Treenode *parent = node->getParent();
	parent->removeChild(node);
	parent->addChild(child);
	double dist= child->getDist() + node->getDist();
	child->addDistance(dist);
	if(evo)
	    evo->addBranchLength(dist);
    }
    else{
	child->makeRoot();
    }
    treenodes.remove(node);
    delete node;   
}
/*
 * MAP inference given a set of weights for leaf nodes
 * if the flag 'fixLeafLabels' is turned on MAP inference is carried out
 * with leaf labels are fixed to the corresponding labels in the hect
 */
double PhyloTree::weightedMAP(OrthoExon &hect, ExonEvo &evo, bool fixLeafLabels){

    int states = evo.getNumStates();
    double k =evo.getPhyloFactor(); //scaling factor 

 start:
    for(list<Treenode*>::iterator node = treenodes.begin(); node != treenodes.end(); node++){
	if((*node)->isLeaf()){
	    // initialization
	    int pos = findIndex((*node)->getSpecies());
	    if(!hect.orthoex[pos]){ //in the case that the exon does not exist, we remove the node
		Treenode* tmp=*node;
		if(node == treenodes.begin()){
		    drop(tmp, &evo);
		    goto start;
		}
		else{
		    node--;
		    drop(tmp, &evo);
		}
	    }
	    else{
		double weight =hect.weights.at(pos);
		if(fixLeafLabels){
		    (*node)->resizeTable(states,-std::numeric_limits<double>::max());
		    int label = hect.labels[pos];
		    (*node)->setTable(label,weight*label);
		}
		else{
		    (*node)->resizeTable(states,0);
		    (*node)->setTable(1,weight);
		} 
		(*node)->bestAssign.clear();
		(*node)->bestAssign.resize(states);
	    }
	}
	else{
	    //recursion for the interior nodes
	    (*node)->resizeTable(states);
	    (*node)->bestAssign.clear();
	    (*node)->bestAssign.resize(states);
	    for(int i=0; i<states; i++){
		double score = 0.0;
		for(list<Treenode*>iterator it = (*node)->children.begin(); it != (*node)->children.end(); it++){
		    double max = -std::numeric_limits<double>::max();
		    int bestAssign = -1;
		    gsl_matrix *P = evo.getSubMatrixLogP(0,(*it)->getDist());
		    for(int j=0; j<states; j++){
			double branch_score = (k*gsl_matrix_get(P,i,j)) + (*it)->getTable(j);
			if(max < branch_score){
			    max = branch_score;
			    bestAssign=j;
			}
		    }
		    score+=max;
		    (*it)->bestAssign[i]=bestAssign;
		}
		(*node)->setTable(i,score);
	    }
	}
    }
    //printRecursionTable();

    // backtracking
    double max = -std::numeric_limits<double>::max();
    if(!treenodes.empty()){
	list< pair<Treenode*,int> > stack;
	int bestAssign = -1;
	Treenode* root = treenodes.back();	
	for(int i=0; i<states; i++){
	    double root_score = (k*evo.getLogPi(i)) + root->getTable(i);
	    if(max < root_score){
		max = root_score;
		bestAssign=i;
	    }
	}
	stack.push_front(make_pair(root,bestAssign));
	while(!stack.empty()){
	    pair<Treenode *,int> p = stack.front();
	    stack.pop_front();
	    Treenode *node=p.first;
	    bestAssign=p.second;
	    if(node->isLeaf()){
		int pos = findIndex(node->getSpecies());
		if(fixLeafLabels && hect.labels[pos] != bestAssign){
		    throw ProjectError("in PhyloTree:weightedMAP. Leaf Labels changed although the flag fixLeafLabels is turned on");
		}
		hect.labels[pos]=bestAssign;
	    }
	    else{
	      for(list<Treenode*>::iterator it=node->children.begin(); it!=node->children.end(); it++){
		    stack.push_front(make_pair(*it, (*it)->bestAssign[bestAssign]));
		    
		}
	    }
	}
    }
    return max;
}
