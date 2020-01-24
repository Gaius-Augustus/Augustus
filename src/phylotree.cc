/*
 * phylotree.cc
 *
 * License: Artistic License, see file LICENSE.TXT or 
 *          https://opensource.org/licenses/artistic-license-1.0
 * 
 * Description: representation of the phylogenetic relationship of the species
 *              in comparative gene prediction
 */

//project includes
#include "phylotree.hh"
#include "properties.hh"

#include <queue>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

#ifdef COMPGENEPRED
#include "parser/parser.h"
#endif

double PhyloTree::phylo_factor=1;

void Treenode::printNode() const {
    cout<<"Species: "<<this->species<<" (idx="<<this->idx <<")"<<" Distance: "<<this->distance;
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

PhyloTree::PhyloTree(string filename){

#ifdef COMPGENEPRED
    filebuf fb;
    fb.open(filename.c_str(),ios::in);
    if (fb.is_open()){
	istream istrm(&fb);
	vector<string> species;
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
	/*
	 * check if all branch lengths >0
	 */
	vector<double> bs;
	getBranchLengths(bs);
	for(int i=0; i < bs.size(); i++){
	    if(bs[i] <= 0)
		throw ProjectError("Branch lengths must be > 0 in " + filename + "."); 
	}
	numSp=species.size();
	/*
	 * if only a subset of the species is sought, drop all leaf nodes
	 * which are not in the given subset
	 */
	string only_species = Constant::speciesfilenames; // by default use only the species from speciesfilenames
	/* old code: read in species list from separate file
	 * only required with mysql dbaccess
	 * if(only_species.empty())
	 *    Properties::assignProperty("/CompPred/only_species", only_species);
	 */
	if(!only_species.empty()){

	    ifstream ifstrm(only_species.c_str());
	    if (ifstrm.is_open()){
		vector<string> keep; // the subset of species to be kept
		char buf[256];
		while(ifstrm.getline(buf,255)){
		    stringstream stm(buf);
		    string s;
		    if(stm >> s)
			keep.push_back(s);
		}
		ifstrm.close();
		for(int i=0; i<species.size(); i++){
		    bool found=false;
		    for(int j=0; j<keep.size(); j++){
			if(species[i] == keep[j]){
			    found=true;
			    break;
			}
		    }
		    if(!found){ //if species name ist not in list, remove leaf
			drop(species[i]);
		    }
		}
		if(species.size() < 2 || species.size() < keep.size())
		    throw ProjectError(only_species + " has the wrong format. correct format:\n\n" + 
				       "hg19\n" + "mm9\n" + "galGal3\n" + "...\n" + 
				       "The species identifiers must be the same as in the phylogenetic tree. At least two species identifiers must be given.");
	    }
	    else
		throw ProjectError("Could not open input file " + only_species);
	    //printNewick("subtree.nwk");
	}
    }
    else
	throw ProjectError("PhyloTree::PhyloTree: Could not open file " + filename);

#else
    throw ProjectError("Comparative gene prediction not possible with this compiled version. Please recompile with flag COMPGENEPRED = true set in common.mk.");
#endif
}

PhyloTree::PhyloTree(const vector<string> &speciesnames, double b){
    Treenode *root = new Treenode;
    for(int i=0;i<speciesnames.size(); i++){
	Treenode *leaf = new Treenode(speciesnames[i],b);
	root->addChild(leaf);
	this->treenodes.push_back(leaf);
    }
    this->treenodes.push_back(root);
    numSp=speciesnames.size();
}

//copy constructor
PhyloTree::PhyloTree(const PhyloTree &other){
    this->numSp=other.numSp;
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
		  Treenode *copy = new Treenode((*it)->getSpecies(), (*it)->getDist(), (*it)->getIdx());
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
}

void PhyloTree::printNewick(string filename) const {

    if(!treenodes.empty()){
	ofstream file;
	file.open(filename.c_str());
	recursiveNWK(file,treenodes.back());
	file.close();
    }
 
}

void PhyloTree::recursiveNWK(ofstream &file, Treenode *node) const {
    if(node->isLeaf()){
	file<<node->getSpecies()<<":"<<node->getDist();
    }
    else{
	file<<"(";
	for(list<Treenode*>::iterator it=node->children.begin(); it !=node->children.end(); it++){
	    if( it != node->children.begin() )
		file<<",";
	    recursiveNWK(file,*it);
	}
	if(node->isRoot()){
	    file<<")"<<node->getSpecies()<<";";
	}
	else{
	    file<<")"<<node->getSpecies()<<":"<<node->getDist();
	}    
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
	    int c = tuple[(*node)->getIdx()];
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
	    //	    cout<<"transition matrix P for omega "<<u<<endl;	
	    for(int i=0; i<states; i++){
		double score = 1.0;
		for(list<Treenode*>::iterator it = (*node)->children.begin(); it != (*node)->children.end(); it++){
		    double sum=0;
		    gsl_matrix *P = evo->getSubMatrixP(u,(*it)->getDist());
		    //cout<<"---------- Transition Matrix for omega nr "<<u<<" and branch lenght = "<<(*it)->getDist()<<"---------"<<endl;
		    //printCodonMatrix(P);
		    for(int j=0; j<states; j++){
		      sum += gsl_matrix_get(P, i, j) * (*it)->getTable(j);
		      //cout<<gsl_matrix_get(P, i, j)<<"\t";
		    }
		    score *= sum; 
		}
		//cout<<endl;
		(*node)->setTable(i, score);	
	    }
	}
    }
    //printRecursionTable();

    //in the root, we take the weighted average over all states
    double tree_score = 0;
    for(int i=0; i<states; i++){
	double ts = (evo->getPi(i) * treenodes.back()->getTable(i));
	tree_score += ts;
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
    for(list<Treenode*>::const_iterator it = treenodes.begin(); it != treenodes.end(); it++){
	cout<<setw(15)<<(*it)->getSpecies();
	for(int i=0; i<(*it)->table.size(); i++){
	    if( (*it)->getTable(i) == - std::numeric_limits<double>::max()){
		cout<<setw(15)<<"-inf";
	    }
	    else{
		cout<<setw(15)<<(*it)->getTable(i);
	    }
	}
	cout<<endl;
    }
    cout<<"+-------------------------------------------------------------------+\n";

}

void PhyloTree::getBranchLengths(vector<double> &branchset) const {

  for(list<Treenode*>::const_iterator node = treenodes.begin(); node != treenodes.end(); node++){
	if(!((*node)->isRoot()))
	    branchset.push_back((*node)->getDist());
    }
}

void PhyloTree::getSpeciesNames(vector<string> &speciesnames) {
    for(list<Treenode*>::const_iterator node = treenodes.begin(); node != treenodes.end(); node++){
        if((*node)->isLeaf()){
            speciesnames.push_back((*node)->getSpecies());
	    // set species idx
	    int idx = speciesnames.size()-1 ;
	    (*node)->setIdx(idx);
	}
    }
}

Treenode *PhyloTree::getLeaf(string species) const {
  for(list<Treenode*>::const_iterator node = treenodes.begin(); node != treenodes.end(); node++){
	if( (*node)->getSpecies() == species)
	    return (*node);
    }
    throw ProjectError("species " + species + " not in tree") ;
}

void PhyloTree::drop(string species){

    Treenode *node=getLeaf(species);
    drop(node,NULL);
}

void PhyloTree::drop(Treenode *node, Evo *evo){
    
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
	numSp--;
    }
}

void PhyloTree::collapse(Treenode *node, Evo *evo){

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
 * with leaf labels fixed to the labels in the vector
 */

double PhyloTree::MAP(vector<int> &labels, vector<double> &weights, Evo *evo, bool fixLeafLabels){

    int states = evo->getNumStates();
    
    for(list<Treenode*>::iterator node = treenodes.begin(); node != treenodes.end(); node++){
	if((*node)->isLeaf()){
	    // initialization
	    int idx = (*node)->getIdx();
	    int c = labels[idx];
	    if(c >= 4 || c < 0)
		throw ProjectError("PhyloTree::MAP(): species with index "+ itoa(idx) + " does not exist in tree");
	    double weight = weights[idx];
	    (*node)->resizeTable(states,-std::numeric_limits<double>::max());
	    if(fixLeafLabels || c >= 2){ // c=2: EC absent, but aligned is always a fixed state
		(c == 1)? (*node)->setTable(c,weight) : (*node)->setTable(c,0);
	    }
	    else{
		for(int i = 0; i < 2; i++)
		    (*node)->setTable(i,i*weight);
	    } 
	    (*node)->bestAssign.clear();
	    (*node)->bestAssign.resize(states);
	}
	else{
	    //recursion for the interior nodes
	    (*node)->resizeTable(states);
	    (*node)->bestAssign.clear();
	    (*node)->bestAssign.resize(states);
	    for(int i=0; i<states; i++){
		double score = 0.0;
		for(list<Treenode*>::iterator it = (*node)->children.begin(); it != (*node)->children.end(); it++){
		    double max = -std::numeric_limits<double>::max();
		    int bestAssign = -1;
		    gsl_matrix *P = (*it)->getLogP();
		    if(!P){ // upon the first call, set pointer to the corresponding probability matrix
			P = evo->getSubMatrixLogP(0,(*it)->getDist());
			(*it)->setLogP(P);
		    }
		    for(int j=0; j<states; j++){
			double branch_score = (PhyloTree::phylo_factor*(gsl_matrix_get(P,i,j))) + (*it)->getTable(j);
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

    double max = -std::numeric_limits<double>::max();
    if(!treenodes.empty()){
	int bestAssign = -1;
	Treenode* root = treenodes.back();	
	for(int i=0; i<states; i++){
	    double root_score = (PhyloTree::phylo_factor*(evo->getLogPi(i))) + root->getTable(i);
	    if(max < root_score){
		max = root_score;
		bestAssign=i;
	    }
	}
	// backtracking to assign leaf nodes to the MAP labels
	//if(!fixLeafLabels)
	MAPbacktrack(labels, root, bestAssign, fixLeafLabels);
    }
    return max;
}

void PhyloTree::MAPbacktrack(vector<int> &labels, Treenode* root, int bestAssign, bool fixLeafLabels){

    list< pair<Treenode*,int> > stack;
    stack.push_front(make_pair(root,bestAssign));
    while(!stack.empty()){
	pair<Treenode *,int> p = stack.front();
	stack.pop_front();
	Treenode *node=p.first;
	bestAssign=p.second;
	if(node->isLeaf()){
	    int idx = node->getIdx();
	    if( (fixLeafLabels || labels[idx] >= 2 ) && bestAssign != labels[idx])
		throw ProjectError("in MAPbacktrack: fixLeafNodes is one but different node labels");
	    labels[idx]=bestAssign;
	}
	else{
	    for(list<Treenode*>::iterator it=node->children.begin(); it!=node->children.end(); it++){
		stack.push_front(make_pair(*it, (*it)->bestAssign[bestAssign]));
	    }
	}
    }
}

/*
 * calculate diversity (sum of branch lengths)
 */
double PhyloTree::getDiversity(){
    double div=0.0;
    for(list<Treenode*>::const_iterator node = treenodes.begin(); node != treenodes.end(); node++){
	if(!((*node)->isRoot()))
	    div+=(*node)->getDist();
    }
    return div;
}

/*
 * scale tree by multiplying each branch length by a factor
 */
void PhyloTree::scaleTree(double factor){
    for(list<Treenode*>::const_iterator node = treenodes.begin(); node != treenodes.end(); node++){
	if(!((*node)->isRoot())){
	    double dist = (*node)->getDist();
	    dist*=factor;
	    (*node)->addDistance(dist);
	}
    }
}

/*
 * prune all leaf nodes of species that are not present as indicated by a bit vector
 * (i-th bit in the vector is 1 if species i is present and 0 if species i is absent)
 */
void PhyloTree::prune(bit_vector &bv, Evo *evo){
 start:
    for(list<Treenode*>::iterator node = treenodes.begin(); node != treenodes.end(); node++){
	if((*node)->isLeaf()){
	    if(!bv[(*node)->getIdx()]){ 
		Treenode* tmp=*node;
		if(node == treenodes.begin()){
		    drop(tmp,evo);
		    goto start;
		}
		else{
		    node--;
		    drop(tmp,evo);
		}
	    }
	}
    }
}

int PhyloTree::fitch(vector<int> &labels, int states){
    
    for(list<Treenode*>::iterator node = treenodes.begin(); node != treenodes.end(); node++){
	if((*node)->isLeaf()){
	    // initialization
	    int idx = (*node)->getIdx();
	    int c = labels[idx];
	    if(c >= states || c < 0)
		throw ProjectError("PhyloTree::fitch(): index "+ itoa(c) + " out of bounds.");

	    (*node)->resizeTable(states, 100000); // any number > 1
	    (*node)->setTable(c,0);
	    (*node)->bestAssign.clear();
	    (*node)->bestAssign.resize(states);
	}
	else{
	    //recursion for the interior nodes
	    (*node)->resizeTable(states);
	    (*node)->bestAssign.clear();
	    (*node)->bestAssign.resize(states);
	    for(int i=0; i<states; i++){
		double score = 0.0;
		for(list<Treenode*>::iterator it = (*node)->children.begin(); it != (*node)->children.end(); it++){
		    double min = std::numeric_limits<double>::max();
		    int bestAssign = -1;
		    for(int j=0; j<states; j++){
			double branch_score = (*it)->getTable(j);
			if(i != j)
			    branch_score++; // count one substitution
			if(min > branch_score){
			    min = branch_score;
			    bestAssign=j;
			}
		    }
		    score+=min;
		    (*it)->bestAssign[i]=bestAssign;
		}
		(*node)->setTable(i,score);
	    }
	}
    }
    double min = std::numeric_limits<double>::max();
    if(!treenodes.empty()){
	Treenode* root = treenodes.back();	
	for(int i=0; i<states; i++){
	    if(min > root->getTable(i))
		min = root->getTable(i);
	}
    }else{
      min = -1;
    }
    return min;
}
