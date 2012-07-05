/**********************************************************************
 * file:    speciesgraph.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  builds a directed acyclic graph from a set of sampled genes.
 *          The underlying auxiliary structure of the graph consists of seven
 *          neutral lines each representing a type of non-coding segment.
 *          In comparative gene prediction for each species an object of
 *          this class is created.
 * authors: Stefanie König
 *
 * date    |   author      |  changes
 * --------|---------------|------------------------------------------
 * 18.06.12| Stefanie König| creation of the file
 **********************************************************************/


#include "speciesgraph.hh"

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
	
    /* the seven lines of neutral nodes, each line represents a non-coding type between to exons
     * intergenetic -> intergenetic region
     * plus -> introns on the plus strand in all three phase 0,1,2
     * minus -> introns on the minus strand in all three phase 0,1,2
     */
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

    printGraph(speciesname + ".dot");
	
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

Node* SpeciesGraph::getPredecessor(Node* node){
    if(node->n_type > minus2){
	return existingNodes[itoa(node->begin) + ":" + itoa( (NodeType)(fromNeutralLine(node)))];
    }
    else{
	throw ProjectError("in SpeciesGraph::getPredecessor: node is not an exon");
    }
    return NULL;
}
Node* SpeciesGraph::getSuccessor(Node* node){
    if(node->n_type > minus2){
	return existingNodes[itoa(node->end) + ":" + itoa( (NodeType)(toNeutralLine(node)))];
    }
    else{
	throw ProjectError("in SpeciesGraph::getSuccessor: node is not an exon");
    }
    return NULL;
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

void SpeciesGraph::printGraph(string filename){

    //creates inputfile for graphviz
    ofstream file;
    file.open((filename).c_str());
  
    file<<"digraph MEAgraph {\n";
    file<<"rankdir=LR;\n";
 
    file<<"\tnode[shape=box];\n";
    Node *pos = head;
    Node *current_path = head;
    while(pos != NULL){

	if(!(pos->n_type == unsampled_exon && pos->label == 0)){
   
	    for(list<Edge>::iterator it=pos->edges.begin(); it!=pos->edges.end(); it++){

		if(!(it->to->n_type == unsampled_exon && it->to->label == 0)){

		    string name1 = "";
		    string name2 = "";
		    StateType type1 = pos->castToStateType();
		    StateType type2 = it->to->castToStateType();

		    if(pos == head){
			name1 += "head";
		    }
		    else if (pos == tail){
			name1 += "tail";
		    }
		    else if(pos->item != NULL)
			name1 += (string)stateTypeIdentifiers[type1] +  "_" + itoa(pos->begin+1) + "_"+ itoa(pos->end+1);
		    else
			name1 += nodeTypeIdentifiers[pos->n_type] + "_" + itoa(pos->begin+1) + "_" + itoa(pos->end+1);

		    if(it->to == head){
			name2 += "head";
		    }
		    else if (it->to == tail){
			name2 += "tail";
		    }
		    else if(it->to->item != NULL)
			name2 += (string)stateTypeIdentifiers[type2] + "_" + itoa(it->to->begin+1) + "_" + itoa(it->to->end+1);
		    else
			name2 += nodeTypeIdentifiers[it->to->n_type] + "_" + itoa(it->to->begin+1) + "_" + itoa(it->to->end+1);

	   
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
				if(it->to->label == 1 && pos == current_path){
				    file<<"color=blue,";
				    current_path = it->to;
				}
				else
				    file<<"color=red,";
			    }
			    else{
				if(it->to->label == 1 && pos == current_path){
				    file<<"weight=100,color=blue,";
				    current_path = it->to;
				}
				else
				    file<<"weight=100,color=red,";
			    }
			}
			else if(it->to->label == 1 && pos == current_path){
			    file<<"color=blue,";
			    current_path = it->to;
			}

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
	}
	pos = pos->topSort_next;
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
    if(node != tail){
	node->topSort_next = temp;
    }
    temp = node;
}

void SpeciesGraph::relax(Node *begin, Node *end){

    Node *next = begin;
  
    //initialize
    while(next != end){
	next = next->topSort_next;
	next->score =  - numeric_limits<double>::max();  // reset node distances
	next->label = 0;                                 // reset path labels  
    }
    begin->score = 0;
    begin->label= 0;

    //relax
    next = begin;
    while(next != end){
	for(list<Edge>::iterator edge = next->edges.begin(); edge != next->edges.end(); edge++){
	    if(next->score + edge->score > edge->to->score){
		// update exon distance
		edge->to->score = next->score + edge->score;
		edge->to->pred = next;
	    }
	}
	next = next->topSort_next;
    } 

    //set node labels
    while(next != begin){
	next->label = 1;
	next = next->pred;
    }
    next->label = 1;

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

    /*if(local_head == local_tail){
      goLeftOnPath(step_size);
      goRightOnPath(step_size);
      }*/

    
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

    cout << "\nSpecies: " << speciesname <<"\n\nbefore local change" << endl;
    double local_score = - getScorePath(move->getHead(), move->getTail());  //old score of local path

    move->addWeights();
 
    relax(move->getHead(), move->getTail());            // relax edges

    move->undoAddWeights();

    cout << "\nafter local change" << endl;
    local_score += getScorePath(move->getHead(), move->getTail());  // new score of local path

    printGraph(speciesname + ".opt.dot"); 
    return local_score;
}

double SpeciesGraph::getScorePath(Node *begin, Node *end){

    double score = 0;

    while(begin != end){	
	Node *next = begin->topSort_next;
	while(next->label == 0){
	    next = next->topSort_next;
	}
	for(list<Edge>::iterator edge = begin->edges.begin(); edge != begin->edges.end(); edge++){
	    if(edge->to == next){
		score += edge->score;
		cout << edge->to->begin << "\t" << edge->to->end << "\t" << score  << endl;
		begin = edge->to;
		break;
	    }
	}   
    }
    return score;
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
