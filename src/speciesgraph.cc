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
     * extend size of NeutralLine in case that an additional Exons has a smaller start position or a 
     * greater end position than any sampled exon
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
    vector<Node*> intergenic;   
    neutralLines.push_back(intergenic);
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

    vector<Node*> T_plus1;                   // intron in phase 1 on forward strand with preceeding 'T' :  ..T|GT... 
    neutralLines.push_back(T_plus1);
    vector<Node*> TA_plus2;                  // intron in phase 2 on forward strand with preceeding 'TA' : .TA|GT...
    neutralLines.push_back(TA_plus2);
    vector<Node*> TG_plus2;                  // intron in phase 2 on forward strand with preceeding 'TG' :  .TG|GT...
    neutralLines.push_back(TG_plus2);
    vector<Node*> T_minus1;                  // intron in phase 1 on reverse strand with preceeding 'T' :  ..T|CT...
    neutralLines.push_back(T_minus1);
    vector<Node*> C_minus1;                  // intron in phase 1 on reverse strand with preceeding 'C' :  ..C|CT...
    neutralLines.push_back(C_minus1);
    vector<Node*> YY_minus0;                 // intron in phase 0 on reverse strand with preceeding 'TC', 'CT' or 'TT' :  .TC|CT..., .TT|CT... or .CT|CT 
    neutralLines.push_back(YY_minus0);
    

	
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
	
    // add all sampled states to graph (states which are sampled in namgene)
#ifdef DEBUG
    cout << "adding sampled states and additional exon candidates" << endl;
#endif
	
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
		    if(!mergedStopcodon(e1,e2)){
			addIntron(e1, e2, it->next);
		    }
		}
	    }
	}
	if(it->next == NULL) // end of transcript
	    gene_end = true;
    }
	
    //add additional exoncandidates
    if(!additionalExons.empty()){
	for(list<ExonCandidate*>::iterator it = additionalExons.begin(); it!=additionalExons.end(); it++){
	    addExon(*it, neutralLines);
	}
    }

#ifdef DEBUG
    cout << "---------------------------------------------------------" << endl;
    cout << "sampled\t\t" << count_sampled << endl;
    cout << "additional\t" << count_additional << endl;
    cout << "overlap\t\t" << count_overlap << endl << endl;
#endif

    //create neutral lines by connecting neutral nodes in the vector
    for(int j=0; j<neutralLines.size(); j++){
	createNeutralLine(neutralLines.at(j));
    }
	
    //add node weight to the edge weight of all outgoing edges
    addWeightToEdge();
	
    //find topological order of nodes in graph and set pointers topSort_next and topSort_pred
    topSort();

    //relax all nodes in topological order and label all nodes with 1 if on max weight path
    relax();

#ifdef DEBUG
    printGraph(speciesname + ".dot");
#endif
    
}

Node* SpeciesGraph::addExon(Status *exon, vector< vector<Node*> > &neutralLines){

    if(!alreadyProcessed(exon)){
#ifdef DEBUG
	count_sampled++;
#endif
	//cout << "sampled_exon\t\t"<< exon->begin << "\t\t" << exon->end << "\t\t" << (string)stateTypeIdentifiers[((State*)exon->item)->type] << "\t"<< endl;
	Node *ex = new Node(exon->begin, exon->end, setScore(exon), exon->item, sampled);
	printSampledExon(ex);
	nodelist.push_back(ex);
	addToHash(ex);
	addNeutralNodes(ex, neutralLines);    
	return ex;
    }
    return getNode(exon);
}

void SpeciesGraph::addExon(ExonCandidate *exon, vector< vector<Node*> > &neutralLines){
#ifdef DEBUG
    count_additional++;
#endif
    if(!alreadyProcessed(exon)){
	//cout << "unsampled_exon\t\t"<< exon->begin << "\t\t" << exon->end << "\t\t" <<(string)stateTypeIdentifiers[exon->getStateType()] << endl;
	Node *ex = new Node(exon->begin, exon->end, ec_score, exon, unsampled_exon);
	nodelist.push_back(ex);
	addToHash(ex);
	addNeutralNodes(ex, neutralLines);
    }
    else{
#ifdef DEBUG
	count_overlap++;
#endif
    }
}

void SpeciesGraph::addNeutralNodes(Node *node,vector< vector<Node*> > &neutralLines){

    list<NodeType> neut_types = fromNeutralLine(node);
    for(list<NodeType>::iterator it = neut_types.begin(); it != neut_types.end(); it++){
	Node *neut_to = new Node(node->begin, node->begin, 0.0, NULL, *it );
	Edge intron(node, false);
	if(!alreadyProcessed(neut_to)){
	    neutralLines.at(*it).at(node->begin-min) = neut_to;
	    neut_to->edges.push_back(intron);
	    nodelist.push_back(neut_to);
	    addToHash(neut_to);
	}
	else{
	    getNode(neut_to)->edges.push_back(intron);
	    delete neut_to;
	}
    }
  
    Node *neut_from = new Node(node->end, node->end, 0.0, NULL, toNeutralLine(node));
    if(!alreadyProcessed(neut_from)){
	neutralLines.at(neut_from->n_type).at(node->end-min) = neut_from;
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

void SpeciesGraph::addIntron(Node* exon1, Node* exon2, Status *intr){


    if( !edgeExists(exon1,exon2) ){
	//cout << "sampled_intron\t\t"<< intr->begin << "\t\t" << intr->end << "\t\t" << (string)stateTypeIdentifiers[((State*)intr->item)->type] << endl;
	Edge in(exon2, false, setScore(intr), intr->item);
	exon1->edges.push_back(in);
    }
}

void SpeciesGraph::printSampledExon(Node *node){
    streambuf *coutbuf = cout.rdbuf(); //save old buf
    cout.rdbuf(sampled_exons->rdbuf()); //redirect std::cout to species file
    cout << getSeqID() << "\tSAMPLED_ECs\texon\t";
    if(strand == plusstrand){
	cout <<  node->begin + getSeqOffset() + 1 << "\t" << node->end + getSeqOffset() + 1;
    }
    else{
	cout << getSeqLength() - node->end + getSeqOffset() << "\t" << getSeqLength() - node->begin + getSeqOffset();
    }
    cout <<"\t" << node->score << "\t.\t.\tName=" << (string)stateTypeIdentifiers[node->castToStateType()] <<"|"<< node->score << "|";
    if (node->n_type == sampled) {
	cout << ((State*)(node->item))->apostprob << endl;
    }
    cout.rdbuf(coutbuf); //reset to standard output again 
}

list<NodeType> SpeciesGraph::fromNeutralLine(Node *node){

   list<NodeType> fromTypes; 

    StateType type = node->castToStateType();
    int frame = mod3(stateReadingFrames[type]);

    if( isFirstExon(type) ){
	fromTypes.push_back(IR);
    }
    else if(isInternalExon(type) || type == terminal){
	NodeType ntype = (NodeType)(mod3(frame - mod3(node->end - node->begin + 1)) + 1);
	fromTypes.push_back(ntype);
	if(ntype == plus1 && !( strncasecmp(sequence + node->begin, "ag", 2) == 0 || strncasecmp(sequence + node->begin, "aa", 2) == 0 || strncasecmp(sequence + node->begin, "ga", 2) == 0) ){
	    fromTypes.push_back(T_plus1);
	}
	if(ntype == plus2){
	    if(!(strncasecmp(sequence + node->begin, "a", 1) == 0 || strncasecmp(sequence + node->begin, "g", 1) == 0) ){
		fromTypes.push_back(TA_plus2);
	    }
	    if(!(strncasecmp(sequence + node->begin, "a", 1) == 0)){
		fromTypes.push_back(TG_plus2);
	    }
	}
    }
    else if(isRInternalExon(type) || type == rinitial){
	NodeType ntype = (NodeType)(mod3(frame + mod3(node->end - node->begin + 1)) + 4);
	fromTypes.push_back(ntype);
	if(ntype == minus1){
	    if( !( strncasecmp(sequence + node->begin, "ca", 2) == 0 || strncasecmp(sequence + node->begin, "ta", 2) == 0 ) ){
		fromTypes.push_back(T_minus1);
	    }
	    if(!( strncasecmp(sequence + node->begin, "ta", 2) == 0) ){
		fromTypes.push_back(C_minus1);
	    }
	}
	if(ntype == minus0 && !( strncasecmp(sequence + node->begin, "a", 1) == 0) ){
	    fromTypes.push_back(YY_minus0);
	}
    }
    else
	throw ProjectError("in SpeciesGraph::fromNeutralLine(): node " + getKey(node)); 
    return fromTypes;
}

NodeType SpeciesGraph::toNeutralLine(Node *node){

    StateType type = node->castToStateType();
    int frame = mod3(stateReadingFrames[type]);

    if( isLastExon(type)){
	return IR;
    }
    else if( isInternalExon(type) || isInitialExon(type) ){
	if(frame == 1 && strncasecmp(sequence + node->end, "t", 1) == 0 ){
	    return T_plus1;
	}
	else if(frame == 2 && strncasecmp(sequence + node->end - 1 , "ta", 2) == 0){
	    return TA_plus2;
	}
	else if(frame == 2 && strncasecmp(sequence + node->end - 1, "tg", 2) == 0){
	    return TG_plus2;
	}
	else{
	    return (NodeType)(frame + 1);
	}
    }
    else if( isRInternalExon(type) || isRTerminalExon(type) ){
	if(frame == 1 && strncasecmp(sequence + node->end, "t", 1) == 0 ){
	    return T_minus1;
	}
	else if(frame == 1 && strncasecmp(sequence + node->end, "c", 1) == 0){
	    return C_minus1;
	}
	else if(frame == 0 && ( strncasecmp(sequence + node->end - 1, "tt", 2) == 0 || strncasecmp(sequence + node->end - 1, "tc", 2) == 0 || strncasecmp(sequence + node->end - 1, "ct", 2) == 0 ) ){
	    return YY_minus0;
	}
	else{
	    return (NodeType)(frame + 4);
	}
    }
    else
	throw ProjectError("in SpeciesGraph::toNeutralLine(): node " +  getKey(node)); 
    return NOT_KNOWN;
}

void SpeciesGraph::printGraph(string filename, Node *begin, Node *end, bool only_sampled){

    //creates inputfile for graphviz
    ofstream file;
    file.open((filename).c_str());
  
    file<<"digraph MEAgraph {\n";
    file<<"rankdir=LR;\n";
 
    file<<"\tnode[shape=box];\n";
    Node *pos = begin;
    Node *current_path = begin;
    while(pos != end){

	if(only_sampled && (pos->n_type == unsampled_exon && pos->label == 0)){
	}
	else{
	    for(list<Edge>::iterator it=pos->edges.begin(); it!=pos->edges.end(); it++){

		if(only_sampled && (it->to->n_type == unsampled_exon && it->to->label == 0)){
		}
		else{
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
			name1 += "ex" + itoa(type1) +  "_" + itoa(pos->begin+1) + "_"+ itoa(pos->end+1);
		    else
			name1 += nodeTypeIdentifiers[pos->n_type] + "_" + itoa(pos->begin+1) + "_" + itoa(pos->end+1);

		    if(it->to == head){
			name2 += "head";
		    }
		    else if (it->to == tail){
			name2 += "tail";
		    }
		    else if(it->to->item != NULL)
			name2 += "ex" + itoa(type2)  + "_" + itoa(it->to->begin+1) + "_" + itoa(it->to->end+1);
		    else
			name2 += nodeTypeIdentifiers[it->to->n_type] + "_" + itoa(it->to->begin+1) + "_" + itoa(it->to->end+1);

	   
		    file<<name1<<"[";
		    if(pos==head)
			file<<"style=filled,label=head";
		    else if(pos->item == NULL)
			file<<"shape=point";
		    else if(pos->n_type == unsampled_exon)
			file<<"style=filled,fillcolor=lightgrey,";
		    file<<"];\n";

		    file<<name2<<"[";
		    if(it->to==tail)
			file<<"style=filled,label=tail";
		    else if(it->to->item == NULL)
			file<<"shape=point";
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

			if(pos == head && it->to->n_type >= IR && it->to->n_type <=YY_minus0)
			    file<<"label=" << nodeTypeIdentifiers[it->to->n_type] << "];\n";
			else if(it->to == tail && pos->n_type >= IR && pos->n_type <= YY_minus0)
			    file<<"label=" << nodeTypeIdentifiers[pos->n_type] << "];\n";
			else if( pos->n_type >= IR && pos->n_type <= YY_minus0  &&  it->to->n_type >= IR && it->to->n_type <= YY_minus0){
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

    for(list<Node*>::iterator node=nodelist.begin(); node!=nodelist.end(); node++){   
	if((*node)->label == 0){
	    dfs(*node);
	}
    }
    // set pointers to the preceding nodes in topSort
    Node *next = head;
    while(next != tail){
	next->topSort_next->topSort_pred = next;
	next = 	next->topSort_next;
    }
}

void SpeciesGraph::dfs(Node *node){

    static Node* next = NULL;

    node->label = 1;
    for(list<Edge>::iterator edge=node->edges.begin(); edge!=node->edges.end(); edge++){
	if(edge->to->label == 0)
	    dfs(edge->to);
    }
    if(node != tail){
	node->topSort_next = next;
    }
    next = node;
}

double SpeciesGraph::relax(Node *begin, Node *end){

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
		// update node distance
		edge->to->score = next->score + edge->score;
		edge->to->pred = next;
	    }
	}
	next = next->topSort_next;
    } 

    //set new path labels
    while(next != begin){
	next->label = 1;
	next = next->pred;
    }
    next->label = 1;
    return end->score;

}

Node* SpeciesGraph::getTopSortPred(Node *node){

 if(node != head){
     do{
	 node = node->topSort_pred;
     }
     while(node->label == 0);
    }
    return node;


} 
Node* SpeciesGraph::getTopSortNext(Node *node){

 if(node != tail){
     do{
	 node = node->topSort_next;
     }
     while(node->label == 0);
 }
 return node;
}

Node* SpeciesGraph::getPredExonOnPath(Node *node, size_t step){

 while(step > 0 && node != head){
     do{
	 node = node->topSort_pred;
     }
     while(node->label == 0 || ( node->n_type >= IR && node->n_type <=YY_minus0 ));
     step--;
    }
    return node;
}

Node* SpeciesGraph::getNextExonOnPath(Node *node, size_t step){

 while(step > 0 && node != tail){
     do{
	 node = node->topSort_next;
     }
     while(node->label == 0 || ( node->n_type >= IR && node->n_type <= YY_minus0 ) );
     step--;
 }
 return node;
}

string SpeciesGraph::getKey(Node *n){

    if(n->item == NULL)
	return (itoa(n->begin) +  ":" + itoa((n->n_type)));
    else
	return (itoa(n->begin) + ":" + itoa(n->end) + ":" + itoa( (int)(n->castToStateType()) )); 
}


double SpeciesGraph::setScore(Status *st){

  double score = AugustusGraph::setScore(st);
  updateMaxWeight(score);
  return score;
}


double SpeciesGraph::localChange(Move *move){

#ifdef DEBUG
    //cout << "\nSpecies: " << speciesname <<"\n\nbefore local change" << endl;
#endif
    double local_score = - getScorePath(move->getHead(), move->getTail());  //old score of local path
    move->addWeights();
    relax(move->getHead(), move->getTail());            // relax edges
    move->undoAddWeights();
#ifdef DEBUG
    //cout << "\nafter local change" << endl;
#endif
    local_score += getScorePath(move->getHead(), move->getTail());  // new score of local path
    return local_score;
}

double SpeciesGraph::getScorePath(Node *begin, Node *end){

    double score = 0;
#ifdef DEBUG
    //cout << begin->begin << "\t" << begin->end << "\t" << score  << endl;
#endif
    while(begin != end){
	Node *next = getTopSortNext(begin);
	for(list<Edge>::iterator edge = begin->edges.begin(); edge != begin->edges.end(); edge++){
	    if(edge->to == next){
		score += edge->score;
#ifdef DEBUG
		//cout << edge->to->begin << "\t" << edge->to->end << "\t" << score  << endl;
#endif
		begin = edge->to;
		break;
	    }
	}   
    }
    return score;
}

void SpeciesGraph::printNode(Node *node){

    if(node->n_type >= IR && node->n_type <= YY_minus0){
	cout << node->begin + getSeqOffset() << "\t" << node->end + getSeqOffset() << "\t" << nodeTypeIdentifiers[node->n_type]<< "\t" << endl;
    }
    else if(node->n_type > YY_minus0){
	if(strand == plusstrand){
	    cout <<  node->begin + getSeqOffset() + 1 << "\t" << node->end + getSeqOffset() + 1;
	}
	else{
	    cout << getSeqLength() - node->end + getSeqOffset() << "\t" << getSeqLength() - node->begin + getSeqOffset();
	}
	cout  << "\t" << stateTypeIdentifiers[node->castToStateType()] << "\t" << endl;
    }
    else if( node->begin == -1)
	cout << "head" << endl;
    else{
	cout << "tail" << endl;
    }
}

void Move::addWeights(){
    for (list<MoveNode>::iterator it = nodes.begin(); it != nodes.end(); it++){
	for (list<Edge>::iterator iter =  it->node->edges.begin(); iter != it->node->edges.end(); iter++){
	    iter->score += it->weight;
	}
    }
    for (list<MoveEdge>::iterator iter = edges.begin(); iter != edges.end(); iter++){
	iter->edge->score += iter->weight;
    }
}

void Move::undoAddWeights(){
    for (list<MoveNode>::iterator it = nodes.begin(); it != nodes.end(); it++){
	for (list<Edge>::iterator iter =  it->node->edges.begin(); iter != it->node->edges.end(); iter++){
	    iter->score -= it->weight;
	}
    }
    for (list<MoveEdge>::iterator iter = edges.begin(); iter != edges.end(); iter++){
	iter->edge->score -= iter->weight;
    }
}

void Move::initLocalHeadandTail(){

    Node *left = NULL;
    Node *right = NULL;

    if(nodes.empty() && edges.empty()){
	throw ProjectError("MoveObject::initLocalHeadandTail(): no nodes and no edges specified in Move");
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
    cout <<  graph->getSpeciesname() << endl;
    if (left->label != right->label)
	throw ProjectError("Move::initLocalHeadandTail(): either all nodes have label 1 or all node have label 0");
    else{
	local_head = graph->getPredExonOnPath(left, step_size);
	local_tail = graph->getNextExonOnPath(right, step_size);
    }
}

