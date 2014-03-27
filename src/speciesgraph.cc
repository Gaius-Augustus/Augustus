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
 * 18.06.12| Stefanie Koenig| creation of the file
 **********************************************************************/

#include "speciesgraph.hh"

void SpeciesGraph::buildGraph(){
    
    vector< vector<Node*> > neutralLines; //represents the seven neutral lines
    int seqlen = getSeqLength();

    /* the seven lines of neutral nodes, each line represents a non-coding type between to exons
     * intergenetic -> intergenetic region
     * plus -> introns on the plus strand in all three phase 0,1,2
     * minus -> introns on the minus strand in all three phase 0,1,2
     */
    vector<Node*> intergenic(seqlen);   
    neutralLines.push_back(intergenic);
    vector<Node*> plus0(seqlen);
    neutralLines.push_back(plus0);
    vector<Node*> plus1(seqlen);
    neutralLines.push_back(plus1);
    vector<Node*> plus2(seqlen);
    neutralLines.push_back(plus2);
    vector<Node*> minus0(seqlen);
    neutralLines.push_back(minus0);
    vector<Node*> minus1(seqlen);
    neutralLines.push_back(minus1);
    vector<Node*> minus2(seqlen);
    neutralLines.push_back(minus2);

    vector<Node*> T_plus1(seqlen);                   // intron in phase 1 on forward strand with preceeding 'T' :  ..T|GT... 
    neutralLines.push_back(T_plus1);
    vector<Node*> TA_plus2(seqlen);                  // intron in phase 2 on forward strand with preceeding 'TA' : .TA|GT...
    neutralLines.push_back(TA_plus2);
    vector<Node*> TG_plus2(seqlen);                  // intron in phase 2 on forward strand with preceeding 'TG' :  .TG|GT...
    neutralLines.push_back(TG_plus2);
    vector<Node*> T_minus1(seqlen);                  // intron in phase 1 on reverse strand with preceeding 'T' :  ..T|CT...
    neutralLines.push_back(T_minus1);
    vector<Node*> C_minus1(seqlen);                  // intron in phase 1 on reverse strand with preceeding 'C' :  ..C|CT...
    neutralLines.push_back(C_minus1);
    vector<Node*> YY_minus0(seqlen);                 // intron in phase 0 on reverse strand with preceeding 'TC', 'CT' or 'TT' :  .TC|CT..., .TT|CT... or .CT|CT 
    neutralLines.push_back(YY_minus0);
    	
    head = new Node(-1,-1); // initialize nodelist of the graph with head and tail
    nodelist.push_back(head);
    tail = new Node(seqlen,seqlen);
    nodelist.push_back(tail);
	
    calculateBaseScores();
	
    // add all sampled exons and introns to the graph
#ifdef DEBUG
    cout << "adding sampled states and additional exon candidates" << endl;
#endif
	
    Node *pred = NULL;

    for(list<Status>::iterator it=statelist->begin(); it!=statelist->end(); it++){
		
	if(it->name == CDS || it->name == utr3 || it->name == utr5){ // add an exon (CDS or UTR)
	    Node *node = addExon(&(*it), neutralLines);
	    if( !pred && !isGeneStart(node))
		addAuxilaryEdge(head,getPredUTRSS(node));
	    pred = node;
	    if(!it->next){ // end of gene
		if( !isGeneEnd(node) )
		    addAuxilaryEdge(getSuccUTRSS(pred),tail);
		pred=NULL;
	    }
	}
	else if(it->name == intron || it->name == utr3Intron || it->name == utr5Intron){ // add an intron and its succeeding exon
	    Node *succ=NULL;
	    if(it->next){
		succ = addExon(it->next,neutralLines);
	    }
	    if(it->name == intron){ // CDS intron
		if( (pred && succ && !mergedStopcodon(pred,succ)) || !(pred && succ) ){
		    addIntron(pred, succ, &(*it));
		}
	    }
	    else{ // UTR intron
		if(pred)
		    pred=getSuccUTRSS(pred);
		if(succ)
		    succ=getPredUTRSS(succ);
		addIntron(pred,succ, &(*it));
	    }
	    pred = succ;
	}
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

    cout << "onlyCompleteGenes="<< onlyCompleteGenes << endl;
    //create neutral lines by connecting neutral nodes in the vector
    for(int j=0; j<neutralLines.size(); j++){
	createNeutralLine(neutralLines.at(j),onlyCompleteGenes);
    }
    	
    //find topological order of nodes in graph and set pointers topSort_next and topSort_pred
    topSort();

    //relax all nodes in topological order and label all nodes with 1 if on max weight path
    relax();
    
    //#ifdef DEBUG
    printGraph(speciesname + ".dot");
    //#endif
    
}

Node* SpeciesGraph::addExon(Status *exon, vector< vector<Node*> > &neutralLines){

    if(!alreadyProcessed(exon)){
	//cout << "sampled\t\t"<< exon->begin << "\t\t" << exon->end << "\t\t" << (string)stateTypeIdentifiers[((State*)exon->item)->type] << "\t"<< endl;
	NodeType ntype = utrExon;
	if(exon->name == CDS){
	    ntype = sampled;
#ifdef DEBUG
	    count_sampled++;
#endif
	}
	Node *ex = new Node(exon->begin, exon->end, setScore(exon), exon->item, ntype);
	printSampledExon(ex);
	nodelist.push_back(ex);
	addToHash(ex);
	connectToPred(ex, neutralLines);
	connectToSucc(ex,neutralLines);
	return ex;
    }
    return getNode(exon);
}

void SpeciesGraph::addExon(ExonCandidate *exon, vector< vector<Node*> > &neutralLines){
#ifdef DEBUG
    count_additional++;
#endif
    if(!alreadyProcessed(exon)){
	//cout << "unsampled\t\t"<< exon->begin << "\t\t" << exon->end << "\t\t" <<(string)stateTypeIdentifiers[exon->getStateType()] << endl;
	Node *ex = new Node(exon->begin, exon->end, ec_score, exon, unsampled_exon);
	nodelist.push_back(ex);
	addToHash(ex);
	connectToPred(ex, neutralLines);
	connectToSucc(ex,neutralLines);
    }
    else{
#ifdef DEBUG
	count_overlap++;
#endif
    }
}

void SpeciesGraph::addAuxilaryEdge(Node *pred, Node *succ){
    if(!edgeExists(pred,succ)){
	Edge edge(succ,true,pred->score);
	pred->edges.push_back(edge);
    }
}

Node* SpeciesGraph::addAuxilaryNode(NodeType type, int pos, vector< vector<Node*> >&neutralLines){
    string key = itoa(pos) +  ":" + itoa(type);
    if(!alreadyProcessed(key)){
	Node *node = new Node(pos, pos, 0.0, NULL, type);
	if(type >=IR && type <= YY_minus0)
	    neutralLines.at(type).at(pos) = node;
	nodelist.push_back(node);
	addToHash(node);
	if(genesWithoutUTRs){
	    if(type == TLstop || type == rTLstart){
		addAuxilaryEdge(node,addAuxilaryNode(IR,pos,neutralLines));
	    }
	    if(type == TLstart || type == rTLstop){
		addAuxilaryEdge(addAuxilaryNode(IR,pos,neutralLines),node);
	    }
	}
	return node;
    }
    else{
	return getNode(key);
    }
}

void SpeciesGraph::connectToPred(Node *node,vector< vector<Node*> > &neutralLines){

    list<NodeType> pred_types = getPredTypes(node);
    for(list<NodeType>::iterator it = pred_types.begin(); it != pred_types.end(); it++){
	int begin = node->begin;
	if( *it == rTLstart || *it == TLstop )
	    begin--;
	Node *pred=addAuxilaryNode(*it,begin,neutralLines);
	addAuxilaryEdge(pred,node);
    }
}

Node* SpeciesGraph::getPredUTRSS(Node *node){

    list<NodeType> pred_types = getPredTypes(node);
    list<NodeType>::iterator it = pred_types.begin();
    int begin = node->begin;
    if( *it == rTLstart || *it == TLstop )
	begin--;
    string key = itoa(begin) +  ":" + itoa(*it);
    if(!alreadyProcessed(key)){
	throw ProjectError("Internal error in SpeciesGraph::getPredUTRSS: preceeding UTR is missing.");
    }
    return getNode(key);
}

Node* SpeciesGraph::getSuccUTRSS(Node *node){
    return node->edges.front().to;
}

void SpeciesGraph::connectToSucc(Node *node,vector< vector<Node*> > &neutralLines){

    NodeType succ_type = getSuccType(node);
    int end = node->end;
    if(  succ_type == TLstart ||  succ_type == rTLstop )
	end++;
    Node *succ = addAuxilaryNode(succ_type, end, neutralLines);
    addAuxilaryEdge(node,succ);
}

void SpeciesGraph::addIntron(Node* pred, Node* succ, Status *intr){

    if(!pred)
	pred = head;
    if(!succ)
	succ = tail;
    if( !edgeExists(pred,succ) ){
	//cout << "sampled_intron\t\t"<< intr->begin << "\t\t" << intr->end << "\t\t" << (string)stateTypeIdentifiers[((State*)intr->item)->type] << endl;
	double intr_score = pred->score;
	if(intr->name == intron) // only CDS introns have a posterior probability
	    intr_score += setScore(intr);
	Edge in(succ, false, intr_score, intr->item);
	pred->edges.push_back(in);
    }
}

void SpeciesGraph::printSampledExon(Node *node){
    streambuf *coutbuf = cout.rdbuf(); //save old buf
    cout.rdbuf(sampled_exons->rdbuf()); //redirect std::cout to species file
    cout << getSeqID() << "\tSAMPLED_ECs\t";
    if(node->n_type == sampled){
	cout << "exon\t";
    }
    else{
	cout << "UTR\t";
    }
    if(strand == plusstrand){
	cout <<  node->begin + getSeqOffset() + 1 << "\t" << node->end + getSeqOffset() + 1;
    }
    else{
	cout << getSeqLength() - node->end + getSeqOffset() << "\t" << getSeqLength() - node->begin + getSeqOffset();
    }
    cout <<"\t" << node->score << "\t.\t.\tName=" << (string)stateTypeIdentifiers[node->castToStateType()] <<"|"<< node->score << "|";
    if (node->n_type == sampled || node->n_type == utrExon) {
	cout << ((State*)(node->item))->apostprob << endl;
    }
    cout.rdbuf(coutbuf); //reset to standard output again 
}

list<NodeType> SpeciesGraph::getPredTypes(Node *node) {

   list<NodeType> predTypes; 

    StateType type = node->castToStateType();
    int frame = mod3(stateReadingFrames[type]);

    if( isFirstExon(type) ){
	if(!utr)
	    predTypes.push_back(IR);
	else if(utr && (type == singleG || isInitialExon(type) ) )
	    predTypes.push_back(TLstart);
	else
	    predTypes.push_back(rTLstop);
    }
    else if(isInternalExon(type) || type == terminal){
	NodeType ntype = (NodeType)(mod3(frame - mod3(node->end - node->begin + 1)) + 1);
	predTypes.push_back(ntype);
	if(ntype == plus1 && !( strncasecmp(sequence + node->begin, "ag", 2) == 0 || strncasecmp(sequence + node->begin, "aa", 2) == 0 || strncasecmp(sequence + node->begin, "ga", 2) == 0) ){
	    predTypes.push_back(T_plus1);
	}
	if(ntype == plus2){
	    if(!(strncasecmp(sequence + node->begin, "a", 1) == 0 || strncasecmp(sequence + node->begin, "g", 1) == 0) ){
		predTypes.push_back(TA_plus2);
	    }
	    if(!(strncasecmp(sequence + node->begin, "a", 1) == 0)){
		predTypes.push_back(TG_plus2);
	    }
	}
    }
    else if(isRInternalExon(type) || type == rinitial){
	NodeType ntype = (NodeType)(mod3(frame + mod3(node->end - node->begin + 1)) + 4);
	predTypes.push_back(ntype);
	if(ntype == minus1){
	    if( !( strncasecmp(sequence + node->begin, "ca", 2) == 0 || strncasecmp(sequence + node->begin, "ta", 2) == 0 ) ){
		predTypes.push_back(T_minus1);
	    }
	    if(!( strncasecmp(sequence + node->begin, "ta", 2) == 0) ){
		predTypes.push_back(C_minus1);
	    }
	}
	if(ntype == minus0 && !( strncasecmp(sequence + node->begin, "a", 1) == 0) ){
	    predTypes.push_back(YY_minus0);
	}
    }
    else if(isFirstUTRExon(type)){
	predTypes.push_back(IR);
    }
    else if(type == utr5internal || type == utr5term){
	predTypes.push_back(utr5intr);
    }
    else if(type == rutr5internal || type == rutr5init){
	predTypes.push_back(rutr5intr);
    }
    else if(type == rutr5single || type == rutr5term){
	predTypes.push_back(rTLstart);
    }
    else if(type == utr3internal || type == utr3term){
	predTypes.push_back(utr3intr);
    }
    else if(type == rutr3internal || type == rutr3init){
	predTypes.push_back(rutr3intr);
    }
    else if(type == utr3single || type == utr3init){
	predTypes.push_back(TLstop);
    }
    else
	throw ProjectError("in SpeciesGraph::getPredTypes(): node " + getKey(node)); 
    return predTypes;
}

NodeType SpeciesGraph::getSuccType(Node *node) {

    StateType type = node->castToStateType();
    int frame = mod3(stateReadingFrames[type]);

    if( isLastExon(type)){
	if(!utr)
	    return IR;
	else if(utr && (type == terminal || type == singleG) )
	    return TLstop;
	else
	    return rTLstart;
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
    else if(isLastUTRExon(type)){
	return IR;
    }
    else if(type == utr5internal || type == utr5init){
	return utr5intr;
    }
    else if(type == rutr5internal || type == rutr5term){
	return rutr5intr;
    }
    else if(type == utr5single || type == utr5term){
	return TLstart;
    }
    else if(type == utr3internal || type == utr3init){
	return utr3intr;
    }
    else if(type == rutr3internal || type == rutr3term){
	return rutr3intr;
    }
    else if(type == rutr3single || type == rutr3init){
	return rTLstop;
    }
    else
	throw ProjectError("in SpeciesGraph::getSuccType(): node " +  getKey(node)); 
    return NOT_KNOWN;
}

bool SpeciesGraph::isGeneStart(Node *exon){ 

    if(!utr)
	return true;
    StateType type = exon->castToStateType();
    if( !isCodingExon(type) && !isFirstUTRExon(type) )
	return false;
    if ( !genesWithoutUTRs && isFirstExon(type) )
	return false;
    return true;
}

bool SpeciesGraph::isGeneEnd(Node *exon){

    if(!utr)
	return true;
    StateType type = exon->castToStateType();
    if( !isCodingExon(type) && !isLastUTRExon(type) )
	return false;
     if ( !genesWithoutUTRs && isLastExon(type) )
	 return false;
     return true;
}

void SpeciesGraph::printGraph(string filename){

    //creates inputfile for graphviz
    ofstream file;
    file.open((filename).c_str());
  
    file<<"digraph MEAgraph {\n";
    file<<"rankdir=LR;\n";
 
    Node *pos = tail;
    map<Node *, int> nodeIDs;
    int IDcount = 0;

    while(pos != NULL){
	
	if(pos->n_type != unsampled_exon){
	    file<<IDcount<<"[" + getDotNodeAttributes(pos) + "];\n";
	    
	    for(list<Edge>::iterator it=pos->edges.begin(); it!=pos->edges.end(); it++){
		if( !(pos == head && it->to == tail) && it->to->n_type != unsampled_exon){
		    std::map<Node*,int>::iterator mit = nodeIDs.find(it->to);
		    if(mit == nodeIDs.end()){
			throw ProjectError("Internal error in SpeciesGraph::printGraph: topological ordering of nodes is not correct.");
		    }
		    file<<IDcount<<"->"<<mit->second<<"[" + getDotEdgeAttributes(pos, &(*it)) + "];\n";
		}
	    }
	    nodeIDs.insert(std::pair<Node*,int>(pos,IDcount));
	    IDcount++;
	}
	pos = pos->topSort_pred;
    }
    file<<"}\n";
    file.close();
}

string SpeciesGraph::getDotNodeAttributes(Node *node){

    string attr = "";

    if(node == head){
	attr += "shape=box,style=filled,peripheries=2,color=red,label=head";
    }
    else if (node == tail){
	attr += "shape=box,style=filled,peripheries=2,color=red,label=tail";
    }
    else if(node->item != NULL){
	if(node->n_type == utrExon)
	    attr += "style=filled,";
	StateType type = node->castToStateType();
	attr = attr + "shape=box,label=" + (string)stateTypeIdentifiers[type] +  "_" + itoa(node->begin+1) + "_" + itoa(node->end+1);
    }
    else{
	if(node->n_type == TLstart || node->n_type == TLstop || node->n_type == rTLstart || node->n_type == rTLstop )
	    attr = attr + "shape=diamond, style=filled, color=yellow, label=" + nodeTypeIdentifiers[node->n_type];
	else
	    attr += "shape=point,";

    }
    return attr;
}

string SpeciesGraph::getDotEdgeAttributes(Node *pred, Edge *edge){

    string attr = "";
    Node *succ=edge->to;

    if(succ->pred == pred && succ->label == 1)
	attr += "color=blue, ";
    else if (pred->item == NULL && succ->item == NULL){
	if(edge->item != NULL)
	    attr += "color=green, ";
	    else
		attr += "color=red, ";
    }
    if( edge->score != 0.0 ){
	attr = attr + "label=" + ftoa(edge->score) + ",";
    }
    else if(pred == head && succ->item == NULL){
	attr = attr + "label=" + nodeTypeIdentifiers[succ->n_type] + ",";
    }
    if(pred->n_type > IR && pred->n_type < utrExon &&  succ->n_type > IR && succ->n_type < utrExon ){
	attr += "rank=same";
    }
    return attr;
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
     while(node->label == 0 || ( node->n_type >= IR && node->n_type < utrExon ));
     step--;
    }
    return node;
}

Node* SpeciesGraph::getNextExonOnPath(Node *node, size_t step){

 while(step > 0 && node != tail){
     do{
	 node = node->topSort_next;
     }
     while(node->label == 0 || ( node->n_type >= IR && node->n_type < utrExon ) );
     step--;
 }
 return node;
}

string SpeciesGraph::getKey(Node *n) {

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

    if(node->n_type >= IR && node->n_type < utrExon){
	cout << node->begin + getSeqOffset() << "\t" << node->end + getSeqOffset() << "\t" << nodeTypeIdentifiers[node->n_type]<< "\t" << endl;
    }
    else if(node->n_type >= utrExon){
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
