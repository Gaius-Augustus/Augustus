/**********************************************************************
 * file:    speciesgraph.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  builds a directed acyclic graph from a set of sampled genes.
 *          The underlying auxiliary structure of the graph consists of seven
 *          neutral lines each representing a type of non-coding segment.
 *          In comparative gene prediction for each species an object of
 *          this class is created.
 * authors: Stefanie Koenig
 *
 * date    |   author       |  changes
 * --------|----------------|------------------------------------------
 * 18.06.12| Stefanie Koenig| creation of the file
 **********************************************************************/

#include "speciesgraph.hh"

double SpeciesGraph::ec_thold = 0.0;
double SpeciesGraph::ic_thold = 0.0;
double SpeciesGraph::maxCostOfExonLoss = 0.0;

void SpeciesGraph::buildGraph(double meanIntrLen){
  
    vector< vector<Node*> > neutralLines; //represents the seven neutral lines
    int seqlen = getSeqLength()+1;

    /*
     * hash that stores all auxiliary nodes 
     * key encodes both position and node type and uses at most 25 of the 32 bits of an int32_t
     */
    unordered_map<int32_t,Node*> auxiliaryNodes; // hash of auxiliary nodes

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

    vector<Node*> ncintr(seqlen);                   // intron between two non-coding exons on forward strand
    neutralLines.push_back(ncintr);
    vector<Node*> rncintr(seqlen);                  // intron between two non-coding exons on reverse strand
    neutralLines.push_back(rncintr);


    	
    head = new Node(-1,-1); // initialize nodelist of the graph with head and tail
    nodelist.push_back(head);
    tail = new Node(seqlen,seqlen);
    nodelist.push_back(tail);
	
    calculateBaseScores();
	
    // add all sampled exons and introns to the graph
#ifdef DEBUG
    cout << "adding sampled states and additional exon candidates" << endl;
#endif

    // add all expicit introns and TLstarts/stops
    Status* pred_state = NULL;
    for(list<Status>::iterator it=statelist->begin(); it!=statelist->end(); it++){
	if( it->isIntron()){ // add an intron
	    if(!mergedStopcodon(pred_state, it->next))
	       addIntron(addRightSS(pred_state, neutralLines, auxiliaryNodes), addLeftSS(it->next,neutralLines,auxiliaryNodes),&(*it));
	}
	else if(isTlstartOrstop(pred_state, &(*it))){ // add translation start or stop
	    addLeftSS(&(*it), neutralLines, auxiliaryNodes);
	}
	pred_state = &(*it);
	if(!it->next)
	    pred_state = NULL;
    }

    // add all exons
    Node *pred = NULL;
    for(list<Status>::iterator it=statelist->begin(); it!=statelist->end(); it++){		
	if(it->isExon()){ 
	    Node *node = addExon(&(*it), neutralLines, auxiliaryNodes);
	    if( !pred && !isGeneStart(node))
		addAuxilaryEdge(head,node);
	    pred = node;
	    if(!it->next){ // end of gene
		if( !isGeneEnd(node) )
		    addAuxilaryEdge(pred,tail);
		pred=NULL;
	    }
	}
    }
#ifdef DEBUG
    int num_sampled = existingNodes.size();
#endif

    //add additional exoncandidates
    if(additionalExons && !additionalExons->empty()){
	for(list<ExonCandidate*>::iterator it = additionalExons->begin(); it!=additionalExons->end();it++){
	    Node *node = addExon(*it, neutralLines, auxiliaryNodes);
	    (*it)->setScore(node->score);
	}
    }

#ifdef DEBUG
    cout << "--------------------------------------------------------" << endl;
    cout << "sampled\t" << num_sampled << endl;
    cout << "additional\t" << existingNodes.size() - num_sampled << endl;
#endif

    //create neutral lines by connecting neutral nodes in the vector
    for(int j=0; j<neutralLines.size(); j++){
	double weight = 0.0;
	if(j>0 && j<13){ // neutralLines representing introns between coding exons
	    weight = log(1.0 - 1.0/meanIntrLen);
	}
	createNeutralLine(neutralLines.at(j),weight,onlyCompleteGenes);
    }

    // add backedges without cycles from nodes with edge to IR to the first IR-node from which you can't reach the node again (adds auxiliary-IR-nodes if exons starts at the same position)
    if (overlapComp)
      addBackEdgesComp();

    //topological sorting of nodelist
    topSort();

    // execute tarjan strongest component algorithm to find cycles in O(N+M)
    if (overlapComp)
      tarjan();

    //relax all nodes in topological order and label all nodes with 1 if on max weight path
    relax();
}

Node* SpeciesGraph::addNode(Status *exon){
    NodeType ntype = utrExon;
    if(exon->name == CDS)
	ntype = sampled;
    double score;
    if(!Constant::logreg){
      score = setScore(exon);
    }else{
      score = ec_thold 
	+ Constant::ex_sc[0] // intercept  
	+ Constant::ex_sc[1] // for not having omega
	+ Constant::ex_sc[2] // for not beeing an OE
	+ Constant::ex_sc[3] * log(exon->getLen())
	+ Constant::ex_sc[4] * exon->getPostProb()
	+ Constant::ex_sc[5] * getAvgBaseProb(exon);
      // if (exon->hasEvidence() && exon->name == CDS)
      //  score += maxCostOfExonLoss;
    }
   
    /*

    if (string("fly") == Properties::getProperty("species")){
	score=ec_thold - 9.9121118 + 7.2057311 * exon->getPostProb() + 2.9993128 * getAvgBaseProb(exon) + 0.3998047 * log(exon->getLen());
	if(exon->hasEvidence() && exon->name == CDS)
	    score+=maxCostOfExonLoss;
    } else if (string("arabidopsis") == Properties::getProperty("species")){
	score = ec_thold 
	    - 3.6803 // intercept
	    - 5.1385 // for not having omega
	    + 0.9453 * log(exon->getLen())
	    + 4.2741 * exon->getPostProb() 
	    + 2.5422 * getAvgBaseProb(exon);
	if (exon->hasEvidence() && exon->name == CDS)
	  score += maxCostOfExonLoss;
    } else if (string("human") == Properties::getProperty("species")){ // need to be more specific if another human alignment is used, length not significant
      score = ec_thold 
	- 4.2444 // intercept
	- 4.8556 // for not having omega
	+ 5.4296 * exon->getPostProb() 
	+ 4.3545 * getAvgBaseProb(exon);
      if (exon->hasEvidence() && exon->name == CDS)
	score += maxCostOfExonLoss;
    } else {
      score = ec_thold
	- 6.2313204 // intercept
	- 3.6918148 // for not having omega
	- 0.3606701 // for not beeing an OE
	+ 0.3235385 * log(exon->getLen())
	+ 5.3554965 * exon->getPostProb()
	+ 4.9943482 * getAvgBaseProb(exon);
      if (exon->hasEvidence() && exon->name == CDS)
	score += maxCostOfExonLoss;
    }

    */


    Node *node = new Node(exon->begin, exon->end, score, exon->item, ntype);
    printSampledGF(exon,score);
    nodelist.push_back(node);
    addToHash(node);
    return node;
}

Node* SpeciesGraph::addNode(ExonCandidate *exon){

  double score;
  if(!Constant::logreg){
    score = ec_score;
  }else{
    score = ec_thold
      + Constant::ex_sc[0] // intercept
      + Constant::ex_sc[1] // for not having omega
      + Constant::ex_sc[2] // for not beeing an OE
      + Constant::ex_sc[3] * log(exon->len())
      + Constant::ex_sc[12]; // for not beeing sampled
  }

  /*
    if (string("fly") == Properties::getProperty("species"))
	score =ec_thold - 9.9121118 + 0.3998047 * log(exon->len());
    else if (string("arabidopsis") == Properties::getProperty("species"))
	score = ec_thold 
	    - 3.6803  // intercept
	    - 5.1385 // for not having omega
	    + 0.9453 * log(exon->len());
    
    else if (string("human") == Properties::getProperty("species"))
      score = ec_thold
	- 4.2444 // intercept
	- 4.8556; // for not having omega 
    
    else
      score = ec_thold
	- 6.2313204 // intercept
        - 3.6918148 // for not having omega
        - 0.3606701 // for not beeing an OE
	+ 0.3235385 * log(exon->len());
  */

    Node *node = new Node(exon->begin, exon->end, score, exon, unsampled_exon);
    nodelist.push_back(node);
    addToHash(node);
    return node;
}

Node* SpeciesGraph::addNode(NodeType type, int pos){
    Node *node = new Node(pos, pos, 0.0, NULL, type);
    nodelist.push_back(node);
    return node;
}

template<class T> Node* SpeciesGraph::addExon(T *exon, vector< vector<Node*> > &neutralLines, unordered_map<int32_t,Node*> &auxiliaryNodes){
    
    if(!alreadyProcessed(exon)){
	Node *node = addNode(exon);
	NodeType pred_type = getPredType(node->castToStateType(), node->begin, node->end);
	int begin = node->begin;
	NodeType succ_type = getSuccType(node->castToStateType());
	int end = node->end + 1;
	/*
	 * connect to the preceeding features
	 */
	if(pred_type > IR){
	    Node *pred = getAuxilaryNode(pred_type,begin, auxiliaryNodes);
	    if(pred)
		addAuxilaryEdge(pred,node);
	    else if(genesWithoutUTRs && (pred_type == TLstart || pred_type == rTLstop))
		addAuxilaryEdge(addAuxNodeToLine(IR,begin,neutralLines),node);
	}
	if(pred_type >= IR && pred_type <= rncintr){ // add auxiliary nodes to the neutral lines
	    list<NodeType> pred_types = getPredTypes(node);
	    for(list<NodeType>::iterator it = pred_types.begin(); it != pred_types.end(); it++){
		Node *pred = addAuxNodeToLine(*it, begin, neutralLines);
		addAuxilaryEdge(pred,node);
	    }
	}

	/*
	 * connect to the succeeding features
	 */
	if(succ_type > IR){
	    Node *succ = getAuxilaryNode(succ_type, end, auxiliaryNodes);
	    if(succ)
		addAuxilaryEdge(node,succ);   
	    else if(genesWithoutUTRs && (succ_type == rTLstart || succ_type == TLstop))
		addAuxilaryEdge(node,addAuxNodeToLine(IR,end,neutralLines));
	}
	if(succ_type >= IR && succ_type <= rncintr){ // add auxiliary nodes to the neutral lines
	    if(succ_type>IR)
		end+=Constant::min_intron_len; // min intron length
	    if(end > getSeqLength())
		end = getSeqLength();
	    list<NodeType> succ_types = getSuccTypes(node);
	    for(list<NodeType>::iterator it = succ_types.begin(); it != succ_types.end(); it++){
		Node *succ = addAuxNodeToLine(*it, end, neutralLines);
		addAuxilaryEdge(node,succ);
	    }
	}
	return node;
    }
    return getNode(exon);
}

void SpeciesGraph::addAuxilaryEdge(Node *pred, Node *succ){
    if(!edgeExists(pred,succ)){
	Edge edge(succ,true,pred->score);
	pred->edges.push_back(edge);
    }
}

Node* SpeciesGraph::addAuxilaryNode(NodeType type, int pos, vector< vector<Node*> >&neutralLines, unordered_map<int32_t,Node*> &auxiliaryNodes){

    int32_t key = (pos << 5) //20 bits
	+ type; // 5 bits
    unordered_map<int32_t,Node*>::iterator it = auxiliaryNodes.find(key);
    if(it == auxiliaryNodes.end()){ // insert new auxiliary node
	Node *node = addNode(type, pos);
	auxiliaryNodes.insert(pair<int32_t,Node*>(key,node));
	if(genesWithoutUTRs){
	    if(type == TLstop || type == rTLstart){
		addAuxilaryEdge(node,addAuxNodeToLine(IR,pos,neutralLines));
	    }
	    if(type == TLstart || type == rTLstop){
		addAuxilaryEdge(addAuxNodeToLine(IR,pos,neutralLines),node);
	    }
	}
	return node;
    }
    else{
	return it->second;
    }
}

Node* SpeciesGraph::getAuxilaryNode(NodeType type, int pos, unordered_map<int32_t,Node*> &auxiliaryNodes) const{

    int32_t key = (pos << 5) //20 bits
	+ type; // 5 bits
    unordered_map<int32_t,Node*>::iterator it = auxiliaryNodes.find(key);
    if(it == auxiliaryNodes.end()){ // insert new auxiliary node
	return NULL;
    }
    else{
	return it->second;
    }
}



Node* SpeciesGraph::addAuxNodeToLine(NodeType type, int pos, vector< vector<Node*> >&neutralLines){

    Node *auxNode = neutralLines.at(type).at(pos);
    if (!auxNode){
	auxNode = new Node(pos, pos, 0.0, NULL, type);
	nodelist.push_back(auxNode);
	neutralLines.at(type).at(pos) = auxNode;
    }
    return auxNode;
}

void SpeciesGraph::addIntron(Node* pred, Node* succ, Status *intr){

    if(!pred)
	pred = head;
    if(!succ)
	succ = tail;
    if( !edgeExists(pred,succ) ){
	//cout << "intron\t\t"<< intr->begin << "\t\t" << intr->end << "\t\t" << (string)stateTypeIdentifiers[((State*)intr->item)->type] << endl;
	double intr_score = 0.0;
	if(intr->name == intron){ // only CDS introns have a posterior probability                    
	  if(!Constant::logreg)
	    intr_score = setScore(intr);
	  else
	    intr_score = ic_thold
	      + Constant::in_sc[0] // intercept
	      + Constant::in_sc[1] * intr->getPostProb()
	      + Constant::in_sc[2] * getAvgBaseProb(intr)
	      + Constant::in_sc[3] * log(intr->getLen());
	  
	  /*

	    if (string("fly") == Properties::getProperty("species")){
		intr_score = ic_thold - 5.64405 + 5.640821 * intr->getPostProb() + 4.740363 * getAvgBaseProb(intr) - 0.155695 * log(intr->getLen());
	    } else if (string("arabidopsis") == Properties::getProperty("species")){
		intr_score = ic_thold 
		    -1.34679
		    + 5.20390 * intr->getPostProb()
		    + 4.08664 * getAvgBaseProb(intr)
		    - 0.70702 * log(intr->getLen());
	    } else if (string("human") == Properties::getProperty("species")){
		intr_score = ic_thold 
		    - 2.03403
		    + 4.89456 * intr->getPostProb()
		    + 3.59555 * getAvgBaseProb(intr)
		    - 0.42217 * log(intr->getLen());
	    }else{
	      intr_score = ic_thold
		- 4.693283
		+ 5.772046 * intr->getPostProb()
		+ 4.170951 * getAvgBaseProb(intr)
		- 0.261357 * log(intr->getLen());
	    }

	  */

	  // if (intr->hasEvidence())
	  // intr_score += maxCostOfExonLoss;
	}

	Edge in(succ, false, intr_score, intr->item);
	pred->edges.push_back(in);
	printSampledGF(intr,intr_score);
    }
}

Node* SpeciesGraph::addLeftSS(Status *exon, vector< vector<Node*> >&neutralLines, unordered_map<int32_t,Node*> &auxiliaryNodes){

    if(!exon)
	return NULL;
    NodeType ntype = getPredType(((State*)exon->item)->type, exon->begin, exon->end);
    int begin = exon->begin;
    return addAuxilaryNode(ntype,begin,neutralLines, auxiliaryNodes);
}

Node* SpeciesGraph::addRightSS(Status *exon, vector< vector<Node*> >&neutralLines, unordered_map<int32_t,Node*> &auxiliaryNodes){

    if(!exon)
	return NULL;
    NodeType ntype = getSuccType(((State*)exon->item)->type);
    int end = exon->end + 1;
    return addAuxilaryNode(ntype, end, neutralLines, auxiliaryNodes);	
}


void SpeciesGraph::printSampledGF(Status *st, double score){
    streambuf *coutbuf = cout.rdbuf(); //save old buf
    cout.rdbuf(sampled_GFs->rdbuf()); //redirect std::cout to species file
    cout << getSeqID() << "\tSAMPLING\t";
    cout << stateNameIdentifiers[st->name] << "\t";
    if(strand == plusstrand){
	cout <<  st->begin + getSeqOffset() + 1 << "\t" << st->end + getSeqOffset() + 1;
    }
    else{
	cout << getSeqLength() - st->end + getSeqOffset() << "\t" << getSeqLength() - st->begin + getSeqOffset();
    }
    cout << "\t" << score;
    StateType type = ((State*)st->item)->type;
    // the gff strand of the exon is the "strand product" of the alignment strand and exon type strand
    // e.g. "-" x "-" = "+"
    cout << "\t" << (((isOnFStrand(type) || type == intron_type) == (strand == plusstrand))? '+' : '-');
    if(!st->isCDS()) // frame
	cout << "\t.";
    else if (isOnFStrand(type))
	cout << "\t" << mod3(3-(st->getFrame() - st->getLen()));
    else
	cout << "\t" << mod3(2-st->getFrame());
    cout << "\tName=" << (string)stateTypeIdentifiers[((State*)st->item)->type] <<";postProb="<< st->getPostProb() << ";avgBaseProb=" <<getAvgBaseProb(st)<< endl;
    cout.rdbuf(coutbuf); //reset to standard output again 
}

NodeType SpeciesGraph::getPredType(StateType type, int begin, int end){

    int frame = mod3(stateReadingFrames[type]);


    if(isNcExon(type)){
	if(type == ncsingle || type == ncinit || type == rncsingle || type == rncterm)
	    return IR;
	else if(type == ncinternal || type == ncterm)
	    return ncintr;
	else if(type == rncinternal || type == rncinit)
	    return rncintr;
    }
    else if( isFirstExon(type) ){
        if(!utr)
            return IR;
        else if(utr && (type == singleG || isInitialExon(type) ) )
            return TLstart;
        else
            return rTLstop;
    }
    else if(isInternalExon(type) || type == terminal){
	return (NodeType)(mod3(frame - mod3(end - begin + 1)) + 1);
    }
    else if(isRInternalExon(type) || type == rinitial){
	return (NodeType)(mod3(frame + mod3(end - begin + 1)) + 4);
    }
    else if(isFirstUTRExon(type)){
	return IR;
    }
    else if(type == utr5internal || type == utr5term){
	return utr5intr;
    }
    else if(type == rutr5internal || type == rutr5init){
	return rutr5intr;
    }
    else if(type == rutr5single || type == rutr5term){
	return rTLstart;
    }
    else if(type == utr3internal || type == utr3term){
	return utr3intr;
    }
    else if(type == rutr3internal || type == rutr3init){
	return rutr3intr;
    }
    else if(type == utr3single || type == utr3init){
	return TLstop;
    }
    else
	throw ProjectError("in SpeciesGraph::getPredTypes()"); 
    return NOT_KNOWN;

}

list<NodeType> SpeciesGraph::getPredTypes(Node *node) {

   list<NodeType> predTypes; 

    StateType type = node->castToStateType();
    NodeType ntype = getPredType(type, node->begin, node->end);

    predTypes.push_back(ntype);
    if(ntype == plus1 && !( strncasecmp(sequence + node->begin, "ag", 2) == 0 || strncasecmp(sequence + node->begin, "aa", 2) == 0 || strncasecmp(sequence + node->begin, "ga", 2) == 0) ){
	predTypes.push_back(T_plus1);
    }
    else if(ntype == plus2){
	if(!(strncasecmp(sequence + node->begin, "a", 1) == 0 || strncasecmp(sequence + node->begin, "g", 1) == 0) ){
	    predTypes.push_back(TA_plus2);
	}
	if(!(strncasecmp(sequence + node->begin, "a", 1) == 0)){
	    predTypes.push_back(TG_plus2);
	}
    }
    else if(ntype == minus1){
	if( !( strncasecmp(sequence + node->begin, "ca", 2) == 0 || strncasecmp(sequence + node->begin, "ta", 2) == 0 ) ){
	    predTypes.push_back(T_minus1);
	}
	if(!( strncasecmp(sequence + node->begin, "ta", 2) == 0) ){
	    predTypes.push_back(C_minus1);
	}
    }
    else if(ntype == minus0 && !( strncasecmp(sequence + node->begin, "a", 1) == 0) ){
	predTypes.push_back(YY_minus0);
    }
    return predTypes;
}

NodeType SpeciesGraph::getSuccType(StateType type){

    int frame = mod3(stateReadingFrames[type]);

    if(isNcExon(type)){
	if(type == ncsingle || type == ncterm || type == rncsingle || type == rncinit)
	    return IR;
	else if(type == ncinternal || type == ncinit)
	    return ncintr;
	else if(type == rncinternal || type == rncterm)
	    return rncintr;
    }
    else if( isLastExon(type)){
	if(!utr)
	    return IR;
	else if(utr && (type == terminal || type == singleG) )
	    return TLstop;
	else
	    return rTLstart;
    }
    else if( isInternalExon(type) || isInitialExon(type) ){
	return ((NodeType)(frame + 1));
    }
    else if( isRInternalExon(type) || isRTerminalExon(type) ){
	return ((NodeType)(frame + 4));
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
	throw ProjectError("in SpeciesGraph::getSuccTypes()"); 
    return NOT_KNOWN;

}

list<NodeType> SpeciesGraph::getSuccTypes(Node *node) {

    list<NodeType> succTypes; 
    StateType type = node->castToStateType();
    NodeType ntype = getSuccType(type);

    if(ntype == plus1 && strncasecmp(sequence + node->end, "t", 1) == 0 ){
	succTypes.push_back(T_plus1);
    }
    else if(ntype == plus2 && strncasecmp(sequence + node->end - 1 , "ta", 2) == 0){
	succTypes.push_back(TA_plus2);
    }
    else if(ntype == plus2 && strncasecmp(sequence + node->end - 1, "tg", 2) == 0){
	succTypes.push_back(TG_plus2);
    }
    if(ntype == minus1 && strncasecmp(sequence + node->end, "t", 1) == 0 ){
	succTypes.push_back(T_minus1);
    }
    else if(ntype == minus1 && strncasecmp(sequence + node->end, "c", 1) == 0){
	succTypes.push_back(C_minus1);
    }
    else if(ntype == minus0 && ( strncasecmp(sequence + node->end - 1, "tt", 2) == 0 || strncasecmp(sequence + node->end - 1, "tc", 2) == 0 || strncasecmp(sequence + node->end - 1, "ct", 2) == 0 ) ){
	succTypes.push_back(YY_minus0);
    }
    else{
	succTypes.push_back(ntype);
    }
    return succTypes;
}

bool SpeciesGraph::isGeneStart(Node *exon){ 

    if(!utr)
	return true;
    StateType type = exon->castToStateType();
    if(isNcExon(type))
	return true;
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
    if(isNcExon(type))
	return true;
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
 
    map<Node *, int> nodeIDs;
    int IDcount = 0;

    for(list<Node*>::reverse_iterator rit=nodelist.rbegin(); rit!=nodelist.rend(); ++rit){
	
	Node *pos = (*rit);
	if(pos->n_type != unsampled_exon){
	    file<<IDcount<<"[" + getDotNodeAttributes(pos) + "];\n";

	    for(list<Edge>::iterator it=pos->edges.begin(); it!=pos->edges.end(); it++){
		if( !(pos == head && it->to == tail) && it->to->n_type != unsampled_exon){
		    std::map<Node*,int>::iterator mit = nodeIDs.find(it->to);
		    /*if(mit == nodeIDs.end()){
			throw ProjectError("Internal error in SpeciesGraph::printGraph: topological ordering of nodes is not correct.");
			}*/
		    if(mit != nodeIDs.end()){
			file<<IDcount<<"->"<<mit->second<<"[" + getDotEdgeAttributes(pos, &(*it)) + "];\n";
		    }
		    else{
			cout << "Warning: " << it->to << "comes after tail in topSort" << endl;
		    }
		}
	    }
	    nodeIDs.insert(std::pair<Node*,int>(pos,IDcount));
	    IDcount++;
	}
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
    if(pred->n_type <= IR && succ->n_type <= IR)
	attr +="style=bold";
    if(pred->n_type > IR && pred->n_type < utrExon &&  succ->n_type > IR && succ->n_type < utrExon ){
	attr += "rank=same";
    }
    return attr;
}

void SpeciesGraph::topSort(){

    list<Node*> unsorted = nodelist;
    nodelist.clear();

    for(list<Node*>::iterator node=unsorted.begin(); node!=unsorted.end(); node++){   
	if((*node)->label == 0){ // not visited, yet
	    dfs(*node);
	}
    }
}

void SpeciesGraph::dfs(Node *node){

    node->label = 1; // mark es visited
    for(list<Edge>::iterator edge=node->edges.begin(); edge!=node->edges.end(); edge++){
	if(edge->to->label == 0) 
	    dfs(edge->to);
    }
    nodelist.push_front(node);
}

double SpeciesGraph::relax(){
  
    //initialize
    for(list<Node*>::iterator it=nodelist.begin(); it!=nodelist.end(); it++){
	(*it)->score =  - numeric_limits<float>::max();   // reset node distances
	(*it)->label = 0;                                 // reset path labels  
	(*it)->pred = NULL;
    }
    head->score = 0;

    //relax
    for(list<Node*>::iterator it=nodelist.begin(); it!=nodelist.end(); it++){
	for(list<Edge>::iterator edge = (*it)->edges.begin(); edge != (*it)->edges.end(); edge++){
	    if((*it)->score + edge->score > edge->to->score){
		// update node distance
		edge->to->score = (*it)->score + edge->score;
		edge->to->pred = (*it);
	    }
	}
    } 

    //set new path labels
    Node *next = tail;
    while(next != head){
	next->label = 1;
	next = next->pred;
    }
    next->label = 1;
    return tail->score;
}

void SpeciesGraph::printCurrentPath(){

    Node *next = tail;
    
    cout << "-------------------------------" << endl;
    while(next != head){
	if(next->isSampled())
	    cout << next << endl;
	else if(next->n_type == IR)
	    cout << "gene end" << endl;
	next = next->pred;
    }
    cout << "-------------------------------" << endl;
}

string SpeciesGraph::getKey(Node *n) {

    if(n->item == NULL)
	throw ProjectError("SpeciesGraph::getKey() called on an auxiliary node.");
    else
	return (itoa(n->begin) + ":" + itoa(n->end) + ":" + itoa( (int)(n->castToStateType()) )); 
}

/*
 * old code: optimize cgp by making small moves
 */ 

/*Node* SpeciesGraph::getTopSortPred(Node *node){

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
*/
