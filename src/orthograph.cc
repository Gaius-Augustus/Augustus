/**********************************************************************
 * file:    orthograph.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  orthologous graphs for comparative gene prediction
 * authors: Stefanie König
 *
 * date    |   author      |  changes
 * --------|---------------|------------------------------------------
 * 23.03.12|Stefanie König | creation of the file
 **********************************************************************/


#include "orthograph.hh"
#include "geneMSA.hh"
#include "mea.hh"

using namespace std;

PhyloTree *OrthoGraph::tree = NULL;
vector<ofstream*> OrthoGraph::filestreams;
size_t OrthoGraph::numSpecies;


OrthoGraph::OrthoGraph(){
    graphs.resize(numSpecies);
    orthoSeqRanges.resize(numSpecies);
}



OrthoGraph::~OrthoGraph(){
    for(int i = 0; i < numSpecies; i++){
	delete graphs[i];
	delete ptrs_to_alltranscripts[i];
	delete orthoSeqRanges[i];
    }
}

void OrthoGraph::outputGenes(){

    Boolean noInFrameStop;

    try {
	noInFrameStop = Properties::getBoolProperty("noInFrameStop");
    } catch (...) {
	noInFrameStop = false;
    }
 
    for (size_t pos = 0; pos < numSpecies; pos++){

	list<Gene> *genes = new list<Gene>;
	
	Node* current = graphs[pos]->tail;
	Node* head =  graphs[pos]->head;
	Node* predcurrent;
	
	Gene *currentGene = new Gene();
	
	// start backtracking
	while(current != head){
	    
	    while(current->item == NULL){ //skip all neutral nodes
		if (current == head){
		    goto end;
		}
		current = graphs[pos]->getTopSortPred(current);  //exon1
	    }
	    State *ex;
	    if(current->n_type == sampled){
		ex = new State(*((State*)(current->item)));
	    }
	    else{
		ex = new State(current->begin, current->end, current->castToStateType());
	    }
	    addExonToGene(currentGene, ex);
	    predcurrent = graphs[pos]->getTopSortPred(current);
	    if(predcurrent->n_type  == IR){ //end of gene
		setGeneProperties(currentGene);
		genes->push_front(*currentGene);
		delete currentGene;
		currentGene = new Gene();
		current = predcurrent;
	    }
	    else{
		while(predcurrent->item == NULL){
		    if(predcurrent == head){
			setGeneProperties(currentGene);
			genes->push_front(*currentGene);
			delete currentGene;
			goto end;
		    }
		    predcurrent = graphs[pos]->getTopSortPred(predcurrent); //exon2
		}
		addIntronToGene(currentGene, predcurrent, current); //add intron exon2->exon1
		current = predcurrent;
	    }  
	}
    end:

	list<Gene> *filteredTranscripts = new list<Gene>;
	filteredTranscripts = Gene::filterGenePrediction(genes, orthoSeqRanges[pos]->sequence, bothstrands, noInFrameStop);
	list<AltGene> *agl = groupTranscriptsToGenes(filteredTranscripts);

	delete genes;
	/*
	 * possibly more filter steps
	 */
	agl->sort();

	static vector<int> geneid(numSpecies, 1); // made this static so gene numbering goes across sequences and is unique
	int transcriptid;

	// shift gene coordinates, set sequence name, gene and transcript names
	for (list<AltGene>::iterator agit = agl->begin(); agit != agl->end(); ++agit){
	    agit->shiftCoordinates(this->orthoSeqRanges[pos]->offset);
	    agit->seqname =  this->orthoSeqRanges[pos]->seqname;
	    agit->id = "g"+itoa(geneid[pos]);
	    agit->sortTranscripts();

	    transcriptid = 1;
	    for(list<Gene*>::iterator it = agit->transcripts.begin();it != agit->transcripts.end(); ++it ) {
		(*it)->seqname = this->orthoSeqRanges[pos]->seqname;
		(*it)->id = "t" + itoa(transcriptid);
		(*it)->geneid = agit->id;
		transcriptid++;
	    }
	    geneid[pos]++;
	}
	// print the genes
	streambuf *coutbuf = cout.rdbuf(); //save old buf
	cout.rdbuf(filestreams[pos]->rdbuf()); //redirect std::cout to species file!

	//print sequence information

	cout << "#----- prediction on sequence range " << this->orthoSeqRanges[pos]->offset  << "-" << this->orthoSeqRanges[pos]->offset + this->orthoSeqRanges[pos]->length << " (length = "
	     << this->orthoSeqRanges[pos]->length << ", name = "
	     << this->orthoSeqRanges[pos]->seqname << ") -----" << endl << "#" << endl;

	if(!agl->empty()){
	    printGeneList(agl, this->orthoSeqRanges[pos], Constant::codSeqOutput, Constant::proteinOutput, false);
	}
	else{
	    cout << "# (none)" << endl;
	}
	cout.rdbuf(coutbuf); //reset to standard output again   
    }
}

void OrthoGraph::initOutputFiles(){

    string outdir;  //direction for output files

    try {
	outdir = Properties::getProperty("/CompPred/outdir");
    } catch (...) {
	outdir = "";
    }

    outdir = expandHome(outdir); //replace "~" by "$HOME"
    filestreams.resize(numSpecies);
    for(size_t pos = 0; pos < numSpecies; pos++){
	string filename = outdir + tree->species[pos] + ".gff";
	ofstream *out = new ofstream(filename.c_str());
	if(out){
	    filestreams[pos] = out;
	    (*out) << PREAMBLE << endl;
	    (*out) << "#\n#----- prediction for species '" << tree->species[pos] << "' -----" << endl << "#" << endl;
	}
    }
}

void OrthoGraph::closeOutputFiles(){

    for(size_t pos = 0; pos < numSpecies; pos++){
	if(filestreams[pos]){
	    if(filestreams[pos]->is_open()){
		filestreams[pos]->close();
		delete filestreams[pos];
	    }
	}
    }
}


void OrthoGraph::optimize(){

    //loop over OrthoExons
    if(!all_orthoex.empty()){
	for(list<OrthoExon>::iterator orthoex = all_orthoex.begin(); orthoex != all_orthoex.end(); orthoex++){
	    vector<MoveObject*> orthomove = majorityRuleMove(&(*orthoex));
	    if(!orthomove.empty()){
		localMove(orthomove);
		//delete MoveObjects
		for(size_t pos = 0; pos < numSpecies; pos++){
		    delete orthomove[pos];
		} 
	    }   
	}
    }
}

void OrthoGraph::localMove(vector<MoveObject*> &orthomove){


    // number of iterations, the local bounds are expanded and the move is repeated,
    // makes sure, that the move is not kept too locally and get's stuck in local optimum
    size_t maxIterations = 1;  

    // the magnitude of the expansion in each iteration. The number of nodes, the local_heads are shifted left
    // and the local_tails are shifted right on the current path.
   

    for(size_t iter = 0; iter < maxIterations; iter++){

	cout << "iteration: " << iter << endl;

	if(iter > 0){

	    //shift local_heads and local_heads
	    for(size_t pos = 0; pos < numSpecies; pos++){
		if(orthomove[pos]){
		    orthomove[pos]->shiftHead();
		    orthomove[pos]->shiftTail();
		    cout << "local_head: " << orthomove[pos]->getHead() << endl;
		    cout << "local_tail: " << orthomove[pos]->getTail() << endl;
		}
	    }
	}
	
	// do the local changes and determine the  difference  between the new and the old local score
	double score = calculateScoreDiff(orthomove);
	if(score >= 0){
	    cout << "overall score has improved, accept move\n" << endl;
	    break;
	}
	else{
	    cout << "score as not improved, undo changes\n" << endl;
	    // no improvement, undo local changes
	    for(size_t pos = 0; pos < numSpecies; pos++){
		if(orthomove[pos]){
		    graphs[pos]->relax(orthomove[pos]->getHead(), orthomove[pos]->getTail());
		}
	    }
	}
    }
}


double OrthoGraph::calculateScoreDiff(vector<MoveObject*> &orthomove){

    double graph_score = 0; //difference between the new and the old scores of the local paths
    double phylo_score = 0; //diference between the new and the old score of the phylogenetic edges

    // determine all OrthoExons in that range
    // TODO: if there are other OrthoExons within that range, relax them simulateneously
  
    list<OrthoExon> local_orthoexons = orthoExInRange(orthomove);
    //addOrthoIntrons(orthomove,local_orthoexons );
 
    //calculate phylo_score
    phylo_score -= pruningAlgor(local_orthoexons);
    cache::printCache(local_orthoexons);

    // do local changes for each graph  
    for(size_t pos = 0; pos < numSpecies; pos++){
	if(orthomove[pos]){
	    graph_score += graphs[pos]->localChange(orthomove[pos]);
	}
    }
    //calculate new phylo_score
    phylo_score += pruningAlgor(local_orthoexons);
    cache::printCache(local_orthoexons);

    cout << "graph_score: " << graph_score << endl;
    cout << "phylo_score: " << phylo_score << endl;

    return (graph_score + phylo_score);
}


double OrthoGraph::pruningAlgor(list<OrthoExon> &orthoex){

    double tree_score = 0;

    for(list<OrthoExon>::iterator ortho = orthoex.begin(); ortho != orthoex.end(); ortho++){

	string labelpattern = getLabelpattern(*ortho);
	if(cache::inHash(labelpattern)){
	    tree_score += cache::getScore(labelpattern);
	}
	else{
	    double score = tree->pruningAlgor(labelpattern);
	    cache::addToHash(labelpattern, score);
	    tree_score += score;
	}
    }
    return tree_score;
}

string OrthoGraph::getLabelpattern(OrthoExon &ex){

    ex.labelpattern = "";
    for (size_t i = 0; i < ex.orthoex.size(); i++){
	if (ex.orthoex.at(i) == NULL){
	    ex.labelpattern += "2";
	}
	else {
	    map<string, Node*>::iterator it = graphs.at(i)->existingNodes.find(graphs.at(i)->getKey(ex.orthoex.at(i)));
	    if (it != graphs.at(i)->existingNodes.end()){
		bool label = it->second->label;
		if (label == 1){
		    ex.labelpattern += "1";
		}
		else if (label == 0){
		    ex.labelpattern += "0";
		}
	    }
	    else
		throw ProjectError("Error in OrthoExon::getKey: exon " + graphs.at(i)->getKey(ex.orthoex.at(i)) +" not in graph!");
	}
    }
    return ex.labelpattern;
}


list<OrthoExon> OrthoGraph::orthoExInRange(vector<MoveObject*> &orthomove){

    list<OrthoExon> local_orthoexons;

    for(list<OrthoExon>::const_iterator ortho = all_orthoex.begin(); ortho != all_orthoex.end(); ortho++){
	for(size_t pos = 0; pos < numSpecies; pos++){
	    if(orthomove[pos]){
		if(ortho->orthoex[pos]){
		    if(ortho->orthoex[pos]->begin >= orthomove[pos]->getHead()->begin && ortho->orthoex[pos]->end <= orthomove[pos]->getTail()->end){
			local_orthoexons.push_back(*ortho);
			break;
		    }
		}
	    }
	}
    }
    return local_orthoexons;
}

void OrthoGraph::addOrthoIntrons(vector<MoveObject*> &orthomove, list<OrthoExon> &local_orthoexons){

    bool flag = true;
    cout << "schaue rechts" << endl;
    for(size_t pos = 0; pos < numSpecies; pos++){
	if(orthomove[pos]){
	    if(!orthomove[pos]->nodesIsEmpty()){
		Node *node = orthomove[pos]->getNodeBack();
		while(flag){
		    //cout << "entering loop" << endl;
		    //cout << "node " << node << endl; 
		    flag = false;
	  
		    for(list<OrthoExon>::const_iterator ortho = local_orthoexons.begin(); ortho != local_orthoexons.end(); ortho++){
			Node *target = graphs[pos]->getNode(ortho->orthoex[pos]);
			for(list<Edge>::iterator it = node->edges.begin(); it!= node->edges.end(); it++){
			    if(it->to == target){
				cout <<"Kante gefunden: "<< *it << endl;
				cout <<"füge Knoten " << target << " zu dem MoveObject hinzu " << endl;
				orthomove[pos]->addNodeBack(target, graphs[pos]->getMaxWeight());
				orthomove[pos]->addEdgeBack(&(*it), graphs[pos]->getMaxWeight());
				flag = true;
				node = target;
				break;
			    }
			}
			if (flag == true){
			    break;
			}
		    }
		}
	    }
	}
    }
    flag = true;
    cout << "schaue links" << endl;
  
    for(size_t pos = 0; pos < numSpecies; pos++){
	if(orthomove[pos]){
	    if(!orthomove[pos]->nodesIsEmpty()){
		Node *node = orthomove[pos]->getNodeFront();
		while(flag){
		    //cout << "entering loop" << endl;
		    //cout << "node " << node << endl; 
		    flag = false;
		    for(list<OrthoExon>::const_iterator ortho = local_orthoexons.begin(); ortho != local_orthoexons.end(); ortho++){
			Node *target = graphs[pos]->getNode(ortho->orthoex[pos]);
			for(list<Edge>::iterator it = target->edges.begin(); it!= target->edges.end(); it++){
			    if(it->to == node){
				cout <<"Kante gefunden: "<< *it << endl;
				cout <<"füge Knoten " << target << " zu dem MoveObject hinzu " << endl;
				orthomove[pos]->addNodeFront(target, graphs[pos]->getMaxWeight());
				orthomove[pos]->addEdgeFront(&(*it), graphs[pos]->getMaxWeight());
				flag = true;
				node = target;
				break;
			    }
			}
		    }
		}
	    }
	}
    }
}

vector<MoveObject*> OrthoGraph::majorityRuleMove(OrthoExon *orthoex){

    vector<MoveObject*> orthomove;
  
    //count number of zeros and ones in labelpattern of an OrthoExon
  
    size_t numOnes = 0;
    size_t numZeros = 0;

  
    for(string::iterator string_it = orthoex->labelpattern.begin(); string_it < orthoex->labelpattern.end(); string_it++){
	if(*string_it == '1'){
	    numOnes++;
	}
	if(*string_it == '0'){
	    numZeros++;
	}
    }
  
    if( numOnes == 0 || numZeros == 0 ){           // nothing has to be done
    }
    else {  
	cout << "Majority Rule Move " << endl;
	cout << *orthoex << "\t" << orthoex->labelpattern << endl;
	orthomove.resize(numSpecies);
	if ( numOnes >= numZeros ){                  // make all zeros to ones	
	    for(size_t pos = 0; pos < numSpecies; pos++){
		if( (orthoex->labelpattern[pos] == '0')  ){
		    MoveObject *move = new MoveObject(graphs[pos], 3);
		    move->addNodeBack( graphs[pos]->getNode(orthoex->orthoex[pos]), graphs[pos]->getMaxWeight() );
		    move->initLocalHeadandTail();
		    orthomove[pos] = move;
		    cout << "local_head: " << orthomove[pos]->getHead() << endl;
		    cout << "local_tail: " << orthomove[pos]->getTail() << endl;
		}
	    }
	}
	else {                                       // make all ones to zeros
	    cout << "make all ones to zeros" << endl;	
	    for(size_t pos = 0; pos < numSpecies; pos++){
		if( (orthoex->labelpattern[pos] == '1')  ){
		    MoveObject *move = new MoveObject(graphs[pos], 3);
		    move->addNodeBack( graphs[pos]->getNode(orthoex->orthoex[pos]), - graphs[pos]->getMaxWeight() );
		    move->initLocalHeadandTail();
		    orthomove[pos] = move;
		    cout << "local_head: " << orthomove[pos]->getHead() << endl;
		    cout << "local_tail: " << orthomove[pos]->getTail() << endl;
		}
	    }
	}
    }
    return orthomove;
}


map<string, Score> cache::labelscore; //stores score of prunning algorithm for each pattern (leaf labelling)

bool cache::inHash(string key){
    return ( labelscore.find(key) != labelscore.end() );
}

void cache::addToHash(string key, double score){
    Score s;
    s.treescore = score;
    labelscore[key] = s;

}

void cache::resetCounter(){

    for(map<string, Score>::iterator it = labelscore.begin(); it != labelscore.end(); it++){
	it->second.count = 0;
    }

}

void cache::printCache(list<OrthoExon> &ortho){

    resetCounter();

    cout << "*************************************************************************" << endl;
    cout << "--- orthologous exons + labelpattern ---" << endl;
    for(list<OrthoExon>::iterator it = ortho.begin(); it != ortho.end(); it++){
	incrementCounter(it->labelpattern);
	//cout << *it << "\t" << it->labelpattern << endl;
    }
    cout << "\n--- cache summary ---" << endl;
    cout.width(4); cout << "";
    cout.width(12); cout << "score";
    cout.width(4); cout << "#" << endl;

    for(map<string, Score>::iterator it = labelscore.begin(); it != labelscore.end(); it++){

	cout.width(4); cout << it->first;
	cout.width(12); cout << it->second.treescore;
	cout.width(4); cout << it->second.count << endl;

    }
    cout << "*************************************************************************" << endl;
}

double cache::getScore(string key){
    return labelscore[key].treescore;
}
void cache::incrementCounter(string key){
    labelscore[key].count++;
}
