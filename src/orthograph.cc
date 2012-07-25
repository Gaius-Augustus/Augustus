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
    ptrs_to_alltranscripts.resize(numSpecies);
}



OrthoGraph::~OrthoGraph(){
    for(int i = 0; i < numSpecies; i++){
	if(orthoSeqRanges[i]){
	    delete graphs[i];
	    delete ptrs_to_alltranscripts[i];
	    delete orthoSeqRanges[i];
	}
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

	if(graphs[pos]){

	    list<Gene> *genes = new list<Gene>;
	
	    Node* current = graphs[pos]->tail;
	    Node* head =  graphs[pos]->head;
	    Node* predcurrent;
	
	    Gene *currentGene = new Gene();
	
	    // convert node labeling of graph into a list of genes
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
	    cout.rdbuf(filestreams[pos]->rdbuf()); //redirect std::cout to species file

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

    //create MoveObjects

    //for now, loop over OrthoExons
    if(!all_orthoex.empty()){
	for(list<OrthoExon>::iterator orthoex = all_orthoex.begin(); orthoex != all_orthoex.end(); orthoex++){
	    vector<MoveObject*> orthomove = majorityRuleMove(*orthoex);
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


    bool retry = false;  //if true, the move is repeated on a 'larger' subgraph

    int maxIterations = 3;  // max number a move can be repeated
    int iter = 0;
    int shift_size = 2;  // number of nodes local_head/local_tail is shifted to the left/right on the current path
   
    do{

	if(retry){

	    cout << "iteration:\t" << iter <<"\tshift_size:\t"<<shift_size <<endl;

	    //shift local_heads and local_tails for each species
	    for(size_t pos = 0; pos < numSpecies; pos++){
		if(orthomove[pos]){
		    orthomove[pos]->shiftHead(shift_size);
		    orthomove[pos]->shiftTail(shift_size);
		    cout << graphs[pos]->getSpeciesname() << endl;
		    cout << "local_head:\t" << orthomove[pos]->getHead() << endl;
		    cout << "local_tail:\t" << orthomove[pos]->getTail() << endl;
		}
	    }
	}

	retry = false;

	double graph_score = 0; //difference between the new and the old scores of the local paths
	double phylo_score = 0; //diference between the new and the old phylogenetic score

	// determine all OrthoExons in that range
	
	list<OrthoExon> local_orthoexons = orthoExInRange(orthomove);
 
	//calculate phylo_score

	phylo_score -= pruningAlgor(local_orthoexons);

	// do local changes for each graph  
	for(size_t pos = 0; pos < numSpecies; pos++){
	    if(orthomove[pos]){
		graph_score += graphs[pos]->localChange(orthomove[pos]);
	    }
	}
	print_change = true;
	//calculate new phylo_score
	phylo_score += pruningAlgor(local_orthoexons);
	print_change = false;

	cout << "-------------------------------------" << endl;
	cout << "graph_score\t" << graph_score << endl;
	cout << "phylo_score\t" << phylo_score << endl;

	if( graph_score == 0 && phylo_score == 0 ){
	    //nothing changed, repeat move on larger subgraph
	    if( iter < maxIterations ){
		cout << "repeat\n" << endl;
		retry = true;
	    }
	}
	
	else if( (graph_score + phylo_score) > 0 ){
	    //score improved, accept move
	    cout << "accept\n" << endl; 
	}
	else{
	    cout << "undo\n" << endl;
	    // no improvement, undo local changes
	    for(size_t pos = 0; pos < numSpecies; pos++){
		if(orthomove[pos]){
		    graphs[pos]->relax(orthomove[pos]->getHead(), orthomove[pos]->getTail());
		}
	    }
	}
	iter++;
    }
    while( retry );
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

    string old_labelpattern = ex.labelpattern;

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
    if(!old_labelpattern.empty() && old_labelpattern != ex.labelpattern && print_change ){
	cout << ex << "\t" + old_labelpattern + " --> " + ex.labelpattern << endl;
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
	
vector<MoveObject*> OrthoGraph::majorityRuleMove(OrthoExon &orthoex){

    vector<MoveObject*> orthomove;

    int shift_size = 5;
  
    //count number of zeros and ones in labelpattern of an OrthoExon
  
    size_t numOnes = 0;
    size_t numZeros = 0;

    //get current label pattern, important to update this !!!!
    string labelpattern = getLabelpattern(orthoex);

  
    for(string::iterator string_it = labelpattern.begin(); string_it < labelpattern.end(); string_it++){
	if(*string_it == '1'){
	    numOnes++;
	}
	if(*string_it == '0'){
	    numZeros++;
	}
    }
  
    if( numOnes == 0 || numZeros == 0 ){ // nothing has to be done
    }
    else {  
	cout << "Majority Rule Move " << endl;
	cout << orthoex << "\t" << orthoex.labelpattern << endl;
	cout << "iteration:\t" << 0  <<"\tshift_size:\t" <<shift_size<<endl;
	orthomove.resize(numSpecies);
	if ( numOnes >= numZeros ){   // make all zeros to ones	
	    for(size_t pos = 0; pos < numSpecies; pos++){
		if( (orthoex.labelpattern[pos] == '0')  ){
		    MoveObject *move = new MoveObject(graphs[pos],shift_size);
		    move->addNodeBack( graphs[pos]->getNode(orthoex.orthoex[pos]), graphs[pos]->getMaxWeight() );
		    move->initLocalHeadandTail();
		    orthomove[pos] = move;
		    cout << "local_head:\t" << orthomove[pos]->getHead() << endl;
		    cout << "local_tail:\t" << orthomove[pos]->getTail() << endl;
		}
	    }
	}
	else {   // make all ones to zeros	
	    for(size_t pos = 0; pos < numSpecies; pos++){
		if( (orthoex.labelpattern[pos] == '1')  ){
		    MoveObject *move = new MoveObject(graphs[pos], shift_size);
		    move->addNodeBack( graphs[pos]->getNode(orthoex.orthoex[pos]), - graphs[pos]->getMaxWeight() );
		    move->initLocalHeadandTail();
		    orthomove[pos] = move;
		    cout << "local_head:\t" << orthomove[pos]->getHead() << endl;
		    cout << "local_tail:\t" << orthomove[pos]->getTail() << endl;
		}
	    }
	}
    }
    return orthomove;
}

map<string, Score> cache::labelscore; //stores score of prunning algorithm for each pattern (leaf labeling)

bool cache::inHash(string labelpattern){
    return ( labelscore.find(labelpattern) != labelscore.end() );
}

void cache::addToHash(string labelpattern, double score){
    Score s;
    s.treescore = score;
    labelscore[labelpattern] = s;

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

double cache::getScore(string labelpattern){
    return labelscore[labelpattern].treescore;
}
void cache::incrementCounter(string labelpattern){
    labelscore[labelpattern].count++;
}
