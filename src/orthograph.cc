/**********************************************************************
 * file:    orthograph.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  orthologous graphs for comparative gene prediction
 * authors: Stefanie Koenig
 *
 * date    |   author       |  changes
 * --------|----------------|------------------------------------------
 * 23.03.12|Stefanie Koenig | creation of the file
 **********************************************************************/


#include "orthograph.hh"
#include "contTimeMC.hh"
#include "mea.hh"

using namespace std;

PhyloTree *OrthoGraph::tree = NULL;

size_t OrthoGraph::numSpecies;


void OrthoGraph::buildGeneList(vector< list<Transcript*>* > &genelist) {
    for (size_t pos = 0; pos < numSpecies; pos++){
	if (graphs[pos]){
	    // delete old genelist
	    if (genelist[pos]){
		delete genelist[pos];
	    }
	    list<Transcript*> *genes = new list<Transcript*>;
	    Node* current = graphs[pos]->tail;
	    Node* head =  graphs[pos]->head;
	    Node* succExon= NULL;
	    Transcript *currentGene = NULL;
	
	    // convert node labeling of graph into a list of genes (backtracking from tail)

	    State *intr = current->pred->getIntron(current);
	    if (intr){
		// this is a bit messy, but the problem is that
		// the type of a non-coding intron is wrong
		// so its hard to figure out if the intron ins between coding or non-coding exons
		if(current->pred->pred->n_type == ncintr || current->pred->pred->n_type == rncintr) // non-coding gene
		    currentGene = new Transcript();
		else
		    currentGene = new Gene();
		addIntronToGene(currentGene, intr);  
	    }
	    while (current != NULL){
		if (current == head && succExon){
		    State *intr = current->getIntron(succExon->pred);
		    if(intr)
			addIntronToGene(currentGene,intr); 
		    setGeneProperties(currentGene);
		    genes->push_front(currentGene);
		    currentGene = NULL;
		}
		if(current->n_type == IR && succExon){
		    setGeneProperties(currentGene);
		    genes->push_front(currentGene);
		    currentGene = NULL;
		    succExon = NULL;
		}
		if(current->n_type >= utrExon){ // add an exon to the current gene
		    State *ex;
		    if(current->n_type == sampled || current->n_type == utrExon){
			ex = new State(*((State*)(current->item)));
		    }
		    else{
			ex = new State(current->begin, current->end, current->castToStateType());
		    }
		    if(!currentGene){
			if(isNcExon(ex->type))
			    currentGene = new Transcript();
			else
			    currentGene = new Gene();
		    }
		    addExonToGene(currentGene, ex);
		    if(succExon){ // if the current exon is not the last, add an intron from the current exon to the succeding exon
			if(current->end+1 < succExon->begin){
			    State *intr = current->edges.begin()->to->getIntron(succExon->pred);
			    if(!intr){ // if no explicit intron exists, convert auxiliary edge to an intron
				intr = new State(current->end+1, succExon->begin-1, getIntronStateType((State*)current->item,(State*)succExon->item));
			    }
			    addIntronToGene(currentGene, intr);   
			}
		    }
		    succExon=current;
		}
		current=current->pred;
	    }
	    genelist[pos] = genes;
	}
    }
}

void OrthoGraph::filterGeneList(vector< list<Transcript*> *> &genelist, vector<int> &geneid){
    
    vector< list<AltGene> *> agls(numSpecies);

    for (size_t pos = 0; pos < numSpecies; pos++){	
	if (genelist[pos]){
	    AnnoSequence *annoseq = graphs[pos]->getAnnoSeq();

	    list<AltGene> *agl = groupTranscriptsToGenes(*genelist[pos]);

	    if(sfcs[pos] && sfcs[pos]->collection->hasHintsFile){
		// compile extrinsic evidence
		for (list<AltGene>::iterator git = agl->begin(); git != agl->end(); git++) {
		    for (list<Transcript*>::iterator trit = git->transcripts.begin(); trit != git->transcripts.end(); trit++) {
			Gene *g = dynamic_cast<Gene*> (*trit);
			if (g)
			    g->compileExtrinsicEvidence(sfcs[pos]->groupList);
		    }
		}
	    }

	    agl->sort();

	    int transcriptid;

	    // shift gene coordinates, set sequence name, gene and transcript names
	    for (list<AltGene>::iterator agit = agl->begin(); agit != agl->end(); ++agit){
		agit->seqname =  annoseq->seqname;
		agit->id = "g"+itoa(geneid[pos]);
		agit->sortTranscripts();

		transcriptid = 1;
		for (list<Transcript*>::iterator it = agit->transcripts.begin();it != agit->transcripts.end(); ++it ) {
		    (*it)->seqname = annoseq->seqname;
		    (*it)->id = "t" + itoa(transcriptid);
		    (*it)->geneid = agit->id;
		    transcriptid++;
		}
		geneid[pos]++;
	    }
	    geneLists[pos]=agl;
	}
    }
}
    

void OrthoGraph::printGenelist(vector<ofstream*> &filestreams){

    for (size_t pos = 0; pos < numSpecies; pos++){	
	if (geneLists[pos]){
	    AnnoSequence *annoseq = graphs[pos]->getAnnoSeq();
	    Strand strand = graphs[pos]->getSeqStrand();

	    bool withEvidence = false;
	    if(sfcs[pos] && sfcs[pos]->collection->hasHintsFile)
		withEvidence = true;

	    list<AltGene> *agl = geneLists[pos];
	    
	    if (strand == minusstrand){
 		agl = reverseGeneList(agl, annoseq->length - 1);
		agl->sort();
	    }

	    // shift gene coordinates
	    for (list<AltGene>::iterator agit = agl->begin(); agit != agl->end(); ++agit)
		agit->shiftCoordinates(annoseq->offset);
	    
	    //print the genes
	    streambuf *coutbuf = cout.rdbuf(); //save old buf
	    cout.rdbuf(filestreams[pos]->rdbuf()); //redirect std::cout to species file

	    //print sequence information

	    cout << "#----- prediction on sequence range " << annoseq->seqname << ":" << annoseq->offset + 1  << "-" << annoseq->offset + annoseq->length << " (" << annoseq->length << "bp) -----" << endl << "#" << endl;

	    if(!agl->empty()){
		if(strand == minusstrand){
		    char *DNA = annoseq->sequence;
		    char *reverseDNA = reverseComplement(DNA);
		    annoseq->sequence = reverseDNA;
		    printGeneList(agl, annoseq, Constant::codSeqOutput, Constant::proteinOutput, withEvidence);
		    annoseq->sequence = DNA;
		    delete [] reverseDNA;
		}
		else{
		    printGeneList(agl, annoseq, Constant::codSeqOutput, Constant::proteinOutput, withEvidence);
		}
	    }
	    else{
		cout << "# (none)" << endl;
	    }
	    cout.rdbuf(coutbuf); //reset to standard output again
	}
    }
}



vector<ofstream*> initOutputFiles(string outdir, string extension){

    vector<ofstream*> filestreams;
    filestreams.resize(OrthoGraph::numSpecies);
    vector<string> species;
    OrthoGraph::tree->getSpeciesNames(species);
    for(size_t pos = 0; pos < OrthoGraph::numSpecies; pos++){
	string filename = outdir + species[pos] + extension + ".gff";
	if(Gene::gff3){
	    filename += "3";	    
	}
	ofstream *out = new ofstream(filename.c_str());
	if(out){
	    filestreams[pos] = out;
	    (*out) << PREAMBLE << endl;
	    (*out) << "#\n#----- prediction for species '" << species[pos] << "' -----" << endl << "#" << endl;
	}	
    }
    return filestreams;
}

void closeOutputFiles(vector<ofstream*> filestreams){

    for(size_t pos = 0; pos < OrthoGraph::numSpecies; pos++){
	if(filestreams[pos]){
	    if(filestreams[pos]->is_open()){
		filestreams[pos]->close();
		delete filestreams[pos];
	    }
	}
    }
}

void OrthoGraph::addScoreSelectivePressure(){

    double a;
    double b;
    try {
	a = Properties::getdoubleProperty("/CompPred/ec_addend");
    } catch (...) {
	a = 0;
    }
    try {
	b = Properties::getdoubleProperty("/CompPred/ec_factor");
    } catch (...) {
	b = 0;
    }
    
    // reward/penalty that each EC receives
    for(size_t pos = 0; pos < numSpecies; pos++){
	if(graphs[pos]){
	    for(list<Node*>::iterator node = graphs[pos]->nodelist.begin(); node != graphs[pos]->nodelist.end(); node++){
		if( (*node)->n_type >= sampled ){
		    for(list<Edge>::iterator edge = (*node)->edges.begin(); edge != (*node)->edges.end(); edge++){
			// default EC-filter on:
			edge->score += a + (-1.396*b);
			}
		}
	    }
	}
    }
    // reward/penalty that only EC receives which are part of an OE
    static bool detection = false;
    if(!all_orthoex.empty()){
	for(list<OrthoExon>::const_iterator it = all_orthoex.begin(); it != all_orthoex.end(); it++){
	    for(size_t pos = 0; pos < numSpecies; pos++){
		if(it->orthoex[pos]){
		    Node* node = it->orthonode[pos];
		    int len =  node->end - node->begin + 1;
 		    double x = b*(0.015 * len + 1.277 * it->getConsScore() + 4.37 * it->getDiversity() -0.92);
		    if (string("honeybee1") == Properties::getProperty("species") && !detection){
			cout << "honeybee1 detected" << endl;
			detection = true;
			double Eomega = it->getEomega(),
			    cons = it->getConsScore(),
			    cont = it->getContainment(),
			    s = it->numExons(),
			    len = it->getAliLen();
			x = b*(1.396 - 6.585 -3.80 * Eomega + 0.860 * Eomega * Eomega -5.147 * cons -0.0084 * cont + 1.446 * s + 1.245 * log(len));
		    }
		    for (list<Edge>::iterator iter =  node->edges.begin(); iter != node->edges.end(); iter++){
			// default EC-filter on:  
		       	iter->score += x;
		    }
		}
	    }
	}
    }
}

double OrthoGraph::globalPathSearch(){

    double score=0;

    for(size_t pos = 0; pos < numSpecies; pos++){
	if(graphs[pos]){
	    score += graphs[pos]->relax();
	}
    }
    return score;
}

double OrthoGraph::dualdecomp(ExonEvo &evo, vector< list<Transcript*> *> &genelist, int gr_ID, int T, vector<double> &c){

    cout << "dual decomposition on gene Range " << gr_ID << endl;
    cout<<"round\titer\tstep_size\tprimal\tdual\t#inconsistencies"<<endl;

    double best_dual = std::numeric_limits<double>::max();   // best dual value so far
    double best_primal = -std::numeric_limits<double>::max(); // best primal value so far
    double initial_primal = std::numeric_limits<double>::max(); // initial primal value (in 0-the iteration)
	
    for(size_t r=0; r < c.size(); r++){ // number of rounds
	
	/*
	 * initialization
	 */
	double old_dual = 0.0; // dual value of the previous iteration
	double delta = 0.0; // step size	
	int v = 0; 	// number of iterations prior to t where the dual value increases

	cout.precision(10);
	for(int t=0; t<T;t++){
	    double path_score = globalPathSearch();
	    int numInconsistent = 0;
	    double hect_score = 0;
	    if(t == 0){
		// initially set labels of hects by a simple majority rule
		// set all labels of a hect to 1, if the majority
		// iff the corresponding labels in the graph are labelled with one/
		// otherwise, set all labels to 0
		hect_score += init(evo,numInconsistent);
	    }
	    else{
		hect_score += treeMAPInf(evo,numInconsistent);
	    }
	    double current_dual = path_score + hect_score;       // dual value of the t-th iteration 
	    best_dual = min(best_dual,current_dual);              // update upper bound
	    if( (t >= 1) && (old_dual < current_dual) )  // update v
		v++;
	    double current_primal = path_score + makeConsistent(evo); // primal value of the t-the iteration
	    if(best_primal < current_primal){
		best_primal = current_primal;
		buildGeneList(genelist); // save new record
	    }
	    if(t == 0 && r == 0){
		initial_primal = current_primal;
		//c*=initial_gap/numInconsistent; // adjust step size parameter to problem size
	    }
	    cout<<r<<"\t"<<t<<"\t"<<delta<<"\t"<<current_primal<<"\t"<<current_dual<<"\t"<<numInconsistent<<endl;
	    
	    if(numInconsistent == 0 || abs(best_dual - best_primal) < 1e-8){ // exact solution is found
		double best_gap = abs(best_dual - best_primal);
		double initial_gap = abs(best_dual - initial_primal);
		double perc_gap = (initial_gap > 0 )? best_gap/initial_gap : 0;
		cout<<"dual decomposition reduced initial duality gap of "<<initial_gap<<" to "<<best_gap<<" (to "<<perc_gap<<"%)"<<endl;
		return best_gap;
	    }
	    
	    // determine new step size
	    delta = getStepSize(c[r],t,v);
	    
	    // updated weights
	    for(list<OrthoExon>::iterator hects = all_orthoex.begin(); hects != all_orthoex.end(); hects++){
		for(size_t pos = 0; pos < numSpecies; pos++){
		    if(hects->orthoex[pos]){
			// get corresponding node in graph
			Node* node = hects->orthonode[pos];
			bool h = node->label;
			bool v = hects->labels[pos];
			if(v != h){  //shared nodes are labelled inconsistently in the two subproblems
			    float weight = delta*(v-h);
			    //update weights
			    node->addWeight(weight);
			    hects->weights[pos] -= weight;     
			}
		    }
		}
	    }
	    old_dual = current_dual;
	}
	// reset weights
	for(list<OrthoExon>::iterator hects = all_orthoex.begin(); hects != all_orthoex.end(); hects++){
	    for(size_t pos = 0; pos < numSpecies; pos++){
		if(hects->orthoex[pos]){
		    Node* node = hects->orthonode[pos];
		    node->addWeight(hects->weights[pos]);
		    hects->weights[pos] = 0;     
		}
	    }
	}
    }
    double best_gap = abs(best_dual - best_primal);
    double initial_gap = abs(best_dual - initial_primal);
    double perc_gap = (initial_gap > 0 )? best_gap/initial_gap : 0;
    cout<<"dual decomposition reduced initial duality gap of "<<initial_gap<<" to "<<best_gap<<" (to "<<perc_gap<<"%)"<<endl;
    return best_gap;
}

/*
 * requirement delta -> 0 for t -> infinity and
 * sum(deltas) = infinity
 * v ist the number of iterations prior to t where the dual value increases
 * the purpose of v is to decrease the step size only if we move in the wrong direction
 */
double OrthoGraph::getStepSize(double c, int t, int v){    
    return c/sqrt(v+1);
}

double OrthoGraph::init(ExonEvo &evo, int &numInconsistent){

    double score = 0;
    double k =evo.getPhyloFactor(); //scaling factor 

    for(list<OrthoExon>::iterator hects = all_orthoex.begin(); hects != all_orthoex.end(); hects++){
	//count number of zeros   
	size_t numOnes = 0;
	size_t numZeros = 0;

	for(int pos=0; pos < hects->orthonode.size(); pos++){
	    if(hects->orthonode[pos])
		(hects->orthonode[pos]->label == 1) ? numOnes++ : numZeros++;
	}
	size_t majority = 0;
	(numOnes >= numZeros) ? majority=1 : majority=0;


	for(int pos=0; pos < hects->orthoex.size(); pos++){
	    if(hects->orthoex[pos])
		(hects->labels[pos] = majority);
	}
	PhyloTree *temp = hects->getTree();
	Evo *evo_base = &evo;
	score += temp->MAP(hects->labels, hects->weights, evo_base, k, true);
	double maxi =  max(numOnes,numZeros);
	numInconsistent += ((numOnes+numZeros) - maxi);
	
    }
    return score;
    
}

double OrthoGraph::treeMAPInf(ExonEvo &evo, int &numInconsistent){

    double score=0;
    double k =evo.getPhyloFactor(); //scaling factor 

    for(list<OrthoExon>::iterator hects = all_orthoex.begin(); hects != all_orthoex.end(); hects++){
	PhyloTree *temp = hects->getTree();
	Evo *evo_base = &evo;
	score += temp->MAP(hects->labels, hects->weights, evo_base, k);
	for(int pos=0; pos < hects->orthonode.size(); pos++){
            if(hects->orthonode[pos] && hects->orthonode[pos]->label != hects->labels[pos])
                numInconsistent++;
        }

    }
    return score;
}

double OrthoGraph::makeConsistent(ExonEvo &evo){

    double score = 0;
    double k =evo.getPhyloFactor(); //scaling factor 

    for(list<OrthoExon>::iterator hects = all_orthoex.begin(); hects != all_orthoex.end(); hects++){
	PhyloTree *temp = hects->getTree();
	Evo *evo_base = &evo;
	vector<int> labels(numSpecies,2);
	for(int pos=0; pos < hects->orthonode.size(); pos++){
	    if(hects->orthonode[pos])
		labels[pos] = hects->orthonode[pos]->label;
	}
	score += temp->MAP(labels, hects->weights, evo_base, k, true);
    }
    return score;
}	

void OrthoGraph::linkToOEs(list<OrthoExon> &oes){

    all_orthoex = oes;
    for(list<OrthoExon>::iterator it = all_orthoex.begin(); it != all_orthoex.end(); it++){
	double oe_score = it->getLogRegScore();
	for(size_t pos = 0; pos < OrthoGraph::numSpecies; pos++){ 
	    if(it->orthoex[pos]==NULL){
		it->labels[pos]=2;
	    }
	    else{
		if(!graphs[pos])
		    throw ProjectError("Internal error OrthoGraph::linkToOEs: graph does not exist.");
		Node* node = graphs[pos]->getNode(it->orthoex[pos]);
		if(!node){
		    throw ProjectError("Internal error OrthoGraph::linkToOEs: EC has no corrpesonding node in OrthoGraph.");
		}
		it->orthonode[pos]=node;
		node->addWeight(oe_score); // add OE score to all outgoing edges
	    }
	    it->setLabelpattern();
	}	
    }  
}

/* old code: optimize cgp by making small moves
void OrthoGraph::optimize(ExonEvo &evo){

    //create MoveObjects

    //for now, loop over OrthoExons and apply majority rule move if possible

    int shift_size=3;     // number of exons local_head/local_tail is shifted to the left/right on the current path 
    if(!all_orthoex.empty()){
	for(list<OrthoExon>::iterator orthoex = all_orthoex.begin(); orthoex != all_orthoex.end(); orthoex++){
	    vector<Move*> orthomove = majorityRuleMove(*orthoex, shift_size);
	    if(!orthomove.empty()){
		localMove(orthomove,evo,shift_size);
		//delete MoveObjects
		for(size_t pos = 0; pos < numSpecies; pos++){
		    delete orthomove[pos];
		} 
	    }   
	}
    }
}

void OrthoGraph::localMove(vector<Move*> &orthomove, ExonEvo &evo, int shift_size){

    int maxIterations=0;  // max number a move can be repeated                                                                                                                                         
    bool retry = false;  //if true, the move is repeated on a 'larger' subgraph
    int iter = 0;  // number of current repetition
   
    do{

	if(retry){

	    cout << "repetition:\t" << iter <<"\tshift_size:\t"<<shift_size <<endl;

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
	phylo_score -= pruningAlgor(local_orthoexons,evo);

	// do local changes for each graph  
	for(size_t pos = 0; pos < numSpecies; pos++){
	    if(orthomove[pos]){
		graph_score += graphs[pos]->localChange(orthomove[pos]);
	    }
	}
	print_change = true;
	//calculate new phylo_score
	phylo_score += pruningAlgor(local_orthoexons, evo);
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
	    // no improvement, undo local changes by relaxing nodes again this time with original weights
	    // TODO: has to be replaced by an undo buffer or something similar in future
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

double OrthoGraph::pruningAlgor(list<OrthoExon> &orthoex, ExonEvo &evo){

    double tree_score = 0;
    Evo *evo_base = &evo;

    for(list<OrthoExon>::iterator ortho = orthoex.begin(); ortho != orthoex.end(); ortho++){

	string labelpattern = getLabelpattern(*ortho);
	if(cache::inHash(labelpattern)){
	    tree_score += cache::getScore(labelpattern);
	}
	else{
	    PhyloTree temp(*tree);
	    double score = temp.pruningAlgor(labelpattern, evo_base);
	    score=score * evo.getPhyloFactor();
	    cache::addToHash(labelpattern, score);
	    tree_score += score;
	}
    }
    return tree_score;
}

string OrthoGraph::getLabelpattern(OrthoExon &ex){

    string labelpattern;
    for (size_t i = 0; i < ex.orthoex.size(); i++){

	if (ex.orthoex.at(i) == NULL){
	    labelpattern += "2";
	}
	else {
	    map<string, Node*>::iterator it = graphs.at(i)->existingNodes.find(graphs.at(i)->getKey(ex.orthoex.at(i)));
	    if (it != graphs.at(i)->existingNodes.end()){
		bool label = it->second->label;
		if (label == 1){
		    labelpattern += "1";
		}
		else if (label == 0){
		    labelpattern += "0";
		}
	    }
	    else
		throw ProjectError("Error in OrthoExon::getKey: exon " + graphs.at(i)->getKey(ex.orthoex.at(i)) +" not in graph!");
	}
    }
    if ( ex.labelpattern != labelpattern){
	if(print_change){
	    cout << ex << "\t" << ex.labelpattern << "-->" << labelpattern << endl;
	    //temp: html output for gBrowse
	    printHTMLgBrowse(ex);
	}
    }

    ex.labelpattern = labelpattern;
    return ex.labelpattern;
}


list<OrthoExon> OrthoGraph::orthoExInRange(vector<Move*> &orthomove){

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
	
vector<Move*> OrthoGraph::majorityRuleMove(OrthoExon &orthoex, int shift_size){

    vector<Move*> orthomove;
  
    //count number of zeros and ones in labelpattern of an OrthoExon
  
    size_t numOnes = 0;
    size_t numZeros = 0;

    //get current label pattern
    string labelpattern = getLabelpattern(orthoex);

  
    for(string::iterator string_it = labelpattern.begin(); string_it < labelpattern.end(); string_it++){
	if(*string_it == '1'){
	    numOnes++;
	}
	if(*string_it == '0'){
	    numZeros++;
	}
    }
    if ( numOnes > numZeros && numZeros > 0 ){                    // make all zeros to ones	
	cout << "Majority Rule Move " << endl;
	cout << orthoex << "\t" << labelpattern << endl;
	cout << "iteration:\t" << 0  << "\tshift_size:\t" << shift_size << endl;
	orthomove.resize(numSpecies);
	for(size_t pos = 0; pos < numSpecies; pos++){
	    if( (labelpattern[pos] == '0')  ){
		Move *move = new Move(graphs[pos], shift_size + 1);
		move->addNodeBack( graphs[pos]->getNode(orthoex.orthoex[pos]), graphs[pos]->getMaxWeight() );
		move->initLocalHeadandTail();
		orthomove[pos] = move;
		cout << "local_head:\t";
		graphs[pos]->printNode(orthomove[pos]->getHead());
		cout << "local_tail:\t";
		graphs[pos]->printNode(orthomove[pos]->getTail());
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


void OrthoGraph::printCache(list<OrthoExon> &ortho){

    cache::resetCounter();

    cout << "*************************************************************************" << endl;
    cout << "--- orthologous exons + labelpattern ---" << endl;
    for(list<OrthoExon>::iterator it = ortho.begin(); it != ortho.end(); it++){
	cache::incrementCounter(getLabelpattern(*it));
	cout << *it << "\t" << getLabelpattern(*it)  << endl;
	printHTMLgBrowse(*it);
    }
    cout << "\n--- cache summary ---" << endl;
    cout.width(4); cout << "";
    cout.width(12); cout << "score";
    cout.width(4); cout << "#" << endl;

    for(map<string, Score>::iterator it = cache::labelscore.begin(); it != cache::labelscore.end(); it++){

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

void OrthoGraph::printHTMLgBrowse(OrthoExon &ex){

    for (int j=0; j<ex.orthoex.size(); j++) {
	if (ex.orthoex.at(j)!=NULL) {

	    int offset = graphs[j]->getSeqOffset();
	    Strand strand = graphs[j]->getSeqStrand();
	    int length = graphs[j]->getSeqLength();
	    char* seqname = graphs[j]->getSeqID();

	    cout << "http://bioinf.uni-greifswald.de/gb2/gbrowse/vergl" << itoa(j+1) <<"/?name=" << seqname  << "%3A";
	    if (strand == plusstrand) {
		cout << ex.orthoex.at(j)->begin + offset + 1 << ".." << ex.orthoex.at(j)->end + offset + 1 << endl;
	    } else {
		cout << length - ex.orthoex.at(j)->end + offset << ".." << length - ex.orthoex.at(j)->begin + offset << endl;
	    }
	}
    }
}
*/

/* 
 * Patrick Balmerths playground
 */
void OrthoGraph::createOrthoGenes(const GeneMSA *geneRange){
    // create all_orthogenes

    // In case the alignment in geneRange is actually not needed and
    // the orthoexons are sufficient, then remove this argument.
}

void OrthoGraph::printOrthoGenes(){
  // ouput all_orthogenes
}
