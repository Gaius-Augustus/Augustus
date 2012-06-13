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
#include "mea.hh"



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

void OrthoGraph::outputGenes(Strand strand){

  Boolean noInFrameStop;
  double minmeanexonintronprob;
  double minexonintronprob;  // lower bound on probabilities of all exons and introns in the coding region

  try {
      noInFrameStop = Properties::getBoolProperty("noInFrameStop");
  } catch (...) {
      noInFrameStop = false;
  }
  try {
    minmeanexonintronprob = Properties::getdoubleProperty("minmeanexonintronprob");
  } catch (...) {
    minmeanexonintronprob = 0.0;
  }
  try {
    minexonintronprob = Properties::getdoubleProperty("minexonintronprob");
  } catch (...) {
    minexonintronprob = 0.0;
  }

  for (size_t pos = 0; pos < numSpecies; pos++){
    /*
     * backtracking
     */
    list<Node*> longest_path; //stores path nodes
    Node *back = this->graphs[pos]->tail; //pointer to tail
    longest_path.push_front(back);
    while(back->pred != NULL){
      longest_path.push_front(back->pred);
      back = back->pred;
    }
    /*cout <<"longest path: "<< endl;
    for(list<Node*>::iterator it = longest_path.begin(); it != longest_path.end(); it++){
      cout << *it << endl;
      }*/
    list<Gene> *genes = new list<Gene>;
    /*
     * convert list of path nodes to list of genes
     */
    getMeaGenelist7(longest_path, genes);
   
    list<Gene> *filteredTranscripts = new list<Gene>;
    filteredTranscripts = Gene::filterGenePrediction(genes, this->orthoSeqRanges[pos]->sequence, strand, noInFrameStop, minmeanexonintronprob, minexonintronprob);
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

  /*
   * for now, loop ofer all OrthoExons,  if # of 0's in labelpattern of an OrthoExon is smaller than # of 1's,
   * all nodes with the labels 0 are made to 1 and vice versa
   * this can be changed later, when it's more clear, what good changes are
   */

  cout << "before local change" << endl;
  cache::printCache(all_orthoex);

  //loop over OrthoExons
  for(list<OrthoExon>::iterator orthoex = all_orthoex.begin(); orthoex != all_orthoex.end(); orthoex++){

    cout << *orthoex << endl;
    cout << orthoex->labelpattern << endl;

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
      cout << "nothing has to be done" << endl;
    }
    else {

      vector<MoveObject*> orthomove(numSpecies);

      if ( numOnes >= numZeros ){                  // make all zeros to ones	
	cout << "make all zeros to ones" << endl;
	for(size_t pos = 0; pos < numSpecies; pos++){
	  if( (orthoex->labelpattern[pos] == '0')  ){
	    MoveObject *move = new MoveObject(graphs[pos]);
	    move->addNode( graphs[pos]->getNode(orthoex->orthoex[pos]), graphs[pos]->getMaxWeight() );
	    move->initLocalHeadandTail(2);
	    orthomove[pos] = move;
	    cout << "local_head: " << orthomove[pos]->getHead() << endl;
	    cout << "local_tail: " << orthomove[pos]->getTail() << endl;
	  }
	}
	localMove(orthomove);
      }
      else {                                       // make all ones to zeros
	cout << "make all ones to zeros" << endl;	
	for(size_t pos = 0; pos < numSpecies; pos++){
	  if( (orthoex->labelpattern[pos] == '1')  ){
	    MoveObject *move = new MoveObject(graphs[pos]);
	    move->addNode( graphs[pos]->getNode(orthoex->orthoex[pos]), - graphs[pos]->getMaxWeight() );
	    move->initLocalHeadandTail(2);
	    orthomove[pos] = move;
	    cout << "local_head: " << orthomove[pos]->getHead() << endl;
	    cout << "local_tail: " << orthomove[pos]->getTail() << endl;
	  }
	}
	localMove(orthomove);
      }
      //delete MoveObjects
      for(size_t pos = 0; pos < numSpecies; pos++){
	delete orthomove[pos];
      } 
    }
  }
}

void OrthoGraph::localMove(vector<MoveObject*> &orthomove){


  // number of iterations, the local bounds are expanded and the move is repeated,
  // makes sure, that the move is not kept too locally and get's stuck in local optimum
  size_t maxIterations = 3;  

  // the magnitude of the expansion in each iteration. The number of nodes, the local_heads are shifted left
  // and the local_tails are shifted right on the current path.
  size_t step_size =1;

  for(size_t iter = 0; iter < maxIterations; iter++){

    cout << "iteration: " << iter << endl;

    if(iter > 0){

      //shift local_heads and local_heads
      for(size_t pos = 0; pos < numSpecies; pos++){
	if(orthomove[pos]){
	  orthomove[pos]->goLeftOnPath(step_size);
	  orthomove[pos]->goRightOnPath(step_size);
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
	  graphs[pos]->undoChange(orthomove[pos]);
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
  addOrthoIntrons(orthomove,local_orthoexons );
 
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

  double score = 0;

  for(list<OrthoExon>::iterator ortho = orthoex.begin(); ortho != orthoex.end(); ortho++){
    score += tree->pruningAlgor(*ortho, *this);
  }
  return score;

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

  
  for(size_t pos = 0; pos < numSpecies; pos++){
    if(orthomove[pos]){
      if(!orthomove[pos]->nodes.empty()){
	Node *begin = orthomove[pos]->nodes.front().node;
	if(begin->n_type > minus2){
	  cout << "MoveNode: " << begin << endl;
	  string key1 = itoa(begin->end) + ":" + itoa( (NodeType)(graphs[pos]->toNeutralLine(begin)));
	  Node* node1 = graphs[pos]->existingNodes[key1];
	  cout << "neutraler Knoten: " << node1 <<  endl;
	  for(list<OrthoExon>::const_iterator ortho = local_orthoexons.begin(); ortho != local_orthoexons.end(); ortho++){
	    Node *end = graphs[pos]->getNode(ortho->orthoex[pos]);
	    string key2 = itoa(end->begin) + ":" + itoa( (NodeType)(graphs[pos]->fromNeutralLine(end)));
	    Node* node2 = graphs[pos]->existingNodes[key2];
	    cout << "OrthoEx: " << end << endl;
	    cout << "neutraler Knoten in or: " << node2 <<  endl;
	    for(list<Edge>::iterator it = node1->edges.begin(); it!= node1->edges.end(); it++){
	      if(it->to == node2){
		cout <<"Kante gefunden: "<< *it << endl;
		orthomove[pos]->addNode(end, graphs[pos]->getMaxWeight());
		orthomove[pos]->addEdge(&(*it), graphs[pos]->getMaxWeight());
		
	      }
	    }
	  }
	}
      }
    }
  }
  for(size_t pos = 0; pos < numSpecies; pos++){
    if(orthomove[pos]){
      if(!orthomove[pos]->nodes.empty()){
	Node *begin = orthomove[pos]->nodes.front().node;
	if(begin->n_type > minus2){
	  cout << "MoveNode: " << begin << endl;
	  string key1 = itoa(begin->begin) + ":" + itoa( (NodeType)(graphs[pos]->fromNeutralLine(begin)));
	  Node* node1 = graphs[pos]->existingNodes[key1];
	  cout << "neutraler Knoten: " << node1 <<  endl;
	  for(list<OrthoExon>::const_iterator ortho = local_orthoexons.begin(); ortho != local_orthoexons.end(); ortho++){
	    Node *end = graphs[pos]->getNode(ortho->orthoex[pos]);
	    string key2 = itoa(end->end) + ":" + itoa( (NodeType)(graphs[pos]->toNeutralLine(end)));
	    Node* node2 = graphs[pos]->existingNodes[key2];
	    cout << "OrthoEx: " << end << endl;
	    cout << "neutraler Knoten: " << node2 <<  endl;
	    for(list<Edge>::iterator it = node2->edges.begin(); it!= node2->edges.end(); it++){
	      if(it->to == node1){
		cout <<"Kante gefunden: "<< *it << endl;
		orthomove[pos]->addNodeFront(end, graphs[pos]->getMaxWeight());
		orthomove[pos]->addEdgeFront(&(*it), graphs[pos]->getMaxWeight());
	      }
	    }
	  }
	}
      }
    }
  }
}
