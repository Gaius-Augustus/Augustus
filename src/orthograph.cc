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
#include "meaPath.hh"
#include "mea.hh"
#include "namgene.hh"



PhyloTree *OrthoGraph::tree = NULL;
vector<ofstream*> OrthoGraph::filestreams;


OrthoGraph::OrthoGraph(RandSeqAccess *rsa){

  //TODO: für GeneMSA Object umschreiben

  NAMGene namgene; // creates and initializes the states
  FeatureCollection extrinsicFeatures; // hints, empty for now, will later read in hints for sequence ranges from database
  SequenceFeatureCollection sfc(&extrinsicFeatures); 
  StateModel::readAllParameters(); // read in the parameter files: species_{igenic,exon,intron,utr}_probs.pbl

  score = 0.0;
  graphs.resize(tree->species.size());
  orthoSeqRanges.resize(tree->species.size());

  //temp:
  all_orthoex = readOrthoExons(Constant::orthoexons);
  vector<string> speciesname = tree->species;
  vector<string> chrName;
  chrName.push_back("chr21");
  chrName.push_back("chr17");
  chrName.push_back("chr31");
  int start[] = {44836600, 31983200, 39789600} ;
  int end[] = {44846200, 31992000, 39798300}; 

  for (int s = 0; s < 3; s++) {
    orthoSeqRanges[s] = rsa->getSeq(speciesname[s], chrName[s], start[s], end[s]);
    if (orthoSeqRanges[s]){
      //cout << seqRange->seqname << "\t" << seqRange->length << "\t" << seqRange->offset << endl;
      //cout << seqRange->sequence << endl;
      /*
       * build list of additional exoncandidates, which are inserted in the graph
       */
      list<Status> additionalExons;
      for(list<OrthoExon>::iterator it = this->all_orthoex.begin(); it !=  this->all_orthoex.end(); it++){
	if(it->orthoex[s] != NULL){
	  it->orthoex[s]->begin -= orthoSeqRanges[s]->offset;
	  it->orthoex[s]->end -= orthoSeqRanges[s]->offset;
	  Status state = Status(CDS, it->orthoex[s]->begin, it->orthoex[s]->end, 0.0, it->orthoex[s]);
	  additionalExons.push_back(state);
	}
      }
      /*for(list<Status>::iterator it = additionalExons.begin(); it != additionalExons.end(); it++){
	cout << it->begin <<"\t"<< it->end <<"\t"<<((State*) it->item)->type << endl;
	}*/

      namgene.doViterbiPiecewise(sfc, orthoSeqRanges[s], minusstrand); // builds graph for each species

      list<Gene> *alltranscripts = namgene.getAllTranscripts();
      if(alltranscripts){
	cout << "building Graph for " << speciesname[s] << endl;
	if(!alltranscripts->empty()){
	  /*
	   * build datastructure for graph representation
	   * @stlist : list of all sampled states
	   */
	  list<Status> stlist;
	  buildDatastructure(alltranscripts, false, stlist);
	  //build graph
	  AugustusGraph *singleGraph = new AugustusGraph(&stlist, orthoSeqRanges[s]->length);
	  singleGraph->buildGraph(additionalExons);
	  // determine initial path and add score to total score
	  MEApath path(singleGraph);
	  this->score += path.findMEApath7();
	  //find correct position in vector and add graph for species to OrthoGraph
	  size_t pos = tree->getVectorPositionSpecies(speciesname[s]);
	  if (pos < this->graphs.size()){
	    this->graphs[pos] = singleGraph;
	  }
	  else{
	    cerr << "species names in Orthograph and OrthoExon don't match" << endl;
	  }
	  storePtrsToAlltranscripts(alltranscripts); //save pointers to transcripts and delete them after gene list is build
	}
      }    
    } else {
      cerr << "random sequence access failed on " << speciesname[s] << ", " << chrName[s] << ", " << start[s] << ", " <<  end[s] << ", " << endl;
    }
  }
  for(list<OrthoExon>::iterator it =this->all_orthoex.begin(); it != this->all_orthoex.end(); it++){
    tree->pruningAlgor(*it, *this);
  }
}



OrthoGraph::~OrthoGraph(){
  for(int i = 0; i < graphs.size(); i++){
    delete graphs[i];
    delete ptrs_to_alltranscripts[i];
    delete orthoSeqRanges[i];
  }
}

void OrthoGraph::optimizeLabels(){

  for( map<string, Score>::iterator it = cache::labelscore.begin(); it != cache::labelscore.end(); it++){
    cout << it->first << "\t" << it->second.treescore << "\t" << it->second.count << endl;
  }
  cout<<"Orthograph score: " << this->score << endl;
  cout<<"Phylo Score: " << cache::getOverallScore() << endl;

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

  for (size_t pos = 0; pos < this->graphs.size(); pos++){
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

    static vector<int> geneid(tree->species.size(), 1); // made this static so gene numbering goes across sequences and is unique
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
  filestreams.resize(tree->species.size());
  for(size_t pos = 0; pos < tree->species.size(); pos++){
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

  for(size_t pos = 0; pos < filestreams.size(); pos++){
    if(filestreams[pos]){
      if(filestreams[pos]->is_open()){
	filestreams[pos]->close();
	delete filestreams[pos];
      }
    }
  }
}
