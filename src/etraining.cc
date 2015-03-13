/**********************************************************************
 * file:    etraining.cc
 * licence: Artistic Licence, see file LICENCE.TXT or
 *          http://www.opensource.org/licenses/artistic-license-2.0.php
 *          
 * descr.:  etraining executable that determines species-specific parameters
 * authors: Mario Stanke, mario@gobics.de
 *
 **********************************************************************/

// project includes
#include "types.hh"
#include "gene.hh"
#include "extrinsicinfo.hh"
#include "genbank.hh"
#include "statemodel.hh"
#include "merkmal.hh"


typedef ProjectError EHMMTrainingError;

/**
 * class EHMMTraining
 *
 * @authors  Emmanouil Stafilarakis, Mario Stanke
 */
class EHMMTraining{
public:
    EHMMTraining();
    ~EHMMTraining( );
    void HMMbuildParameters(int gcIdx, AnnoSequence* AnnoSeqList, 
			 FeatureCollection &extrinsicFeatures, BaseCount *bc = NULL);
    void CRFbuildParameters(AnnoSequence* AnnoSeqList, FeatureCollection &extrinsicFeatures);
    void deleteModelParameterFiles(const char* suffix = NULL);
    void printParameters(const char* suffix = NULL);
    void saveConfig(const char* file = 0);
    static AnnoSequence *split2SingleGeneSeqs(const AnnoSequence *multiannoseq);

private:
    GBProcessor* gbank;
    StateModel** states;
    int          statecount;
};


/*
 * Set the global variables according to command line options or configuration file
 */
void setParameters();


/*
 * main
 */
int main( int argc, char* argv[] ){
    string configfile;
    string filename;
    int verbosity = 1;

    LLDouble::setOutputPrecision(3);

    try{
	// determination of the configuration file
	if( argc <= 1 ){
	    cout << GREETING << endl << endl;
	    cout << HELP_USAGE_ETRAINING << endl << endl;
	    exit(0);
	}
	
	try {
	    Properties::init( argc, argv );
	} catch (ProjectError e) {
	    cout << e.getMessage() << endl << endl;
	    exit(1);
	}

	Constant::dss_gc_allowed = true; // default value for training!
	Constant::init();
	BaseCount::init();  // replaces train.init()
	GeneticCode::init();
	if (Properties::hasProperty("translation_table")) // check whether a different set of stop codons was specified
	    GeneticCode::chooseTranslationTable(Properties::getIntProperty("translation_table"));

	EHMMTraining train;
	setParameters();

	StateModel::init();
    
	AnnoSequence  *trainannoseq, //contains the annotated training sequences (multi genes on both strands possible)
	    *singletrainannoseq; // contains single-gene annotations on the forward strand only
	GBProcessor traingbank(Properties::getProperty("queryfile"));
	if (traingbank.fileType() != genbank) 
	    throw ProjectError("Input file not in genbank format.");
	try {
	  verbosity = Properties::getIntProperty("/augustus/verbosity");
	} catch (...) {}
	trainannoseq = traingbank.getAnnoSequenceList();
	singletrainannoseq = EHMMTraining::split2SingleGeneSeqs(trainannoseq);

	if (!singletrainannoseq){
	  throw ProjectError("Could not properly read Annotation.");
	}
	//AnnoSequence::deleteSequence(trainannoseq); // for some reason this causes a "double free or corruption" error

 	train.deleteModelParameterFiles();
	
	// TODO: read in the hints for training
	FeatureCollection extrinsicFeatures;
	

	BaseCount *bc = new BaseCount();
	ContentDecomposition cd;
	if (verbosity)
	    examineBaseCountOfGeneSeq(singletrainannoseq); // for parameter settings
	

	// 1. HMM-train all parameters in each GC class separately
	if (verbosity)
	  cout << "HMM-training the parameters..." << endl;
	if (Constant::useCRFtraining)
	  train.deleteModelParameterFiles(".HMM");
	for (int i=0; i< cd.n; i++) {
	    *bc = cd.getBaseCount(i);
	    if (verbosity)
	      cout << "i= " << i  << " bc= " << *bc << endl;
	    train.HMMbuildParameters(i, singletrainannoseq, extrinsicFeatures, bc);
	}
	
	if (!Constant::useCRFtraining){
	  if (verbosity)
	    cout << "Storing parameters to file..." << endl;
	  train.printParameters();
	} else {
	  if (verbosity)
	    cout << "Storing intermediate parameters for reference to files with .HMM suffix ..." << endl;
	  train.printParameters(".HMM");
	}
	
	// 2. CRF train the CRF parameter subset (both GC dependend and independent)
	if (Constant::useCRFtraining){
	  train.CRFbuildParameters(singletrainannoseq, extrinsicFeatures);
	  if (verbosity)
	    cout << "Overwriting parameters with the ones from CRF-training..." << endl;
	  train.deleteModelParameterFiles();
	  train.printParameters();
	}
	AnnoSequence::deleteSequence(singletrainannoseq);
    } catch( ProjectError& err ){
        cerr << "\n" <<  argv[0] << ": ERROR\n\t" << err.getMessage( ) << "\n\n";
        return 1;
    }
    return 0;
}

/*--- EHMMTraining methods -------------------------------------------*/

/*
 * EHMMTraining constructor/destructor
 */
EHMMTraining::EHMMTraining( ) : gbank(NULL) {
    statecount = Properties::getIntProperty( "/EHMMTraining/statecount" );
    states = new StateModel*[statecount];
    
    for( int i = 0; i < statecount; i++){
        char statestring[256];
        sprintf( statestring, "/EHMMTraining/state%02d", i );
        states[i] = StateModel::newStateModelPtr( Properties::getProperty( statestring ) );
	if (states[i]==NULL)
            throw EHMMTrainingError("void HMMTraining::buildTheta( ) : "
				    "couldn't initialize one of the models");
    }
}

EHMMTraining::~EHMMTraining( ){
    if( gbank )
        delete gbank;
    for( int i = 0; i < statecount; i++ )
        delete states[i];
    delete[] states;
}

/*
 * HMMbuildParameters
 * Use simple local maximum likelihood estimation of all the probabilistic models.
 */
void EHMMTraining::HMMbuildParameters(int gcIdx, AnnoSequence* AnnoSeqList, 
				   FeatureCollection &extrinsicFeatures, BaseCount *bc) {
    AnnoSequence *annoseq;
    AnnoSequence *first = (AnnoSequence*)0;

    first = AnnoSeqList;
    // weight the genes according to their base composition
    annoseq = first;
    while (annoseq) {
	if (Constant::decomp_num_steps > 1) {
	    annoseq->bc.normalize();
	    annoseq->setWeight(BaseCount::weight(*bc, annoseq->bc));
	    //cout << "weight(" << *bc << ", " << annoseq->bc << ") = " << BaseCount::weight(*bc, annoseq->bc) << endl;
	} else {
	    annoseq->setWeight(1);
	}
	annoseq = annoseq->next;
    }
    
    // use all (but weighted) genes for HMM-training
    for( int i = 0; i < statecount; i++ )
        states[i]->buildModel( first, gcIdx+1 );
    
    StateModel::storeGCPars(gcIdx);
}

/*
 * CRFbuildParameters
 * Starting from the HMM parameters use CRF training to iteratively (re)estimate some parameters.
 */
void EHMMTraining::CRFbuildParameters(AnnoSequence* AnnoSeqList, FeatureCollection &extrinsicFeatures) {
  Parameters* parameters = new Parameters(); // holds the parameters that are trained by the CRF

  // make vector from linked list of training genes
  vector<AnnoSequence*> AnnoSeqVec;
  for (AnnoSequence* annoseq = AnnoSeqList; annoseq; annoseq = annoseq->next)
    AnnoSeqVec.push_back(annoseq);
  cout << "Using a set of " << AnnoSeqVec.size() << " sequences for CRF-training" << endl;

  /*
   * check whether there is an evaluation set of genes and get read the annotated sequences
   */
  string evalFileName;
  try {
    evalFileName = expandHome(Properties::getProperty("evalset"));
  } catch  (...) {}
  if (evalFileName.length() > 0){
    GBProcessor gbank(evalFileName);
    vector<AnnoSequence*> evalAnnoSeqList;
    for (AnnoSequence* as = gbank.getAnnoSequenceList(); as; as = as->next)
      evalAnnoSeqList.push_back(as);
    cout << "using a set of " << evalAnnoSeqList.size() << " sequences for CRF-evaluation" << endl;
    CRF::setEvalAnnoSeqs(evalAnnoSeqList);
  }

  // register all parameters with the submodels
  for( int i = 0; i < statecount; i++ )
    states[i]->registerPars( parameters );
 
  vector<double> startWeights = parameters->getWeights();
  cout << "Parameters before CRF training" << endl;
  parameters->print(-1, false);
  CRF::setPrintPars(statecount, states);
  CRF::train(parameters, AnnoSeqVec, extrinsicFeatures);
  cout << "Parameters after CRF training" << endl;
  parameters->print(-1, false);
  vector<double> endWeights = parameters->getWeights();
  CRF::compareWeights(parameters, startWeights, endWeights);
}


void EHMMTraining::deleteModelParameterFiles(const char* suffix){
  BaseCount bc;
  for (int i = 0; i < statecount; i++)
    states[i]->printProbabilities( -1, &bc, suffix );
}

void EHMMTraining::printParameters(const char* suffix){
  ContentDecomposition cd;
  BaseCount bc;
  for (int idx=0; idx < cd.n; idx++) {
    bc = cd.getBaseCount(idx);
    for (int i = 0; i < statecount; i++) {
        states[i]->printProbabilities(idx, &bc, suffix);
    }
  }
}

void EHMMTraining::saveConfig(const char* file) {
    if (!file) {
        file = Properties::getProperty( "/HMMTraining/savefile" );
    }
    ofstream ofstrm( file );

    if (ofstrm) {
        ofstrm.close();
    }
}

/**
 * split2SingleGeneSeqs
 *
 * The input multiannoseq can be a list of annotated sequences where the sequences contain multiple
 * genes on either strand.
 * For training, I need single gene sequences on the forward strand.
 * this function splits each multi-gene sequence into several single-gene sequences and computes the 
 * reverse complement of the sequence and annotation if the gene is on the reverse strand.
 * Prerequisite: The genes are in order and non-overlapping.
 * 
 */
AnnoSequence* EHMMTraining::split2SingleGeneSeqs(const AnnoSequence *multiannoseq) {
  AnnoSequence *singleGeneSeqs=NULL, *singleGeneSeqsTail=NULL;
  AnnoSequence *newAS;
  Gene *g, *newG, *nextgene;
  int pieceIndex, pieceBegin, pieceEnd, windowBegin, windowEnd, geneEnd, geneBegin, nextBegin;
  int padding = 10000;
  bool multiGene = false;
  while (multiannoseq) {
    if (multiannoseq->anno) {
      /*
       * seq        : -----------------------------------------------------------------------------------
       * genes      :      11  1111   11                          22222            33333     3333
       * piecebounds: |                          |                         |                            |
       *                                     pieceBegin                 pieceEnd
       * singleseqs :   -------------------                  --------------    ----------------------
       *                                                     |             |
       *                                                 windowBegin windowEnd
       */

      pieceIndex = 0;
      pieceBegin = 0;
      g = dynamic_cast<Gene*> (multiannoseq->anno->genes);
      multiGene = (g && g->next);
      int skipped = 0;
      while (g) { // loop over all genes of a certain anno sequence
	pieceIndex++;
	geneEnd = (g->transend >0 )? g->transend : g->codingend;
	geneBegin = (g->transstart >0 )? g->transstart : g->codingstart;
	if (geneEnd<0)
	  throw ProjectError(string("Encountered gene without any exons: ")+ g->geneid);
	nextgene = g;
	pieceEnd =  -1;
	// search for the next gene that is completely to the right of g
	skipped--;
	while (nextgene && (pieceEnd <= geneEnd)){
	  nextgene = (Gene*) nextgene->next;
	  skipped++;
	  // --------**** g ****------------***** g->next ******---------
	  //                          X
	  // pieceEnd is in the middle of the two neighboring genes
	  if (nextgene) {
	    nextBegin = (nextgene->transstart >0 )? nextgene->transstart : nextgene->codingstart;
	    pieceEnd = (geneEnd + nextBegin)/2;
	  } else 
	    pieceEnd = multiannoseq->length-1; // just take the rest
	}

	windowBegin = pieceBegin;
	if (windowBegin < geneBegin - padding)
	  windowBegin = geneBegin - padding;
	windowEnd = pieceEnd;
	if (windowEnd > geneEnd + padding)
	  windowEnd = geneEnd + padding;

	// make a new AnnoSequence containing just g and the sequence from windowBegin to windowEnd
	newAS = new AnnoSequence();
	// copy the sequence interval
	newAS->sequence = newstrcpy(multiannoseq->sequence + windowBegin, windowEnd-windowBegin+1);
	newAS->length = windowEnd-windowBegin+1;
	char *newSeqName = new char[strlen(multiannoseq->sequence)+7];
	if (multiGene)
	  sprintf(newSeqName, "%s.p%d", multiannoseq->seqname, pieceIndex);
	else 
	  strcpy(newSeqName, multiannoseq->seqname);
	newAS->seqname = newSeqName;
	// copy the gene
	newG = new Gene(*g);
	newG->seqname =newSeqName;
	newG->next = NULL;
	// shift the coordinates to the new origin
	newG->shiftCoordinates(-windowBegin);
	newG->length = newAS->length;
	newAS->anno = new Annotation();
	newAS->anno->appendGene(newG);
	if (g->strand == minusstrand) { //reverse complement the sequence and the annotation
	  AnnoSequence *temp = newAS;
	  newAS = temp->getReverseComplement();
	  delete temp;
	}
	newAS->bc = BaseCount(newAS->sequence);
	//newAS->printGFF();
	// append newAS to the result AnnoSequence singleGeneSeqs
	if (singleGeneSeqsTail) {// list nonempty
	  singleGeneSeqsTail->next = newAS;
	  singleGeneSeqsTail = newAS;
	} else {
	  singleGeneSeqs = singleGeneSeqsTail = newAS;
	}

	// the next interval starts right after the end of this interval
	pieceBegin = pieceEnd + 1;
	       
	g = nextgene;
      }
      if (skipped > 0){
	cerr << "Skipped " << skipped << " gene(s) because their transcribed regions were overlapping in sequence "
	     << multiannoseq->seqname << endl;
      }

    }
    multiannoseq = multiannoseq->next;
  }
  return singleGeneSeqs;
}

/*--------------------------------------------------------------------*/


/*
 * setParameters
 */
void setParameters(){
    try {
	Gene::print_start = Properties::getBoolProperty("start");
    } catch (...) {}
    try {
	Gene::print_stop = Properties::getBoolProperty("stop");
    } catch (...) {}
    try {
	Gene::print_tss = Properties::getBoolProperty("tss");
    } catch (...) {}
    try {
	Gene::print_tts = Properties::getBoolProperty("tts");
    } catch (...) {}
    try {
	Gene::print_introns = Properties::getBoolProperty("introns");
    } catch (...) {}
    try {
	Gene::gff3 = Properties::getBoolProperty("gff3");
    } catch (...) {}
    try {
      Gene::stopCodonExcludedFromCDS = Properties::getBoolProperty("stopCodonExcludedFromCDS");
    } catch (...) {}
    try {
	Gene::print_cds = Properties::getBoolProperty("cds");
    } catch (...) {} 
    try {
	Gene::print_exonnames = Properties::getBoolProperty("exonnames");
    } catch (...) {} 
    try {
	Gene::print_utr = Properties::getBoolProperty("print_utr");
    } catch (...) {}
}
