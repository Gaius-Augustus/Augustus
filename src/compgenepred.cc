/**********************************************************************
 * file: compgenepred.cc licence: Artistic Licence, see file
 * LICENCE.TXT or
 * http://www.opensource.org/licenses/artistic-license.php descr.:
 * comparative gene prediction on multiple species authors: Mario
 * Stanke, Alexander Gebauer, Stefanie Koenig
 *
 * date    |   author       |  changes
 * --------|----------------|------------------------------------------
 * 07.03.12| Mario Stanke   | creation of the file
 * 06.10.12|Stefanie Koenig | construction of an OrthoGraph from a set of
 *         |                | orthologous sequences + integration of
 *         |                | the optimization method
 ***********************************************************************/

#include "compgenepred.hh"
#include "orthograph.hh"
#include "mea.hh"
#include "genomicMSA.hh"
#include "geneMSA.hh"
#include "orthoexon.hh"
#include "namgene.hh"

#include "contTimeMC.hh"
#include <gsl/gsl_matrix.h>
#include <ctime>
#include <sys/stat.h>

CompGenePred::CompGenePred() : tree(Constant::treefile) {

    vector<string> speciesNames;
    tree.getSpeciesNames(speciesNames);
    
    
    if (Constant::Constant::dbaccess.empty()) { // give priority to database in case both exist
        rsa = new MemSeqAccess(speciesNames);
    } else {
	string dbaccess = Constant::dbaccess;
	if (dbaccess.find(',') != string::npos){ // assuming mysql access
	    cout << "# assuming a MySQL database" << endl;
#ifdef AMYSQL
	    rsa = new MysqlAccess(speciesNames);
#else
	    throw ProjectError("Database access not possible with this compiled version. Please recompile with flag MYSQL.");
#endif

	}
	else if(dbaccess.find('.') != string::npos){ // assuming sqlite access
	    cout << "# assuming an SQLite database" << endl;
	    if(Constant::speciesfilenames.empty())
		throw ProjectError("Missing parameter speciesfilenames.");
#ifdef SQLITE
	    rsa = new SQLiteAccess(dbaccess.c_str(), speciesNames);
#else
	    throw ProjectError("Database access not possible with this compiled version. Please recompile with flag SQLITE.");
#endif
	    
	}
	else{
	    throw ProjectError("cannot determine the type of database.\n\
        for MySQL access, pass a comma separated list of connection details\n\
        --dbaccess=dbname,host,user,passwd\n\
        for SQLite access, pass the database file (file extension has to be .db)\n\
        --dbaccess=dbname.db\n");
	}
    }
}

void CompGenePred::start(){

    // read in alignment, determine orthologous sequence fragments
    
    ExonEvo evo;
    vector<double> branchset;
    tree.getBranchLengths(branchset);
    evo.setBranchLengths(branchset);
    evo.computeLogPmatrices();
    OrthoGraph::tree = &tree;
    GeneMSA::setTree(&tree);
    OrthoGraph::numSpecies = OrthoGraph::tree->numSpecies();
    Boolean noprediction = false;

#ifdef DEBUG
    cout << "-------------------------------\nparameters phylogenetic model\n-------------------------------" << endl;
    cout << "rate exon loss:\t" << evo.getMu() << endl;
    cout << "rate exon gain:\t" << evo.getLambda() << endl;
    cout << "phylo factor:\t" << evo.getPhyloFactor() <<  "\n-------------------------------" << endl;
#endif

   
    bool dualdecomp;
    try {
	dualdecomp = Properties::getBoolProperty("/CompPred/dualdecomp");
    } catch (...) {
	dualdecomp = true;
    }
    int maxIterations; // maximum number of dual decomposition iterations
    try {
	maxIterations = Properties::getIntProperty("/CompPred/maxIterations");
    } catch (...) {
	maxIterations = 100;
    }
    if(maxIterations <= 0){
	cerr << "Warning: /CompPred/maxIterations was set to "<<maxIterations<<". At least one iteration must be made." << endl;
	maxIterations =1;
    }
    double dd_factor; // parameter of the dual decomposition step size function 
    try {
	dd_factor = Properties::getdoubleProperty("/CompPred/dd_factor");
    } catch (...) {
	dd_factor = 15;
    }
    try {
	noprediction = Properties::getBoolProperty("noprediction");
    } catch (...) {}
    if(Constant::alternatives_from_evidence){
	cerr << "Warning: The option 'alternatives_from_evidence' is only available in the single species mode. Turned it off." << endl;
	Constant::alternatives_from_evidence=false;
    }
    bool genesWithoutUTRs;
    try {
	genesWithoutUTRs = Properties::getBoolProperty("/CompPred/genesWithoutUTRs");
    } catch (...) {
	genesWithoutUTRs = true;
    }
    bool onlyCompleteGenes = false;
    const char* genemodelValue = Properties::hasProperty("genemodel") ? Properties::getProperty("genemodel") : "partial";
    if(strcmp(genemodelValue, "complete") == 0){
	onlyCompleteGenes = true;
    }
    else if(strcmp(genemodelValue, "partial") != 0){
	throw ProjectError("in cgp mode only the options --genemodel=partial and --genemodel=complete are implemented.");
    }
    if(onlyCompleteGenes && Constant::utr_option_on)
	genesWithoutUTRs = false;
    bool featureScoreHects;
    try {
        featureScoreHects = Properties::getBoolProperty("/CompPred/featureScoreHects");
    } catch (...) {
        featureScoreHects = false;
    }
    bool conservationTrack;
    try {
        conservationTrack = Properties::getBoolProperty("/CompPred/conservation");
    } catch (...) {
        conservationTrack = false;
    }

    string outdir;  //direction for output files                                                                                                                                                                    
    try {
        outdir = Properties::getProperty("/CompPred/outdir");
    } catch (...) {
        outdir = "";
    }
    if(outdir != ""){
        if(outdir[outdir.size()-1] != '/')
            outdir += '/';                   // append slash if neccessary                                                                                                                                           
        outdir = expandHome(outdir);         // replace "~" by "$HOME"                                                                                                                                               
        // Does the directory actually exist?                                                                                                                                                                       
        struct stat buffer;
        if( stat(outdir.c_str(), &buffer) == -1 || !S_ISDIR(buffer.st_mode)){
            int err = mkdir(outdir.c_str(),0755);
            if(err != 0)
                throw ProjectError("could not create directory " + outdir);
        }
    }
    

    //initialize output files of initial gene prediction and optimized gene prediction
    vector<ofstream*> baseGenes = initOutputFiles(outdir,".mea"); // equivalent to MEA prediction
    vector<int> base_geneid(OrthoGraph::numSpecies, 1); // gene numbering
    vector<ofstream*> initGenes;
    vector<int> init_geneid(OrthoGraph::numSpecies, 1);
    if(featureScoreHects)
	initGenes = initOutputFiles(outdir,".init"); // score added to all orthologous exons and penalty added to all non orthologous exons, then global path search repeated
    vector<ofstream*> optGenes = initOutputFiles(outdir,".cgp");  //optimized gene prediction by applying majority rule move
    vector<int> opt_geneid(OrthoGraph::numSpecies, 1);
    vector<ofstream*> sampledExons = initOutputFiles(outdir,".sampled_ECs");

    BaseCount::init();
    PP::initConstants();
    NAMGene namgene; // creates and initializes the states
    StateModel::readAllParameters(); // read in the parameter files: species_{igenic,exon,intron,utr}_probs.pbl

    // temporary tests of codon rate matrix stuff (Mario)
    double *pi = ExonModel::getCodonUsage();
    CodonEvo codonevo;
    codonevo.setKappa(4.0);
    codonevo.setPi(pi);
    vector<double> b; // all branch lengths occuring in the tree
    b.push_back(.5); // 50% codon substitutions between D.mel and D.pseudoo., 40% between human and mouse
    codonevo.setBranchLengths(b, 10);
    // codonevo.printBranchLengths();
    
    codonevo.setOmegas(20);
    codonevo.setPrior(0.5);
    /*cout << "Omegas, for which substitution matrices are stored:" << endl;
      codonevo.printOmegas();*/
    codonevo.computeLogPmatrices();
    
    // gsl_matrix *P = codonevo.getSubMatrixLogP(0.3, 0.25);
    // printCodonMatrix(P);
    GeneMSA::setCodonEvo(&codonevo);
    vector<string> speciesNames;
    OrthoGraph::tree->getSpeciesNames(speciesNames);
    PhyloTree::setRSA(rsa);
    if (0) {
	vector<string> seqtuple(5,"");
	seqtuple[0] = "aaaaaaaaaaaa";
	seqtuple[1] = "aagaataataat";
	seqtuple[2] = "------------";
	seqtuple[3] = "------------";
	seqtuple[4] = "------------";
	int subst = 0;
	double omega = codonevo.estOmegaOnSeqTuple(seqtuple, &tree, subst);
	for (int i=0; i<seqtuple.size(); i++)
	    cout << seqtuple[i]<< endl;
	cout << "omega= " << omega << " subst= " << subst << endl;
	exit(1);
    }
    GenomicMSA msa(rsa);
    msa.readAlignment(Constant::alnfile);  // reads the alignment
    // rsa->printStats();
    // msa.compactify(); // Mario: commented out as this excludes paths through the alignment graph 
                         //(trivial mergers of neighboring alignments)
    msa.findGeneRanges(); // nontrivial summary of alignments
  
    //msa.printAlignment("");
    //exit(0);

    GeneMSA::openOutputFiles(outdir);
    while (GeneMSA *geneRange = msa.getNextGene()) {
	cout << "processing next gene range:" << endl;
	geneRange->printStats();
        OrthoGraph orthograph;
	vector<AnnoSequence> seqRanges(speciesNames.size());
        for (int s = 0; s < speciesNames.size(); s++) {
            string seqID = geneRange->getSeqID(s);
            if (!seqID.empty()) {
		int start = geneRange->getStart(s); // start, end refer to plus strand
		int end = geneRange->getEnd(s);
		AnnoSequence *as = rsa->getSeq(speciesNames[s], seqID, start, end, geneRange->getStrand(s));
                if (!as) {
                    cerr << "random sequence access failed on " << speciesNames[s] << ", " << seqID << ", " 
			 << start << ", " << end << ", " << endl;
                    break;
                } else {
		    seqRanges[s] = *as; // DNA seqs will be reused when omega is computed AND gene lists are processed for output  

		    list<Gene> *alltranscripts = NULL;

		    if (!noprediction){
			SequenceFeatureCollection* sfc = rsa->getFeatures(speciesNames[s],seqID,start,end,geneRange->getStrand(s));
			sfc->prepare(as, false, rsa->withEvidence(speciesNames[s]));
			namgene.doViterbiPiecewise(*sfc, as, bothstrands); // sampling
			alltranscripts = namgene.getAllTranscripts();
			orthograph.sfcs[s] = sfc;
		    } else {
			// turn whole sequence to lowercase characters
			for (unsigned pos = 0; pos < as->length; pos++)
			    as->sequence[pos] = tolower(as->sequence[pos]);
		    }
		    // this is needed for IntronModel::dssProb in GenomicMSA::createExonCands
                    namgene.getPrepareModels(as->sequence, as->length); 
		    
		    // identifies exon candidates in the sequence for species s
                    geneRange->createExonCands(s, as->sequence);
		    list<ExonCandidate*> *additionalExons = geneRange->getExonCands(s);
		    
		     if (alltranscripts){
			// cout << "building Graph for " << speciesNames[s] << endl;
			// build datastructure for graph representation
			// @stlist : list of all sampled states
			list<Status> stlist;
			if(!alltranscripts->empty()){
                            buildStatusList(alltranscripts, Constant::utr_option_on, stlist);
                        }
                        // build graph
                        orthograph.graphs[s] = new SpeciesGraph(&stlist, as, additionalExons, speciesNames[s], 
								geneRange->getStrand(s), genesWithoutUTRs, onlyCompleteGenes, sampledExons[s]);
                        orthograph.graphs[s]->buildGraph();
                        //orthograph.graphs[s]->printGraph(outdir + speciesNames[s] + "." + itoa(GeneMSA::geneRangeID) + ".dot");
			
			//save pointers to transcripts and delete them after gene list is build
                        orthograph.ptrs_to_alltranscripts[s] = alltranscripts;
		    }
                }
            }
        }
	geneRange->printGeneRanges();	    
	if (Constant::exoncands) // by default, ECs are not printed
	    geneRange->printExonCands();
	geneRange->createOrthoExons();
	try { // Kathrin Middendorf's playground
	    if (Properties::getBoolProperty("/CompPred/compSigScoring"))
		geneRange->comparativeSignalScoring(); 
	} catch (...) {}
	
	//geneRange->computeOmegas(seqRanges); // omega and number of substitutions is stored as OrthoExon attribute
	if(conservationTrack)
	    geneRange->printConsScore(seqRanges, outdir);

	if (!noprediction){
	    list<OrthoExon> hects = geneRange->getOrthoExons();
	    orthograph.linkToOEs(hects); // link ECs in HECTs to nodes in orthograph	    

	    orthograph.outputGenes(baseGenes,base_geneid);
	    
	    if(featureScoreHects){
		//add score for selective pressure of orthoexons
		orthograph.addScoreSelectivePressure();
		//determine initial path
		orthograph.globalPathSearch();
		orthograph.outputGenes(initGenes,init_geneid);
	    }	    
	    if(!orthograph.all_orthoex.empty()){
		if (dualdecomp){ // optimization via dual decomposition
		    vector< list<Gene> *> genelist(OrthoGraph::numSpecies);
		    orthograph.dualdecomp(evo,genelist,GeneMSA::geneRangeID-1,maxIterations, dd_factor);
		    orthograph.filterGeneList(genelist,optGenes,opt_geneid);
		} else { // optimization by making small changes (moves)
		    orthograph.pruningAlgor(evo);
		    orthograph.printCache();
		    orthograph.optimize(evo);
		    // transfer max weight paths to genes + filter + ouput
		    orthograph.outputGenes(optGenes, opt_geneid);
		}
	    }else{
		orthograph.outputGenes(optGenes, opt_geneid);
	    }
	    //geneRange->printOrthoExons(); //TODO: two copies of list<OrthoExon> (class geneRange and class OrthoGraph) -> replace one copy by a pointer to the other 
	    for(list<OrthoExon>::iterator it=orthograph.all_orthoex.begin(); it != orthograph.all_orthoex.end(); it++){
		geneRange->printSingleOrthoExon(*it,true);
	    }
	}
	seqRanges.clear(); // delete sequences
	delete geneRange;
    }

    GeneMSA::closeOutputFiles();
    if(featureScoreHects)
	closeOutputFiles(initGenes);
    closeOutputFiles(baseGenes);
    closeOutputFiles(optGenes);

}
