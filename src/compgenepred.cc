/*
 * compgenepred.cc
 *
 * License: Artistic License, see file LICENSE.TXT or 
 *          https://opensource.org/licenses/artistic-license-1.0
 */
  
#include "compgenepred.hh"
#include "orthograph.hh"
#include "mea.hh"
#include "genomicMSA.hh"
#include "geneMSA.hh"
#include "orthoexon.hh"
#include "namgene.hh"
#include "contTimeMC.hh"
#include "liftover.hh"
#include "intronmodel.hh"
#include "train_logReg_param.hh"

#include <gsl/gsl_matrix.h>
#include <ctime>
#include <sys/stat.h>


/*
*   added by Giovanna Migliorelli 10.10.2019 
*   signatures for functions in charge of PRINTING MSA associated with some OE
*/
void ec2Chr(int ecStart, int ecEnd, int offset, int chrLen, Strand strand, int& chrStart, int& chrEnd, bool referToPlus);
void printMSA_printEC_helper(ostream& os, GeneMSA* geneRange, int sp, vector<string>& speciesNames, RandSeqAccess* rsa, int chrStart, int chrEnd, AnnoSequence* as, bool expandUnaligned);
void printMSA_printEC(ostream& os, vector<string>& speciesNames, RandSeqAccess* rsa, GeneMSA* geneRange, int sp, ExonCandidate* ec, bool downcase, bool expandUnaligned);
void printMSA(ostream& os, vector<string>& speciesNames, RandSeqAccess* rsa, GeneMSA* geneRange, OrthoExon& oe, bool extraInfo, bool expandUnaligned);


/*
*   added by Giovanna Migliorelli 11.10.2019 
*   signatures for functions in charge of RETRIEVING MSA associated with some OE as a vector of strings
*/
bool getMSA_getEC_helper(vector<string>& msa, GeneMSA* geneRange, int sp, vector<string>& speciesNames, RandSeqAccess* rsa, int chrStart, int chrEnd, AnnoSequence* as, bool includeUnaligned);
void getMSA_getEC(vector<string>& msa, vector<string>& speciesNames, RandSeqAccess* rsa, GeneMSA* geneRange, int sp, ExonCandidate* ec, bool downcase, bool includeUnaligned);
void getMSA(vector<string>& msa, vector<string>& speciesNames, RandSeqAccess* rsa, GeneMSA* geneRange, OrthoExon& oe, bool downcase, bool includeUnaligned);


CompGenePred::CompGenePred() : tree(Constant::treefile) {

    rsa = NULL;
    vector<string> speciesNames;
    tree.getSpeciesNames(speciesNames);
    cout<<"-------Speciesnames:--------"<<endl;
    for(int j=0; j<speciesNames.size(); j++){
      cout<<"species "<<j<<"\t"<<speciesNames[j]<<endl;
    }

    if (Properties::hasProperty("trainFeatureFile")) {
      trainFeature_from_file();
      train_OEscore_params(speciesNames.size());
    } else {

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
}

void CompGenePred::start(){


  // check if training mode is turned on and logreg features shall be read from file
  
  if(!Properties::hasProperty("trainFeatureFile")){

    // read in alignment, determine orthologous sequence fragments

    int num_phylo_states; // number of states in ExonEvo model
    try {
	num_phylo_states = Properties::getIntProperty("/CompPred/phylo_model");
    } catch (...) {
	num_phylo_states = 2;
    }
    if( num_phylo_states > 4 || num_phylo_states < 2)
	throw ProjectError("/CompPred/phylo_model must be 2,3 or 4");
    ExonEvo evo(num_phylo_states);
    vector<double> branchset;
    tree.getBranchLengths(branchset);
    evo.setBranchLengths(branchset);
    //evo.printBranchLengths();
    evo.computeLogPmatrices();
    OrthoGraph::tree = &tree;
    GeneMSA::setTree(&tree);
    OrthoGraph::numSpecies = OrthoGraph::tree->numSpecies();
    vector<string> speciesNames;
    OrthoGraph::tree->getSpeciesNames(speciesNames);    
    Boolean noprediction = false;
    Boolean useLocusTrees = false; // reestimate tree for each gene range


#ifdef DEBUG
    cout << "-------------------------------\nparameters phylogenetic model\n-------------------------------" << endl;
    cout << "rate exon loss:\t" << evo.getMu() << endl;
    cout << "rate exon gain:\t" << evo.getLambda() << endl;
#endif

    Constant::temperature = 3;
    Properties::assignProperty("temperature", Constant::temperature);
    if (Constant::temperature < 0){
	Constant::temperature = 0;
    }
    if (Constant::temperature > 7){
	Constant::temperature = 7;
    }
    try {
	PhyloTree::phylo_factor  = Properties::getdoubleProperty("/CompPred/phylo_factor");
    } catch (...) {
	PhyloTree::phylo_factor = 1;
    }
    if(PhyloTree::phylo_factor <= 0.0){
	throw ProjectError("/CompPred/phylo_factor must to be real positive number.");
    }
    
    // used for mapping an existing annotation to the other genomes
    // each exon/intron that is supported by a hint (e.g. exon/intron from the annotation) gets this
    // additional score to make sure that it is transferred to the other genomes and not the
    // other way round. It has to be at least has high as the maximum cost of an exon gain or loss event
    SpeciesGraph::maxCostOfExonLoss = - log((evo.getMu()+evo.getLambda())*evo.minBranchLength())*PhyloTree::phylo_factor*100;

    double ctree_scaling_factor = 1; // scaling factor to scale branch lengths in codon tree to one codon substitution per time unit
    try {
	ctree_scaling_factor = Properties::getdoubleProperty("/CompPred/scale_codontree");
    } catch (...) {
	ctree_scaling_factor = 1;
    }
    if(ctree_scaling_factor <= 0.0){
	cerr << "No negative scaling factor allowed. /CompPred/scale_codontree must be a positive real number. Will use =1." << endl;	
	ctree_scaling_factor=1;
    }
    int maxIterations; // maximum number of dual decomposition iterations in each round
    try {
	maxIterations = Properties::getIntProperty("/CompPred/maxIterations");
    } catch (...) {
	maxIterations = 500;
    }
    if(maxIterations <= 0){
	cerr << "Warning: /CompPred/maxIterations must be a pos. Will use =500." << endl;
	maxIterations = 500;
    }
    int rounds; // number of dual decomposition rounds
    try {
	rounds = Properties::getIntProperty("/CompPred/dd_rounds");
    } catch (...) {
	rounds = 5;
    }
    if(rounds <= 0){
	cerr << "Warning: /CompPred/rounds was set to "<<rounds<<". At least one round must be made." << endl;
	rounds =1;
    }
    const char* dd_step_rule = Properties::hasProperty("/CompPred/dd_step_rule") ? Properties::getProperty("/CompPred/dd_step_rule") : "mixed";
    OrthoGraph::setStepRule(dd_step_rule);
    
    string dd_param_s; 
    // parameter that defines the step size in dual decomposition.
    // If a range is given, all values in that range are tried until convergence is achieved
    vector<double> dd_factors; 
    try {
        dd_param_s = Properties::getProperty("/CompPred/dd_factor");
    } catch (...) {
        dd_param_s = "1-4"; // default: 1st round "polyak", 2nd-5th round "square_root" with c=1,2,3,4
    }
    if(OrthoGraph::step_rule == polyak) // only one round possible for this step size rule
	dd_param_s = "1";
    try{
	size_t i = dd_param_s.find('-');
	if(i != std::string::npos){
	    double start;
	    if( !(stringstream(dd_param_s.substr(0,i)) >> start))
		throw ProjectError("Cannot read interval start.");
	    double end;
	    if( !(stringstream(dd_param_s.substr(i+1, string::npos)) >> end))
		throw ProjectError("Cannot read interval end.");
	    if(start > end)
		throw ProjectError("Interval start greater than interval end.");
	    if(OrthoGraph::step_rule == mixed){ // 1st round: polyak, 2nd round and onwards: square_root
		dd_factors.push_back(1);        // dummy, not required in first round
		for (int i=0; i < rounds-1; i++)
                    dd_factors.push_back(start+i*(end-start)/(rounds-2));
	    }
	    else{
		for (int i=0; i < rounds; i++)
		    dd_factors.push_back(start+i*(end-start)/(rounds-1));	
	    }
	}
	else{
	    double pos;	    
	    if( !(stringstream(dd_param_s) >> pos))
		throw ProjectError("Is not numeric.");
	    dd_factors.push_back(pos);
	}
    } catch (ProjectError &e) {
	throw ProjectError("Format error parsing parameter --/CompPred/dd_factor=" + dd_param_s +".\n" + e.getMessage());
    }
    bool onlySampling = false;
    try {
	noprediction = Properties::getBoolProperty("noprediction");
    } catch (...) {}
    if(Properties::hasProperty("referenceFile")){
      onlySampling = true;
      cout << "# AUGUSTUS is running in training mode. No prediction will be done!" << endl;
      try {
	Constant::refSpecies = Properties::getProperty("refSpecies");
	if(Properties::hasProperty("param_outfile")){
	  cout << "# Using file " << Properties::getProperty("param_outfile") << " to store logReg parameters." << endl;
	}else{
	  cout << "# No outfile for logReg parameters specified. Writing parameters to " << Constant::configPath <<  "/cgp/log_reg_parameters_trained.cfg" << endl;
	}
      } catch (ProjectError &e) {
	throw ProjectError("For parameter training a reference species must be specified. Use --refSpecies=<SPECIES> and note, that <SPECIES> must be identical to one of the species names provided in the alignment and tree files.");
      }
    }
    try {
	useLocusTrees = Properties::getBoolProperty("locustree");
    } catch (...) {}
    

    if(Constant::alternatives_from_evidence){
	cerr << "Warning: The option 'alternatives-from-evidence' is only available in single species mode. Turned it off." << endl
	     << "Rerun with --alternatives-from-evidence=0 to remove this warning." << endl;
	Constant::alternatives_from_evidence=false;
    }
    // optional config file should contain feature scores from a logistic regression (otherwise defaults are used)
    if (Properties::hasProperty(EXTERNAL_KEY)) {
	string optCfgFile = expandHome(Properties::getProperty(EXTERNAL_KEY));
	cout << "# Optional config file " << optCfgFile << " is used." << endl;
    }
    bool genesWithoutUTRs;
    try {
	genesWithoutUTRs = Properties::getBoolProperty("/CompPred/genesWithoutUTRs");
    } catch (...) {
	genesWithoutUTRs = true;
    }
    bool overlapComp;
    try {
      overlapComp = Properties::getBoolProperty("/CompPred/overlapcomp");
    } catch (...) {
      overlapComp = false;
    }
    bool onlyCompleteGenes = false;
    const char* genemodelValue = Properties::hasProperty("genemodel") ? Properties::getProperty("genemodel") : "partial";
    if(strcmp(genemodelValue, "complete") == 0){
	onlyCompleteGenes = true;
    }
    else if(strcmp(genemodelValue, "partial") != 0 && strcmp(genemodelValue, "bacterium") != 0){
	throw ProjectError("in cgp mode only the options --genemodel=partial, --genemodel=bacterium and --genemodel=complete are implemented.");
    }
    if(onlyCompleteGenes && Constant::utr_option_on)
	genesWithoutUTRs = false;
    bool conservation;
    try {
        conservation = Properties::getBoolProperty("/CompPred/conservation");
    } catch (...) {
      if(Constant::logreg)
	conservation = true;
      else
        conservation = false;
    }
    double thold;
    try {
	thold = Properties::getdoubleProperty("/CompPred/ec_thold");
    } catch (...) {
	thold = 0.0;
    }
    SpeciesGraph::setECThold(thold);
    try {
	thold = Properties::getdoubleProperty("/CompPred/ic_thold");
    } catch (...) {
	thold = 0.0;
    }
    SpeciesGraph::setICThold(thold);
   
    string outdir;  //direction for output files 
    try {
        outdir = Properties::getProperty("/CompPred/outdir");
    } catch (...) {
        outdir = "";
    }
    double mil_factor; // mean intron length factor (>=1), the higher the less are long introns penalized 
    try {
        mil_factor = Properties::getdoubleProperty("/CompPred/mil_factor");
    } catch (...) {
        mil_factor = 1.0; // a value of 100 roughly corresponds to not penalizing long introns
    }
    double meanIntrLen = -1.0; // initialized later
    /*
     * by default only likely exon candidates (the ones from sampling) are lifted over to
     *  the other genomes. If this flag is turned on ALL exon candidates are lifted over
     */
    bool liftover_all_ECs;
    try {
        liftover_all_ECs = Properties::getBoolProperty("/CompPred/liftover_all_ECs");
    } catch (...) {
	liftover_all_ECs = false;
    }

    bool use_omega;
    try{ 
      use_omega = Properties::getBoolProperty("/CompPred/omega");
    } catch (...){
      if(Constant::logreg && Constant::ex_sc[6] != 0)
	use_omega = true;
      else
	use_omega = false;
    }
    
    //initialize output files of initial gene prediction and optimized gene prediction
    vector<ofstream*> baseGenes, optGenes, sampledGFs;
    if (Constant::printMEA)
	baseGenes = initOutputFiles(outdir,".mea"); // equivalent to MEA prediction
    vector<int> base_geneid(OrthoGraph::numSpecies, 1); // gene numbering
    optGenes = initOutputFiles(outdir,".cgp");  //optimized gene prediction by applying majority rule move
    vector<int> opt_geneid(OrthoGraph::numSpecies, 1);
    if (Constant::printSampled)
	sampledGFs = initOutputFiles(outdir,".sampled_GFs"); // prints sampled exons/introns and their posterior probs to file

    bool printCodons;
    try {
      printCodons =  Properties::getBoolProperty("/CompPred/printOrthoExonAli");
    } catch(...){
      printCodons = 0;
    }
    ofstream codonAli;   // prints codon alignments of all orthoexons in maf format
    if(printCodons){
      codonAli.open(outdir + "orthoexons_codonAlignment.maf");
    }

    BaseCount::init();
    PP::initConstants();
    NAMGene namgene; // creates and initializes the states
    StateModel::readAllParameters(); // read in the parameter files: species_{igenic,exon,intron,utr}_probs.pbl

    int k; // number of omega values for which rate matrices are stored
    try {
      k = Properties::getIntProperty("/CompPred/num_omega");
    } catch(...){
      k = 20;
    }

    // initializing codon rate matricies, for exon evolution see code above (evo)
    PhyloTree ctree(tree); // codon tree
    ctree.scaleTree(ctree_scaling_factor); // scale branch lengths to codon substitutions 
    vector<double> ct_branchset;
    ctree.getBranchLengths(ct_branchset);
    double *pi = ExonModel::getCodonUsage();
    CodonEvo codonevo;
    codonevo.setKappa(4.0);
    codonevo.setPi(pi);
    codonevo.setBranchLengths(ct_branchset, 25);
    //codonevo.printBranchLengths();
    codonevo.setOmegas(k);
    // TODO: different prior for coding and noncoding
    codonevo.setPrior(0.5);
    if(Constant::useAArates){
      codonevo.setAAPostProbs();
    }
    //cout << "Omegas, for which substitution matrices are stored:" << endl;
    //codonevo.printOmegas();
    codonevo.computeLogPmatrices();
    
    // gsl_matrix *P = codonevo.getSubMatrixLogP(0.3, 0.25);
    // printCodonMatrix(P);
    GeneMSA::setCodonEvo(&codonevo);
    GenomicMSA msa(rsa);
    msa.readAlignment(Constant::alnfile);  // reads the alignment
    // msa.printAlignment("");    
    // rsa->printStats();
    // msa.compactify(); // Mario: commented out as this excludes paths through the alignment graph 
                         //(trivial mergers of neighboring alignments)
    msa.findGeneRanges(); // nontrivial summary of alignments

    GeneMSA::openOutputFiles(outdir);
    while (GeneMSA *geneRange = msa.getNextGene()) {
	cout << "processing next gene range:" << endl;
	geneRange->printStats();
	
	if (useLocusTrees){ // Charlotte Janas playground, off by default
	    geneRange->constructTree();
	}

        OrthoGraph orthograph;
	vector<AnnoSequence*> seqRanges(speciesNames.size());
	vector<map<int_fast64_t,ExonCandidate*> > exoncands(speciesNames.size()); // EC hash: collection of all ECs of all species

	// retrieval of sequences and sampling of gene structures
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
		    seqRanges[s] = as; // DNA seqs will be reused when omega is computed AND gene lists are processed for output  

		    list<Transcript*> *transcripts = NULL;

		    if (!noprediction){
		      SequenceFeatureCollection* sfc = rsa->getFeatures(speciesNames[s],seqID,start,end,geneRange->getStrand(s));
		      sfc->prepare(as, false, rsa->withEvidence(speciesNames[s]));
		      namgene.doViterbiPiecewise(*sfc, as, bothstrands); // sampling
		      transcripts = namgene.getAllTranscripts();
		      orthograph.sfcs[s] = sfc;
		      orthograph.ptrs_to_alltranscripts[s] = transcripts;
		    } else {
		      // turn whole sequence to lowercase characters
		      for (unsigned pos = 0; pos < as->length; pos++)
			as->sequence[pos] = tolower(as->sequence[pos]);
		    }
		    // insert sampled exons into the EC hash
		    if (transcripts){		    
			for (list<Transcript*>::iterator geneit = transcripts->begin(); geneit != transcripts->end(); geneit++) {
			    if ((*geneit)->isCoding()){ // noncoding comparative prediction not (yet) implemented
				Gene *g = dynamic_cast<Gene*> (*geneit);
				if (g && orthograph.sfcs[s])
				    g->compileExtrinsicEvidence(orthograph.sfcs[s]->groupList);
				State *st = (*geneit)->exons;
				while (st) {
				    // include framemod into type
				    st->includeFrameModIntoType();
				    ExonCandidate *ec = new ExonCandidate(toExonType(stateTypeIdentifiers[st->type]),st->begin,st->end);
				    int_fast64_t key = ec->getKey();
				    map<int_fast64_t, ExonCandidate*>::iterator ecit;
				    ecit = exoncands[s].find(key);
				    if (ecit == exoncands[s].end()){ // insert new EC
					exoncands[s].insert(pair<int_fast64_t, ExonCandidate*>(key,ec));
				    } else {
					delete ec;
				    }
				    st = st->next;
				}
			    }
			}
		    }
		}
	    }
	}
	// liftover of sampled exons to other species
	vector<int> offsets = geneRange->getOffsets();
	LiftOver lo(geneRange->getAlignment(), offsets);
	map<int_fast64_t, list<pair<int,ExonCandidate*> > > alignedECs; // hash of aligned ECs
	// liftover of sampled ECs from genome to alignment space
	lo.projectToAli(exoncands,alignedECs);
	
	/*
	 * liftover of sampled ECs from alignment to genome space
	 * - creates new ECs e_j that are sampled in a subset of species i != j
	 * - marks e_j as "absent", if both boundaries of some EC e_i, i != j are aligned to j,
	 *   but the exon signals (e.g. splice sites, open reading frame) are missing
	 */
	lo.projectToGenome(alignedECs, seqRanges, exoncands, true);

	// create additional ECs for each species and  insert them into exoncands and alignedECs
	vector<map<int_fast64_t,ExonCandidate*> > addECs(speciesNames.size()); // new ECs that need to be mapped to the alignment
	for (int s = 0; s < speciesNames.size(); s++) {
	    if (seqRanges[s]) {
	        AnnoSequence *as = seqRanges[s];
	        try {
	            // this is needed for IntronModel::dssProb in GenomicMSA::createExoncands
	            namgene.getPrepareModels(as->sequence, as->length);
	            // identifies exon candidates in the sequence for species s
	            geneRange->createExonCands(s, as->sequence, exoncands[s], addECs[s]);
	        } catch (ProjectError &e) {
		    cerr << "CGP error when creating additional ECs for " << speciesNames[s] << endl
			 << e.getMessage();
		    throw e;
	        }
	    }
	}

	// liftover of additional ECs from genome to alignment space
	lo.projectToAli(addECs,alignedECs);
	addECs.clear(); // not needed anymore
	
	// liftover of additional ECs from alignment to genomes space 
	lo.projectToGenome(alignedECs, seqRanges, exoncands, liftover_all_ECs);
	
	geneRange->setExonCands(exoncands);
	exoncands.clear(); // not needed anymore, exoncands are now stored as a vector of lists of ECs in geneRange

	// create HECTS
	list<OrthoExon> hects;  // list of ortholog exons found in a gene Range
	geneRange->createOrthoExons(hects, alignedECs, &evo);

	if(meanIntrLen<0.0)
	    meanIntrLen = mil_factor * IntronModel::getMeanIntrLen(); // initialize mean intron length

	// build graph from sampled gene structures and additional ECs
        for (int s = 0; s < speciesNames.size(); s++) {	    
	    if (orthograph.ptrs_to_alltranscripts[s]){
		list<Transcript*> *alltranscripts = orthograph.ptrs_to_alltranscripts[s];
		cout << "building Graph for " << speciesNames[s] << endl;
		// build datastructure for graph representation
		// @stlist : list of all sampled states
		list<Status> stlist;
		if (!alltranscripts->empty()){
		    buildStatusList(*alltranscripts, Constant::utr_option_on, stlist);
		}
		// build graph
		ofstream *gf = NULL;
		if (Constant::printSampled)
		    gf = sampledGFs[s];
		orthograph.graphs[s] = new SpeciesGraph(&stlist, seqRanges[s], geneRange->getExonCands(s), speciesNames[s], 
							geneRange->getStrand(s), genesWithoutUTRs, onlyCompleteGenes, gf, overlapComp);
		orthograph.graphs[s]->buildGraph(meanIntrLen);
		//orthograph.graphs[s]->printGraph(outdir + speciesNames[s] + "." + itoa(GeneMSA::geneRangeID) + ".dot");
		
	    }
	}
 
	geneRange->printGeneRanges();
	if (Constant::exoncands) // by default, ECs are not printed
	    geneRange->printExonCands();
	try { // Kathrin Middendorf's playground
	    if (Properties::getBoolProperty("/CompPred/compSigScoring"))
		geneRange->comparativeSignalScoring(hects); 
	} catch (...) {}
       
	if(use_omega){
	    geneRange->computeOmegasEff(hects, seqRanges, &ctree, &codonAli); // omega and number of substitutions is stored as OrthoExon attribute
	    // calculates an omega for every single codon alignment and prints wiggle trac for ever reading frame and species combination that exists in an ortho exon
	    //geneRange->printOmegaForCodon(outdir);
	    //inefficient omega calculation, only use for debugging purpose 
	    //geneRange->computeOmegas(hects, seqRanges, &ctree);
	}
	if (conservation)
	    geneRange->calcConsScore(hects, seqRanges, outdir);

	if(!noprediction && !onlySampling){
	    orthograph.linkToOEs(hects); // link ECs in HECTs to nodes in orthograph	    
	    orthograph.globalPathSearch();
	    if (Constant::printMEA)
		orthograph.outputGenes(baseGenes,base_geneid);
	    	    
	    if(!hects.empty()){
	    	// optimization via dual decomposition
		vector< list<Transcript*> *> genelist(OrthoGraph::numSpecies);
		orthograph.dualdecomp(hects,evo,genelist,GeneMSA::geneRangeID-1,maxIterations, dd_factors);
		orthograph.filterGeneList(genelist,opt_geneid);
		orthograph.createOrthoGenes(geneRange);
		orthograph.printOrthoGenes();
		orthograph.printGenelist(optGenes);

	    }else{
		orthograph.outputGenes(optGenes, opt_geneid);
	    }
	}
	if(Constant::printOEs)
	    geneRange->printOrthoExons(hects);
	    
	// store hect features globally for training
	if(Properties::hasProperty("referenceFile")){
	  cout << "collect sample features" << endl;
	  int speciesID = find(speciesNames.begin(), speciesNames.end(), Constant::refSpecies) - speciesNames.begin();
	  if(speciesID >= speciesNames.size()){
	    throw ProjectError("Species " + Constant::refSpecies + " not found. Use one of the names specified in the alignment file as a reference!");
	  }else{
	    geneRange->collect_features(speciesID, &hects, orthograph.graphs[speciesID]);
	  }
	}

	// delete sequences
	for (int i=0; i<seqRanges.size(); i++) {
		delete seqRanges[i];
	}
	// delete geneRange
	delete geneRange;
    }

    GeneMSA::closeOutputFiles();
    if (Constant::printMEA)
	closeOutputFiles(baseGenes);
    closeOutputFiles(optGenes);
    if (Constant::printSampled)
	closeOutputFiles(sampledGFs);

    // delete all trees                                           
    for(unordered_map< bit_vector, PhyloTree*, boost::hash<bit_vector>>::iterator topit = GeneMSA::topologies.begin(); topit !=GeneMSA::topologies.end(); topit++){
	delete topit->second;
    }                                                                                                                                                              
    GeneMSA::topologies.clear(); 
  
    if(Properties::hasProperty("referenceFile")){
      // initialize training of log reg parameters
      train_OEscore_params(speciesNames.size());
    }
  } 
}



/*
*   added by Giovanna Migliorelli 10.10.2019 
*   functions to print the MSA associated with some ortho exon
*/

// conversion from ec to chr coordinate sys (output coordinates are invariably referred to the positive strand if referToPlus = true)
void ec2Chr(int ecStart, int ecEnd, int offset, int chrLen, Strand strand, int& chrStart, int& chrEnd, bool referToPlus){
	if(!referToPlus || strand == Strand::plusstrand){
		chrStart = offset + ecStart;
		chrEnd = offset + ecEnd;
	}
	else{
		chrStart = chrLen - 1 - offset - ecEnd;
		chrEnd = chrLen - 1 - offset - ecStart;
	}
}

// helper function to print MSA
void printMSA_printEC_helper(ostream& os, GeneMSA* geneRange, int sp, vector<string>& speciesNames, RandSeqAccess* rsa, int chrStart, int chrEnd, AnnoSequence* as, bool expandUnaligned){

	os << speciesNames[sp] << "\t\t";
	
	vector<int> pattern;
	
	int res = geneRange->getAlignment()->rows[sp]->getSeqInfo(chrStart, chrEnd, pattern);
	
	if(res>-1){
		if(pattern.empty()){
			os << " erroneous condition due to empty pattern: " << chrStart << " " << chrEnd << " " << res;
			os << endl;
		}
		else{

			for(int pos = 0, i, k=0;k<pattern.size();k+=3){
				if(expandUnaligned){
					for(i=0;i<pattern[k];++i){
						os << "X";
						++pos;
					}
				}
				else{	
					if(pattern[k]>0){			
						os << "X" << pattern[k];
						pos += pattern[k];
					}
				}
						
				for(i=0;i<pattern[k+2];++i){
					os << "-";
				}
				
				for(i=0;i<pattern[k+1];++i){
					os << as->sequence[pos];
					++pos;
				}
			}
		}
	}
	else 
		os << " erroneous condition";

	os << endl;
}

// helper function to print MSA
void printMSA_printEC(ostream& os, vector<string>& speciesNames, RandSeqAccess* rsa, GeneMSA* geneRange, int sp, ExonCandidate* ec, bool downcase, bool expandUnaligned){

	int chrStart, chrEnd, chrStartStranded, chrEndStranded;
	string seqID = geneRange->getSeqID(sp);
	Strand strand;
	AnnoSequence* as = NULL;
	

	// coordinate conversion from exon cand space to chromosome space (we want coo referred to the + strand because getSeq from RandSeqAccess requires so)
	ec2Chr(ec->begin, ec->end, geneRange->getOffsets()[sp], rsa->getChrLen(sp, seqID), geneRange->getStrand(sp), chrStart, chrEnd, true);
	
	if(geneRange->getStrand(sp) == plusstrand){
		if(isPlusExon(ec->type))
			strand = plusstrand;
		else 
			strand = minusstrand;
	}
	else{
		if(isPlusExon(ec->type))
			strand = minusstrand;
		else 
			strand = plusstrand;
	}

	as = rsa->getSeq(speciesNames[sp], seqID, chrStart, chrEnd, strand);
	
	if(as){
		if(downcase){
			for(int k=0;k<as->length;++k){
				switch(as->sequence[k]){
					case 'A':
						as->sequence[k] = 'a';
					break;
					case 'C':
						as->sequence[k] = 'c';
					break;
					case 'G':
						as->sequence[k] = 'g';
					break;
					case 'T':
						as->sequence[k] = 't';
					break;
				}
			}
		}

		// coordinate conversion from exon cand space to chromosome space (in this case coordinates are referred to the strand the ec lies on
		// because all methods dealing with the alignment usually require so)		
		ec2Chr(ec->begin, ec->end, geneRange->getOffsets()[sp], rsa->getChrLen(sp, seqID), geneRange->getStrand(sp), chrStartStranded, chrEndStranded, false);
		printMSA_printEC_helper(os, geneRange, sp, speciesNames, rsa, chrStartStranded, chrEndStranded, as, expandUnaligned);
		delete as;
	}
}

// main function to print MSA
void printMSA(ostream& os, vector<string>& speciesNames, RandSeqAccess* rsa, GeneMSA* geneRange, OrthoExon& oe, bool extraInfo, bool expandUnaligned){
	if(rsa==NULL || geneRange==NULL)
		return;

	int lenMod3;
	ExonType etype;
	
	for(int sp=0;sp<oe.orthoex.size();++sp){
		if(oe.orthoex[sp] != NULL){
			etype = oe.orthoex[sp]->getExonType();
			lenMod3 = (oe.orthoex[sp]->getEnd() - oe.orthoex[sp]->getStart() + 1) % 3;
			break;
		}
	}

	if(extraInfo)
		os << "Printing OE (type: " << stateExonTypeIdentifiers[etype] << " lenMod3: " << lenMod3 << ")" << endl;
	
	for(int sp=0;sp<oe.orthoex.size();++sp){
		if(oe.orthoex[sp] != NULL){
			printMSA_printEC(os, speciesNames, rsa, geneRange, sp, oe.orthoex[sp], true, true);
		}
	}	
	os << endl;
}


/*
*   added by Giovanna Migliorelli 11.10.2019 
*   alternative version of printMSA which allows to get the MSA as a vector of strings instead of printing it
*/

// helper function to get MSA : from the ec sequence and frags list the alignment is derived (it includes unaligned parts optionally)
bool getMSA_getEC_helper(vector<string>& msa, GeneMSA* geneRange, int sp, vector<string>& speciesNames, RandSeqAccess* rsa, int chrStart, int chrEnd, AnnoSequence* as, bool includeUnaligned){

	vector<int> pattern;
	
	int res = geneRange->getAlignment()->rows[sp]->getSeqInfo(chrStart, chrEnd, pattern);
	
	if(res<0 || pattern.empty())
		return false;


	int alilen = 0;
	for(int k=0;k<pattern.size();++k){
		if(includeUnaligned || k%3 != 0)
			alilen += pattern[k];
	}

	msa[sp].resize(alilen);
	
	for(int alipos = 0, pos = 0, i, k=0;k<pattern.size();k+=3){

		if(includeUnaligned){
			for(i=0;i<pattern[k];++i){
				msa[sp][alipos] = 'X';
				++pos;
				++alipos;
			}
		}
						
		for(i=0;i<pattern[k+2];++i){
			msa[sp][alipos] = '-';
			++alipos;			
		}
		
		for(i=0;i<pattern[k+1];++i){
			msa[sp][alipos] = as->sequence[pos];
			++pos;
			++alipos;
		}
	}

	return true;	
}

// helper function to get MSA : retrieves ec sequence
void getMSA_getEC(vector<string>& msa, vector<string>& speciesNames, RandSeqAccess* rsa, GeneMSA* geneRange, int sp, ExonCandidate* ec, bool downcase, bool includeUnaligned){

	int chrStart, chrEnd, chrStartStranded, chrEndStranded;
	string seqID = geneRange->getSeqID(sp);
	Strand strand;
	AnnoSequence* as = NULL;
	

	// coordinate conversion from exon cand space to chromosome space (we want coo referred to the + strand because getSeq from RandSeqAccess requires so)
	ec2Chr(ec->begin, ec->end, geneRange->getOffsets()[sp], rsa->getChrLen(sp, seqID), geneRange->getStrand(sp), chrStart, chrEnd, true);
	
	if(geneRange->getStrand(sp) == plusstrand){
		if(isPlusExon(ec->type))
			strand = plusstrand;
		else 
			strand = minusstrand;
	}
	else{
		if(isPlusExon(ec->type))
			strand = minusstrand;
		else 
			strand = plusstrand;
	}

	as = rsa->getSeq(speciesNames[sp], seqID, chrStart, chrEnd, strand);
	
	if(as){
		if(downcase){
			for(int k=0;k<as->length;++k){
				switch(as->sequence[k]){
					case 'A':
						as->sequence[k] = 'a';
					break;
					case 'C':
						as->sequence[k] = 'c';
					break;
					case 'G':
						as->sequence[k] = 'g';
					break;
					case 'T':
						as->sequence[k] = 't';
					break;
				}
			}
		}

		// coordinate conversion from exon cand space to chromosome space (in this case coordinates are referred to the strand the ec lies on
		// because all methods dealing with the alignment usually require so)
		
		ec2Chr(ec->begin, ec->end, geneRange->getOffsets()[sp], rsa->getChrLen(sp, seqID), geneRange->getStrand(sp), chrStartStranded, chrEndStranded, false);
		getMSA_getEC_helper(msa, geneRange, sp, speciesNames, rsa, chrStartStranded, chrEndStranded, as, includeUnaligned);
		delete as;
	}
}

// main function to get the MSA form some OE
void getMSA(vector<string>& msa, vector<string>& speciesNames, RandSeqAccess* rsa, GeneMSA* geneRange, OrthoExon& oe, bool downcase, bool includeUnaligned){
	if(rsa==NULL || geneRange==NULL)
		return;

	msa.resize(speciesNames.size());

	for(int sp=0;sp<oe.orthoex.size();++sp){
		if(oe.orthoex[sp] != NULL){
			getMSA_getEC(msa, speciesNames, rsa, geneRange, sp, oe.orthoex[sp], downcase, includeUnaligned);
		}
	}	
}

