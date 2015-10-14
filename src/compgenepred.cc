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
#include "liftover.hh"
#include "intronmodel.hh"

#include <gsl/gsl_matrix.h>
#include <ctime>
#include <sys/stat.h>

CompGenePred::CompGenePred() : tree(Constant::treefile) {

    vector<string> speciesNames;
    tree.getSpeciesNames(speciesNames);
    cout<<"-------Speciesnames:--------"<<endl;
    for(int j=0; j<speciesNames.size(); j++){
      cout<<"species "<<j<<"\t"<<speciesNames[j]<<endl;
    }

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
    //evo.printBranchLengths();
    evo.computeLogPmatrices();
    OrthoGraph::tree = &tree;
    GeneMSA::setTree(&tree);
    OrthoGraph::numSpecies = OrthoGraph::tree->numSpecies();
    Boolean noprediction = false;
    Boolean useLocusTrees = false; // reestimate tree for each gene range

    // used for mapping an existing annotation to the other genomes
    // each CDSexon/intron that is supported by a hint (e.g. CDSexon/intron from the annotation) gets this
    // additional score to make sure that it is transferred to the other genomes and not the
    // other way round. It has to be at least has high as the maximum cost of an exon gain or loss event
    SpeciesGraph::maxCostOfExonLoss = - log((evo.getMu()+evo.getLambda())*evo.minBranchLength())*evo.getPhyloFactor();

#ifdef DEBUG
    cout << "-------------------------------\nparameters phylogenetic model\n-------------------------------" << endl;
    cout << "rate exon loss:\t" << evo.getMu() << endl;
    cout << "rate exon gain:\t" << evo.getLambda() << endl;
    cout << "phylo factor:\t" << evo.getPhyloFactor() <<  "\n-------------------------------" << endl;
#endif

    Constant::temperature = 3;
    Properties::assignProperty("temperature", Constant::temperature);
    if (Constant::temperature < 0){
	Constant::temperature = 0;
    }
    if (Constant::temperature > 7){
	Constant::temperature = 7;
    }
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
	cerr << "Warning: /CompPred/maxIterations was set to "<<maxIterations<<". At least one iteration must be made." << endl;
	maxIterations =1;
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
    string dd_param_s; 
    // parameter that defines the step size in dual decomposition.
    // If a range is given, all values in that range are tried until convergence is achieved
    vector<double> dd_factors; 
    try {
        dd_param_s = Properties::getProperty("/CompPred/dd_factor");
    } catch (...) {
        dd_param_s = "";
    }
    try{
	if(!dd_param_s.empty()){
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
		for (int i=0; i < rounds; i++){
		    dd_factors.push_back(start+i*(end-start)/(rounds-1));
		}	
	    }
	    else{
		double pos;	    
		if( !(stringstream(dd_param_s) >> pos))
		    throw ProjectError("Is not numeric.");
		dd_factors.push_back(pos);
	    }
	}
	else{
	    dd_factors.push_back(20); // by default, only do one round with dd parameter set to 20.
	}
    } catch (ProjectError e) {
	throw ProjectError("Format error parsing parameter --/CompPred/dd_factor=" + dd_param_s +".\n" + e.getMessage());
    }
    try {
	noprediction = Properties::getBoolProperty("noprediction");
    } catch (...) {}
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
    bool featureScoreHects;
    try {
        featureScoreHects = Properties::getBoolProperty("/CompPred/featureScoreHects");
    } catch (...) {
        featureScoreHects = false;
    }
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
    
    //initialize output files of initial gene prediction and optimized gene prediction
    vector<ofstream*> baseGenes = initOutputFiles(outdir,".mea"); // equivalent to MEA prediction
    vector<int> base_geneid(OrthoGraph::numSpecies, 1); // gene numbering
    vector<ofstream*> initGenes;
    vector<int> init_geneid(OrthoGraph::numSpecies, 1);
    if(featureScoreHects)
	initGenes = initOutputFiles(outdir,".init"); // score added to all orthologous exons and penalty added to all non orthologous exons, then global path search repeated
    vector<ofstream*> optGenes = initOutputFiles(outdir,".cgp");  //optimized gene prediction by applying majority rule move
    vector<int> opt_geneid(OrthoGraph::numSpecies, 1);
    vector<ofstream*> sampledGFs = initOutputFiles(outdir,".sampled_GFs"); // prints sampled exons/introns and their posterior probs to file

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
	lo.projectToAli(exoncands,alignedECs);
	lo.projectToGenome(alignedECs, seqRanges, exoncands); // inserts ECs that are mappable to other species both into alignedECs and exoncands

	// create additional ECs for each species and  insert them into exoncands and alignedECs
	vector<map<int_fast64_t,ExonCandidate*> > addECs(speciesNames.size()); // new ECs that need to be mapped to the alignment
        for (int s = 0; s < speciesNames.size(); s++) {
	    if (seqRanges[s]) {
		AnnoSequence *as = seqRanges[s];		
		// this is needed for IntronModel::dssProb in GenomicMSA::createExoncands
		namgene.getPrepareModels(as->sequence, as->length); 
		// identifies exon candidates in the sequence for species s
		geneRange->createExonCands(s, as->sequence, exoncands[s], addECs[s]);
	    }
	}
	exoncands.clear(); // not needed anymore, exoncands are now stored as a vector of lists of ECs in geneRange

	// map additional ECs to alignment space and insert them into alignedECs
	lo.projectToAli(addECs,alignedECs);
	addECs.clear(); // not needed anymore

	// create HECTS
	geneRange->createOrthoExons(alignedECs, &evo);

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
		orthograph.graphs[s] = new SpeciesGraph(&stlist, seqRanges[s], geneRange->getExonCands(s), speciesNames[s], 
							geneRange->getStrand(s), genesWithoutUTRs, onlyCompleteGenes, sampledGFs[s], overlapComp);
		orthograph.graphs[s]->buildGraph(meanIntrLen);
		//orthograph.graphs[s]->printGraph(outdir + speciesNames[s] + "." + itoa(GeneMSA::geneRangeID) + ".dot");
		
	    }
	}
 
	geneRange->printGeneRanges();
	if (Constant::exoncands) // by default, ECs are not printed
	    geneRange->printExonCands();
	try { // Kathrin Middendorf's playground
	    if (Properties::getBoolProperty("/CompPred/compSigScoring"))
		geneRange->comparativeSignalScoring(); 
	} catch (...) {}
	
	bool use_omega;
	try{ 
	    use_omega = Properties::getBoolProperty("/CompPred/omega");
	} catch (...){
	  if(Constant::logreg)
	    use_omega = true;
	  else
	    use_omega = false;
	}
	if(use_omega)
	    geneRange->computeOmegasEff(seqRanges, &ctree, &codonAli); // omega and number of substitutions is stored as OrthoExon attribute
	    //inefficient omega calculation, only use for debugging purpose 
	    //geneRange->computeOmegas(seqRanges, &ctree);
	
	if (conservation)
	    geneRange->calcConsScore(seqRanges, outdir);

	if (noprediction){
	    geneRange->printOrthoExons();
	}
	else{
	  list<OrthoExon> hects = geneRange->getOrthoExons();
	    orthograph.linkToOEs(hects); // link ECs in HECTs to nodes in orthograph	    
	    orthograph.globalPathSearch();
	    orthograph.outputGenes(baseGenes,base_geneid);
	    if(featureScoreHects){
		//add score for selective pressure of orthoexons
		orthograph.addScoreSelectivePressure();
		//determine initial path
		orthograph.globalPathSearch();
		orthograph.outputGenes(initGenes,init_geneid);
	    }	    
	    if(!orthograph.all_orthoex.empty()){
		// optimization via dual decomposition
		vector< list<Transcript*> *> genelist(OrthoGraph::numSpecies);
		orthograph.dualdecomp(evo,genelist,GeneMSA::geneRangeID-1,maxIterations, dd_factors);
		orthograph.filterGeneList(genelist,opt_geneid);
		orthograph.createOrthoGenes(geneRange);
		orthograph.printOrthoGenes();
		orthograph.printGenelist(optGenes);

	    }else{
		orthograph.outputGenes(optGenes, opt_geneid);
	    }
	    //geneRange->printOrthoExons(); //TODO: two copies of list<OrthoExon> (class geneRange and class OrthoGraph) -> replace one copy by a pointer to the other 
	    for(list<OrthoExon>::iterator it=orthograph.all_orthoex.begin(); it != orthograph.all_orthoex.end(); it++){
		geneRange->printSingleOrthoExon(*it,true);
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
    if(featureScoreHects)
	closeOutputFiles(initGenes);
    closeOutputFiles(baseGenes);
    closeOutputFiles(optGenes);
    closeOutputFiles(sampledGFs);

    // delete all trees                                           
    for(unordered_map< bit_vector, PhyloTree*, boost::hash<bit_vector>>::iterator topit = GeneMSA::topologies.begin(); topit !=GeneMSA::topologies.end(); topit++){
	delete topit->second;
    }                                                                                                                                                              
    GeneMSA::topologies.clear(); 

}
