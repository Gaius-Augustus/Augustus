/*
 * compgenepred.cc
 *
 * License: Artistic License, see file LICENSE.TXT or 
 *          https://opensource.org/licenses/artistic-license-1.0
 */
  
#include "compgenepred.hh"
#include "orthograph.hh"
#include "mea.hh"
#include "orthoexon.hh"
#include "contTimeMC.hh"
#include "liftover.hh"
#include "intronmodel.hh"
#include "train_logReg_param.hh"

#include <gsl/gsl_matrix.h>
#include <ctime>
#include <sys/stat.h>

/*
*   added by Giovanna Migliorelli 14.05.2020 
*   code responsible for serialization/deserialization
*/
#ifdef TESTING
#include "tokenizer.hh"

void serializeAlignment(string archname, Alignment* ali){

    cout << "Serializing alignment into " << archname << " archive" << endl;

    {
        std::ofstream ofs(archname.c_str());
        boost::archive::text_oarchive oa(ofs);
        oa << ali;
    }
}

void deserializeAlignment(string archname, Alignment*& ali){

    cout << "Deserializing alignment from " << archname << " archive" << endl;

    {
        std::ifstream ifs(archname.c_str());
        boost::archive::text_iarchive ia(ifs);
        ia >> ali;
    }
}

void serializeGeneMSA(string archname, GeneMSA* geneRange){

    cout << "Serializing gene range into " << archname << " archive" << endl;

    {
        std::ofstream ofs(archname.c_str());
        boost::archive::text_oarchive oa(ofs);
        oa << geneRange;
    }
}

void deserializeGeneMSA(string archname, GeneMSA*& geneRange){

    cout << "Deserializing gene range from " << archname << " archive" << endl;

    {
        std::ifstream ifs(archname.c_str());
        boost::archive::text_iarchive ia(ifs);
        ia >> geneRange;
    }
}
#endif

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
#ifdef M_MYSQL
	    rsa = new MysqlAccess(speciesNames);
#else
	    throw ProjectError("Database access not possible with this compiled version. Please recompile with flag MYSQL = true added to common.mk.");
#endif

	}
	else if(dbaccess.find('.') != string::npos){ // assuming sqlite access
	    cout << "# assuming an SQLite database" << endl;
	    if(Constant::speciesfilenames.empty())
		throw ProjectError("Missing parameter speciesfilenames.");
#ifdef M_SQLITE
	    rsa = new SQLiteAccess(dbaccess.c_str(), speciesNames);
#else
	    throw ProjectError("Database access not possible with this compiled version. Please recompile with flag SQLITE = true added to common.mk.");
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
    // added to choose between different modes : predict or run test
    #ifdef TESTING
    string testMode; 
    try {
        testMode = Properties::getProperty("/Testing/testMode");
    } catch (...) {
        testMode = "none"; 
    }

    if(testMode=="run"){
        cout << "Running test..." << endl;
        runPredictionOrTest();
        postprocTest();
    }
    else{
        cout << "Running prediction..." << endl; // (if you want to switch to test mode, set /Testing/testMode=run at the command line)
        runPredictionOrTest();
    }
    #else
    cout << "Running prediction..." << endl; // (if you want to switch to test mode, rebuild Augustus after setting TRAINING to true as explained in longrunning-test-CGP.md)
    runPredictionOrTest();
    #endif
}



void CompGenePred::runPredictionOrTest(){

    #ifdef TESTING
    string testMode; 
    try {
        testMode = Properties::getProperty("/Testing/testMode");
    } catch (...) {
        testMode = "predict"; 
    }
    #endif
   
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

    evo.getRateMatrices();
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
    unordered_map<string,int> ref_class; // reference classification (1/0) of CDS, introns for training
    try {
	noprediction = Properties::getBoolProperty("noprediction");
    } catch (...) {}
    if(Properties::hasProperty("referenceFile")){
	onlySampling = true;
	cout << "# AUGUSTUS is running in training mode. No prediction will be done!" << endl;
	try {
	    Constant::refSpecies = Properties::getProperty("refSpecies");
	    if (Properties::hasProperty("param_outfile")){
		cout << "# Using file " << Properties::getProperty("param_outfile") << " to store logReg parameters." << endl;
	    } else {
		cout << "# No outfile for logReg parameters specified. Writing parameters to " << Constant::configPath <<  "/cgp/log_reg_parameters_trained.cfg" << endl;
	    }
	    reference_from_file(&ref_class);
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

    #ifdef TESTING
    string outdir, outdirRun, outdirPrepare;  //direction for output files 
    
    if(testMode=="run"){
        try {
            outdir = Properties::getProperty("/CompPred/outdir");
            expandDir(outdir);
        } catch (...) {
            if(outdir.empty())
            throw ProjectError("Missing parameter /CompPred/outdir.");
        }

        outdirRun = outdir;
        size_t pos = outdir.find("run");
        if(pos == string::npos){
            throw ProjectError("/CompPred/outdir parameter should end with the string \"run\".");
        }
        outdirPrepare = outdir.substr(0, outdir.length() - 4) + "prepare/";
    }
    else{   // testMode = predict
        try {
            outdir = Properties::getProperty("/CompPred/outdir");
        } catch (...) {
            outdir = "";
        }
    }
    #else
    string outdir;  //direction for output files 
    try {
        outdir = Properties::getProperty("/CompPred/outdir");
    } catch (...) {
        outdir = "";
    }
    #endif







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

    codonevo.getRateMatrices();
    codonevo.computeLogPmatrices();
    
    // gsl_matrix *P = codonevo.getSubMatrixLogP(0.3, 0.25);
    // printCodonMatrix(P);
    GeneMSA::setCodonEvo(&codonevo);
    GenomicMSA msa(rsa);

    #ifdef TESTING
    if(testMode!="run")
        msa.readAlignment(Constant::alnfile);  // reads the alignment only if testMode != "run"
    #else
        msa.readAlignment(Constant::alnfile);  // reads the alignment
    #endif
    
    // msa.printAlignment("");    
    // rsa->printStats();
    // msa.compactify(); // Mario: commented out as this excludes paths through the alignment graph 
                         //(trivial mergers of neighboring alignments)

    #ifdef TESTING
    if(testMode=="run"){        
	    string dbdir;
        try {
            dbdir = Properties::getProperty("/Testing/workingDir");
            expandDir(dbdir);
        } catch (...) {
            if(outdir.empty())
            throw ProjectError("Missing parameter /Testing/workingDir.");
        }

        dbdir += "names/";
        
        msa.readNameDB(dbdir); // default was ../examples/cgp12way/names/ but currently an exception is raised if this parameter is not passed, temporarily added to solve a problem with string serialization 
    }
    else 
        msa.findGeneRanges(); // nontrivial summary of alignments, disabled in case of testing since deserialization does the job 
    #else
    msa.findGeneRanges(); // nontrivial summary of alignments, disabled in case of testing since deserialization does the job 
    #endif

    GeneMSA::openOutputFiles(outdir);

    int numGeneRange = 1;

    #ifdef TESTING
    int n = 0;
    struct stat buffer;  
    string filename;
    #endif
    
    #ifdef TESTING
    if(testMode=="run"){        
        while(true){
            filename = outdirPrepare + "generange_" + to_string(1+n) + ".bed";   // interspecies
            if(stat(filename.c_str(), &buffer) != 0)
                break;
            ++n;
        }
    }
    #endif

    #ifdef TESTING
    Alignment* ali = NULL;
    vector<list<tuple<string,int,int> > > grlist(n);        // contains intervals for all generanges within the current chunk (the file has been created by prepareTest)
    vector<list<tuple<string,int,int> > > mergedlist(speciesNames.size());    // contains what remains of intervals for all generange within the current chunk after having merged original ones (the file has been created after prepareTest outside Augustus) 
    #endif
    
    #ifdef TESTING
    if(testMode=="run"){        
        // read bed containing original intervals for each gene range and bed containing the same inervals after merging by species  : required for conversion of alignment
        while(true){
            filename = outdirPrepare + "generange_" + to_string(numGeneRange) + ".bed";   // interspecies

            if(stat(filename.c_str(), &buffer) != 0)
                break;
                
            if(!readInterval(filename, grlist[numGeneRange-1]))
                cout << "File " << filename << " absent : cannot recover gene range " << numGeneRange << endl;
            cout << "read gene range num " << numGeneRange << " " << grlist[numGeneRange-1].size() << endl;
            int ii = 0;
            for(list<tuple<string,int,int> >::iterator it=grlist[numGeneRange-1].begin();it!=grlist[numGeneRange-1].end();++it, ++ii)
                cout << speciesNames[ii] << " " << get<0>(*it) << " " << get<1>(*it) << " " << get<2>(*it) << endl;
            ++numGeneRange;       
        }
        
        for(int s=0;s<speciesNames.size();++s){
            filename = outdirPrepare + speciesNames[s] + ".bed";   //removed here .MERGED                // intraspecies        
            if(!readInterval(filename, mergedlist[s]))
                cout << "File " << filename << " absent : cannot recover merged interval list for generange " << numGeneRange << " for species " << speciesNames[s] << endl;
        }

        numGeneRange = 1;
    }
    #endif
    
    while(true){
    GeneMSA *geneRange = NULL;
    
    #ifdef TESTING
    if(testMode=="run"){  

        // files containing serialized data for gene ranges are no longer named after the reference interval within the gr 
        filename = outdirPrepare + "generange_" + to_string(numGeneRange); // to_string(itGR->first) + "_" + to_string(itGR->second);	

        if(stat(filename.c_str(), &buffer) != 0){
            cout << "File " << filename << " absent" << endl;
            break;
        }
        else{		    
            // can be made static or unique for all geneRanges 
            deserializeAlignment(filename, ali);
        
            if(ali==NULL){
                ++numGeneRange;
                continue;
            }

            // restore string names we couldn't serialize, reset sequence lengths
            for(int r=0;r<ali->rows.size();++r){
                if(ali->rows[r]){
                    ali->rows[r]->seqID = msa.seqIDarhive2seqIDConversion(r, ali->rows[r]->seqIDarchive);   // old names! as before minimizing fasta
                    // postponed : rsa->setLength(r, ali->rows[r]->seqID, ali->rows[r]->chrLen);
                }
            }

            // extract new interval the current gr falls within new fasta (eg 200-400 completely included within 0-1000, uses 0-1000)
            list<tuple<string, int, int> >::iterator itGR = grlist[numGeneRange-1].begin();
            for (int s = 0; s < speciesNames.size(); s++, ++itGR) {
                // cout << "Extracting interval " << speciesNames[s] << " " << get<0>(*itGR) << " " << get<1>(*itGR) << " " << get<2>(*itGR) << endl;
                int newStart = -1, newEnd = -1;
                // search the larger interval the original interval has been merged into
                for(list<tuple<string,int,int> >::iterator itMERGED=mergedlist[s].begin();itMERGED!=mergedlist[s].end();++itMERGED){
                    if(get<0>(*itGR)==get<0>(*itMERGED) && get<1>(*itGR)>=get<1>(*itMERGED) &&  get<2>(*itGR)<= get<2>(*itMERGED)){// both BED
                        // cout << "Interval detected " << get<0>(*itMERGED) << " " << get<1>(*itMERGED)  << " " << get<2>(*itMERGED) << endl;
                        // maximal containing interval found
                        newStart = get<1>(*itMERGED);
                        newEnd = get<2>(*itMERGED) - 1; // BED no longer need, we close the interval
                        break;
                    }
                }
                if(newStart>-1 && newEnd>-1)
                    ali->convertAlignment(s, newStart, newEnd);
                else{
                    // cout << "ERROR " << newStart << " " << newEnd << " " << get<0>(*itGR) << " " << get<1>(*itGR) << " " << get<2>(*itGR) << endl;
                }
            }

            for(int r=0;r<ali->rows.size();++r){
                if(ali->rows[r]){
                    rsa->setLength(r, ali->rows[r]->seqID, ali->rows[r]->chrLen);
                    // cout << "SETLEN POSTPONED " << ali->rows[r]->chrLen << endl;
                }
            }

            // create a gene range from deserialized alignment and run predixtion over it
            geneRange = new GeneMSA(rsa, ali);
        }
    }
    else
    {
        geneRange = msa.getNextGene();
        if(geneRange == NULL)
            break;
    }    
    #else
    geneRange = msa.getNextGene();
    if(geneRange == NULL)
        break;
	#endif
    
    cout << "processing gene range number " << numGeneRange << endl;
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
	if (Constant::printOEs)
	    geneRange->printOrthoExons(hects);
	    
	// store hect features globally for training
	if (Properties::hasProperty("referenceFile")){
	    cout << "collecting sample features" << endl;
	    int speciesID = find(speciesNames.begin(), speciesNames.end(), Constant::refSpecies) - speciesNames.begin();
	    if (speciesID >= speciesNames.size()){
		throw ProjectError("Species " + Constant::refSpecies + " not found. Use one of the names specified in the alignment file as a reference!");
	    }
	    
	    if (!Constant::printExonCandsMSA){
		geneRange->collect_features(speciesID, &hects, orthograph.graphs[speciesID]);
	    } else{// training set generation for codon evolution model
                geneRange->getAllOEMsas(speciesID, &hects, &ref_class, seqRanges);
	    }
	}

	// delete sequences
	for (int i=0; i<seqRanges.size(); i++) {
            delete seqRanges[i];
	}
	// delete geneRange
	delete geneRange;
    ++numGeneRange;
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
  
    if (Properties::hasProperty("referenceFile")){
	// initialize training of log reg parameters
        if (! Constant::printExonCandsMSA)
            train_OEscore_params(speciesNames.size());
    }
  } 
}



#ifdef TESTING
/*
*   added by Giovanna Migliorelli 14.05.2020 
*   functions for moving from shifted (small db) to original coordinate system (full size sequences)
*/

enum SplitMode {open, closed};

bool splitString(string line, char splitAt, SplitMode modeleft, SplitMode moderight, string& left, string& right){
    // e.g. s = "aaa-bbb" splitAt = '-'
    size_t pos = line.find(splitAt);
    if(pos == string::npos)
        return false;
    if(modeleft==open)
        left = line.substr(0, pos);             // e.g. left = "aaa"
    else
        left = line.substr(0, pos+1);           // e.g. left = "aaa-"
    
    if(moderight==open)
        right = line.substr(pos+1);             // e.g. right = "bbb"
    else
        right = line.substr(pos);               // e.g. left = "-bbb"

    return true;
}

bool extractData(string line, int& shift, string& seqID, string& newline){
    seqID = "";
    shift = -1;
    
    // e.g. #----- prediction on sequence range chr1_892852_2501000:551407-640869 (89463bp) -----
    string strA, strB, strC, left, right;
	size_t pos;
    vector<string> tokens;

    pos = line.find(':');
    if(pos == string::npos)
        return false;           // the line should contain ':'
    strA = line.substr(0, pos); // strA = "#----- prediction on sequence range chr1_892852_2501000"
    strB = line.substr(pos+1);  // strB = "551407-640869 (89463bp) -----"

    Tokenizer s(strA, " ");
	while(s.NextTokenBlank())
		tokens.push_back(s.GetToken());
    
    if(tokens.empty())
        return false;       // should never be here 
    strA = tokens.back();   // strA = "chr1_892852_2501000"

    tokens.clear();
    Tokenizer s1(strA, "_");
	while(s1.NextTokenBlank())
		tokens.push_back(s1.GetToken());
    
    if(tokens.size()<3)
        return false;
    shift = atoi(tokens[tokens.size()-2].c_str());  // because the seq name can contain underscores as well
    seqID = tokens[0];
    for(int i=1;i<tokens.size()-2;++i)
        seqID += "_" + tokens[i];
    tokens.clear();

    splitString(strB, '(', open, closed, left, right);

    strC = right;
    splitString(left, '-', open, open, left, right);

    newline = "#----- prediction on sequence range " + seqID + ":" + to_string(shift + atoi(left.c_str())) + "-" + to_string(shift + atoi(right.c_str())) + " " + strC;
    return true;
}

bool shiftFeature(string line, string seqID, int shift, string& newline){
    vector<string> tokens;
    Tokenizer s(line, "\t");
	while(s.NextTokenBlank())
		tokens.push_back(s.GetToken());
    if(tokens.size()!=9)
        return false;   // not valid GFF
    newline = seqID + "\tAUGUSTUS\t" + tokens[2] + "\t" + to_string(shift + stoi(tokens[3].c_str())) + "\t" + to_string(shift + stoi(tokens[4].c_str())) + "\t" + tokens[5] + "\t" + tokens[6] + "\t" + tokens[7] + "\t" + tokens[8];
    return true;
}

/*
*   added on 14.05.2020 
*   reads a GFF containing the prediction obtained over minimal FASTAs and shift coordinates back to full length genomes  
*/
bool shiftGFF(string filenameIn, string filenameOut){
    ifstream ifs(filenameIn);
    if(!ifs.is_open()){
        cout << "shiftGFF failed to open " << filenameIn << endl;
        return false;
    }
    
    ofstream ofs(filenameOut);
    if(!ofs.is_open()){
        cout << "shiftGFF failed to open " << filenameOut << endl;
        ifs.close();
        return false;
    }

    string line = "", newline, seqID;
    int shift = -1;
    bool pending = false;   // to print real empty lines without printing any additional empty line at the end of the file 

    while(ifs){
        if(pending)
            ofs << line << endl;

        getline(ifs,line);
        if(line.find("prediction on sequence range") != string::npos){
            extractData(line, shift, seqID, newline);
            ofs << newline << endl;
            line = "";
            pending = false;
        }
        else{
            vector<string> tokens;
            Tokenizer s(line, "\t");
            while(s.NextTokenBlank())
                tokens.push_back(s.GetToken());
                    
            if(tokens.size() == 9 && tokens[1] == "AUGUSTUS"){
                // when here, both shift and seqID should be well defined (optionally handle errors)
                shiftFeature(line, seqID, shift, newline);
                ofs << newline << endl;
                line = "";
                pending = false;
            }
            else{
                // ofs << line << endl;
                pending = true;
            }
        }
    }
    
    return true;
}

/*
*   added on 14.05.2020 
*   functions for merging/reading/writing intervals for all generanges within the current chunk (BED format)  
*/
bool sortInterval(const tuple<string,int,int>& a, const tuple<string,int,int>& b){
    if(get<0>(a) == get<0>(b))
        return get<1>(a) <= get<1>(b);
    return get<0>(a)<=get<0>(b);
}

void mergeInterval(list<tuple<string,int,int> >& intervals){
    if(intervals.empty())
        return;
    
    stack<tuple<string,int,int> > s;
    tuple<string,int,int> intv;
    intervals.sort(sortInterval);
    intv = intervals.front();
    intervals.pop_front();
    
    s.push(intv);
    while(!intervals.empty()){
        intv = intervals.front();
        intervals.pop_front();

        if(get<0>(intv) != get<0>(s.top()))
            s.push(intv);   // different seqID
        else{
            if(get<1>(intv)>=get<2>(s.top()) || get<1>(intv)<=get<1>(s.top()))
                s.push(intv);   // no overlap
            else if(get<2>(intv)>get<2>(s.top()))
                get<2>(s.top()) = get<2>(intv); 
        }
    }

    while(!s.empty()){
        intervals.push_front(s.top());
        s.pop();
    }
}

void mergeIntervals(vector<string>& speciesNames, vector<list<tuple<string,int,int> > >& intervals){
    for(int s=0;s<speciesNames.size();++s){
        intervals[s].sort();
        intervals[s].unique();
        mergeInterval(intervals[s]);
    }
}

void writeIntervals(string dirname, vector<string>& speciesNames, vector<list<tuple<string,int,int> > >& intervals){
    string filename;
    for(int s=0;s<speciesNames.size();++s){
        filename = dirname + speciesNames[s] + ".bed";
        ofstream ofs(filename);
        if(!ofs.is_open())
            continue;
        for(list<tuple<string,int,int> >::iterator it=intervals[s].begin();it!=intervals[s].end();++it)
            ofs << get<0>(*it) << "\t" << get<1>(*it) << "\t" << get<2>(*it) << endl;
        ofs.close();
    }
}

bool CompGenePred::readInterval(string filename, list<tuple<string,int,int> >& grlist){
	ifstream ifs(filename);
	if(!ifs.is_open()){
		cout << "failed to open " << filename << endl;
		return false;
	}

    string line;
	vector<string> tokens;

	while(ifs){
		getline(ifs, line);

		if(line.empty())
            continue;
        Tokenizer s(line, "\t");
		while(s.NextTokenBlank())
			tokens.push_back(s.GetToken());

		grlist.push_back(make_tuple(tokens[0], stoi(tokens[1].c_str()), stoi(tokens[2].c_str())));
        
        tokens.clear();
        line = "";
    }

	ifs.close();
	return true;
}

void CompGenePred::postprocTest(){
    string outdir, filenameIn, filenameOut;
    map<string,string> speciesNames;
    
    try {
        outdir = Properties::getProperty("/CompPred/outdir");
        expandDir(outdir);
    } catch (...) {
        outdir = "";
    }

    speciesNames = getFileNames(Constant::speciesfilenames);

    for(map<string,string>::iterator it=speciesNames.begin();it!=speciesNames.end();++it){
        filenameIn = outdir + it->first + ".cgp.SHIFTED.gff";
        filenameOut = outdir + it->first + ".cgp.gff";

        cout << "Shifting GFFs " << filenameIn << " " << filenameOut << endl;
        shiftGFF(filenameIn, filenameOut);
    }
}

#endif
