
/**********************************************************************
 * file:    properties.cc
 * licence: Artistic Licence, see file LICENCE.TXT or
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:
 * authors: Stafilarakis, Mario Stanke (mario@gobics.de)
 *
 * date    |   author      |  changes
 * --------|---------------|------------------------------------------
 * 26.09.01| Stafilarakis  | creations of the class
 * 04.08.03| Stanke        | changed operator Boolean( )
 * 01.08.05| Stanke        | added new species
 * 13.08.05| Mario Stanke  | code species independent: species is prefix of config file name
 * 04.04.07| Mario Stanke  | added parameters
 **********************************************************************/

#include "properties.hh"

// standard C/C++ includes
#include <cstring>     // for tolower
#include <fstream>     // for ifstream
#include <sys/stat.h>  // for stat
#include <iomanip>     // for setfill,setw
#include <cstdlib>     // for getenv
#include <algorithm>   // for for_each
#include <iostream>
#include <unistd.h>
#include <limits.h>

void makelow(char& c) {
    c = tolower(c);
}
void downcase(string& s) {
    for_each(s.begin(), s.end(), makelow);
}


map<string, string> Properties::properties;
const char* Properties::parameternames[NUMPARNAMES]= 
{
"allow_hinted_splicesites",
"alnfile",
"alternatives-from-evidence",
"alternatives-from-sampling",
"/augustus/verbosity",
"/BaseCount/weighingType",
"/BaseCount/weightMatrixFile",
"bridge_genicpart_bonus",
"canCauseAltSplice",
"capthresh",
"cds",
CGP_CONFIG_KEY,
CGP_PARS_KEY,
"checkExAcc",
"codingseq",
"codonAlignmentFile",
"complete_genes",
"/CompPred/assmotifqthresh",
"/CompPred/assqthresh",
"/CompPred/dssqthresh",
"/CompPred/compSigScoring",
"/CompPred/conservation",
"/CompPred/covPen",
"/CompPred/ec_score",
"/CompPred/exon_gain",
"/CompPred/exon_loss",
"/CompPred/ali_error",
"/CompPred/computeNumSubs",
"/CompPred/maxIterations",
"/CompPred/mil_factor",
"/CompPred/ec_factor",
"/CompPred/ec_addend",
"/CompPred/ec_thold",
"/CompPred/ic_thold",
"/CompPred/genesWithoutUTRs",
"/CompPred/liftover_all_ECs",
"/CompPred/logreg",
"/CompPred/maxCov",
"/CompPred/omega",
"/CompPred/num_bins",
"/CompPred/num_omega",
"/CompPred/num_features",
"/CompPred/scale_codontree",
"/CompPred/oeExtensionWidth",
"/CompPred/onlySampledECs",
"/CompPred/only_species",
"/CompPred/outdir_orthoexons",
"/CompPred/outdir",
"/CompPred/printOEcodonAli",
"/CompPred/printConservationWig",
"/CompPred/phylo_factor",
"/CompPred/phylo_model",
"/CompPred/dd_factor",
"/CompPred/dd_rounds",
"/CompPred/dd_step_rule",
"/CompPred/dualdecomp",
"/CompPred/overlapcomp",
"/CompPred/lambda",
"/CompPred/alpha",
"/CompPred/configFile",
"/CompPred/parsFile",
"/Constant/almost_identical_maxdiff",
"/Constant/amberprob",
"/Constant/ass_end",
"/Constant/ass_maxbinsize",
"/Constant/ass_start",
"/Constant/ass_upwindow_size",
"/Constant/decomp_num_at",
"/Constant/decomp_num_gc",
"/Constant/decomp_num_steps",
"/Constant/dss_end",
"/Constant/dss_maxbinsize",
"/Constant/dss_start",
"/Constant/gc_range_max",
"/Constant/gc_range_min",
"/Constant/init_coding_len",
"/Constant/intterm_coding_len",
"/Constant/max_contra_supp_ratio",
"/Constant/min_coding_len",
"/Constant/ochreprob",
"/Constant/opalprob",
"/Constant/probNinCoding",
"/Constant/subopt_transcript_threshold",
"/Constant/tis_maxbinsize",
"/Constant/trans_init_window",
"/Constant/tss_upwindow_size",
"contentmodels",
"CRF",
"CRF_N",
"CRF_TRAIN",
"CRFtrainCDS",
"CRFtrainIntron",
"CRFtrainIgenic",
"CRFtrainSS",
"CRFtrainUTR",
"CRFtrainTIS",
"dbaccess",
"dbhints",
"/EHMMTraining/state00",
"/EHMMTraining/state01",
"/EHMMTraining/state02",
"/EHMMTraining/state03",
"/EHMMTraining/statecount",
"/EHMMTraining/trainfile",
"emiprobs",
"errfile",
"evalset",
"exoncands",
"/ExonModel/etorder",
"/ExonModel/etpseudocount",
"/ExonModel/exonlengthD",
"/ExonModel/infile",
"/ExonModel/k",
"/ExonModel/lenboostE",
"/ExonModel/lenboostL",
"/ExonModel/maxexonlength",
"/ExonModel/minexonlength",
"/ExonModel/minPatSum",
"/ExonModel/minwindowcount",
"/ExonModel/outfile",
"/ExonModel/patpseudocount",
"/ExonModel/slope_of_bandwidth",
"/ExonModel/tisalpha",
"/ExonModel/tis_motif_memory",
"/ExonModel/tis_motif_radius",
"/ExonModel/verbosity",
"exonnames",
EXTRFILE_KEY,
"GCwinsize",
"GD_stepsize",
"/genbank/verbosity",
"gff3",
HINTSFILE_KEY,
"/HMMTraining/savefile",
"/IGenicModel/infile",
"/IGenicModel/k",
"/IGenicModel/outfile",
"/IGenicModel/patpseudocount",
"/IGenicModel/verbosity",
"/IntronModel/allow_dss_consensus_gc",
"/IntronModel/ass_motif_memory",
"/IntronModel/ass_motif_radius",
"/IntronModel/asspseudocount",
"/IntronModel/d",
"/IntronModel/dssneighborfactor",
"/IntronModel/dsspseudocount",
"/IntronModel/infile",
"/IntronModel/k",
"/IntronModel/minwindowcount",
"/IntronModel/non_ag_ass_prob",
"/IntronModel/non_gt_dss_prob",
"/IntronModel/outfile",
"/IntronModel/patpseudocount",
"/IntronModel/sf_with_motif",
"/IntronModel/slope_of_bandwidth",
"/IntronModel/splicefile",
"/IntronModel/ssalpha",
"/IntronModel/verbosity",
"introns",
"keep_viterbi",
"label_flip_prob",
"learning_rate",
"lossweight", // temp
"locustree",
"maxDNAPieceSize",
"maxOvlp",
"maxtracks",
"max_sgd_epoch",
"mea",
"mea_evaluation",
"/MeaPrediction/no_compatible_edges",
"/MeaPrediction/alpha_E",
"/MeaPrediction/alpha_I",
"/MeaPrediction/x0_E",
"/MeaPrediction/x0_I",
"/MeaPrediction/x1_E",
"/MeaPrediction/x1_I",
"/MeaPrediction/y0_E",
"/MeaPrediction/y0_I",
"/MeaPrediction/i1_E",
"/MeaPrediction/i1_I",
"/MeaPrediction/i2_E",
"/MeaPrediction/i2_I",
"/MeaPrediction/j1_E",
"/MeaPrediction/j1_I",
"/MeaPrediction/j2_E",
"/MeaPrediction/j2_I",
"/MeaPrediction/weight_base",
"/MeaPrediction/r_be",
"/MeaPrediction/r_bi",
"/MeaPrediction/weight_exon",
"/MeaPrediction/weight_gene",
"/MeaPrediction/weight_utr",
"minexonintronprob",
"min_intron_len",
"minmeanexonintronprob",
"mutationRate",
"noInFrameStop",
"noprediction",
EXTERNAL_KEY, // optCfgFile
"printOEs",
"outfile",
"param_outfile",
"predictionEnd",
"predictionStart",
"print_blocks",
"print_utr",
"progress",
"protein",
"/ProteinModel/allow_truncated",
"/ProteinModel/exhaustive_substates",
"/ProteinModel/block_threshold_spec",
"/ProteinModel/block_threshold_sens",
"/ProteinModel/blockpart_threshold_spec",
"/ProteinModel/blockpart_threshold_sens",
"/ProteinModel/global_factor_threshold",
"/ProteinModel/absolute_malus_threshold",
"/ProteinModel/invalid_score",
"/ProteinModel/weight",
"proteinprofile",
"referenceFile",
"refSpecies",
"rescaleBoni",
"rLogReg",
"sample",
"scorediffweight", // temp
SINGLESTRAND_KEY,
"softmasking",
"speciesfilenames",
"start",
"stop",
"stopCodonExcludedFromCDS",
"strand",
"temperature",
"tieIgenicIntron",
"trainFeatureFile",
"translation_table",
"treefile",
"truncateMaskedUTRs",
"tss",
"tts",
"uniqueCDS",
"uniqueGeneId",
"useAminoAcidRates",
"useNonCodingModel",
"use_sgd",
"useTFprob",
UTR_KEY,
"/UtrModel/d_polya_cleavage_max",
"/UtrModel/d_polya_cleavage_min",
"/UtrModel/d_polyasig_cleavage",
"/UtrModel/d_tss_tata_max",
"/UtrModel/d_tss_tata_min",
"/UtrModel/exonlengthD",
"/UtrModel/infile",
"/UtrModel/k",
"/UtrModel/max3singlelength",
"/UtrModel/max3termlength",
"/UtrModel/maxexonlength",
"/UtrModel/minwindowcount",
"/UtrModel/outfile",
"/UtrModel/patpseudocount",
"/UtrModel/prob_polya",
"/UtrModel/slope_of_bandwidth",
"/UtrModel/tata_end",
"/UtrModel/tata_pseudocount",
"/UtrModel/tata_start",
"/UtrModel/tss_end",
"/UtrModel/tss_start",
"/UtrModel/tssup_k",
"/UtrModel/tssup_patpseudocount",
"/UtrModel/tts_motif_memory",
"/UtrModel/utr3patternweight",
"/UtrModel/utr5patternweight",
"/UtrModel/utr3prepatternweight",
"/UtrModel/utr5prepatternweight",
"/UtrModel/verbosity"};


void Properties::readFile( string filename, int fileTypeNr ) throw( PropertiesError ){
    ifstream strm;
    strm.open( filename.c_str() );
    if (!strm) {
	cerr << "Could not find the config file " << filename << "." << endl;
	int len = filename.length();
	string::size_type pos = filename.rfind ("/", len-1);
	if (pos != string::npos)
	    filename = filename.substr(pos+1, len-pos-1);
	if (filename != "") {
	    cerr << " Looking for " << filename << " in the configuration directory instead ...";
	    strm.clear();
	    filename = properties[CFGPATH_KEY]+filename;
	    strm.open(filename.c_str());
	    if (!strm)  cerr << " not";
	    cerr << " found." << endl;
	}
	if (!strm)
	    throw PropertiesError( "Properties::readFile: "
				   "Could not open the this file!" );
    }
    if(fileTypeNr == 0){
      while( strm )
	readLine( strm );
    }else if(fileTypeNr == 1){
      readCGPconfigFile(strm);
    }else if(fileTypeNr == 2){
      readCGPparsFile(strm);
    }
    strm.close();
}

void Properties::init( int argc, char* argv[] ){
    string name;
    int arg;

/*
 * The rightmost parameter without '--' is the input file, unless
 * the parameter 'queryfile' is specified
 *
 */
    arg = argc;
    while( arg > 1 ){
        arg--;
	string argstr(argv[arg]);
	if (argstr == "--version")
	    throw HelpException(GREETING);
	if (argstr == "--paramlist") {
	    ostringstream strm;
	    strm << "Possible parameter names are:" << endl;
	    for (int i=0; i<NUMPARNAMES; i++)
		strm << parameternames[i] << endl;
	    throw HelpException(strm.str());
	}
	if (argstr == "--help")
	    throw HelpException(HELP_USAGE);
	if (argstr.substr(0,2) == "--") {
	    string::size_type pos = argstr.find('=');
	    name = argstr.substr(0,pos).substr(2);
	    if( name == GENEMODEL_KEY ||
		name == NONCODING_KEY ||
		name == SINGLESTRAND_KEY ||
		name == SPECIES_KEY ||
		name == EXTERNAL_KEY ||
		name == CFGPATH_KEY ||
		name == ALN_KEY ||
		name == TREE_KEY ||
		name == DB_KEY ||
		name == SEQ_KEY ||
		name == CODONALN_KEY ||
		name == REF_EXON_KEY ||
		name == CGP_CONFIG_KEY ||
		name == CGP_PARS_KEY)
	    {
		if (pos >= argstr.length()-1)
		    throw ProjectError(string("Wrong argument format for ") +  name + ". Use: --argument=value");
		properties[name] = argstr.substr(pos+1);
	    }
	} else if (!hasProperty(INPUTFILE_KEY))
	    properties[INPUTFILE_KEY] = argv[arg];
	else
	    throw ProjectError(string("Error: 2 query files given: ") + properties[INPUTFILE_KEY] + " and " + argv[arg] +"."
			       + "\nparameter names must start with '--'");
    }

    Properties::assignProperty(TREE_KEY, Constant::treefile);
    Properties::assignProperty(SEQ_KEY, Constant::speciesfilenames);
    Properties::assignProperty(DB_KEY, Constant::dbaccess);
    Properties::assignProperty(ALN_KEY, Constant::alnfile);
    Properties::assignProperty(REF_EXON_KEY, Constant::referenceFile);
    Properties::assignProperty(CGP_CONFIG_KEY, Constant::cgpConfigFile);
    Properties::assignProperty(CGP_PARS_KEY, Constant::cgpParsFile);
    if(Constant::codonalnfile.empty()){
      if (!Constant::alnfile.empty() && !Constant::treefile.empty() && (!Constant::speciesfilenames.empty() || !Constant::dbaccess.empty())){
	Constant::MultSpeciesMode = true;
      } else if (!(Constant::alnfile.empty() && Constant::speciesfilenames.empty() && Constant::dbaccess.empty()) ){
	throw ProjectError("In comparative gene prediction mode you must specify parameters alnfile, treefile\n\
        and one of the following combinations of parameters\n\n		\
        - dbaccess (retrieving genomes from a MySQL db)\n		\
        - speciesfilenames and dbaccess (retrieving genomes from an SQLite db)\n \
        - speciesfilenames (retrieving genomes from flat files)\n\n	\
        In single species mode specify none of these parameters.\n");
      }
    }
    // set configPath variable
    string& configPath = properties[CFGPATH_KEY]; // first priority: command line
    if (configPath == "") {
	char *dir = getenv(CFGPATH_KEY); // second priority: environment variable
	if (dir){
	    configPath = string(dir);
	} else { // third priority: relative to the path of the executable
	    configPath = findLocationOfSelfBinary();
	}
    }
    if (configPath[configPath.size()-1] != '/')
	configPath += '/';  // append slash if neccessary

    // does the directory actually exist?
    struct stat buffer;
    if( stat(configPath.c_str(), &buffer) == -1 || !S_ISDIR(buffer.st_mode))
	throw ProjectError(configPath + " is not a directory. Could not locate directory " CFGPATH_KEY ".");
    properties[CFGPATH_KEY] = configPath;

    string optCfgFile = "";
    if (hasProperty(EXTERNAL_KEY)) {
	optCfgFile = expandHome(properties[EXTERNAL_KEY]);
	if (!hasProperty(SPECIES_KEY))
	// read in species and logReg parameters from extra config file
	    readFile(optCfgFile);
    } else if (Constant::MultSpeciesMode) {
	optCfgFile = expandHome(configPath + "cgp/log_reg_parameters_default.cfg");
	try {
	    readFile(optCfgFile);
	} catch (PropertiesError e){
	  //throw ProjectError("Default configuration file with parameters from logistic regression " + optCfgFile + " not found!");
	}
    }
  
    
    // determine species
    if (!hasProperty(SPECIES_KEY))
	throw ProjectError("No species specified. Type \"augustus --species=help\" to see available species.");
    string& speciesValue = properties[SPECIES_KEY];
    string speciesParFileName = speciesValue + "_parameters.cfg";

    if (speciesValue == "help")
	throw HelpException(SPECIES_LIST);

    // check query filename
    if (hasProperty(INPUTFILE_KEY)){// expand the ~ if existent
	string filename =  properties[INPUTFILE_KEY];
	filename = expandHome(filename);
	properties[INPUTFILE_KEY] = filename;
    }

    // read in file with species specific parameters
    string& speciesDir = properties[SPECIESDIR_KEY];
    speciesDir = SPECIES_SUBDIR + speciesValue + "/";
    try {
	readFile(configPath + speciesDir + speciesParFileName);
    } catch (PropertiesError e){
	throw ProjectError("Species-specific configuration files not found in " + configPath + SPECIES_SUBDIR + ". Type \"augustus --species=help\" to see available species.");
    }


    if (Constant::MultSpeciesMode || hasProperty(REF_EXON_KEY)){
      string cgp_file = "";
      Constant::lr_features.push_back(new LRfeatureGroup(0, "intercept", 0, 1)); // exon intercept
      Constant::lr_features.push_back(new LRfeatureGroup(20, "intron_intercept", 0, 1)); // intron intercept
      if (hasProperty(CGP_CONFIG_KEY))
	cgp_file = expandHome(properties[CGP_CONFIG_KEY]);
      else     
	cgp_file = expandHome(configPath + "cgp/logReg_feature.cfg");
      readFile(cgp_file, 1);
      
      if(!hasProperty(REF_EXON_KEY)){
	if (hasProperty(CGP_PARS_KEY))
	  cgp_file = expandHome(properties[CGP_PARS_KEY]);
	else
	  cgp_file = expandHome(configPath + "cgp/logReg_feature.pars");
	readFile(cgp_file, 2);
      }

      struct sortFeature {
	bool operator() (LRfeatureGroup *i,LRfeatureGroup *j) { return (i->id < j->id);}
      } featureSort;

      sort(Constant::lr_features.begin(), Constant::lr_features.end(), featureSort);
    }
    

// if mult-species mode is turned on, try to read cgp config file (not needed anymore)
    /*if(Constant::MultSpeciesMode){
	string cgpParFileName = speciesValue + "_parameters.cgp.cfg";
	string savedSpecies = speciesValue;
	try {
	    readFile(configPath + speciesDir + cgpParFileName);
	} catch (PropertiesError e){
	    cerr <<"Resuming without "<< cgpParFileName << endl;
	}
	speciesValue = savedSpecies;
	}*/

    // read in extra config file (again), but do not change species
    if (optCfgFile != "") {
	string savedSpecies = speciesValue;
	readFile(optCfgFile);
	speciesValue = savedSpecies;
    }

    // read in the other parameters, command line parameters have higher priority
    // than those in the config files
    arg = argc;
    while ( arg > 1 ){
        arg--;
 	string argstr(argv[arg]);
	if (argstr.substr(0,2) != "--")
	    continue;
	argstr.erase(0,2);
	string::size_type pos = argstr.find('=');
	name = argstr.substr(0,pos);
	if (name == GENEMODEL_KEY || name == NONCODING_KEY || name == SINGLESTRAND_KEY ||
	    name == SPECIES_KEY || name == CFGPATH_KEY ||
	    name == EXTERNAL_KEY || name == ALN_KEY ||
	    name == TREE_KEY || name == DB_KEY || name == SEQ_KEY || name == CODONALN_KEY || name == REF_EXON_KEY || name == CGP_CONFIG_KEY || name == CGP_PARS_KEY)
	    continue;
	if (pos == string::npos)
	    throw PropertiesError(string("'=' missing for parameter: ") + name);
	
	// check that name is actually the name of a parameter
	for (int i=0; i<NUMPARNAMES; i++)
	    if (name == parameternames[i]) {
		properties[name] = argstr.substr(pos+1);
		goto ENDWHILE;
	    }
	throw PropertiesError("Unknown parameter: \"" + name + "\". Type \"augustus\" for help.");
      ENDWHILE: {}
    }

    // singlestrand mode or not?
    bool singleStrand = hasProperty(SINGLESTRAND_KEY) && getBoolProperty(SINGLESTRAND_KEY);
    string strandCFGName = singleStrand? "singlestrand" : "shadow";

    // genemodel {partial, complete, atleastone, exactlyone, intronless, bacterium}
    string genemodelValue = hasProperty(GENEMODEL_KEY) ? properties[GENEMODEL_KEY] : "partial";
    if (genemodelValue != "partial" && genemodelValue != "complete" && genemodelValue != "atleastone" &&
	genemodelValue != "exactlyone" && genemodelValue != "intronless" && genemodelValue != "bacterium")
        throw ProjectError(string("Unknown value for parameter ") + GENEMODEL_KEY + ": " + genemodelValue);
    string transfileValue = "trans_" + strandCFGName + "_" + genemodelValue;

    bool utr_option_on = false;
    if (hasProperty(UTR_KEY)){
	try {
	    utr_option_on = getBoolProperty(UTR_KEY);
	} catch (...) {
	    throw ProjectError("Unknown option for parameter UTR. Use --UTR=on or --UTR=off.");
	}
    }
    bool nc_option_on = false;
    if (hasProperty(NONCODING_KEY)){
	try {
	    nc_option_on = getBoolProperty(NONCODING_KEY);
	} catch (...) {
	    throw ProjectError("Unknown option for parameter " NONCODING_KEY ". Use --"
			       NONCODING_KEY "=on or --" NONCODING_KEY "=off.");
	}
    }
    if (nc_option_on) {
	if (singleStrand || !(genemodelValue == "partial" || genemodelValue == "complete"))
	    throw ProjectError("Noncoding model (--" NONCODING_KEY "=1) only implemented with shadow and --genemodel=partial or --genemodel=complete.");
	if (!utr_option_on){
	    cerr << "Noncoding model (--" NONCODING_KEY "=1) only implemented with --UTR=on. Will turn UTR prediction on..." << endl;
	    utr_option_on = true;
	    properties[UTR_KEY] = "on";
	}
    }

    if (utr_option_on) {
	if (singleStrand || !(genemodelValue == "partial" || genemodelValue == "complete"))
	    throw ProjectError("UTR only implemented with shadow and partial or complete.");
	transfileValue += "_utr";
    } 
    if (nc_option_on)
	transfileValue += "_nc";

    transfileValue += ".pbl";
    properties[TRANSFILE_KEY] = transfileValue;

    // determine state config file
    string stateCFGfilename = "states" + ("_" + strandCFGName);
    if (genemodelValue == "atleastone" || genemodelValue == "exactlyone")
      stateCFGfilename += "_2igenic";
    else if (genemodelValue == "intronless"){
      stateCFGfilename += "_intronless";
    } else if (genemodelValue == "bacterium"){
      stateCFGfilename += "_bacterium";
      Constant::overlapmode = true;
    } else if (utr_option_on) {
	stateCFGfilename += "_utr";
	if (nc_option_on)
	    stateCFGfilename += "_nc";
    }

    stateCFGfilename += ".cfg";

    // read in config file with states
    readFile(configPath + MODEL_SUBDIR + stateCFGfilename);

    // reread a few command-line parameters that can override parameters in the model files
    arg = argc;
    while( arg > 1 ){
        arg--;
 	string argstr(argv[arg]);
	if (argstr.substr(0,2) != "--")
	    continue;
	argstr.erase(0,2);
	string::size_type pos = argstr.find('=');
	name = argstr.substr(0,pos);
	if (name == "/EHMMTraining/statecount" || name == "/EHMMTraining/state00" ||name == "/EHMMTraining/state01" ||
	    name == "/EHMMTraining/state02" || name == "/EHMMTraining/state03"){
	    properties[name] = argstr.substr(pos+1);
	}
    }

    // expand the extrinsicCfgFile filename in case it is specified
    if (hasProperty(EXTRFILE_KEY))
	properties[EXTRFILE_KEY] = expandHome(properties[EXTRFILE_KEY]);

    Properties::assignProperty("softmasking", Constant::softmasking);
    Properties::assignProperty("dbhints", Constant::dbhints);
    if(!Constant::MultSpeciesMode && Constant::dbhints){
	Constant::dbhints=false;
	cerr << "Warning: the option dbhints is only available in comparative gene prediction. Turned it off." << endl;
    }

    // expand hintsfilename
    if (hasProperty(HINTSFILE_KEY) || Constant::softmasking || Constant::dbhints){
	if (hasProperty(HINTSFILE_KEY)){
	    string& filename = properties[HINTSFILE_KEY];
	    filename = expandHome(filename);
	    properties[EXTRINSIC_KEY] = "true";
	}
	if (!hasProperty(EXTRFILE_KEY)) {
	    properties[EXTRFILE_KEY] = configPath + EXTRINSIC_SUBDIR + "extrinsic.cfg";
	    if (Constant::MultSpeciesMode)
		properties[EXTRFILE_KEY] = configPath + EXTRINSIC_SUBDIR + "extrinsic-cgp.cfg";
#ifdef DEBUG
	    cerr << "# No extrinsicCfgFile given. Taking default file: " << properties[EXTRFILE_KEY] << endl;
#endif
	}
    }
}

Boolean Properties::hasProperty( string name) {
    return properties.count(name)>0;
}


void Properties::addProperty( string name, string value ) {
    properties[name] = value;
}


Integer Properties::getIntProperty( string name ) {
    Integer val;
    istringstream strm(getProperty(name));
    if( !(strm >> val) )
        throw ValueFormatError(name, strm.str(), "Integer");
    return val;
}

Double Properties::getDoubleProperty( string name ){
    Double val;
    istringstream strm(getProperty(name));
    if ( !(strm >> val) )
        throw ValueFormatError(name, strm.str(), "LLDouble");
    return val;
}

double Properties::getdoubleProperty( string name ) {
    double val;
    istringstream strm(getProperty(name));
    if ( !(strm >> val) )
        throw ValueFormatError(name, strm.str(), "double");
    return val;
}

Boolean Properties::getBoolProperty( string name ) {
    string str(getProperty(name));
    downcase(str);
    if( str == "true" || str=="1" || str == "t"  || str=="on" || str == "yes" || str == "y")
        return true;
    if( str == "false" || str=="0" || str=="f" || str=="off" || str == "no" || str == "n")
        return false;
    throw ValueFormatError(name, str, "Boolean");
}

const char* Properties::getProperty( string name ) {
    if (!hasProperty(name))
	throw KeyNotFoundError(name);
    return properties[name].c_str();
}

const char* Properties::getProperty(string name, int index) {
    ostringstream strm;
    strm << name << setfill('0') << setw(2) << index;
    return getProperty(strm.str());
}

void Properties::readLine( istream& strm ) {
    string line, name, value;
    getline( strm, line );

    istringstream istrm( line.c_str() );
    istrm >> name >> value;
    if (name =="include"){
	/*
	 * include another file by recursively calling readfile
	 * 
	 */
	readFile(properties[CFGPATH_KEY] + value);
    } else
      if( name != "" && value != "" && name[0] != '#')
	properties[name] = value;
}

void Properties::readCGPconfigFile(ifstream& strm){
  
  while( strm ){
    string line, feature_nr, feature_name, if_std, num_bins;
    getline( strm, line );
    if(line[0] == '#')
      continue;
    istringstream istrm( line.c_str() );
    istrm >> feature_nr >> feature_name >> if_std >> num_bins;
    if(feature_nr != "" && feature_name != "" && if_std != "" && num_bins != ""){
      int nb = stoi(num_bins);
      bool std;
      if(nb == 2)
	throw PropertiesError("Properties::readCGPconfigFile: "
                              "number of bins must be either 1 or at least 3.");
      if(if_std == "1")
	std = true;
      else if(if_std == "0")
	std = false;
      else
	throw PropertiesError("Properties::readCGPconfigFile: "
			      "CGP config file format error: standartization values need to be 1 or 0.");
      LRfeatureGroup* f;
      if(nb == 1)
	f = new LRfeatureGroup(stoi(feature_nr), feature_name, std, nb);
      else
	f = new LRbinnedFeatureGroup (stoi(feature_nr), feature_name, std, nb);
      Constant::lr_features.push_back(f);
      
    }else{
      if(feature_nr[0] != '#' && feature_nr != "")
	throw PropertiesError("Properties::readCGPconfigFile: "
			      "CGP config file has wrong format!");
    }
  }
}

void Properties::readCGPparsFile(ifstream& strm){
  
  while( strm ){
    string line, feature_nr, feature_name, num_bins, binning_bounds, weight_list;
    getline( strm, line );
    if(line[0] == '#')
      continue;
    istringstream istrm( line.c_str() );
    istrm >> feature_nr >> feature_name >> num_bins >> binning_bounds >> weight_list;
    if(feature_nr != "" && feature_name != "" && num_bins != "" && binning_bounds != "" && weight_list != ""){
      int id = stoi(feature_nr);
      int nb = stoi(num_bins);
      for(int i=0; i<Constant::lr_features.size(); i++){
	if(Constant::lr_features[i]->id == id){
	  if(Constant::lr_features[i]->descript != feature_name || Constant::lr_features[i]->num_bins < nb)
	    throw PropertiesError("Properties::readCGPparsFile: "
				  "Error in readCGPparsFile: number of bins and/or feature description not consistent with the CGP config file.");
	  
	  if(nb == 1){
	    int str_length;
	    if(weight_list[weight_list.size()-1] == ',')
	      str_length = weight_list.size()-1;
	    else
	      str_length = weight_list.size();
	    Constant::lr_features[i]->setWeight(stod(weight_list.substr(0, str_length)));
	  }else{
	    LRbinnedFeatureGroup* bf;
	    try{
	      bf = static_cast<LRbinnedFeatureGroup*>(Constant::lr_features[i]);
	    }catch(...){
              throw PropertiesError("Properties::readCGPparsFile: "
                                    "number of bins does not agree with the number of bins in the config file.");
            }
	    if(bf->num_bins != nb){
	      bf->num_bins = nb;
	      bf->num_features = nb-1;
	      bf->bb.resize(nb-1, 0);
	      bf->weights.resize(nb-1, 0);
	    }
	    vector<double> items;
	    istringstream b(binning_bounds);
	    string buffer;
	    while(getline(b, buffer, ',')){
	      items.push_back(stod(buffer));
	    }
	    if(items.size() != bf->bb.size())
	      throw PropertiesError("Properties::readCGPparsFile: "
				    "number of bins does not agree with the number of provided binning boundaries.");
	    bf->bb = items;
	    items.clear();
	    istringstream w(weight_list);
	    while(getline(w, buffer, ',')){
	      items.push_back(stod(buffer));
	    }
	    if(items.size() != bf->weights.size())
              throw PropertiesError("Properties::readCGPparsFile: "
                                    "number of bins does not agree with the number of provided weights.");
	    bf->weights = items;
	    break;
	  }
	}
      }
    }else{
      if(feature_nr[0] != '#' && feature_nr != "")
	throw PropertiesError("Properties::readCGPparsFile: "
			      "CGP parameter file has wrong format! Make sure the first character in each line is not a white space.");
    }
  }  
}

LRfeatureGroup* Properties::findFeatureGroup(int id){
  
  for(int i = 0; i < Constant::lr_features.size(); i++){
    if(Constant::lr_features[i]->id == id)
      return Constant::lr_features[i];
  }
  return NULL;
}

double Properties::calculate_feature_score(Traits* t){

  double score = 0;
  for(int i = 0; i < Constant::lr_features.size(); i++)
    score += calculate_feature_score(Constant::lr_features[i], t);
  return score;
}

double Properties::calculate_feature_score(LRfeatureGroup* lr_feature, Traits* t){

  int id = lr_feature->id;
  bool hasTrait = true;
  double x = getFeature(id, t, &hasTrait);
  
  if(hasTrait){
    if(lr_feature->num_bins==1) // not binned features
      return (x * lr_feature->weight);
    else{        // binned features
      double score = 0;
      LRbinnedFeatureGroup* lr_bf = static_cast<LRbinnedFeatureGroup*>(lr_feature);
      vector<double> bf = getBinnedFeature(lr_bf, x);
      if(bf.size() != lr_bf->weights.size())
	throw ProjectError("Error in Properties::calculate_feature_score: Size of binned feature vector is not equal to size of weight vector!");
      for(int i = 0; i < bf.size(); i++)
	score += bf[i] * lr_bf->weights[i];
      return score;
    }
  }else{
    return 0;
  }
}


// for prediction devide into single and ortho exon feature sets
double Properties::calculate_single_feature_score(Traits* t){
  double score = 0;
  for(int i = 0; i < Constant::lr_features.size(); i++)
    if( Constant::lr_features[i]->id <= 5 )
      score += calculate_feature_score(Constant::lr_features[i], t);
  return score;
}

double Properties::calculate_OE_feature_score(Traits* t){
  double score = 0;
  for(int i = 0; i < Constant::lr_features.size(); i++)
    if( Constant::lr_features[i]->id > 5 && Constant::lr_features[i]->id < 20 )
      score += calculate_feature_score(Constant::lr_features[i], t);
  return score;
}

double Properties::calculate_intron_feature_score(Traits* t){
  double score = 0;
  for(int i = 0; i < Constant::lr_features.size(); i++)
    if( Constant::lr_features[i]->id >= 20 )
      score += calculate_feature_score(Constant::lr_features[i], t);
  return score;
}



double Properties::getFeature(int id, Traits *t, bool *hasTrait){
  if(id <= 5)
    return getSingleFeature(id,t, hasTrait);
  else if(id < 20)
    return getOEfeature(id,static_cast<OEtraits*>(t), hasTrait);
  else
    return getIntronFeature(id,t, hasTrait);
}

vector<double> Properties::getBinnedFeature(LRbinnedFeatureGroup *lr_bf, double x){

  int curr_bin = 0;
  int num_features = lr_bf->num_features;
  vector<double> f;
  // find bin x is in
  for(int j = 0; j < num_features; j++)
    if(x <= lr_bf->bb[j]){
      curr_bin = j;
      break;
    }
  if(x > lr_bf->bb.back())
    curr_bin = num_features;
  
  for(int j = 0; j < num_features; j++){	    
    if( (j == 0 && curr_bin == 0) || (j == num_features - 1 && curr_bin == num_features) )
      f.push_back(1);
    else if(curr_bin == j+1 && j >= 0 && j < num_features - 1)
      f.push_back((lr_bf->bb[j+1] - x) / (lr_bf->bb[j+1] - lr_bf->bb[j]));
    else if(curr_bin == j && j > 0 && j <= num_features - 1)
      f.push_back((x - lr_bf->bb[j-1]) / (lr_bf->bb[j] - lr_bf->bb[j-1]));
    else
      f.push_back(0);
  }
  return f;
}

double Properties::getSingleFeature(int id, Traits *t, bool *hasTrait){
  bool ht = true;
  double x = 0;
  switch(id){
  case 0: x = 1;                           // intercept
    break;
  case 1: x = log( t->getLength() );       // log exon length 
    break;
  case 2:                                  // posterior probability
    if(t->hasPostProb())          
      x = t->getPostProb(); 
    else
      ht = false;
    break;
  case 3:                                  // mean base probability
    if(t->hasMeanBaseProb())
      x = t->getMeanBaseProb();
    else
      ht = false;
    break;
  case 4: x = ( t->isSampled() ) ? 1 : 0;  // is not sampled
    break;
  case 5: x = ( t->isOE() ) ? 0 : 1;       // is not an OE
    break;
  default:
    throw ProjectError("feature number " + to_string(id) + " does not exist!");
    break;
  }
  if(hasTrait)
    *hasTrait = ht;
  return x;
}


double Properties::getOEfeature(int id, OEtraits* t, bool *hasTrait){

  bool ht = true;
  double x = 0;
  switch(id){
  case 6: x = ( t->hasOmega() ) ? 0 : 1;   // not having omega
    break;
  case 7:                                  // omega
    if(t->hasOmega())
      x = log(t->getOmega());
    else
      ht = false;
    break;
  case 8:                                  // variance of omega
    if(t->hasVarOmega())
      x = t->getVarOmega();
    else
      ht = false;
    break;
  case 9:                                  // conservation
    if(t->hasConservation())
      x = t->getConsScore();
    else
      ht = false;
    break;
  case 10:                                  // containment
    if(t->hasContainment())
      x = t->getContainment();
    else
      ht = false;
    break;
  case 11:                                 // diversity
    if(t->hasDiversity())
      x = t->getDiversity();
    else
      ht = false;
    break;
  case 12: x = t->getNumExons();           // number of exons involved in OE
    break;
  case 13:                                 // conservation * diversity
    if(t->hasConservation() && t->hasDiversity())
      x = t->getConsScore() * t->getDiversity();
    else
      ht = false;
    break;
  case 14:                                 // omega * diversity
    if(t->hasOmega() && t->hasDiversity())
      x = t->getOmega() * t->getDiversity();
    else
      ht = false;
    break;
  case 15:                                 // boundary omega
    if(t->hasBoundaryOmegas()){
      double left_bs =  log( t->getLeftIntOmega() ) - log( t->getLeftExtOmega() );
      double right_bs = log( t->getRightIntOmega() ) - log( t->getRightExtOmega() );
      x = ( left_bs * exp( Constant::lambda * left_bs ) + right_bs * exp( Constant::lambda * right_bs ) )
	/ ( exp( Constant::lambda * left_bs ) + exp( Constant::lambda * right_bs ) );
    } else
      ht = false;
    break;
  case 16:                                // coding posterior probability
    if(t->hasCodingPostProb())
      x = t->getCodingPostProb();
    else
      ht = false;
    break;
  case 17:  // score using inverse sigmoid value of coding posterior probability
    if(t->hasCodingPostProb()){
      if(t->getCodingPostProb() > 0 && t->getCodingPostProb() < 1)
	x = -log(1/t->getCodingPostProb() - 1);
      else if(t->getCodingPostProb() == 0)
	x = -50;
      else
	x= 50;
    }else
      ht = false;
    break;
  default:
    throw ProjectError("feature number " + to_string(id) + " does not exist!");
    break;
  }
 if(hasTrait)
    *hasTrait = ht;
  return x;
}

double Properties::getIntronFeature(int id, Traits* t, bool *hasTrait){

  bool ht = true;
  double x = 0;
  switch(id){
  case 20: x = 1;                           // intercept
    break;
  case 21: x = log( t->getLength() );       // log exon length 
    break;
  case 22:                                  // posterior probability
    if(t->hasPostProb())          
      x = t->getPostProb(); 
    else
      ht = false;
    break;
  case 23:                                  // mean base probability
    if(t->hasMeanBaseProb())
      x = t->getMeanBaseProb();
    else
      ht = false;  
    break;
  default:
    throw ProjectError("feature number " + to_string(id) + " does not exist!");
    break;
  }
 if(hasTrait)
    *hasTrait = ht;
  return x;
}


static char* get_self(void){
    char *path = NULL;
    size_t allocated = 256;

    while (1) {
	ssize_t pos = 0;
	if (!(path = (char*) malloc(allocated)))
	    abort();
	if ((pos = readlink( "/proc/self/exe", path, allocated - 1 )) != -1)
	    path[pos] = '\0';
	else {
	    free(path);
	    return NULL;
	}
	if (pos < allocated - 1)
	    break;
	free(path);
	allocated *= 2;
    }
    return path;
}


string findLocationOfSelfBinary(){
    string self;

#ifdef __APPLE__
    char path[16384];
    uint32_t size = static_cast<uint32_t>(sizeof(path));
    if (_NSGetExecutablePath(path, &size) == 0){
	self = path;
    } else {
	throw ProjectError("Could not determine binary path. Buffer too small?");
	// need to program workaround with new/free.\n";
    }
#else // LINUX
    char *path = get_self();
    ssize_t pos = 0;
    if (path)
	pos = strlen(path);
    if (pos > 0){
	self = string(path);
	free(path);
	pos = self.find_last_of("/");
	if (pos>0)
	    pos = self.find_last_of("/", pos-1);
	if (pos >= 0)
	    self.resize(pos);
	self += "/config";
    } else {
	free(path);
	throw ProjectError("/proc/self/exe not found.\nPlease specify environment variable or parameter " CFGPATH_KEY ".");
    }
#endif // WINDOWS not supported
    return self;
}

