/*
 * properties.cc
 *
 * License: Artistic License, see file LICENSE.TXT or 
 *          https://opensource.org/licenses/artistic-license-1.0
 */

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
"/CompPred/num_omega",
"/CompPred/num_features",
"/CompPred/scale_codontree",
"/CompPred/oeExtensionWidth",
"/CompPred/onlySampledECs",
"/CompPred/only_species",
"/CompPred/outdir_orthoexons",
"/CompPred/outdir",
"/CompPred/printOrthoExonAli", // codon alignment (Gerischer)
"/CompPred/printExonCandsMSA", // nucleotide alignment (Mertsch)
"/CompPred/printConservationWig",
"/CompPred/phylo_factor",
"/CompPred/phylo_model",
"/CompPred/dd_factor",
"/CompPred/dd_rounds",
"/CompPred/dd_step_rule",
"/CompPred/dualdecomp",
"/CompPred/overlapcomp",
"/CompPred/lambda",
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
"max_sgd_inter",
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
"noInFrameStop",
"noprediction",
EXTERNAL_KEY, // optCfgFile
"printEvidence",
"printHints",
"printMEA",
"printOEs",
"printSampled",
"printGeneRangesBED",
"printGeneRangesGFF",
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
"/UtrModel/verbosity"
#ifdef TESTING
, "/Testing/testMode"
, "/Testing/workingDir"
#endif
};


void Properties::readFile( string filename ){
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
    while( strm )
	readLine( strm );
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
		name == REF_EXON_KEY)
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

    // check whether multi-species mode is turned on
    Properties::assignProperty(TREE_KEY, Constant::treefile);
    Properties::assignProperty(SEQ_KEY, Constant::speciesfilenames);
    Properties::assignProperty(DB_KEY, Constant::dbaccess);
    Properties::assignProperty(ALN_KEY, Constant::alnfile);
    Properties::assignProperty(REF_EXON_KEY, Constant::referenceFile);
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
	configPath += '/';  // append slash if necessary

    // does the directory actually exist?
    struct stat buffer;
    if( stat(configPath.c_str(), &buffer) == -1 || !S_ISDIR(buffer.st_mode))
	throw ProjectError(configPath + " is not a directory. Could not locate directory " CFGPATH_KEY ".");
    properties[CFGPATH_KEY] = configPath;

    // determine species
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
	} catch (PropertiesError &e){
	    throw ProjectError("Default configuration file with parameters from logistic regression " + optCfgFile + " not found!");
	}
    }

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
    } catch (PropertiesError &e){
	throw ProjectError("Species-specific configuration files not found in " + configPath + SPECIES_SUBDIR + ". Type \"augustus --species=help\" to see available species.");
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
	    name == TREE_KEY || name == DB_KEY || name == SEQ_KEY || name == CODONALN_KEY || name == REF_EXON_KEY)
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

    Constant::softmasking_explicitly_requested = Properties::assignProperty("softmasking", Constant::softmasking);
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
