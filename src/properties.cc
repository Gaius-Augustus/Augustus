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
json Properties::allowedParameters;


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
	if (argstr == "--help")
	    throw HelpException(HELP_USAGE);
	if (argstr.substr(0,2) == "--") {
	    string::size_type pos = argstr.find('=');
	    name = argstr.substr(0,pos).substr(2);
	    if( name == CFGPATH_KEY) {
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

    // read in the other parameters, command line parameters have higher priority
    // than those in the config files
    string paramsFilePath = configPath + "parameters/aug_cmdln_parameters.json";
    ifstream ifile;
    ifile.open(paramsFilePath);
    if (!ifile)
        throw ProjectError("Could not locate command line parameters file: " + paramsFilePath + ".");
    std::ifstream i(paramsFilePath);
    i >> allowedParameters;
    std::cout << "Number of Params: " << allowedParameters.size() << std::endl;

    arg = argc;
    bool miss = false;
    while (arg > 1 && !miss)
    {
        arg--;
        string argstr(argv[arg]);
        if (argstr == "--paramlist")
        {
            ostringstream strm;
            strm << "Possible parameter names are:" << endl;
            for (auto &el : allowedParameters.items())
            {
                // exclude marked parameters
                if (!el.value()["exclude_apps"].is_array() || (el.value()["exclude_apps"].is_array() && !hasValue(el.value()["exclude_apps"], "augustus")))
                {
                    strm << el.value()["name"].get<string>() << endl;
                    if (!el.value()["type"].is_null())
                        strm << "  Type: " << el.value()["type"].get<string>() << endl;
                    if (!el.value()["default_value"].is_null())
                        strm << "  Default value: " << el.value()["default_value"].dump() << endl;
                    if (el.value()["possible_values"].is_array())
                        strm << "  Possible value(s): " << el.value()["possible_values"].dump() << endl;
                    if (!el.value()["usage"].is_null())
                        strm << "  Usage: " << el.value()["usage"].get<string>() << endl;
                    if (el.value()["description"].is_object())
                    {
                        strm << "  Description: " << endl;
                        for (auto &desc : el.value()["description"].items())
                        {
                            strm << "    " << desc.key() << ": " << desc.value().get<string>() << endl;
                        }
                    }
                    else if (!el.value()["description"].is_null() && el.value()["description"].dump() != "\"\"")
                    {
                        strm << "  Description: " << el.value()["description"].get<string>() << endl;
                    }
                    strm << endl;
                }
            }
            throw HelpException(strm.str());
        }
        if (argstr.substr(0, 2) != "--")
            continue;
        argstr.erase(0, 2);
        string::size_type pos = argstr.find('=');
        name = argstr.substr(0, pos);
        if (name == CFGPATH_KEY)
            continue;
        if (pos == string::npos)
            throw PropertiesError(string("'=' missing for parameter: ") + name);

        // check that name is actually the name of a parameter
        bool found = false;
        for (auto &el : allowedParameters.items())
        {
            if (el.value()["name"] == name)
            {
                if (!el.value()["exclude_apps"].is_array()
                        || (el.value()["exclude_apps"].is_array() && !hasValue(el.value()["exclude_apps"], "augustus"))) {
                    properties[name] = argstr.substr(pos + 1);
                    found = true;
                    break;
                }
            }
        }
        if (!found)
            miss = true;
    }
    if (miss)
        throw PropertiesError("Unknown parameter: \"" + name + "\". Type \"augustus\" for help.");

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


	// TODO: check!!!
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

    // set default value if property is not set
    setDefaultValues();
}

Boolean Properties::hasProperty( string name) {
    return properties.count(name)>0;
}


void Properties::addProperty( string name, string value ) {
    properties[name] = value;
}


Integer Properties::getIntProperty( string name ) {
    Integer val;
    if (!isDefinedType("int", name)) {
        throw SpecifiedTypeError(name, "int");
    }
    istringstream strm(getProperty(name));
    if( !(strm >> val) )
        throw ValueFormatError(name, strm.str(), "Integer");
    return val;
}

Double Properties::getDoubleProperty( string name ){
    Double val;
    if (!isDefinedType("float", name)) {
        throw SpecifiedTypeError(name, "float");
    }
    istringstream strm(getProperty(name));
    if ( !(strm >> val) )
        throw ValueFormatError(name, strm.str(), "LLDouble");
    return val;
}

double Properties::getdoubleProperty( string name ) {
    double val;
    if (!isDefinedType("float", name)) {
        throw SpecifiedTypeError(name, "float");
    }
    istringstream strm(getProperty(name));
    if ( !(strm >> val) )
        throw ValueFormatError(name, strm.str(), "double");
    return val;
}

Boolean Properties::getBoolProperty( string name ) {
    string str(getProperty(name));
    if (!isDefinedType("bool", name)) {
        throw SpecifiedTypeError(name, "bool");
    }
    downcase(str);
    if( str == "true" || str=="1" || str == "t"  || str=="on" || str == "yes" || str == "y")
        return true;
    if( str == "false" || str=="0" || str=="f" || str=="off" || str == "no" || str == "n")
        return false;
    throw ValueFormatError(name, str, "Boolean");
}

const char *Properties::getProperty(string name) {
    if (!hasProperty(name))
        throw KeyNotFoundError(name);

    // check value according to json config file
    const string& value = properties[name];
    if (!isPossibleValue(value, name)) {
        throw ValueError(name, value);
    }
    return value.c_str();
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

bool Properties::hasValue(const json &list, const string value) {
    return std::find(list.begin(), list.end(), value) != list.end();
}

/*
 * Checks if the given parameter has been defined with the given type
 * in configuration file.
 * Possible types are: 'string', 'int', 'float', 'bool', 'list<string>'
 */
bool Properties::isDefinedType(const string typeName, const string paramName) {
    for (auto &el : allowedParameters.items())
    {
        if (el.value()["name"] == paramName)
        {
            if (el.value()["type"].is_null()) {
                //TODO: treat not specified types?
                //cerr << "Warning: The parameter " << paramName << " has no specified type in config file." << endl;
                return true;
            }
            if (!el.value()["type"].is_null() && el.value()["type"] != typeName)
                return false;
            if (!el.value()["type"].is_null() && el.value()["type"] == typeName)
                return true;
        }
    }
    cerr << "Warning: The parameter " << paramName << " is not specified in config file." << endl;
    return true;
}

/*
 * Checks if the passed value is included in the list of specified possible values.
 * For this the value 'possible_values' must be set in the json configuration date.
 */
bool Properties::isPossibleValue(const string value, const string paramName) {
    for (auto &el : allowedParameters.items())
    {
        if (el.value()["name"] == paramName)
        {
            if (el.value()["possible_values"].is_array()) {
                return hasValue(el.value()["possible_values"], value);
            }
        }
    }
    return true;
}

void Properties::setDefaultValues()
{
    for (auto &el : allowedParameters.items())
    {
        // exclude marked parameters
        if (!el.value()["exclude_apps"].is_array() || (el.value()["exclude_apps"].is_array() && !hasValue(el.value()["exclude_apps"], "augustus")))
        {
            if (!el.value()["default_value"].is_null() && !hasProperty(el.value()["name"].get<string>()))
            {
                properties[el.value()["name"].get<string>(), el.value()["default_value"].get<string>()];
            }
        }
    }
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
