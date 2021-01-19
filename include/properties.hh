/*
 * properties.hh
 *
 * License: Artistic License, see file LICENSE.TXT or 
 *          https://opensource.org/licenses/artistic-license-1.0
 * 
  * Description: Declaration of the Properties class.
 */

#ifndef _PROPERTIES_HH
#define _PROPERTIES_HH

// project includes
#include "types.hh"

// standard C/C++ includes
#include <map>

#ifdef __APPLE__
#include <mach-o/dyld.h>   // for _NSGetExecutablePath
#endif

#ifdef TESTING
#define NUMPARNAMES 279   // +testMode, workingDir command line options
#else
#define NUMPARNAMES 277
#endif

#define GENEMODEL_KEY "genemodel"
#define NONCODING_KEY "nc"
#define SINGLESTRAND_KEY "singlestrand"
#define SPECIES_KEY "species"
#define CFGPATH_KEY "AUGUSTUS_CONFIG_PATH"
#define INPUTFILE_KEY "queryfile"
#define SPECIESDIR_KEY "speciesdir"
#define UTR_KEY "UTR"
#define TRANSFILE_KEY "transfile"
#define EXTRINSIC_KEY "extrinsic"
#define EXTRFILE_KEY "extrinsicCfgFile"
#define EXTERNAL_KEY "optCfgFile"
#define HINTSFILE_KEY "hintsfile"
#define ALN_KEY "alnfile"
#define TREE_KEY "treefile"
#define DB_KEY "dbaccess"
#define SEQ_KEY "speciesfilenames"
#define CODONALN_KEY "codonAlignmentFile"
#define REF_EXON_KEY "referenceFile"

#define OVLPLENFILE "ovlp_len.pbl"

/**
 * @brief The base exception class for Properties
 *
 * @author Emmanouil Stafilarakis
 */
struct PropertiesError : public ProjectError{
    /**
     * Constructor
     *
     * @param msg The describing error message.
     */
    PropertiesError( string msg ) : ProjectError( msg ) { }
};

/**
 * @author Emmanouil Stafilarakis
 */
struct KeyNotFoundError : public PropertiesError{
    KeyNotFoundError( string key ) :
	PropertiesError("Properties::getProperty(): no such key \"" + key + "\".") 
    {}
};

/**
 * @author Emmanouil Stafilarakis
 */
struct ValueFormatError : public PropertiesError{
    ValueFormatError( string key, string value, string type ) :
	PropertiesError( "Properties::getProperty(): cannot convert value \"" 
			 + value + "\" of key \"" + key + "\" into " + type + "." )
    {}
};

/**
 * @author Emmanouil Stafilarakis
 */
struct HelpException {
    string message;
    HelpException(string s) : message(s) {}
}; // no error

/**
 * @brief The Properties class.
 *
 * @details This class is provided to pass several properties to all
 *      objects in the running programm. <br>
 *
 *      The properties can be given to the Properties object by
 *      the programm with the <i>addProperty</i> method. <br>
 *
 *      The Properties object can also parse a file with the
 *      <i>readFile</i> method and also parse the commandline arguments
 *      with the <i>init</i> method. <br>
 *
 *      If an adding property already exists, the old value will be
 *      overwriten with the new value. Therefore, if the commandline
 *      arguments should be prefered, the method <i>init</i> should be
 *      called after the method <i>readFile</i>! <br>
 *
 *      Example:
 *      <pre><pre style="background-color:white; color:blue; padding:10px;">
        int main( int argc, char* argv[] ){
            string     configfile;
            Properties props;
            ....
            try{
                props.readFile( configfile ); // first parse the config file
                props.init( argc, argv );     // then the commandline options
                ....
            }
            catch( PropertiesError& err ){
                ....
            }
        </pre></pre>
 *
 *      A class can get the properties with the <i>get<b>TYPE</b>Property</i>
 *      methods, where TYPE is one of Int, Double, Bool or String.  <br>
 *
 *      The format of the properties file is line based. Each line contains
 *      the property name and the property value, separated by whitespaces.
 *      Comments begin with a '#' and end an the end of line. <br>
 *
 *      Example: Listing of a properties file
 *      <pre><pre style="background-color:white; color:blue; padding:10px;">
        # For the GBProcessor object
        #---------------------------
        /GBProcessor/file               fly.txt     # The db file.

        # For the EHMMTraining class
        #----------------------------
        /EHMMTraining/statecount        3
        /EHMMTraining/state00           igenicmodel.so
        /EHMMTraining/state01           exonmodel.so
        /EHMMTraining/state02           intronmodel.so

        # For the IGenicModel class
        #---------------------------
        /IGenicModel/outfile            igenic_out.pbl
        /IGenicModel/c                  1
        /IGenicModel/k                  5

        # For the ExonModel class
        #---------------------------
        /ExonModel/outfile              exon_out.pbl
        /ExonModel/c                    1
        /ExonModel/k                    5

        # For the IntronModel class
        #---------------------------
        /IntronModel/outfile            intron_out.pbl
        /IntronModel/c                  1
        /IntronModel/d                  350
        /IntronModel/e                  20
        /IntronModel/k                  5
        /IntronModel/h                  5
        </pre></pre>
 *
 * @author Emmanouil Stafilarakis
 */
class Properties{
    public:
        /**
         *
         *
         * @author Emmanouil Stafilarakis
         * @version 0.1
         */

        /**
         * Read all properties from the given file.
         *
         * @param   file The file with the programm properties
         * @exception PropertiesError
         */
        static void        readFile    ( string file );
        /**
         * @memo    Parse the commandline arguments.
         *
         * @doc     The format of the commandline arguments should be
         *          <b>--/CLASS/PROPERTY=VALUE</b>.<br>So, the class CLASS
         *          can set the value VALUE to its property PROPERTY.
         *
         * @param   argc    The count of the arguments.
         * @param   argv    A NULL terminated field with the arguments.
         */
        static void        init        ( int argc, char* argv[] );
        /**
         *
         */
        static Boolean hasProperty( string name);
    
        /**
         * @doc     Add a new property to the Properties object.
         *
         * @param   name    The name of the property.
         * @param   value   The string value of the property.
         */
        static void        addProperty         ( string name, string value );
        /**
         *
         */
        static Integer     getIntProperty      ( string name );
 	
        /**
         *
         */
        static Double      getDoubleProperty   ( string name );
        /**
         *
         */
        static double      getdoubleProperty   ( string name );
        /**
         *
         */
        static Boolean     getBoolProperty     ( string name );
       /**
         *
         */
        static const char* getProperty   ( string name );
        /**
         *
         */
        static const char* getProperty(string name, int index);
        /**
         *
         */

    
    static string getConfigFilename(string name) {
	if (name != "")
	    return properties[CFGPATH_KEY] + getProperty(name);
	else
	    return properties[CFGPATH_KEY];
    }

    // assign value only if key is present
    static bool assignProperty(string name, Integer& target) {
	try {
	    target = getIntProperty(name);
 	} catch (KeyNotFoundError &e) { return false; }
        return true;
    }
    static bool assignProperty(string name, unsigned& target) {
	try {
	    target = (unsigned) getIntProperty(name);
 	} catch (KeyNotFoundError &e) { return false; }
        return true;
    }
    static bool assignProperty(string name, Double& target) {
	try {
	    target = getDoubleProperty(name);
 	} catch (KeyNotFoundError &e) { return false; }
        return true;
    }
    static bool assignProperty(string name, double& target) {
	try {
	    target = getdoubleProperty(name);
	} catch (KeyNotFoundError &e) { return false; }
        return true;
    }
    static bool assignProperty(string name, Boolean& target) {
	try {
	    target = getBoolProperty(name);
	} catch (KeyNotFoundError &e) { return false; }
        return true;
    }
    static bool assignProperty(string name, string& target) {
	try {
	    target = getProperty(name);
	} catch (KeyNotFoundError &e) { return false; }
        return true;
    }
    static bool assignProperty(string name, const char*& target) {
	try {
	    target = getProperty(name);
	} catch (KeyNotFoundError &e) { return false; }
        return true;
    }
    

    private:
	Properties() {}  // do not construct objects!
        /**
         *
         */
	static void            readLine    (istream& strm );
    private:
	static map<string, string> properties;
	static const char* parameternames[NUMPARNAMES];
};

/*
 * Find the path to the executable, e.g. /opt/augustus/bin/augustus, /home/mario/augustus/bin/augustus
 * code contributed by Bastien Chevreux
 */
string findLocationOfSelfBinary();

#endif  //  _PROPERTIES_HH
