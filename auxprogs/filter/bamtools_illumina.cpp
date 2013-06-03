// ***************************************************************************
// bamtools_illumina.cpp 
// Tonatiuh Pena-Centeno
// ---------------------------------------------------------------------------
// Last modified: 2-June-2013
// ---------------------------------------------------------------------------
// Filter out paired-end reads (so far only Illumina technology is considered) 
// ***************************************************************************

#include "bamtools_illumina.h"
#include "bamtools_version.h"
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <utils/bamtools_options.h>
#include <utils/bamtools_utilities.h>
using namespace BamTools;

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <string>
using namespace std;

// --------------------------------------------------------------------------
// general ResolveTool constants
// --------------------------------------------------------------------------

static const int      NUM_MODELS = 8;
static const string   READ_GROUP_TAG = "RG";
static const double   DEFAULT_COVERAGE = 0.95;
static const double   DEFAULT_CONFIDENCE_INTERVAL = 0.9973;
static const uint16_t DEFAULT_MIN_MAPQUALITY = 1;
static const double   DEFAULT_UNUSEDMODEL_THRESHOLD = 0.1;
// --------------------------------------------------------------------------
// stats file constants
// --------------------------------------------------------------------------

// basic char/string constants
static const char COMMENT_CHAR     = '#';
static const char OPEN_BRACE_CHAR  = '[';
static const char CLOSE_BRACE_CHAR = ']';
static const char EQUAL_CHAR       = '=';
static const char TAB_CHAR         = '\t';

static const string WHITESPACE_CHARS = " \t\n";
static const string TRUE_KEYWORD     = "true";
static const string FALSE_KEYWORD    = "false";

// field counts
static const size_t NUM_OPTIONS_FIELDS    = 2;
static const size_t NUM_READGROUPS_FIELDS = 7;

// header strings
static const string INPUT_TOKEN      = "[Input]";
static const string OPTIONS_TOKEN    = "[Options]";
static const string READGROUPS_TOKEN = "[ReadGroups]";

// option keywords
static const string OPTION_COVERAGE   = "Coverage";
static const string OPTION_CONFIDENCEINTERVAL   = "ConfidenceInterval";
static const string OPTION_MINIMUMMAPQUALITY    = "MinimumMapQuality";
static const string OPTION_UNUSEDMODELTHRESHOLD = "UnusedModelThreshold";
static const string OPTION_FORCEMARKREADGROUPS  = "ForceMarkReadGroups";

// other string constants
static const string RG_FIELD_DESCRIPTION =
    "#<name> <medianFL> <minFL> <maxFL> <topModelID> <nextTopModelID> <isAmbiguous?>";

// --------------------------------------------------------------------------
// unique readname file constants
// --------------------------------------------------------------------------

static const string READNAME_FILE_SUFFIX = ".uniq_names.txt";
static const string DEFAULT_READNAME_FILE = "bt_resolve_TEMP" + READNAME_FILE_SUFFIX;

// --------------------------------------------------------------------------
// ModelType implementation

struct ModelType {

    // data members
    uint16_t ID;
    vector<int32_t> FragmentLengths;

    // ctor
    ModelType(const uint16_t id)
        : ID(id)
    {
        // preallocate space for 10K fragments per model type
        FragmentLengths.reserve(10000);
    }

    // convenience access to internal fragment lengths vector
    vector<int32_t>::iterator begin(void) { return FragmentLengths.begin(); }
    vector<int32_t>::const_iterator begin(void) const { return FragmentLengths.begin(); }
    void clear(void) { FragmentLengths.clear(); }
    vector<int32_t>::iterator end(void) { return FragmentLengths.end(); }
    vector<int32_t>::const_iterator end(void) const { return FragmentLengths.end(); }
    void push_back(const int32_t& x) { FragmentLengths.push_back(x); }
    size_t size(void) const { return FragmentLengths.size(); }

    // constants
    static const uint16_t DUMMY_ID;
};

const uint16_t ModelType::DUMMY_ID = 100;

bool operator>(const ModelType& lhs, const ModelType& rhs) {
    return lhs.size() > rhs.size();
}

uint16_t CalculateModelType(const BamAlignment& al) {

    // localize alignment's mate positions & orientations for convenience
    const int32_t m1_begin = ( al.IsFirstMate() ? al.Position : al.MatePosition );
    const int32_t m2_begin = ( al.IsFirstMate() ? al.MatePosition : al.Position );
    const bool m1_isReverseStrand = ( al.IsFirstMate() ? al.IsReverseStrand() : al.IsMateReverseStrand() );
    const bool m2_isReverseStrand = ( al.IsFirstMate() ? al.IsMateReverseStrand() : al.IsReverseStrand() );

    // determine 'model type'
    if ( m1_begin < m2_begin ) {
        if ( !m1_isReverseStrand && !m2_isReverseStrand ) return 0; // ID: 1
        if ( !m1_isReverseStrand &&  m2_isReverseStrand ) return 1; // ID: 2
        if (  m1_isReverseStrand && !m2_isReverseStrand ) return 2; // ID: 3
        if (  m1_isReverseStrand &&  m2_isReverseStrand ) return 3; // ID: 4
    } else {
        if ( !m2_isReverseStrand && !m1_isReverseStrand ) return 4; // ID: 5
        if ( !m2_isReverseStrand &&  m1_isReverseStrand ) return 5; // ID: 6
        if (  m2_isReverseStrand && !m1_isReverseStrand ) return 6; // ID: 7
        if (  m2_isReverseStrand &&  m1_isReverseStrand ) return 7; // ID: 8
    }

    // unknown model
    return ModelType::DUMMY_ID;
}

// --------------------------------------------------------------------------
// ReadGroupResolver implementation

struct ReadGroupResolver {

    // data members
    int32_t MinFragmentLength;
    int32_t MedianFragmentLength;
    int32_t MaxFragmentLength;
    uint16_t TopModelId;
    uint16_t NextTopModelId;
    bool IsAmbiguous;
    bool HasData;
    vector<ModelType> Models;
    map<string, bool> ReadNames;
    // ctor
    ReadGroupResolver(void);

    // resolving methods
    bool IsValidInsertSize(const BamAlignment& al) const;
    bool IsValidOrientation(const BamAlignment& al) const;

    // select 2 best models based on observed data
    void PrintModels(const string& readGroupName);
    void DetermineTopModels(const string& readGroupName);

    // static settings
    static double Coverage;
    static double ConfidenceInterval;
    static double UnusedModelThreshold;
    static void SetCoverage(const double& coverage);
    static void SetConfidenceInterval(const double& ci);
    static void SetUnusedModelThreshold(const double& umt);
};

double ReadGroupResolver::Coverage   = DEFAULT_COVERAGE;
double ReadGroupResolver::ConfidenceInterval   = DEFAULT_CONFIDENCE_INTERVAL;
double ReadGroupResolver::UnusedModelThreshold = DEFAULT_UNUSEDMODEL_THRESHOLD;

ReadGroupResolver::ReadGroupResolver(void)
    : MinFragmentLength(0)
    , MedianFragmentLength(0)
    , MaxFragmentLength(0)
    , TopModelId(ModelType::DUMMY_ID)
    , NextTopModelId(ModelType::DUMMY_ID)
    , IsAmbiguous(false)
    , HasData(false)
{
    // pre-allocate space for 8 models
    Models.reserve(NUM_MODELS);
    for ( uint16_t i = 0; i < NUM_MODELS; ++i )
        Models.push_back( ModelType(i+1) );
}

bool ReadGroupResolver::IsValidInsertSize(const BamAlignment& al) const {  
    const int32_t absInsertSize = abs(al.InsertSize);
    return ( absInsertSize >= MinFragmentLength &&
             absInsertSize <= MaxFragmentLength );
}

bool ReadGroupResolver::IsValidOrientation(const BamAlignment& al) const {
    const uint16_t currentModelId = CalculateModelType(al) + 1; // convert model type (array index) to ID number
    return ( currentModelId == TopModelId || currentModelId == NextTopModelId );
}

void ReadGroupResolver::PrintModels(const string& readGroupName) {

    // sort models (from most common to least common)
    sort( Models.begin(), Models.end(), std::greater<ModelType>() );
	cout << "Read group: " << readGroupName << endl;
	cout << "Found orientations:" << Models.size() << endl;
    for (int it=0;it<Models.size();it++) {
	  cout << "Config " << Models[it].ID << ": " << Models[it].size() << endl; 
	}
    const unsigned int activeModelCountSum = Models[0].size() + Models[1].size();
    const unsigned int unusedModelCountSum = Models[2].size() + Models[3].size() +
                                             Models[4].size() + Models[5].size() +
                                             Models[6].size() + Models[7].size();    
	cout << "activeModelCountSum: " << activeModelCountSum << endl;
	cout << "unusedModelCountSum: " << unusedModelCountSum << endl;
	cout << "activeModelCountSum+unusedModelCountSum = wc -l file_uniq_names.txt" << endl;
}

void ReadGroupResolver::DetermineTopModels(const string& readGroupName) {

    // sort models (from most common to least common)
    sort( Models.begin(), Models.end(), std::greater<ModelType>() );

    // store top 2 models for later
    TopModelId     = Models[0].ID;
    NextTopModelId = Models[1].ID;

    // make sure that the 2 most common models are some threshold more common
    // than the remaining models
    const unsigned int activeModelCountSum = Models[0].size() + Models[1].size();
    if ( activeModelCountSum == 0 ) return; // skip if no data in this read group
    const unsigned int unusedModelCountSum = Models[2].size() + Models[3].size() +
                                             Models[4].size() + Models[5].size() +
                                             Models[6].size() + Models[7].size();    
    const double unusedPercentage = (double)unusedModelCountSum / (double)activeModelCountSum;
    if ( unusedPercentage > UnusedModelThreshold ) {
        cerr << "WARNING: " << readGroupName << " does not have clearly defined 'top models'" << endl
             << "         The fraction of alignments in bottom 6 models (" << unusedPercentage
             << ") exceeds threshold: " << UnusedModelThreshold << endl;
        IsAmbiguous = true;
    }

    // emit a warning if the best alignment models are non-standard
    const bool isModel1Top = (TopModelId == 1) || (NextTopModelId == 1);
    const bool isModel2Top = (TopModelId == 2) || (NextTopModelId == 2);
    const bool isModel4Top = (TopModelId == 4) || (NextTopModelId == 4);
    const bool isModel5Top = (TopModelId == 5) || (NextTopModelId == 5);
    const bool isModel6Top = (TopModelId == 6) || (NextTopModelId == 6);
    const bool isModel8Top = (TopModelId == 8) || (NextTopModelId == 8);

    bool isMatePair  = ( isModel4Top && isModel5Top ? true : false );
    bool isPairedEnd = ( isModel2Top && isModel6Top ? true : false );
    bool isSolidPair = ( isModel1Top && isModel8Top ? true : false );

    if ( !isMatePair && !isPairedEnd && !isSolidPair ) {
        cerr << "WARNING: Found a non-standard alignment model configuration. " << endl
             << "         Using alignment models " << TopModelId << " & " << NextTopModelId
             << endl;
    }

    // store only the fragments from the best alignment models, then sort
    vector<int32_t> fragments;
    fragments.reserve( Models[0].size() + Models[1].size() );
    fragments.insert( fragments.end(), Models[0].begin(), Models[0].end() );
    fragments.insert( fragments.end(), Models[1].begin(), Models[1].end() );
    sort ( fragments.begin(), fragments.end() );

    // clear out Model fragment data, not needed anymore
    Models.clear();

    // skip if no fragments found for this read group
    if ( fragments.empty() ) {
        HasData = false;
        return;
    } else
        HasData = true;

    // calculate & store the min,median, & max fragment lengths
    const unsigned int numFragmentLengths = fragments.size();
    const double halfNonConfidenceInterval = (1.0 - ReadGroupResolver::ConfidenceInterval)/2.0;
    const unsigned int minIndex    = (unsigned int)(numFragmentLengths * halfNonConfidenceInterval);
    const unsigned int medianIndex = (unsigned int)(numFragmentLengths * 0.5);
    const unsigned int maxIndex    = (unsigned int)(numFragmentLengths * (1.0-halfNonConfidenceInterval));

    MinFragmentLength    = fragments[minIndex];
    MedianFragmentLength = fragments[medianIndex];
    MaxFragmentLength    = fragments[maxIndex];
}

void ReadGroupResolver::SetCoverage(const double& coverage) {
    Coverage = coverage;
}

void ReadGroupResolver::SetConfidenceInterval(const double& ci) {
    ConfidenceInterval = ci;
}

void ReadGroupResolver::SetUnusedModelThreshold(const double& umt) {
    UnusedModelThreshold = umt;
}

// --------------------------------------------------------------------------
// IlluminaSettings implementation

struct IlluminaTool::IlluminaSettings {

    // modes
    bool IsClean;
    bool IsMakeStats;

    // bool IsMarkPairs;
    // bool IsTwoPass;

    // I/O flags
    bool HasInputBamFile;
    bool HasOutputBamFile;
    bool HasStatsFile;
    bool IsForceCompression;

    // resolve option flags
  	bool HasCoverage;
    bool HasConfidenceInterval;
    bool HasForceMarkReadGroups;
    bool HasMinimumMapQuality;
    bool HasUnusedModelThreshold;

    // I/O filenames
    string InputBamFilename;
    string OutputBamFilename;
    string StatsFilename;
    string ReadNamesFilename; //  ** N.B. - Only used internally, not set from cmdline **

    // resolve options
 	double	 Coverage;
    double   ConfidenceInterval;
    uint16_t MinimumMapQuality;
    double   UnusedModelThreshold;

    // constructor
    IlluminaSettings(void)
	    : IsMakeStats(false)
		// ///////////////////////////
		, IsClean(false)
		// ///////////////////////////

        // , IsMarkPairs(false)
        // , IsTwoPass(false)

        , HasInputBamFile(false)
        , HasOutputBamFile(false)
        , HasStatsFile(false)
        , IsForceCompression(false)
        , HasCoverage(false)
        , HasConfidenceInterval(false)
        , HasForceMarkReadGroups(false)
        , HasMinimumMapQuality(false)
        , HasUnusedModelThreshold(false)
        , InputBamFilename(Options::StandardIn())
        , OutputBamFilename(Options::StandardOut())
        , StatsFilename("")
        , ReadNamesFilename(DEFAULT_READNAME_FILE)
        , Coverage(DEFAULT_COVERAGE)
        , ConfidenceInterval(DEFAULT_CONFIDENCE_INTERVAL)
        , MinimumMapQuality(DEFAULT_MIN_MAPQUALITY)
        , UnusedModelThreshold(DEFAULT_UNUSEDMODEL_THRESHOLD)
    { }
};

// --------------------------------------------------------------------------
// ReadNamesFileReader implementation

struct IlluminaTool::ReadNamesFileReader {

    // ctor & dtor
    ReadNamesFileReader(void) { }
    ~ReadNamesFileReader(void) { Close(); }

    // main reader interface
    public:
        void Close(void);
        bool Open(const string& filename);
        bool Read(map<string, ReadGroupResolver>& readGroups);

    // data members
    private:
        ifstream m_stream;
};

void IlluminaTool::ReadNamesFileReader::Close(void) {
    if ( m_stream.is_open() )
        m_stream.close();
}

bool IlluminaTool::ReadNamesFileReader::Open(const string& filename) {

    // make sure stream is fresh
    Close();

    // attempt to open filename, return status
    m_stream.open(filename.c_str(), ifstream::in);
    return m_stream.good();
}

bool IlluminaTool::ReadNamesFileReader::Read(map<string, ReadGroupResolver>& readGroups) {

    // up-front sanity check
    if ( !m_stream.is_open() ) return false;

    // parse read names file
    string line;
    vector<string> fields;
    map<string, ReadGroupResolver>::iterator rgIter;
    map<string, ReadGroupResolver>::iterator rgEnd = readGroups.end();
    while ( getline(m_stream, line) ) {

        // skip if empty line
        if ( line.empty() ) continue;

        // split line on '\t'
        fields = Utilities::Split(line, TAB_CHAR);
        if ( fields.size() != 2 ) continue;

        // look up resolver for read group
        rgIter = readGroups.find( fields[0] );
        if ( rgIter == rgEnd ) return false;
        ReadGroupResolver& resolver = (*rgIter).second;

        // store read name with resolver
        resolver.ReadNames.insert( make_pair<string,bool>(fields[1], true) ) ;
    }

    // if here, return success
    return true;
}

// --------------------------------------------------------------------------
// ReadNamesFileWriter implementation

struct IlluminaTool::ReadNamesFileWriter {

    // ctor & dtor
    ReadNamesFileWriter(void) { }
    ~ReadNamesFileWriter(void) { Close(); }

    // main reader interface
    public:
        void Close(void);
        bool Open(const string& filename);
        void Write(const string& readGroupName, const string& readName);

    // data members
    private:
        ofstream m_stream;
};

void IlluminaTool::ReadNamesFileWriter::Close(void) {
    if ( m_stream.is_open() )
        m_stream.close();
}

bool IlluminaTool::ReadNamesFileWriter::Open(const string& filename) {

    // make sure stream is fresh
    Close();

    // attempt to open filename, return status
    m_stream.open(filename.c_str(), ofstream::out);
    return m_stream.good();
}

void IlluminaTool::ReadNamesFileWriter::Write(const string& readGroupName,
                                             const string& readName)
{
    m_stream << readGroupName << TAB_CHAR << readName << endl;
}

// --------------------------------------------------------------------------
// StatsFileReader implementation

struct IlluminaTool::StatsFileReader {

    // ctor & dtor
    public:
        StatsFileReader(void) { }
        ~StatsFileReader(void) { Close(); }

    // main reader interface
    public:
        void Close(void);
        bool Open(const string& filename);
        bool Read(IlluminaTool::IlluminaSettings* settings,
                  map<string, ReadGroupResolver>& readGroups);

    // internal methods
    private:
        bool IsComment(const string& line) const;
        bool IsWhitespace(const string& line) const;
        bool ParseInputLine(const string& line);
        bool ParseOptionLine(const string& line, IlluminaTool::IlluminaSettings* settings);
        bool ParseReadGroupLine(const string& line, map<string, ReadGroupResolver>& readGroups);
        string SkipCommentsAndWhitespace(void);

    // data members
    private:
        ifstream m_stream;

        enum State { None = 0
                   , InInput
                   , InOptions
                   , InReadGroups };
};

void IlluminaTool::StatsFileReader::Close(void) {
    if ( m_stream.is_open() )
        m_stream.close();
}

bool IlluminaTool::StatsFileReader::IsComment(const string& line) const {
    assert( !line.empty() );
    return ( line.at(0) == COMMENT_CHAR );
}

bool IlluminaTool::StatsFileReader::IsWhitespace(const string& line) const {
    if ( line.empty() )
        return true;
    return ( isspace(line.at(0)) );
}

bool IlluminaTool::StatsFileReader::Open(const string& filename) {

    // make sure stream is fresh
    Close();

    // attempt to open filename, return status
    m_stream.open(filename.c_str(), ifstream::in);
    return m_stream.good();
}

bool IlluminaTool::StatsFileReader::ParseInputLine(const string& /*line*/) {
    // input lines are ignored (for now at least), tool will use input from command line
    return true;
}

bool IlluminaTool::StatsFileReader::ParseOptionLine(const string& line,
                                                   IlluminaTool::IlluminaSettings* settings)
{
    // split line into option, value
    vector<string> fields = Utilities::Split(line, EQUAL_CHAR);
    if ( fields.size() != NUM_OPTIONS_FIELDS )
        return false;
    const string& option = fields.at(0);
    stringstream value(fields.at(1));

    // -----------------------------------
    // handle option based on keyword

    // Coverage
    if ( option == OPTION_COVERAGE ) {
        value >> settings->Coverage;
        settings->HasCoverage = true;
        return true;
    }

    // ConfidenceInterval
    if ( option == OPTION_CONFIDENCEINTERVAL ) {
        value >> settings->ConfidenceInterval;
        settings->HasConfidenceInterval = true;
        return true;
    }

    // ForceMarkReadGroups
    if ( option == OPTION_FORCEMARKREADGROUPS ) {
        value >> settings->HasForceMarkReadGroups;
        return true;
    }

    // MinimumMapQuality
    if ( option == OPTION_MINIMUMMAPQUALITY ) {
        value >> settings->MinimumMapQuality;
        settings->HasMinimumMapQuality = true;
        return true;
    }

    // UnusedModelThreshold
    if ( option == OPTION_UNUSEDMODELTHRESHOLD ) {
        value >> settings->UnusedModelThreshold;
        settings->HasUnusedModelThreshold = true;
        return true;
    }

    // otherwise unknown option
    cerr << "filter ERROR - unrecognized option: " << option << " in stats file" << endl;
    return false;
}

bool IlluminaTool::StatsFileReader::ParseReadGroupLine(const string& line,
                                                      map<string, ReadGroupResolver>& readGroups)
{
    // split read group data in to fields
    vector<string> fields = Utilities::Split(line, WHITESPACE_CHARS);
    if ( fields.size() != NUM_READGROUPS_FIELDS ) return false;

    // retrieve RG name
    const string& name = fields.at(0);

    // populate RG's 'resolver' data
    ReadGroupResolver resolver;

    stringstream dataStream;
    dataStream.str(fields.at(1));
    dataStream >> resolver.MedianFragmentLength;
    dataStream.clear();

    dataStream.str(fields.at(2));
    dataStream >> resolver.MinFragmentLength;
    dataStream.clear();

    dataStream.str(fields.at(3));
    dataStream >> resolver.MaxFragmentLength;
    dataStream.clear();

    dataStream.str(fields.at(4));
    dataStream >> resolver.TopModelId;
    dataStream.clear();

    dataStream.str(fields.at(5));
    dataStream >> resolver.NextTopModelId;
    dataStream.clear();

    resolver.IsAmbiguous = ( fields.at(6) == TRUE_KEYWORD );

    // store RG entry and return success
    readGroups.insert( make_pair<string, ReadGroupResolver>(name, resolver) );
    return true;
}

bool IlluminaTool::StatsFileReader::Read(IlluminaTool::IlluminaSettings* settings,
                                        map<string, ReadGroupResolver>& readGroups)
{
    // up-front sanity checks
    if ( !m_stream.is_open() || settings == 0 )
        return false;

    // clear out read group data
    readGroups.clear();

    // initialize state
    State currentState = StatsFileReader::None;

    // read stats file
    string line = SkipCommentsAndWhitespace();
    while ( !line.empty() ) {

        bool foundError = false;

        // switch state on keyword found
        if ( Utilities::StartsWith(line, INPUT_TOKEN) )
            currentState = StatsFileReader::InInput;
        else if ( Utilities::StartsWith(line, OPTIONS_TOKEN) )
            currentState = StatsFileReader::InOptions;
        else if ( Utilities::StartsWith(line, READGROUPS_TOKEN) )
            currentState = StatsFileReader::InReadGroups;

        // otherwise parse data line, depending on state
        else {
            if ( currentState == StatsFileReader::InInput )
                foundError = !ParseInputLine(line);
            else if ( currentState == StatsFileReader::InOptions )
                foundError = !ParseOptionLine(line, settings);
            else if ( currentState == StatsFileReader::InReadGroups )
                foundError = !ParseReadGroupLine(line, readGroups);
            else
                foundError = true;
        }

        // break out if error found
        if ( foundError )
            return false;

        // get next line
        line = SkipCommentsAndWhitespace();
    }

    // if here, return success
    return true;
}

string IlluminaTool::StatsFileReader::SkipCommentsAndWhitespace(void) {
    string line;
    do {
        if ( m_stream.eof() )
            return string();
        getline(m_stream, line);
    } while ( IsWhitespace(line) || IsComment(line) );
    return line;
}

// --------------------------------------------------------------------------
// StatsFileReader implementation

struct IlluminaTool::StatsFileWriter {

    // ctor & dtor
    public:
        StatsFileWriter(void) { }
        ~StatsFileWriter(void) { Close(); }

    // main reader interface
    public:
        void Close(void);
        bool Open(const string& filename);
        bool Write(IlluminaTool::IlluminaSettings* settings,
                   const map<string, ReadGroupResolver>& readGroups);

    // internal methods
    private:
        void WriteHeader(void);
        void WriteInput(IlluminaTool::IlluminaSettings* settings);
        void WriteOptions(IlluminaTool::IlluminaSettings* settings);
        void WriteReadGroups(const map<string, ReadGroupResolver>& readGroups);

    // data members
    private:
        ofstream m_stream;
};

void IlluminaTool::StatsFileWriter::Close(void) {
    if ( m_stream.is_open() )
        m_stream.close();
}

bool IlluminaTool::StatsFileWriter::Open(const string& filename) {

    // make sure stream is fresh
    Close();

    // attempt to open filename, return status
    m_stream.open(filename.c_str(), ofstream::out);
    return m_stream.good();
}

bool IlluminaTool::StatsFileWriter::Write(IlluminaTool::IlluminaSettings* settings,
                                         const map<string, ReadGroupResolver>& readGroups)
{
    // return failure if file not open
    if ( !m_stream.is_open() )
        return false;

    // write stats file elements
    WriteHeader();
    WriteInput(settings);
    WriteOptions(settings);
    WriteReadGroups(readGroups);

    // return success
    return true;
}

void IlluminaTool::StatsFileWriter::WriteHeader(void) {

    // stringify current bamtools version
    stringstream versionStream("");
    versionStream << "v"
                  << BAMTOOLS_VERSION_MAJOR << "."
                  << BAMTOOLS_VERSION_MINOR << "."
                  << BAMTOOLS_VERSION_BUILD;

    // # bamtools resolve (vX.Y.Z)
    // \n

    m_stream << COMMENT_CHAR << " bamtools resolve (" << versionStream.str() << ")" << endl
             << endl;
}

void IlluminaTool::StatsFileWriter::WriteInput(IlluminaTool::IlluminaSettings* settings) {

    // [Input]
    // filename
    // \n

    m_stream << INPUT_TOKEN << endl
             << settings->InputBamFilename << endl
             << endl;
}

void IlluminaTool::StatsFileWriter::WriteOptions(IlluminaTool::IlluminaSettings* settings) {

    // [Options]
    // ConfidenceInterval=<double>
    // ForceMarkReadGroups=<true|false>
    // MinimumMapQuality=<uint16_t>
    // UnusedModelThreshold=<double>
    // \n

    m_stream << OPTIONS_TOKEN << endl
             << OPTION_COVERAGE   << EQUAL_CHAR << settings->Coverage << endl
             << OPTION_CONFIDENCEINTERVAL   << EQUAL_CHAR << settings->ConfidenceInterval << endl
             << OPTION_FORCEMARKREADGROUPS  << EQUAL_CHAR << boolalpha << settings->HasForceMarkReadGroups << endl
             << OPTION_MINIMUMMAPQUALITY    << EQUAL_CHAR << settings->MinimumMapQuality << endl
             << OPTION_UNUSEDMODELTHRESHOLD << EQUAL_CHAR << settings->UnusedModelThreshold << endl
             << endl;
}

void IlluminaTool::StatsFileWriter::WriteReadGroups(const map<string, ReadGroupResolver>& readGroups) {

    // [ReadGroups]
    // #<name> <medianFL> <minFL> <maxFL> <topModelID> <nextTopModelID> <isAmbiguous?>
    m_stream << READGROUPS_TOKEN << endl
             << RG_FIELD_DESCRIPTION << endl;

    // iterate over read groups
    map<string, ReadGroupResolver>::const_iterator rgIter = readGroups.begin();
    map<string, ReadGroupResolver>::const_iterator rgEnd  = readGroups.end();
    for ( ; rgIter != rgEnd; ++rgIter ) {
        const string& name = (*rgIter).first;
        const ReadGroupResolver& resolver = (*rgIter).second;

        // skip if read group has no data
        if ( !resolver.HasData )
            continue;

        // write read group data
        m_stream << name << TAB_CHAR
                 << resolver.MedianFragmentLength << TAB_CHAR
                 << resolver.MinFragmentLength << TAB_CHAR
                 << resolver.MaxFragmentLength << TAB_CHAR
                 << resolver.TopModelId << TAB_CHAR
                 << resolver.NextTopModelId << TAB_CHAR
                 << boolalpha << resolver.IsAmbiguous
                 << endl;
    }

    // extra newline at end
    m_stream << endl;
}

// --------------------------------------------------------------------------
// IlluminaToolPrivate implementation

struct IlluminaTool::IlluminaToolPrivate {

    // ctor & dtor
    public:
        IlluminaToolPrivate(IlluminaTool::IlluminaSettings* settings)
            : m_settings(settings)
        { }
        ~IlluminaToolPrivate(void) { }

    // 'public' interface
    public:
        bool Run(void);

    // internal methods
    private:
        bool CheckSettings(vector<string>& errors);
		///////////////////////////
        bool Clean(void);
		///////////////////////////
        bool MakeStats(void);
        void ParseHeader(const SamHeader& header);
        bool ReadStatsFile(void);
        void ResolveAlignment(BamAlignment& al);
        bool ResolvePairs(void);
        bool WriteStatsFile(void);

    // data members
    private:
        IlluminaTool::IlluminaSettings* m_settings;
        map<string, ReadGroupResolver> m_readGroups;
};

bool IlluminaTool::IlluminaToolPrivate::CheckSettings(vector<string>& errors) {

    // ensure clean slate
    errors.clear();

    // if MakeStats mode
    if ( m_settings->IsMakeStats ) {

        // error if output BAM options supplied
        if ( m_settings->HasOutputBamFile )
            errors.push_back("Cannot use -out (output BAM file) in -makeStats mode.");

        // make sure required stats file supplied
        if ( !m_settings->HasStatsFile )
            errors.push_back("Ouptut stats filename required for -makeStats mode. Please specify one using -stats option.");

    }

	/////////////////////////////////////////////////////////////////////////
    // if Clean mode
    else if ( m_settings->IsClean ) {
        // // error if output BAM options supplied
        // if ( m_settings->HasOutputBamFile )
        //     errors.push_back("Cannot use -out (output BAM file) in -clean mode.");

        // make sure required stats file supplied
        if ( !m_settings->HasStatsFile )
            errors.push_back("Ouptut stats filename required for -clean mode. Please specify one using -stats option.");
    }
	/////////////////////////////////////////////////////////////////////////

    // no mode selected
    else
        errors.push_back("No <mode> specified. Please select -clean. Type filter illumina --help for more info.");

    // boundary checks on values
    if ( m_settings->HasCoverage ) {
        if ( m_settings->Coverage < 0.0 || m_settings->Coverage > 1.0 )
            errors.push_back("Invalid Coverage. Must be between 0 and 1");
    }
    if ( m_settings->HasConfidenceInterval ) {
        if ( m_settings->ConfidenceInterval < 0.0 || m_settings->ConfidenceInterval > 1.0 )
            errors.push_back("Invalid confidence interval. Must be between 0 and 1");
    }
    if ( m_settings->HasMinimumMapQuality ) {
        if ( m_settings->MinimumMapQuality >= 256 )
            errors.push_back("Invalid minimum map quality. Must be between 0 and 255");
    }
    if ( m_settings->HasUnusedModelThreshold ) {
        if ( m_settings->UnusedModelThreshold < 0.0 || m_settings->UnusedModelThreshold > 1.0 )
            errors.push_back("Invalid unused model threshold. Must be between 0 and 1");
    }

    // return success if no errors found
    return ( errors.empty() );
}

bool IlluminaTool::IlluminaToolPrivate::MakeStats(void) {

    // pull resolver settings from command-line settings
    ReadGroupResolver::SetConfidenceInterval(m_settings->ConfidenceInterval);
    ReadGroupResolver::SetUnusedModelThreshold(m_settings->UnusedModelThreshold);

    // open our BAM reader
    BamReader bamReader;
    if ( !bamReader.Open(m_settings->InputBamFilename) ) {
        cerr << "filter ERROR: could not open input BAM file: "
             << m_settings->InputBamFilename << endl;
        return false;
    }

    // retrieve header & parse for read groups
    const SamHeader& header = bamReader.GetHeader();
    ParseHeader(header);

    // open ReadNamesFileWriter
    IlluminaTool::ReadNamesFileWriter readNamesWriter;
    if ( !readNamesWriter.Open(m_settings->ReadNamesFilename) ) {
        cerr << "filter ERROR: could not open (temp) output read names file: "
             << m_settings->ReadNamesFilename << endl;
        bamReader.Close();
        return false;
    }

    // read through BAM file
    BamAlignment al;
    string readGroup("");
    map<string, ReadGroupResolver>::iterator rgIter;
    map<string, bool>::iterator readNameIter;
    while ( bamReader.GetNextAlignmentCore(al) ) {

        // skip if alignment is not paired, mapped, nor mate is mapped
        if ( !al.IsPaired() || !al.IsMapped() || !al.IsMateMapped() )
            continue;

        // skip if alignment & mate not on same reference sequence
        if ( al.RefID != al.MateRefID ) continue;

        // flesh out the char data, so we can retrieve its read group ID
        al.BuildCharData();

        // get read group from alignment (OK if empty)
        readGroup.clear();
        al.GetTag(READ_GROUP_TAG, readGroup);

        // look up resolver for read group
        rgIter = m_readGroups.find(readGroup);
        if ( rgIter == m_readGroups.end() )  {
            cerr << "bamtools ERROR - unable to calculate stats, unknown read group encountered: "
                 << readGroup << endl;
            bamReader.Close();
            return false;
        }
        ReadGroupResolver& resolver = (*rgIter).second;

        // determine unique-ness of current alignment
        const bool isCurrentMateUnique = ( al.MapQuality >= m_settings->MinimumMapQuality );

        // look up read name
        readNameIter = resolver.ReadNames.find(al.Name);

        // if read name found (current alignment's mate already parsed)
        if ( readNameIter != resolver.ReadNames.end() ) {

            // if both unique mates are unique, store read name & insert size for later
            const bool isStoredMateUnique  = (*readNameIter).second;
            if ( isCurrentMateUnique && isStoredMateUnique ) {

                // save read name in temp file as candidates for later pair marking
                readNamesWriter.Write(readGroup, al.Name);

                // determine model type & store fragment length for stats calculation
                const uint16_t currentModelType = CalculateModelType(al);
                assert( currentModelType != ModelType::DUMMY_ID );
                resolver.Models[currentModelType].push_back( abs(al.InsertSize) );
            }

            // unique or not, remove read name from map
            resolver.ReadNames.erase(readNameIter);
        }

        // if read name not found, store new entry
        else resolver.ReadNames.insert( make_pair<string, bool>(al.Name, isCurrentMateUnique) );
    }

    // close files
    readNamesWriter.Close();
    bamReader.Close();

    // iterate back through read groups
    map<string, ReadGroupResolver>::iterator rgEnd  = m_readGroups.end();
    for ( rgIter = m_readGroups.begin(); rgIter != rgEnd; ++rgIter ) {
        const string& name = (*rgIter).first;
        ReadGroupResolver& resolver = (*rgIter).second;

        // calculate acceptable orientation & insert sizes for this read group
        resolver.DetermineTopModels(name);
		
        // clear out left over read names
        // (these have mates that did not pass filters or were already removed as non-unique)
        resolver.ReadNames.clear();
    }

    // if we get here, return success
    return true;
}

// --------------------------------------------------------------------------
// Implementation of cleaning filter

// In SAM, sumDandIOperations = $qBaseInsert+$tBaseInsert
uint32_t sumDandIOperations(vector<CigarOp> cigar, string printFlag)
{
	int cigarSize = cigar.size();
	uint32_t sumDandI = 0;
 	uint32_t cigarLength;
  	char cigarType;

	// Scanning through all CIGAR operations
	for (int it=0; it<cigarSize; it++)
		{
		// Length is number of bases
		cigarLength = cigar.at(it).Length;	
		// Type means operations, i.e. MINDSHP 
		cigarType = cigar.at(it).Type;  

		if (cigarType == 'D' || cigarType == 'I')
			{
		  	// Sum and print operations
			sumDandI = cigarLength + sumDandI;	
			if (!printFlag.compare("print"))
			  {cout << cigarLength << ";" << cigarType << ",";}
			}	
		} // end for

	return sumDandI;
}


uint32_t sumMandIOperations(vector<CigarOp> cigar, string printFlag)
{
	int cigarSize = cigar.size();
	uint32_t sumMandI = 0;
 	uint32_t cigarLength;
  	char cigarType;

	// Scanning through all CIGAR operations
	for (int it=0; it<cigarSize; it++)
		{
		// Length is number of bases
		cigarLength = cigar.at(it).Length;	
		// Type means operations, i.e. MINDSHP 
		cigarType = cigar.at(it).Type;  

		if (cigarType == 'M' || cigarType == 'I')
			{
		  	// Sum and print operations
			sumMandI = cigarLength + sumMandI;	
			if (!printFlag.compare("print"))
			  {cout << cigarLength << ";" << cigarType << ",";}
			}	
		} // end for

	return sumMandI;
}


/////////////////////////////////////////////////////////////////////////////////

bool IlluminaTool::IlluminaToolPrivate::Clean(void) {

    // pull resolver settings from command-line settings
    ReadGroupResolver::SetCoverage(m_settings->Coverage);
    ReadGroupResolver::SetConfidenceInterval(m_settings->ConfidenceInterval);
    ReadGroupResolver::SetUnusedModelThreshold(m_settings->UnusedModelThreshold);

    // open our BAM reader
    BamReader bamReader;
    if ( !bamReader.Open(m_settings->InputBamFilename) ) {
        cerr << "filter ERROR: could not open input BAM file: "
             << m_settings->InputBamFilename << endl;
        return false;
    }

    // retrieve header & parse for read groups
    const SamHeader& header = bamReader.GetHeader();
    const RefVector& references = bamReader.GetReferenceData();
    ParseHeader(header);

    // open BamWriter
    BamWriter writer;
    if ( !writer.Open(m_settings->OutputBamFilename, header, references) ) {
        cerr << "filter ERROR: could not open "
             << m_settings->OutputBamFilename << " for writing." << endl;
        bamReader.Close();
        return false;
    }


    // open ReadNamesFileWriter
    IlluminaTool::ReadNamesFileWriter readNamesWriter;
    if ( !readNamesWriter.Open(m_settings->ReadNamesFilename) ) {
        cerr << "filter ERROR: could not open (temp) output read names file: "
             << m_settings->ReadNamesFilename << endl;
        bamReader.Close();
        return false;
    }

    // read through BAM file
    BamAlignment al;
    string readGroup("");
    map<string, ReadGroupResolver>::iterator rgIter;
    map<string, bool>::iterator readNameIter;
	int readNameCount;
	int notPaired=0, notMapped=0, notMateMapped=0, mateNotOnSameRef=0;

	vector<CigarOp> cigar;
	char cigarType;
	uint32_t cigarLength;
	uint32_t qLength;
	uint32_t RefID;
	uint32_t editDistance;
	int sumMandI, sumDandI; 
  	float coverage, percId;
  	int minCover = 100*(m_settings->Coverage), minId = 92; 
  	int outPercId = 0, outCover = 0;
	int baseInstert = 0, insertLimit = 10;

	cout << "################################" << endl;
	cout << "minCover=" << minCover << endl;
	cout << "percId=" << minId << endl;
	cout << "################################" << endl;
	// For processing percentage identity (PI) and coverage (CO) 
	std::stringstream field;
	string currentPercId, currentCoverage; 
	float curPercId, curCoverage;
	int illumina;

 nextAlignment:
    while ( bamReader.GetNextAlignmentCore(al) ) {

        // flesh out the char data, so we can retrieve its read group ID
        al.BuildCharData();

        // skip if alignment is not paired, mapped, nor mate is mapped
		if (!al.IsPaired()) {notPaired++; goto nextAlignment;};
	   	if (!al.IsMapped()) {notMapped++; goto nextAlignment;}; 
		if (!al.IsMateMapped()) {notMateMapped++; goto nextAlignment;};

        // skip if alignment & mate not on same reference sequence
        if ( al.RefID != al.MateRefID ) {mateNotOnSameRef++; continue;}

		// Augustus-related filters
		qLength = al.Length; 
		al.GetTag("NM", editDistance);

		// Filtering by percentage identity
		percId = (float)100*(qLength-editDistance)/qLength;
		if (percId < minId) { outPercId++; goto nextAlignment; }
		field << percId;		
		al.AddTag("pi", "Z", field.str());
 		field.str("");

		// Filtering by coverage
		sumMandI = sumMandIOperations(al.CigarData, "no");
		coverage = (float)100*sumMandI/qLength;
		if (coverage < minCover) { outCover++; goto nextAlignment; }
		field << coverage;		
		al.AddTag("co", "Z", field.str());
 		field.str("");

  		// // Intron gap filter
  		// baseInsert = sumDandIOperations(al.CigarData, "no");
  		// if (noIntrons && baseInsert > insertLimit) { outIntrons++; goto nextAlignment; }

        // get read group from alignment (OK if empty)
        readGroup.clear();
        al.GetTag(READ_GROUP_TAG, readGroup);

        // look up resolver for read group
        rgIter = m_readGroups.find(readGroup);
        if ( rgIter == m_readGroups.end() )  {
            cerr << "filter ERROR - unable to calculate stats, unknown read group encountered: "
                 << readGroup << endl;
            bamReader.Close();
            return false;
        }
        ReadGroupResolver& resolver = (*rgIter).second;

		// Retrieving PERCID and COVERAGE values stored in alignment
		al.GetTag("pi", currentPercId);
		curPercId = atof(currentPercId.c_str());
		al.GetTag("co", currentCoverage);
  		curCoverage = atof(currentCoverage.c_str());

        // determine unique-ness of current alignment
        const bool isCurrentMateUnique = ( al.MapQuality >= m_settings->MinimumMapQuality );
		const float currentMateScore = (curCoverage+curPercId)/100;

        // look up read name
        readNameIter = resolver.ReadNames.find(al.Name);
		readNameCount = resolver.ReadNames.count(al.Name);
		if (readNameCount > 2){cout << "al.Name aligned mult. times!" << endl;}
		
        // if read name found (current alignment's mate already parsed)
        if ( readNameIter != resolver.ReadNames.end() ) {

            // if both unique mates are unique, store read name & insert size for later
            const bool isStoredMateUnique  = (*readNameIter).second;
            if ( isCurrentMateUnique && isStoredMateUnique ) {

                // save read name in temp file as candidates for later pair marking
                readNamesWriter.Write(readGroup, al.Name);

                // determine model type & store fragment length for stats calculation
                const uint16_t currentModelType = CalculateModelType(al);
                assert( currentModelType != ModelType::DUMMY_ID );
                resolver.Models[currentModelType].push_back( abs(al.InsertSize) );

				// Saving alignments that made it through filter and are unique
				writer.SaveAlignment(al);

            } 

            // unique or not, remove read name from map
            resolver.ReadNames.erase(readNameIter);
        }

        // if read name not found, store new entry
        else {
				resolver.ReadNames.insert( make_pair<string, bool>(al.Name, isCurrentMateUnique) );
			}	

		// ///////////////////////////////////////////////////////////////
        // // if read name found (current alignment's mate already parsed)
        // if ( readNameIter != resolver.ReadNames.end() ) {

        //     // if both unique mates are unique, store read name & insert size for later
        //     const bool isStoredMateUnique  = (*readNameIter).second;
        //     if ( isCurrentMateUnique && isStoredMateUnique ) {

        //         // save read name in temp file as candidates for later pair marking
        //         readNamesWriter.Write(readGroup, al.Name);

        //         // determine model type & store fragment length for stats calculation
        //         const uint16_t currentModelType = CalculateModelType(al);
        //         assert( currentModelType != ModelType::DUMMY_ID );
        //         resolver.Models[currentModelType].push_back( abs(al.InsertSize) );

		// 		// Saving alignments that made it through filter and are unique
		// 		writer.SaveAlignment(al);

        //     }

        //     // unique or not, remove read name from map
        //     resolver.ReadNames.erase(readNameIter);
        // }

        // // if read name not found, store new entry
        // else resolver.ReadNames.insert( make_pair<string, bool>(al.Name, isCurrentMateUnique) );
		// ///////////////////////////////////////////////////////////////

    } // end while

    // close files
    readNamesWriter.Close();
    bamReader.Close();
	writer.Close();

	// Print stats
	cout << "------------------------------------------------------------" << endl;
	cout << "Summary" << endl;
	cout << "notPaired = " << notPaired << endl;
	cout << "notMapped = " << notMapped << endl;
	cout << "notMateMapped = " << notMateMapped << endl;
	cout << "mateNotOnSameRef = " << mateNotOnSameRef << endl;
	cout << "outPercId = " << outPercId << endl;
	cout << "outCover = " << outCover << endl;
	cout << "illumina = " << illumina << endl;
	cout << "------------------------------------------------------------" << endl;

    // iterate back through read groups
    map<string, ReadGroupResolver>::iterator rgEnd  = m_readGroups.end();
    for ( rgIter = m_readGroups.begin(); rgIter != rgEnd; ++rgIter ) {
        const string& name = (*rgIter).first;
        ReadGroupResolver& resolver = (*rgIter).second;

        // calculate acceptable orientation & insert sizes for this read group
		resolver.PrintModels(name);
        resolver.DetermineTopModels(name);

        // clear out left over read names
        // (these have mates that did not pass filters or were already removed as non-unique)
        resolver.ReadNames.clear();
    }

    // if we get here, return success
    return true;
}


/////////////////////////////////////////////////////////////////////////////////

void IlluminaTool::IlluminaToolPrivate::ParseHeader(const SamHeader& header) {

    // iterate over header read groups, creating a 'resolver' for each
    SamReadGroupConstIterator rgIter = header.ReadGroups.ConstBegin();
    SamReadGroupConstIterator rgEnd  = header.ReadGroups.ConstEnd();
    for ( ; rgIter != rgEnd; ++rgIter ) {
        const SamReadGroup& rg = (*rgIter);
        m_readGroups.insert( make_pair<string, ReadGroupResolver>(rg.ID, ReadGroupResolver()) );
    }
}

bool IlluminaTool::IlluminaToolPrivate::ReadStatsFile(void) {

    // skip if no filename provided
    if ( m_settings->StatsFilename.empty() )
        return false;

    // attempt to open stats file
    IlluminaTool::StatsFileReader statsReader;
    if ( !statsReader.Open(m_settings->StatsFilename) ) {
        cerr << "filter ERROR - could not open stats file: "
             << m_settings->StatsFilename << " for reading" << endl;
        return false;
    }

    // attempt to read stats data
    if ( !statsReader.Read(m_settings, m_readGroups) ) {
        cerr << "filter ERROR - could not parse stats file: "
             << m_settings->StatsFilename << " for data" << endl;
        return false;
    }

    // return success
    return true;
}

void IlluminaTool::IlluminaToolPrivate::ResolveAlignment(BamAlignment& al) {

    // clear proper-pair flag
    al.SetIsProperPair(false);

    // quit check if alignment is not from paired-end read
    if ( !al.IsPaired() ) return;

    // quit check if either alignment or its mate are unmapped
    if ( !al.IsMapped() || !al.IsMateMapped() ) return;

    // quit check if alignment & its mate are on differenct references
    if ( al.RefID != al.MateRefID ) return;

    // quit check if map quality less than cutoff
    if ( al.MapQuality < m_settings->MinimumMapQuality ) return;

    // get read group from alignment
    // empty string if not found, this is OK - we handle empty read group case
    string readGroupName("");
    al.GetTag(READ_GROUP_TAG, readGroupName);

    // look up read group's 'resolver'
    map<string, ReadGroupResolver>::iterator rgIter = m_readGroups.find(readGroupName);
    if ( rgIter == m_readGroups.end() ) {
        cerr << "filter ERROR - read group found that was not in header: "
             << readGroupName << endl;
        exit(1);
    }
    const ReadGroupResolver& resolver = (*rgIter).second;

    // quit check if pairs are not in proper orientation (can differ for each RG)
    if ( !resolver.IsValidOrientation(al) ) return;

    // quit check if pairs are not within "reasonable" distance (can differ for each RG)
    if ( !resolver.IsValidInsertSize(al) ) return;

    // quit check if alignment is not a "candidate proper pair"
    map<string, bool>::const_iterator readNameIter;
    readNameIter = resolver.ReadNames.find(al.Name);
    if ( readNameIter == resolver.ReadNames.end() )
        return;

    // if we get here, alignment is OK - set 'proper pair' flag
    al.SetIsProperPair(true);
}

bool IlluminaTool::IlluminaToolPrivate::ResolvePairs(void) {

    // open file containing read names of candidate proper pairs
    IlluminaTool::ReadNamesFileReader readNamesReader;
    if ( !readNamesReader.Open(m_settings->ReadNamesFilename) ) {
        cerr << "filter ERROR: could not open (temp) inputput read names file: "
             << m_settings->ReadNamesFilename << endl;
        return false;
    }

    // parse read names (matching with corresponding read groups)
    if ( !readNamesReader.Read(m_readGroups) ) {
        cerr << "filter ERROR: could not read candidate read names from file: "
             << m_settings->ReadNamesFilename << endl;
        readNamesReader.Close();
        return false;
    }

    // close read name file reader & delete temp file
    readNamesReader.Close();
    if ( remove(m_settings->ReadNamesFilename.c_str()) != 0 ) {
        cerr << "filter WARNING: could not delete temp file: "
             << m_settings->ReadNamesFilename << endl;
    }

    // open our BAM reader
    BamReader reader;
    if ( !reader.Open(m_settings->InputBamFilename) ) {
        cerr << "filter ERROR: could not open input BAM file: "
             << m_settings->InputBamFilename << endl;
        return false;
    }

    // retrieve header & reference dictionary info
    const SamHeader& header = reader.GetHeader();
    const RefVector& references = reader.GetReferenceData();

    // determine compression mode for BamWriter
    bool writeUncompressed = ( m_settings->OutputBamFilename == Options::StandardOut() &&
                               !m_settings->IsForceCompression );
    BamWriter::CompressionMode compressionMode = BamWriter::Compressed;
    if ( writeUncompressed ) compressionMode = BamWriter::Uncompressed;

    // open BamWriter
    BamWriter writer;
    writer.SetCompressionMode(compressionMode);
    if ( !writer.Open(m_settings->OutputBamFilename, header, references) ) {
        cerr << "filter ERROR: could not open "
             << m_settings->OutputBamFilename << " for writing." << endl;
        reader.Close();
        return false;
    }

    // plow through alignments, setting/clearing 'proper pair' flag
    // and writing to new output BAM file
    BamAlignment al;
    while ( reader.GetNextAlignment(al) ) {
        ResolveAlignment(al);
        writer.SaveAlignment(al);
    }

    // clean up & return success
    reader.Close();
    writer.Close();
    return true;
}

bool IlluminaTool::IlluminaToolPrivate::Run(void) {

    // verify that command line settings are acceptable
    vector<string> errors;
    if ( !CheckSettings(errors) ) {
        cerr << "filter ERROR - invalid settings: " << endl;
        vector<string>::const_iterator errorIter = errors.begin();
        vector<string>::const_iterator errorEnd  = errors.end();
        for ( ; errorIter != errorEnd; ++errorIter )
            cerr << (*errorIter) << endl;
        return false;
    }

    // initialize read group map with default (empty name) read group
    m_readGroups.insert( make_pair<string, ReadGroupResolver>("", ReadGroupResolver()) );

    // init readname filename
    // uses (adjusted) stats filename if provided (req'd for makeStats, markPairs modes; optional for twoPass)
    // else keep default filename
    if ( m_settings->HasStatsFile )
        m_settings->ReadNamesFilename = m_settings->StatsFilename + READNAME_FILE_SUFFIX;

    // -makeStats mode
    if ( m_settings->IsMakeStats ) {

        // generate stats data
        if ( !MakeStats() ) {
            cerr << "filter ERROR - could not generate stats" << endl;
            return false;
        }

        // write stats to file
        if ( !WriteStatsFile() ) {
            cerr << "filter ERROR - could not write stats file: "
                 << m_settings->StatsFilename << endl;
            return false;
        }
    }

	/////////////////////////////////////////////////////////////
    // -clean mode
    if ( m_settings->IsClean ) {

        // generate stats data
        if ( !Clean() ) {
            cerr << "filter clean ERROR - could not filter the file" << endl;
            return false;
        }

        // write stats to file
        if ( !WriteStatsFile() ) {
            cerr << "bamtools clean ERROR - could not write stats file: "
                 << m_settings->StatsFilename << endl;
            return false;
        }
    }
	/////////////////////////////////////////////////////////////

    return true;
}

bool IlluminaTool::IlluminaToolPrivate::WriteStatsFile(void) {

    // skip if no filename provided
    if ( m_settings->StatsFilename.empty() )
        return false;

    // attempt to open stats file
    IlluminaTool::StatsFileWriter statsWriter;
    if ( !statsWriter.Open(m_settings->StatsFilename) ) {
        cerr << "filter ERROR - could not open stats file: "
             << m_settings->StatsFilename << " for writing" << endl;
        return false;
    }

    // attempt to write stats data
    if ( !statsWriter.Write(m_settings, m_readGroups) ) {
        cerr << "filter ERROR - could not write stats file: "
             << m_settings->StatsFilename << " for data" << endl;
        return false;
    }

    // return success
    return true;
}

// --------------------------------------------------------------------------
// IlluminaTool implementation

IlluminaTool::IlluminaTool(void)
    : AbstractTool()
    , m_settings(new IlluminaSettings)
    , m_impl(0)
{
    // set description texts
    const string programDescription = "It will clean (filter) paired-end alignments from illumina technology";
    const string programUsage = "<mode> [options] [-in <filename>] [-out <filename>] [-stats <filename>]";
    const string inputBamDescription = "the input BAM file(s)";
    const string outputBamDescription = "the output BAM file";
    const string statsFileDescription = "input/output stats file, depending on selected mode (see below). "
            "This file is human-readable, storing fragment length data generated per read group, as well as "
            "the options used to configure the -makeStats mode";
	//////////////////////////////////////////////////////////
    const string cleanDescription = "will clean paired-end alignments.";
	//////////////////////////////////////////////////////////
    const string makeStatsDescription = "generates a fragment-length stats file from the input BAM. "
            "Data is written to file specified using the -stats option. "
            "MarkPairs Mode Settings are DISABLED.";
    const string coverageDescription = "minimum amount of bases that the read must cover over the region it aligns; "
            "expressed as a percentage. Default 0.95";
    const string confidenceIntervalDescription = "confidence interval. Set min/max fragment lengths such that we capture "
            "this fraction of pairs";
    const string unusedModelThresholdDescription = "unused model threshold. The resolve tool considers 8 possible orientation models "
            "for pairs. The top 2 are selected for later use when actually marking alignments. This value determines the "
            "cutoff for marking a read group as ambiguous. Meaning that if the ratio of the number of alignments from bottom 6 models "
            "to the top 2 is greater than this threshold, then the read group is flagged as ambiguous. By default, NO alignments "
            "from ambiguous read groups will be marked as proper pairs.";
    const string forceMarkDescription = "forces all read groups to be marked according to their top 2 'orientation models'. "
            "When generating stats, the 2 (out of 8 possible) models with the most observations are chosen as the top models for each read group. "
            "If the remaining 6 models account for more than some threshold ([default=10%], see -umt), then the read group is marked as ambiguous. "
            "The default behavior is that for an ambiguous read group, NONE of its alignments are marked as proper-pairs. "
            "By setting this option, a read group's ambiguity flag will be ignored, and all of its alignments will be compared to the top 2 models.";

    // set program details
    Options::SetProgramInfo("filter illumina", programDescription, programUsage);

    // set up I/O options
    OptionGroup* IO_Opts = Options::CreateOptionGroup("Input & Output");
    Options::AddValueOption("-in",  "BAM filename", inputBamDescription, "",
                            m_settings->HasInputBamFile, m_settings->InputBamFilename,
                            IO_Opts, Options::StandardIn());
    Options::AddValueOption("-out", "BAM filename", outputBamDescription, "",
                            m_settings->HasOutputBamFile, m_settings->OutputBamFilename,
                            IO_Opts, Options::StandardOut());
    Options::AddValueOption("-stats", "STATS filename", statsFileDescription, "",
                            m_settings->HasStatsFile, m_settings->StatsFilename, IO_Opts);
    OptionGroup* ModeOpts = Options::CreateOptionGroup("Filtering Modes (must select ONE of the following)");
	//////////////////////////////////////////////////////////////
    Options::AddOption("-clean", cleanDescription, m_settings->IsClean, ModeOpts);
	//////////////////////////////////////////////////////////////
    Options::AddOption("-makeStats", makeStatsDescription, m_settings->IsMakeStats, ModeOpts);
    // OptionGroup* GeneralOpts = Options::CreateOptionGroup("General Resolve Options (available in all modes)");
    // Options::AddValueOption("-minMQ", "unsigned short", minMapQualDescription, "",
    //                         m_settings->HasMinimumMapQuality, m_settings->MinimumMapQuality, GeneralOpts);

    // OptionGroup* MakeStatsOpts = Options::CreateOptionGroup("MakeStats Mode Options (disabled in -markPairs mode)");
    // Options::AddValueOption("-ci", "double", confidenceIntervalDescription, "",
    //                         m_settings->HasConfidenceInterval, m_settings->ConfidenceInterval, MakeStatsOpts);
    // Options::AddValueOption("-umt", "double", unusedModelThresholdDescription, "",
    //                         m_settings->HasUnusedModelThreshold, m_settings->UnusedModelThreshold, MakeStatsOpts);

	////////////////////////////////////////////////////
    OptionGroup* CleanOpts = Options::CreateOptionGroup("-clean Mode Options");
    Options::AddValueOption("-minCover", "double", coverageDescription, "",
                            m_settings->HasCoverage, m_settings->Coverage, CleanOpts);
    Options::AddValueOption("-ci", "double", confidenceIntervalDescription, "",
                            m_settings->HasConfidenceInterval, m_settings->ConfidenceInterval, CleanOpts);
    Options::AddValueOption("-umt", "double", unusedModelThresholdDescription, "",
                            m_settings->HasUnusedModelThreshold, m_settings->UnusedModelThreshold, CleanOpts);
	////////////////////////////////////////////////////

}

IlluminaTool::~IlluminaTool(void) {

    delete m_settings;
    m_settings = 0;

    delete m_impl;
    m_impl = 0;
}

int IlluminaTool::Help(void) {
    Options::DisplayHelp();
    return 0;
}

int IlluminaTool::Run(int argc, char* argv[]) {

    // parse command line arguments
    Options::Parse(argc, argv, 1);

    // initialize IlluminaTool
    m_impl = new IlluminaToolPrivate(m_settings);

    // run IlluminaTool, return success/failure
    if ( m_impl->Run() )
        return 0;
    else
        return 1;
}
