/*
 * randseqaccess.hh
 *
 * License: Artistic License, see file LICENSE.TXT or 
 *          https://opensource.org/licenses/artistic-license-1.0
 */

#ifndef _RANDSEQACCESS
#define _RANDSEQACCESS

// project includes
#include "gene.hh"
#include "types.hh"
#include "extrinsicinfo.hh"

#include <map>
#include <vector>
#include <cstring>

#ifdef M_MYSQL
#include <mysql++.h>
#endif

#ifdef M_SQLITE
#include "sqliteDB.hh"
#endif

/**
 * @brief SpeciesCollection holds all extrinsic evidence given for the species.
 * @details It consists of a set of group specific FeatureCollections and
 * a default FeatureCollection. Species for which no extrinsic evidence
 * is given, make use of the default collection (identical to ab initio
 * gene prediction, no bonus/malus).
 * Subsets of the species with the same extrinsic config, i.e. same feature table
 * in the extrinsicCfgFile, share one group specific FeatureCollection.
 * 
 * @author Stefanie Koenig
 */
class SpeciesCollection{
public:
    FeatureCollection* getFeatureCollection(string speciesname);
    int getGroupID(string speciesname);
    void addSpeciesToGroup(string skey, int groupID);
    bool withEvidence(string speciesname){return getGroupID(speciesname)>0;}
    // reading in the extrinsicCfgFile and hintsFile
    void readGFFFile(const char* filename); 
    void readExtrinsicCFGFile(vector<string> &speciesNames);
private:
    map<int,FeatureCollection> speciesColl; // maps the group number to a FeatureCollection
    map<string,int> groupIDs; // maps the speciesname to the group number
    FeatureCollection defaultColl; // default FeatureColleciton
    static int groupCount; // number of groups
};

/**
 * @brief abstract class for quick access to an arbitrary sequence segment in genomes
 * needed for comparative gene finding
 * @details random acces to sequence data, e.g. get me chr1:1000-2000 from species 'human'
 * 
 * @author Mario Stanke
 */
class RandSeqAccess {
public:
    int getNumSpecies() {return numSpecies;}
    void setLength(int idx, string chrName, int len);
    int getChrLen(int idx, string chrName);
    void setSpeciesNames(vector<string> speciesNames);
    string getSname(size_t idx) {return speciesNames[idx];}
    int getMaxSnameLen(); // for neat indentation into right column
    int getIdx(string speciesname);
    void printStats();
    bool withEvidence(string speciesname) {return extrinsicFeatures.withEvidence(speciesname);}
    virtual AnnoSequence* getSeq(string speciesname, string chrName, int start, int end, Strand strand) =  0;
    AnnoSequence* getSeq(size_t speciesIdx, string chrName, int start, int end, Strand strand) {
	return getSeq(getSname(speciesIdx), chrName, start, end, strand);
    }
    virtual SequenceFeatureCollection* getFeatures(string speciesname, string chrName, int start, int end, Strand strand) = 0;  
    virtual ~RandSeqAccess() {}
protected:
    RandSeqAccess() {};
    int numSpecies;
    vector<map<string,int> > chrLen;
    vector<string> speciesNames;
    map<string, size_t> speciesIndex; // to quickly access the index for a given species name
    SpeciesCollection extrinsicFeatures; // all hints
};

/**
 * @brief Achieve random access by simply storing all genomes in memory and 
 * then retrieving the required substrings when desired. This may need a lot of RAM.
 * 
 * @author Mario Stanke
 */
class MemSeqAccess : public RandSeqAccess {
public:
    MemSeqAccess(vector<string> s);
    ~MemSeqAccess(){} // TODO: delete DNA sequences from 'sequences' map
    AnnoSequence* getSeq(string speciesname, string chrName, int start, int end, Strand strand);
    SequenceFeatureCollection* getFeatures(string speciesname, string chrName, int start, int end, Strand strand);
    void open(){}
private:
    map<string,string> filenames;
    map<string,char*> sequences;  //keys: speciesname:chrName values: dna sequence
};


/*
 * read an input file of format:
 * human        <TAB> /dir/to/genome/genome.fa
 * Mus musculus <TAB> /dir/to/genome/mouse.fa
 * to a map
 */
map<string,string> getFileNames (string listfile);

/**
 * @brief Random access to sequence segments through a database.
 * @details The sequences must be stored in a database.
 * 
 * @author Stefanie Koenig
 */
class DbSeqAccess : public RandSeqAccess {
public:
    virtual AnnoSequence* getSeq(string speciesname, string chrName, int start, int end, Strand strand)=0;
    virtual SequenceFeatureCollection* getFeatures(string speciesname, string chrName, int start, int end, Strand strand)=0;  
    virtual ~DbSeqAccess() {}

protected:
    DbSeqAccess(vector<string> s = vector<string>());
    string dbaccess;

};

#ifdef M_MYSQL
class MysqlAccess : public DbSeqAccess {
public:
    MysqlAccess(vector<string> s = vector<string>()) : DbSeqAccess(s){
	open();
    }
    ~MysqlAccess() {}
    AnnoSequence* getSeq(string speciesname, string chrName, int start, int end, Strand strand);
    // the following function is for the BGI-style database
    AnnoSequence* getSeq2(string speciesname, string chrName, int start, int end, Strand strand);
    SequenceFeatureCollection* getFeatures(string speciesname, string chrName, int start, int end, Strand strand);  
    void open();
    int split_dbaccess();
    void connect_db(ostream& out=cout);
    template<class T>  
    AnnoSequence* getNextDBSequence(string charName, int start, int end, vector<T>& asm_query_region);
    // template<class T>
    // AnnoSequence* getDBSequenceList(string charName,int start,int end,vector<T>& asm_query_region);
    template<class T>
    int get_region_coord(int seq_region_id, int start, int end, vector<T>& asm_query_region);

private:
    mysqlpp::Connection con;
    vector<string> db_information;
};
#endif // M_MYSQL

#ifdef M_SQLITE

/**
 * @brief Random access to sequence segments through a database.
 * @details The sequences must be stored in a database.
 * 
 * @author Stefanie Koenig
 */
class SQLiteAccess : public DbSeqAccess {
public:
    SQLiteAccess(const char* f, vector<string> s = vector<string>()) : DbSeqAccess(s), db(f) {
	filenames = getFileNames (Constant::speciesfilenames);
    }
    ~SQLiteAccess() {}
    AnnoSequence* getSeq(string speciesname, string chrName, int start, int end, Strand strand);
    SequenceFeatureCollection* getFeatures(string speciesname, string chrName, int start, int end, Strand strand);
private:
    SQLiteDB db;
    map<string,string> filenames;
};
#endif // M_SQLITE

#endif  // _RANDSEQACCESS
