/**********************************************************************
 * file:    genome.hh
 * license: Artistic License, see file LICENSE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  
 * authors: Stefanie Koenig, stefaniekoenig@ymail.com
 *
 **********************************************************************/

#ifndef _GENOME_HH
#define _GENOME_HH

#include <map>
#include <vector>

//project includes
#include "gene.hh"

#ifdef M_SQLITE
#include "sqliteDB.hh"
#endif

#ifdef BOOST
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

struct VertexProperty
{
    std::string name;
};

typedef boost::adjacency_list <boost::setS, boost::vecS, boost::undirectedS, VertexProperty> Graph;
#endif

// Forward declarations
class GeneFeature;
struct GeneInfo;
struct GeneInfoCollection;

/*
 * class Genome:
 * a genome is basically a list of genes
 */    
class Genome {

public:
    Genome(std::string _name, int _idx) : 
	name(_name),
	idx(_idx)
    {}
    ~Genome() {}
    void destroyGeneList();
    void destroyHintList();

    // get and set functions
    void setTmpDir(std::string tmpdir);
    std::string getSeqName(int seqID) const;
    int getSeqID(std::string seqname) const;
    int getIdx() const {return idx;}
    std::string getName(){return name;}

    // read and write functions
    void parse(std::string genefile, std::string hintsfile, std::string dbfile);
    void parseGTF(std::string gtffilename);                          // reads a gene file in gtf format
    void parseExtrinsicGFF(std::string gfffilename);                 // reads a hints file in gff format
    void insertHint(std::string seqname, long int start, long int end, Strand strand, std::string esource, int mult, int frame, std::string f_type);
#ifdef M_SQLITE
    void getDbHints(SQLiteDB &db);
#endif
    //void printGFF(std::string outdir, std::vector<Genome> &genomes, bool detailed=false); // output a gene in gtf format with
    void printBed();
    void readBed(Genome &other);
    void writeGeneFeature(GeneFeature *gf, std::ofstream &of, bool unmapped=false) const;
    void writeGene(Gene *g, std::ofstream &of) const;
    void writeTxLine(Gene *g, std::ofstream &of) const;
    void writeTLStart(Gene *g, std::ofstream &of) const;
    void writeTLEnd(Gene *g, std::ofstream &of) const;

    // system call to halLiftover (external program)
    void liftOverTo(Genome& other, std::string halfile, std::string halLiftover_exec, std::string halParam);

    // mapping of homologous gene features
    void write_hgm_gff(std::vector<Genome> &genomes, std::string outdir, bool detailed=false, bool unmapped=false);
    void mapGeneFeatures(std::vector<Genome> &genomes, std::string outdir, bool detailed=false, bool unmapped=false);
    void print_hgm_info(GeneFeature *g, std::ofstream &of) const;
    void print_hgm_unaligned(GeneFeature *g, std::ofstream &of) const;

    // hash functions
    void insertPos(int seqID, long int pos); // insert start/end positions of gene features in mappedPos
    std::vector< std::list< uint_fast64_t > > findMappedPos(int seqID, long int pos);
    void insertSeqInt(GeneFeature* gf); // insert gene features into gfHash
    void insertSeqInt(GeneFeature* gf, int seqID); // insert hints into gfHash
    void insertSeqInts(Gene *g); // inserts all gene features of g into gfHash and mappedPos
    std::list<GeneFeature*> findSeqInt(uint_fast64_t key, int seqID, Strand strand); // retrieve gene features from gfHash

    // static functions
    static int getNumGenomes(){return no_genomes;}
    static void setNumGenomes(int n){no_genomes=n;}

private:
    std::string name;
    std::string tmpdir;
    int idx; // index of the genome
    std::map<std::string, int> seqnames; // maps seqnames to seqIDs
    std::map<int, std::string> seqIDs  ; // maps seqIDs to seqnames
    std::list<GeneFeature*> hints;
    /*
     * stores all start/end positions of gene features and their homologous positions in the other genomes
     * the key encodes sequence name and position
     * the value is a vector of lists of corrpesonding positions, e.g
     * ["chr11:3133206"]->[{"scaffold278:54826", "C313481:244", "C426955:147702"}, - ,  {"chr16:4507589"}]
     *                                                 ^                           ^            ^
     *               list of homologous positions in   |                           |            |
     *                                              genome 0                    genome 1     genome 2       
     *                                             3 homologous               no homologous  1 homologous
     *                                               position                   positions    position
     */
    std::map<uint_fast64_t, std::vector< std::list< uint_fast64_t > > > mappedPos;
    /*
     * hash to quickly access gene features
     * the first key encodes the sequence name
     * the second key encodes start,len,feature type and strand, e.g.
     * ["chr11"]->[["3142308:500:exon:+"]->{exon1,exon2,...}]
     */
    std::map<int,std::map<uint_fast64_t, std::list<GeneFeature*> > > gfHash;

    static int no_genomes; // number of genomes

public:
    std::list<Gene*> genes;
};

// print a list with homologous transcript IDs, e.g.                                                                                    
// # 0     dana                                                                                                                         
// # 1     dere                                                                                                                         
// # 2     dgri                                                                                                                         
// # 3     dmel                                                                                                                         
// # 4     dmoj                                                                                                                         
// # 5     dper                                                                                                                         
// (0, jg4139.t1), (0, jg4140.t1), (1, jg7797.t1), (2, jg3247.t1), (4, jg6720.t1), (5, jg313.t1)                                        
// (1, jg14269.t1), (3, jg89.t1) (5, jg290.t1)                                                                                          
// ...  
#ifdef BOOST
void printHomGeneList(std::string outfile, std::vector<Genome> &genomes);
#endif

struct GeneInfo{

public:
    GeneInfo(Gene *_gene, int _numMatchingCs, int _numMatchingIs, int _numMatchingEs, bool _frameshift) :
	gene(_gene),
	numMatchingCs(_numMatchingCs),
	numMatchingIs(_numMatchingIs),
	numMatchingEs(_numMatchingEs),
	frameshift(_frameshift)
    {}
    ~GeneInfo(){}

    Gene *gene;
    int numMatchingCs;
    int numMatchingIs;
    int numMatchingEs;
    bool frameshift;
};

struct GeneInfoCollection{

 public:
    GeneInfoCollection(int no_genomes){
	ginfo.resize(no_genomes);
	mappedStatsC.resize(no_genomes);
	mappedStatsI.resize(no_genomes);
	mappedStatsE.resize(no_genomes);
	extrinStatsC.resize(no_genomes+1);
	extrinStatsI.resize(no_genomes+1);
	extrinStatsE.resize(no_genomes+1);
    }
    ~GeneInfoCollection(){}
    void createCollection(GeneFeature *g);
    void printDetailedStats(Gene *g, std::ofstream &of);
    std::vector<std::map<std::string,GeneInfo> > ginfo;
    std::vector<int> mappedStatsC;   // number of CDS with exact homologs in at least k other genomes                                                                    
    std::vector<int> mappedStatsI;   // number of Intr ...                                                                                                               
    std::vector<int> mappedStatsE;   // number of Exons ...                                                                                                               
    std::vector<int> extrinStatsC;   // number of CDS supported by evidence in at least k other genomes                                                                  
    std::vector<int> extrinStatsI;   // number of Intr ...                                                                                                               
    std::vector<int> extrinStatsE;   // number of Exons ... 
};

#endif   //  _GENOME_HH
