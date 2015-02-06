 /**********************************************************************
 * file:    genome.hh
 * licence: Artistic Licence, see file LICENCE.TXT or 
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

// Forward declarations
class GeneFeature;
struct GeneInfo;

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

    // get and set functions
    void setTmpDir(std::string tmpdir);
    std::string getSeqName(int seqID) const;
    int getSeqID(std::string seqname) const;
    int getIdx() const {return idx;}
    std::string getName(){return name;}

    // read and write functions
    void parseGTF(std::string gtffilename);                          // reads a gene file in gtf format
    void parseExtrinsicGFF(std::string gfffilename);                 // reads a hints file in gff format
    void printGFF(std::string outdir, std::vector<Genome> &genomes); // output a gene in gtf format with
    void printBed();
    void readBed(Genome &other);
    void writeGeneFeature(GeneFeature *gf, std::ofstream &of) const;
    void writeGene(Gene *g, std::ofstream &of) const;
    void writeTLStart(Gene *g, std::ofstream &of) const;
    void writeTLEnd(Gene *g, std::ofstream &of) const;

    // system call to halLiftover (external program)
    void liftOverTo(Genome& other, std::string halfile, std::string halLiftover_exec, std::string halParam);

    // mapping of homologous gene features
    void mapGeneFeatures(std::vector<Genome> &genomes);
    void printDetailed(GeneFeature *g, std::ofstream &of) const;

    // hash functions
    void insertPos(int seqID, long int pos, Strand strand); // insert start/end positions of gene features in mappedPos
    std::vector< std::list< uint_fast64_t > > findMappedPos(int seqID, long int pos, Strand strand);
    void insertSeqInt(GeneFeature* gf); // insert gene features into gfHash
    std::list<GeneFeature*> findSeqInt(uint_fast64_t key, int seqID); // retrieve gene features from gfHash

    // static functions
    static int getNumGenomes(){return no_genomes;}
    static void setNumGenomes(int n){no_genomes=n;}

private:
    std::string name;
    std::string tmpdir;
    int idx; // index of the genome
    std::map<std::string, int> seqnames; // maps seqnames to seqIDs
    std::map<int, std::string> seqIDs  ; // maps seqIDs to seqnames
    std::list<Gene*> genes;
    /*
     * stores all start/end positions of gene features and their homologous positions in the other genomes
     * the key encodes sequence name,position and strand
     * the value is a vector of lists of corrpesonding positions, e.g
     * ["chr11:3133206:-"]->[{"scaffold278:54826:+", "C313481:244:-", "C426955:147702:-"}, - ,  {"chr16:4507589:+"}]
     *                                                         ^                           ^            ^
     *                       list of homologous positions in   |                           |            |
     *                                                      genome 0                    genome 1     genome 2       
     *                                                     3 homologous               no homologous  1 homologous
     *                                                       position                   positions    position
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
};



struct GeneInfo{

public:
    GeneInfo(Gene *_gene, int _numMatchingEs, int _numMatchingIs, bool _frameshift) :
	gene(_gene),
	numMatchingEs(_numMatchingEs),
	numMatchingIs(_numMatchingIs),
	frameshift(_frameshift)
    {}
    ~GeneInfo(){}

    Gene *gene;
    int numMatchingEs;
    int numMatchingIs;
    bool frameshift;
};

#endif   //  _GENOME_HH
