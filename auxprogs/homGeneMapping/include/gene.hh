 /**********************************************************************
 * file:    gene.hh
 * license: Artistic License, see file LICENSE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  
 * authors: Stefanie Koenig, stefaniekoenig@ymail.com
 *
 **********************************************************************/

#ifndef _GENE_HH
#define _GENE_HH

#include <string>
#include <list> 
#include <vector>

// Forward declarations
class GeneFeature;
class Gene;

#define NUM_STRAND_TYPES 2
#define NUM_FEATURE_TYPES 6

enum Strand{unknown=-1,plusstrand,minusstrand};
extern std::string strandIdentifiers[NUM_STRAND_TYPES];

enum FeatureType{unkown=-1, CDS, intron, exon, UTR, start, stop};
extern std::string featureTypeIdentifiers[NUM_FEATURE_TYPES];

/*
 * class GeneFeature:
 * a gene feature is a part of a gene, e.g. CDS exon, intron or UTR exon (not implemented yet),
 * with a start and a end position.
 * For introns the start position is start-1, e.g. the end position of the preceding exon
 * and the end position is end+1, e.g. the start position of the succeeding exon
 * It is also possible that gene features do not belong to a gene, in this case
 * they purely represent extrinsic evidence (either intron or CDS)
 */
class GeneFeature {

public:
    GeneFeature(FeatureType _type, long int _start, long int _end, Strand _strand, int _frame=-1, double _score=0.0) : 
	type(_type),
	start(_start),
	len(_end-_start+1),
	strand(_strand),
	frame(_frame),
	score(_score),
	mult(0),
	gene(NULL)
    {}
    ~GeneFeature() {}

    Strand getStrand() const {return strand;}
    std::string getGeneID() const;
    std::string getTxID() const;
    int getSeqID() const;
    std::string getSource() const;
    void setGene(Gene* gene);
    Gene* getGene();
    long int getStart() const {return start;}
    int getLen() const {return len;}
    void setLen(int l) {len = l;}
    long int getEnd() const {return (start+len-1);}
    int getFrame() const {return frame;}
    std::string writeFrame() const;
    double getScore() const {return score;}
    void setEvidence(std::string e);
    void setMult(int m){mult+=m;}
    void setStart(long int s) {start=s;}
    std::string getEvidence() const;
    int getMult() const {return mult;}
    FeatureType getFeatureType() const {return type;}
    int lenMod3() const {return ((len) % 3);}
    bool hasEvidence(std::string e) const;
    bool hasEvidence() const {return !extrinsic.empty();}
    bool isExon() const {return (type == exon);}
    bool isCDS() const {return (type == CDS);}
    bool isIntron() const {return (type == intron);}
    bool isUTR() const {return (type == UTR);}
    bool isMapped() const {return (type == exon || type == CDS || type == intron);}
    bool isType(FeatureType t) const {return (type == t);}
    bool isPartofGene() const {return gene;} // if false, gene feature purely represent extrinsic evidence
    bool sameStrand(Strand other);
    bool sameFrame(int other);
    void appendHomolog(GeneFeature *gf, int idx) {homologs.push_back(std::pair<int,GeneFeature*>(idx, gf));}

private:
    FeatureType type;      // CDS exon, intron or UTR exon
    long int start;
    int len;               
    Strand strand;         // 'unknown' if both strands are possible, e.g. for GeneFeatures
                           // that purely represent extrinsic Evidence and are
                           // not part of a Gene
    int frame;             // -1 if gene feature has no frame
    double score;
    std::vector< std::string > extrinsic; // list of sources of extrinsic info, e.g. 'M' (manual) or 'E' (EST).
                           // empty, if gene feature is not supported by extrinsic evidence, 
    int mult;
    Gene *gene;            // pointer to the gene the feature belongs to
    /*
     * homologous gene features, e.g. gene features that map
     * to the same start and end position in the alignment.
     * the first value in the pair is the index of the genome, e.g.
     * (0,gf1), (0,gf2), (2,gf3), (3,gf4)
     * means that gene feature has 4 homologs:
     * 2 in genome 0 (gf1 and gf2), 1 in genome 2 (gf3) and 1 in genome 3 (gf4)
     */
public:
    std::list<std::pair<int,GeneFeature*> >homologs;
};

bool compareGFs(GeneFeature *a, GeneFeature *b);

/*
 * class Gene:
 * a gene is basically a list of gene features
 */
class Gene {

public:
    Gene(std::string _geneID, std::string _txID, int _seqID, Strand _strand, std::string _source) :
	geneID(_geneID),
	txID(_txID),
	seqID(_seqID),
	strand(_strand),
	source(_source),
	tlStart(-1),
	tlEnd(-1)
    {}
    ~Gene(){
	for(std::list<GeneFeature*>::iterator it=features.begin(); it!=features.end(); it++)
	    delete *it;
    }
    // get and set functions
    void appendFeature(GeneFeature* f){features.push_back(f);}
    Strand getStrand() const {return strand;}
    std::string getGeneID() const {return geneID;}
    std::string getTxID() const {return txID;}
    int getSeqID() const {return seqID;}
    std::string getSource() const {return source;}
    int numGFs() const {return features.size();}
    int numGFs(FeatureType t) const;
    bool hasFeatures() const {return !features.empty();}
    void sortGFs() {features.sort(compareGFs);}
    void includeStopInCDS();
    void setTLstart(long int s){tlStart=s;}
    void setTLend(long int e){tlEnd=e;}
    long int getTLstart() const {return tlStart;}
    long int getTLend() const {return tlEnd;}
    GeneFeature* getFirstGF(FeatureType t);
    GeneFeature* getLastGF(FeatureType t);
    long int getTXstart() {return getFirstGF(exon)->getStart();}
    long int getTXend() {return getLastGF(exon)->getEnd();}
    void insertExons();
    void insertMissingGFs();
    void appendHomolog(Gene *g, int idx) {homologs.push_back(std::pair<int,Gene*>(idx, g));}


private:
    std::string geneID;
    std::string txID;
    int seqID;
    Strand strand;
    std::string source;
    long int tlStart; // translation start
    long int tlEnd;   // translation end
public:
    /*
     * homologous genes:
     * two genes are homologous if all their gene features are
     * homologs and they have no additional gene features.
     * the first value in the pair is the index of the genome, e.g.
     * (0,g2), (2,g3), (3,g4)
     * means that gene feature has 3 homologs in genome 0,2 and 3, respectively.
     */
    std::list<std::pair<int,Gene*> >homologs;
    std::list<GeneFeature*> features;
};

Strand getStrand(std::string token);
int getFrame(std::string token);
FeatureType getType(std::string token);
    
#endif   //  _GENE_HH
