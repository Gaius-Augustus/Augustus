 /**********************************************************************
 * file:    gene.hh
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  
 * authors: Stefanie Koenig, stefaniekoenig@ymail.com
 *
 **********************************************************************/

#ifndef _GENE_HH
#define _GENE_HH

#include <string>
#include <list> 

// Forward declarations
class GeneFeature;
class Gene;

#define NUM_STRAND_TYPES 2
#define NUM_FEATURE_TYPES 2

enum Strand{unknown=-1,plusstrand,minusstrand};
extern std::string strandIdentifiers[NUM_STRAND_TYPES];

enum FeatureType{CDS=0, intron};
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
	extrinsic(""),
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
    long int getEnd() const {return (start+len-1);}
    int getFrame() const {return frame;}
    std::string writeFrame() const;
    double getScore() const {return score;}
    void setEvidence(std::string e){ if(extrinsic.empty()){extrinsic=e;}}
    std::string getEvidence() const {return extrinsic;}
    FeatureType getFeatureType() const {return type;}
    int lenMod3() const {return ((len) % 3);}
    bool hasEvidence() const {return !extrinsic.empty();}
    bool isExon() const {return (type == CDS);}
    bool isIntron() const {return (type == intron);}
    bool isPartofGene() const {return gene;} // if false, gene feature purely represent extrinsic evidence
    bool sameStrand(Strand other);
    bool sameFrame(int other);

    std::list< std::pair<int,GeneFeature*> >getHomologs() const {return homologs;}
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
    std::string extrinsic; // source of extrinsic info, e.g. 'M' (manual) or 'E' (EST).
                           // empty, if gene feature is not supported by extrinsic evidence, 
    Gene *gene;            // pointer to the gene the feature belongs to
    /*
     * homologous gene features, e.g. gene features that map
     * to the same start and end position in the alignment.
     * the first value in the pair is the index of the genome, e.g.
     * (0,gf1), (0,gf2), (2,gf3), (3,gf4)
     * means that gene feature has 4 homologs:
     * 2 in genome 0 (gf1 and gf2), 1 in genome 2 (gf3) and 1 in genome 3 (gf4)
     */
    std::list<std::pair<int,GeneFeature*> >homologs;
};

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
    std::list<GeneFeature*> getFeatureList() const {return features;}
    int numGFs() const {return features.size();}
    void setTLstart(long int s) {tlStart = s;}
    void setTLend(long int e) {tlEnd = e;}
    long int getTLstart() const {return tlStart;}
    long int getTLend() const {return tlEnd;}
    long int getStart() const; // gene start
    long int getEnd() const;   // gene end
    bool hasFeatures() const {return !features.empty();}

    std::list< std::pair<int,Gene*> >getHomologs() const {return homologs;}
    void appendHomolog(Gene *g, int idx) {homologs.push_back(std::pair<int,Gene*>(idx, g));}

private:
    std::string geneID;
    std::string txID;
    int seqID;
    Strand strand;
    std::string source;
    long int tlStart; // translation start
    long int tlEnd;   // translation end
    std::list<GeneFeature*> features;

    /*
     * homologous genes:
     * two genes are homologous if all their gene features are
     * homologs and they have no additional gene features.
     * the first value in the pair is the index of the genome, e.g.
     * (0,g2), (2,g3), (3,g4)
     * means that gene feature has 3 homologs in genome 0,2 and 3, respectively.
     */
    std::list<std::pair<int,Gene*> >homologs;
};

Strand getStrand(std::string token);
int getFrame(std::string token);
    
#endif   //  _GENE_HH
