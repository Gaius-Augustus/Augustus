/**********************************************************************
 * file:    orthoexon.hh
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  maintains orthologous exons for comparative gene prediction
 * authors: Stefanie König
 *
 *********************************************************************/

#ifndef _ORTHOEXON_HH
#define _ORTHOEXON_HH

#include "exoncand.hh"
#include "projectio.hh"
#include "types.hh"
#include "phylotree.hh"

#include <vector>
#include <string>
#include <list>

extern const char*  phyleticPatternIdentifiers[6];

//forward declarations:
class Node;

class OrthoExon {
public:
    OrthoExon(int_fast64_t k, size_t numSpecies);
  ~OrthoExon() {};

    // get and and set functions
    StateType getStateType() const; // all exon candidates agree in type
    int numExons() const;
    bit_vector getBV() const {return bv;}
    vector<int> getRFC(vector<int> offsets) const;
    PhyloTree* getTree() const {return tree;}
    int getAliStart() const {return (key>>22);} // start position of HECT in alignment
    int getAliLen() const {int aliStart=getAliStart(); int n=key-(aliStart<<22); return (n>>7);} // length of HECT + 1
    int getAliEnd() const {return getAliStart() + getAliLen();}
    int getStartInWindow(int s) const {return firstAlignedPos[s];}
    int getEndInWindow(int s) const {return lastAlignedPos[s];}
    bool exonExists(int pos) const; // returns true if OE has a candidate exon at position pos
    bool isUnaligned(int i) const {return labels[i] == 3;} // true, if species i is not aligned
    bool isAbsent(int i) const {return labels[i] == 2;}    // true, if species i is aligned, but ECs is absent
    void setPresent(bit_vector v);
    void setAbsent(bit_vector v);
    
    void setOmega(vector<double>* llo, CodonEvo* codonevo, bool oeStart);
    void storeOmega(double currOmega, double currVarOmega);
    void setSubst(int subs, bool oeStart);
      
    void setTree(PhyloTree* t);
    void setBV(bit_vector b){bv = b;}
    string getPhyleticPattern() const; // phyletic pattern: for an explanation, see .cc file
    void setPhyleticPattern(map<int, list<int> > &pp_init, map<int, list<int> > &pp_opt);
    vector<int> getRFC(vector<int> offsets);
    double getLogRegScore();

    vector<ExonCandidate*> orthoex;
    vector<Node*> orthonode; //corresponding nodes in the graph
    vector<double> weights;
    /*
     * labels:
     * 0 - EC present, but not predicted as exon
     * 1 - EC present and predicted as exon
     * 2 - EC absent, but alignment present
     * 3 - alignment not present
     */
    vector<int> labels;
    int ID;
    vector<int> firstAlignedPos;
    vector<int> lastAlignedPos;
    OEtraits oeTraits;
private:
    int_fast64_t key; // key encodes all of: aliStart aliEnd type lenMod3
  list<vector<double> > loglikOmegaStarts;
  int intervalCount;   
  bit_vector bv; //  stores in one bit for each species its absence/presence (0/1)
    PhyloTree *tree; // corresponding tree topology of an OE
};

/*
 * read and write functions for orthologous exons
 * TODO: - allow orthologous exons to be on different strands
 *       - substract offset; on reverse strand, start/end positions have to be made relative to the
 *         reverse complement
 */
// old code:
//std::list<OrthoExon> readOrthoExons(std::string filename); //reads list of orthologous exons from a file
//void writeOrthoExons(const std::list<OrthoExon> &all_orthoex);


std::ostream& operator<<(std::ostream& ostrm, const OrthoExon &ex_tuple);
std::istream& operator>>(std::istream& istrm, OrthoExon& ex_tuple);

#endif
