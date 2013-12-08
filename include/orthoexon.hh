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
#include "orthograph.hh"
#include "types.hh"

#include <vector>
#include <string>
#include <list>

//forward declarations:
class Node;

class OrthoExon {
public:
    OrthoExon(int_fast64_t k);
    ~OrthoExon() {};
    StateType getStateType() const; // all exon candidates agree in type
    int numExons() const;
    double getOmega() const { return omega; }
    double getSubst() const { return subst; }
    double getConsScore() const {return cons;}
    double getDiversity() const {return diversity;}
    size_t getContainment() const {return containment;}
    int getAliStart() const {return (key>>22);} // start position of HECT in alignment
    int getAliLen() const {int aliStart=getAliStart(); int n=key-(aliStart<<22); return (n>>7);} // length of HECT + 1
    int getAliEnd() const {return getAliStart() + getAliLen();}
    bool exonExists(int pos) const; // returns true if OE has a candidate exon at position pos
    void setOmega(double o){omega=o;}
    void setSubst(int s){ subst=s;}
    void setConsScore(double c){cons=c;}
    void setDiversity(double d){diversity=d;}
    void setContainment(int c) { containment = c; }
    void setLabelpattern(); // updates the labelpattern
    string getStoredLabelpattern() {return labelpattern;} // retrieving labelpattern without updating
    string getCurrentLabelpattern() { setLabelpattern(); return labelpattern;} //retrieving labelpattern with updating

    vector<ExonCandidate*> orthoex;
    vector<Node*> orthonode; //corresponding nodes in the graph
    vector<double> weights;
    vector<int> labels;
    //TODO: instead of an attribute write a function getLabelpattern() which returns the current
    //label pattern. This is the safer way and guarantees to always have the current label pattern.
    std::string labelpattern;
    int ID;
    
private:
    int_fast64_t key; // key encodes all of: aliStart aliEnd type lenMod3
    double omega;
    int subst;
    double cons; // conservation score
    double diversity; // sum of branch lengths of the subtree induced by the OrthoExon (measure of phylogenetic diversity)
    size_t containment; // how many bases overhang on average has the largest OrthoExon that includes this one in the same frame
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
