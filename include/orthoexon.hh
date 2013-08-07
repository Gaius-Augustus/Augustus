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
class ExonCandidate;
class Node;

class OrthoExon {
public:
    OrthoExon();
    ~OrthoExon() {};
    //copy with permutation of vector entries
    OrthoExon(const OrthoExon& other, const std::vector<size_t> &permutation);
    StateType getStateType() const; // all exon candidates agree in type
    std::vector<ExonCandidate*> orthoex;
    std::vector<Node*> orthonode; //corresponding nodes in the graph
    std::vector<double> weights;
    std::vector<int> labels;
    //TODO: instead of an attribute write a function getLabelpattern() which returns the current
    //label pattern. This is the safer way and guarantees to always have the current label pattern.
    std::string labelpattern;
    int ID;
    double getOmega() const { return omega; }
    double getSubst() const { return subst; }
    void setOmega(double o){omega=o;}
    void setSubst(int s){subst=s;}
    
private:
    double omega;
    int subst;
};

/*
 * read and write functions for orthologous exons
 * TODO: - allow orthologous exons to be on different strands
 *       - substract offset; on reverse strand, start/end positions have to be made relative to the
 *         reverse complement
 */
std::list<OrthoExon> readOrthoExons(std::string filename); //reads list of orthologous exons from a file
void writeOrthoExons(const std::list<OrthoExon> &all_orthoex);


std::ostream& operator<<(std::ostream& ostrm, const OrthoExon &ex_tuple);
std::istream& operator>>(std::istream& istrm, OrthoExon& ex_tuple);

#endif
