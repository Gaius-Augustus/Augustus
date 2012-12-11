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

#include <vector>
#include <string>
#include <list>

//forward declarations:
class ExonCandidate;

class OrthoExon{

public:
    OrthoExon();
    ~OrthoExon() {};
    //copy with permutation of vector entries
    OrthoExon(const OrthoExon& other, const std::vector<size_t> &permutation);

    std::vector<ExonCandidate*> orthoex;
    //TODO: instead of an attribute write a function getLabelpattern() which returns the current
    //label pattern. This is the safer way and guarantees to always have the current label pattern.
    std::string labelpattern;
    int ID;
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
