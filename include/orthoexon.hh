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

#include "orthograph.hh"
#include <vector>
#include <string>

//forward declarations:
class ExonCandidate;

class OrthoExon{

public:
    OrthoExon();
    ~OrthoExon() {};
    //copy with permutation of vector entries
    OrthoExon(const OrthoExon& other, const vector<size_t> &permutation);

    vector<ExonCandidate*> orthoex;
    //TODO: instead of an attribute write a function getLabelpattern() which returns the current
    //label pattern
    string labelpattern;
};

/*
 * read and write functions for orthologous exons
 * TODO: allow orthologous exons to be on different strands
 */
list<OrthoExon> readOrthoExons(string filename); //reads list of orthologous exons from a file
void writeOrthoExons(const list<OrthoExon> &all_orthoex);


ostream& operator<<(ostream& ostrm, const OrthoExon &ex_tuple);
istream& operator>>(istream& istrm, OrthoExon& ex_tuple);

#endif
