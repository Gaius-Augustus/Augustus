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
#include <map>
#include <vector>
#include <string>

//forward declarations:
class ExonCandidate;

class OrthoExon{

public:
    OrthoExon();
    ~OrthoExon();
    //copy constructor
    OrthoExon(const OrthoExon& other);
    //copy with permutation of vector entries
    OrthoExon(const OrthoExon& other, const vector<size_t> &permutation);

    vector<ExonCandidate*> orthoex;
    string labelpattern;  //changes dynamically and has to be updated after every optimization step
};

/*
 * read and write functions for orthologous exons
 */
list<OrthoExon> readOrthoExons(string filename); //read list of orthologous exons from a file
void writeOrthoExons(const list<OrthoExon> &all_orthoex);


ostream& operator<<(ostream& ostrm, const OrthoExon &ex_tuple);
istream& operator>>(istream& istrm, OrthoExon& ex_tuple);

#endif
