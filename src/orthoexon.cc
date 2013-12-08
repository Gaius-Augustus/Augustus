/**********************************************************************
 * file:    orthoexon.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  maintains orthologous exons for comparative gene prediction
 * authors: Stefanie König
 *
 * date    |   author      |  changes
 * --------|---------------|------------------------------------------
 * 09.03.12|Stefanie König | creation of the file
 **********************************************************************/

#include "orthoexon.hh"

#include <fstream>
#include <iostream>


OrthoExon::OrthoExon(int_fast64_t k) : key(k), omega(-1.0) , subst(-1), cons(-1.0), diversity(-1.0) {
    orthoex.resize(OrthoGraph::numSpecies);
    orthonode.resize(OrthoGraph::numSpecies);
    weights.resize(OrthoGraph::numSpecies,0);
    labels.resize(OrthoGraph::numSpecies);
}

StateType OrthoExon::getStateType() const{
    for (int s = 0; s < orthoex.size(); s++)
	if (orthoex[s])
	    return orthoex[s]->getStateType();
    return TYPE_UNKNOWN;
}

int OrthoExon::numExons() const{
    int k=0;
    for (size_t s = 0; s < orthoex.size(); s++)
	if (orthoex[s])
	    k++;
    return k;
}

bool OrthoExon::exonExists(int pos) const{
    if(pos >= orthoex.size() || pos < 0)
	throw ProjectError("Internal error in OrthoExon::exonExists: out of bound access");
    if(orthoex[pos])
	return true;
    return false;	    
}

void OrthoExon::setLabelpattern(){
    labelpattern.clear();
    for(size_t pos = 0; pos < orthonode.size(); pos++){
	if(orthonode[pos]){
	    labelpattern+=itoa(orthonode[pos]->label);
	}
	else{
	    labelpattern+="-";
	}
    }
}


// old code:
/*list<OrthoExon> readOrthoExons(string filename){

    list<OrthoExon> all_orthoex;

    ifstream istrm; 
    istrm.open(filename.c_str(), ifstream::in);
    if (istrm) {
	int nspecies;
	string species;
	vector<size_t> permutation;
	istrm >> goto_line_after( "[SPECIES]");
	istrm >> comment >> nspecies;
	if (nspecies != OrthoGraph::numSpecies){
	    throw ProjectError("readOrthoExons: number of species in " + filename + 
			       " is " + itoa(nspecies) + ". Number of species in treefile is " + itoa(OrthoGraph::numSpecies));
	}
	istrm >> comment;
	for (int i = 0; i < nspecies; i++){
	    istrm >> species;
	    int pos = OrthoGraph::tree->findIndex(species);
	    if (pos == OrthoGraph::numSpecies){
		throw ProjectError("readOrthoExons: species name in " + filename + 
				   " is not a species name in treefile.");
	    }
	    permutation.push_back(pos);
	}
	vector<string> chr(nspecies);
	while(istrm){
	    istrm >> goto_line_after( "[CHR]") >> comment;
	    for (int i = 0; i < nspecies; i++){
		istrm >> chr[permutation[i]];
	    } 
	    cout << endl;
	    istrm >> goto_line_after( "[ORTHOEX]");
	    istrm >> comment;
	    while( istrm >> comment >> ws, istrm && istrm.peek() != '[' ){
		OrthoExon ex_tuple;
		istrm >> ex_tuple;
		all_orthoex.push_back(OrthoExon(ex_tuple, permutation));
	    }
	} 
    }
    else
	throw ProjectError("readOrthoExons: Could not open this file!");

    return all_orthoex;
    }

void writeOrthoExons(const list<OrthoExon> &all_orthoex){
    cout << "# orthologous exons\n" << "#\n" <<"[SPECIES]\n" << "# number of species" << endl;
    cout << OrthoGraph::numSpecies << endl;
    vector<string> species;
    OrthoGraph::tree->getSpeciesNames(species);
    cout << "# species names" << endl;
    for (size_t i = 0; i < OrthoGraph::numSpecies; i++){
	cout << species[i] << "\t";
    }
    cout << endl;
    cout << "#[ORTHOEX]" << endl;
    for(list<OrthoExon>::const_iterator it = all_orthoex.begin(); it != all_orthoex.end(); it++){
	cout << *it << endl;
    }
}
*/
ostream& operator<<(ostream& ostrm, const OrthoExon &oe){
    ostrm << stateTypeIdentifiers[oe.getStateType()];
    for (int s = 0; s < oe.orthoex.size(); s++){
	ExonCandidate *ec = oe.orthoex[s];
	if (ec)
	    ostrm << "\t" << ec->begin+1 << "\t" << ec->end - ec->begin + 1;
	else
	    ostrm << "\t" << 0 << "\t" << 0 << "\t";
    }
    return ostrm;
}

istream& operator>>(istream& istrm, OrthoExon& ex_tuple){

    string exontype;
    long int begin, length;

    istrm >> exontype;
    for (int i = 0; i < OrthoGraph::numSpecies; i++){
	istrm >> begin >> length;
	if (begin != 0 && length != 0){
	    ExonCandidate *exoncand = new ExonCandidate(toExonType(exontype.c_str()), begin-1, begin+length-2);
	    ex_tuple.orthoex[i] = exoncand;
	}
    }
    return istrm;
}
