/**********************************************************************
 * file:    geneticcode.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  
 * authors: Mario Stanke, mario@gobics.de
 *
 * date    |   author      |  changes 
 * --------|---------------|------------------------------------------
 * 15.09.02| Mario Stanke  | creation of the class
 * 13.06.07| Oliver Keller | added alternative translation tables
 **********************************************************************/

#include "geneticcode.hh"

// standard C/C++ includes
#include <iostream>
#include <iomanip>  // for setw

const char GeneticCode::aa_symbols_with_stop[NUM_AA+2] = "*GDERKNQSTAVLIFYWHMCP";
const char* const GeneticCode::aa_symbols = aa_symbols_with_stop + 1;
const char* const GeneticCode::aa_names[NUM_AA] = {"GLYCINE","ASPARTIC ACID","GLUTAMIC ACID","ARGININE","LYSINE","ASPARAGINE",
						   "GLUTAMINE","SERINE","THREONINE","ALANINE","VALINE","LEUCINE",
						   "ISOLEUCINE","PHENYLALANINE","TYROSINE","TRYPTOPHAN","HISTIDINE",
						   "METHIONINE","CYSTEINE","PROLINE"};

/*
 * codons are numbered in lexicographic order: aaa=0, aac=1, aag=2, aat=3, aca=4, ...
 * for c = 0, .. , 63  map[c] is the index of the amino acid c is coding for
 * or -1 if c is not coding for an amino acid
 */

int GeneticCode::map[64] = { 4, 5, 4, 5, // aa.
			     8, 8, 8, 8, // ac.
			     3, 7, 3, 7, // ag.
		            12,12,17,12, // at.
			     6,16, 6,16, // ca.
			    19,19,19,19, // cc.
			     3, 3, 3, 3, // cg.
			    11,11,11,11, // ct.
			     2, 1, 2, 1, // ga.
			     9, 9, 9, 9, // gc.
			     0, 0, 0, 0, // gg.
			    10,10,10,10, // gt.
			    -1,14,-1,14, // ta.
			     7, 7, 7, 7, // tc.
			    -1,18,15,18, // tg.
			    11,13,11,13  // tt.
};

const char * const GeneticCode::TranslationTables[24] = 
{"", 
 "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF",
 "KNKNTTTT*S*SMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
 "KNKNTTTTRSRSMIMIQHQHPPPPRRRRTTTTEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
 "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
 "KNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
 "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVQYQYSSSS*CWCLFLF",
 "", "",
 "NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
 "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSCCWCLFLF",
 "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF",
 "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLSLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF",
 "KNKNTTTTGSGSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
 "NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVYY*YSSSSWCWCLFLF",
 "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YQYSSSS*CWCLFLF",
 "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YLYSSSS*CWCLFLF",
 "", "", "", "",
 "NNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
 "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YLY*SSS*CWCLFLF",
 "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWC*FLF"};

int **GeneticCode::syncodons = 0;
int *GeneticCode::codonsOfAA = 0;


void GeneticCode::reverseMap(){
    if (syncodons != 0) 
	return;    // has already been done before
    syncodons = new int*[20];
    codonsOfAA = new int[20];
    int codnumber;
    for (int aa=0; aa<20; aa++) {
	codnumber=0;
	for (int c=0; c<64; c++) {
	    if (map[c]==aa) 
		codnumber++;
	}
	codonsOfAA[aa] = codnumber;
	syncodons[aa] = new int[codnumber];
    }
    // fill syncodons
    for (int aa=0; aa<20; aa++) {
	codnumber=0;
	for (int c=0; c<64; c++) {
	    if (map[c]==aa) { 
		syncodons[aa][codnumber] = c;
		codnumber++;
	    }
	}
    }
}

void GeneticCode::printReverseGeneticMap(){
    Seq2Int s2i(3);    
    for (int aa = 0; aa<20; aa++) {
	cout << setw(12) << aa_names[aa] << " ("<< aa_symbols[aa] << ") : "; 
	for (int cidx=0; cidx<codonsOfAA[aa]; cidx++) {
	    cout << s2i.INV(syncodons[aa][cidx]);
	    if (cidx < codonsOfAA[aa]-1)
		cout << ", ";
	}
	cout << endl;
    }
}

void GeneticCode::chooseTranslationTable(int n) {
    if (n >= 24 || n < 0) n=1;
    if (string(TranslationTables[n]) == "") n=1;
    for (int codon=0; codon<64; codon++) {
	int aa = get_aa_from_symbol(string(TranslationTables[n]).at(codon));
	if (aa != map[codon]) {
	    // DEBUG MESSAGE
	    cout << "# Warning: Using nonstandard genetic code: " 
		 << Seq2Int(3).INV(codon)
		 << " coding for " << (aa < 0 ? "STOP" : GeneticCode::aa_names[aa])
		 << " instead of " << (map[codon] < 0 ? "STOP" : GeneticCode::aa_names[map[codon]])
		 << ".\n";
	}
	map[codon]=aa;
    }
}

char GeneticCode::translate(const char* t) {
    try {
	return translate(Seq2Int(3)(t));
    } catch (InvalidNucleotideError) {
	return 'X';
    }
}
char GeneticCode::revtranslate(const char* t) {
    try {
	return translate(Seq2Int(3).rc(t));
    } catch (InvalidNucleotideError) {
	return 'X';
    }
}
	
/*
 * containsInFrameStopcodon:
 * Check whether the interval from begin to end contains a stop codon.
 * If frame is specified, it is as in the gff format.
 * If frame is not specified (-1), then all reading frames are checked.
 */
bool GeneticCode::containsInFrameStopcodon(const char* dna, int begin, int end, bool forward, int frame) {
    bool hasStop = false;
    int i;
    if (forward) {
	if ( frame >= 0 && frame <=2) {
	    for (i=begin + frame; i < end-2 && !hasStop; i+=3)
		if (isStopcodon(dna+i))
		    hasStop = true;
	} else {
	    for (i=begin; i < end-2 && !hasStop; i++)
		if (isStopcodon(dna+i))
		    hasStop = true;
	}
    } else { // reverse
	if ( frame >= 0 && frame <=2) {
	    for (i=end - 3 - frame; i >= begin && !hasStop; i-=3)
		if (isRCStopcodon(dna+i))
		    hasStop = true;
	} else {
	    for (i=begin; i < end-2 && !hasStop; i++)
		if (isRCStopcodon(dna+i))
		    hasStop = true;
	}
    }
    return hasStop;
}

  
