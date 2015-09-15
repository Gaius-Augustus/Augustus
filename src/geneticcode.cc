/**********************************************************************
 * file:    geneticcode.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  
 * authors: Mario Stanke, mario.stanke@uni-greifswald.de
 *
 * date    |   author      |  changes 
 * --------|---------------|------------------------------------------
 * 15.09.02| Mario Stanke  | creation of the class
 * 13.06.07| Oliver Keller | added alternative translation tables
 * 14.02.13| Mario Stanke  | added alternative start codon tables
 **********************************************************************/

#include "geneticcode.hh"

// standard C/C++ includes
#include <iostream>
#include <iomanip>  // for setw

const char GeneticCode::aa_symbols_with_stop[NUM_AA+2] = "*GDERKNQSTAVLIFYWHMCP";
// {A,C,T}TG are start codons in the default translation table 1
bool GeneticCode::start_codons[64];
// by default, ATG has probability 1, the rest 0
Double GeneticCode::start_codon_probs[64] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
					      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
int GeneticCode::translationtable = 1;
int GeneticCode::numStartCodons = 3;
const char* const GeneticCode::aa_symbols = aa_symbols_with_stop + 1;
const char* const GeneticCode::aa_names[NUM_AA] = {"GLYCINE","ASPARTIC ACID","GLUTAMIC ACID","ARGININE","LYSINE","ASPARAGINE",
						   "GLUTAMINE","SERINE","THREONINE","ALANINE","VALINE","LEUCINE",
						   "ISOLEUCINE","PHENYLALANINE","TYROSINE","TRYPTOPHAN","HISTIDINE",
						   "METHIONINE","CYSTEINE","PROLINE"};

Seq2Int GeneticCode::codon(3);

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

const char * const GeneticCode::TranslationTables[NUM_TRANSTABS+1] = 
{"", 
 "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF", // 1
 "KNKNTTTT*S*SMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF", // 2
 "KNKNTTTTRSRSMIMIQHQHPPPPRRRRTTTTEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF", // 3
 "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF", // 4
 "KNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF", // 5
 "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVQYQYSSSS*CWCLFLF", // 6
 "", "",
 "NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF", // 9
 "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSCCWCLFLF", // 10
 "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF", // 11
 "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLSLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF", // 12
 "KNKNTTTTGSGSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF", // 13
 "NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVYY*YSSSSWCWCLFLF", // 14
 "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YQYSSSS*CWCLFLF", // 15
 "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YLYSSSS*CWCLFLF", // 16
 "", "", "", "",
 "NNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF", // 21
 "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YLY*SSS*CWCLFLF", // 22
 "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWC*FLF", // 23
 "KNKNTTTTSSKSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"};// 24
//AAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTT
//AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT
//ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT

const char * const GeneticCode::StartCodons[NUM_TRANSTABS+1] = 
{"", 
 "--------------M---------------M-------------------------------M-", // 1
 "------------MMMM------------------------------M-----------------", // 2
 "------------M-M-------------------------------------------------", // 3
 "------------MMMM--------------M---------------M-------------M-M-", // 4
 "------------MMMM------------------------------M---------------M-", // 5
 "--------------M-------------------------------------------------", // 6
 "", "",
 "--------------M-------------------------------M-----------------", // 9
 "--------------M-------------------------------------------------", // 10
 "------------MMMM--------------M---------------M---------------M-", // 11
 "--------------M---------------M---------------------------------", // 12
 "------------M-M-------------------------------M---------------M-", // 13
 "--------------M-------------------------------------------------", // 14
 "--------------M-------------------------------------------------", // 15
 "--------------M-------------------------------------------------", // 16
 "", "", "", "",
 "--------------M-------------------------------M-----------------", // 21
 "--------------M-------------------------------------------------", // 22
 "--------------MM------------------------------M-----------------", // 23
 "--------------M---------------M---------------M---------------M-"}; // 24

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
    if (n > NUM_TRANSTABS || n < 0) n=1; // default translation table is 1
    if (string(TranslationTables[n]) == "") n=1;
    numStartCodons = 0;
    for (int c=0; c<64; c++) {
	int aa = get_aa_from_symbol(string(TranslationTables[n]).at(c));
	if (string(StartCodons[n]).at(c) != '-'){
	    numStartCodons++;
	    start_codons[c] = true;
	} else {
	    start_codons[c] = false;
	}
	if (aa != map[c]) {
	    // DEBUG MESSAGE
	    cout << "# Warning: Using nonstandard genetic code: " 
		 << Seq2Int(3).INV(c)
		 << " coding for " << (aa < 0 ? "STOP" : GeneticCode::aa_names[aa])
		 << " instead of " << (map[c] < 0 ? "STOP" : GeneticCode::aa_names[map[c]])
		 << ".\n";
	}
	map[c] = aa;
    }
    reverseMap();
    translationtable = n;
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
 * Also stop codons at the end are considered inframe stop codons.
 * If frame is specified, it is as in the gff format.
 * If frame is not specified (-1), then all reading frames are checked.
 */
bool GeneticCode::containsInFrameStopcodon(const char* dna, int begin, int end, bool forward, int frame) {
    bool hasStop = false;
    int i;
    if (forward) {
	if ( frame >= 0 && frame <=2) {
	    for (i=begin + frame; i <= end-2 && !hasStop; i+=3)
		if (isStopcodon(dna+i))
		    hasStop = true;
	} else {
	    for (i=begin; i <= end-2 && !hasStop; i++)
		if (isStopcodon(dna+i))
		    hasStop = true;
	}
    } else { // reverse
	if ( frame >= 0 && frame <=2) {
	    for (i=end - 2 - frame; i >= begin && !hasStop; i-=3)
		if (isRCStopcodon(dna+i))
		    hasStop = true;
	} else {
	    for (i=begin; i <= end-2 && !hasStop; i++)
		if (isRCStopcodon(dna+i))
		    hasStop = true;
	}
    }
    return hasStop;
}

void GeneticCode::printStartCodons(){
    cout << "# admissible start codons and their probabilities: ";
    bool first = true;
    for (int c=0;c<64;c++){
	if (isStartcodon(c)){
	    if (!first)
		cout << ", ";
	    cout <<  Seq2Int(3).INV(c) << "(" << startCodonProb(c) << ")";
	    first = false;
	}
    }
    cout << endl;
}

void GeneticCode::trainStartCodonProbs(int startcounts[]){
    int sum = 0;
    bool first = true;
    if (Constant::augustus_verbosity)
	cout << "start codon frequencies: ";
    for (int c=0; c<64; c++){
	if (startcounts[c]>0){
	    if (Constant::augustus_verbosity){
		if (!first)
		    cout << ", ";
		cout << Seq2Int(3).INV(c) << "(" << startcounts[c] << ")";
	    }
	    sum += startcounts[c];
	    first = false;
	}
    }
    if (Constant::augustus_verbosity)
	cout << endl;

    // estimate probabilities as relative frequencies
    if (sum > 0){
	for (int c=0; c<64; c++){
	    if (start_codons[c]){
		start_codon_probs[c] = Double(startcounts[c]) / sum;
	    }
	}
    } // otherwise leave default (only ATG)
    if (Constant::augustus_verbosity)
	printStartCodons();
}

void GeneticCode::writeStart(ofstream &out){
    int n=0;
    for (int c=0; c<64; c++)
	if (start_codons[c] && start_codon_probs[c] > 0.0)
	    n++;
    out << "# number of start codons:" << endl << n << endl;;
    out << "# start codons and their probabilities" << endl;
    for (int c=0; c<64; c++){
	if (start_codons[c] && start_codon_probs[c] > 0.0){
	    out << Seq2Int(3).INV(c) << "\t" << start_codon_probs[c] << endl;
	}
    }
}

void GeneticCode::readStart(ifstream &in){
    int n;
    string cod;
    Double prob;
    if (!in)
	throw ProjectError("Error when reading in start codon probabilities.");
    in >> comment >> n >> comment;
    for (int i=0; i<n; i++){
	in >> cod >> prob;
	if (prob < 0.0)
	    throw ProjectError("Start codon probability is negative.");
	try {
	    int pn = Seq2Int(3)(cod.c_str());
	    if (pn < 0 || pn > 63)
		throw ProjectError();
	    if (start_codons[pn]){
		start_codon_probs[pn] = prob;
	    } else {
		cerr << cod << " is not a start codon in the chosen translation table "
		     << translationtable << ". Ignoring it." << endl;
	    }
	} catch (...) {
	    throw ProjectError(string("Invalid start codon ")+cod);
	}
    }
    printStartCodons();
}

/*
 * find a longest ORF in a given DNA sequence
 * At the 3' end it is delineated by a stop codon or the sequence end (then complete3prime is false).
 * At the 5' end it is delineated by the most upstream start codon (preferred) or the sequence end
 * (if none). complete5prime is true iff the ORF starts with a start codon and there is an upstream
 * in-frame stop codon.
 */

ORF GeneticCode::longestORF(const char* dna){
    ORF orf, longestOrf;
    unsigned n = strlen(dna);
    Double startPthresh = 0.01; // minimal start codon probability for ORF finding
    for (int si=0; si<=1; si++){ // the more elegant loop does not compile at UCSC:  for (Strand s : {plusstrand, minusstrand})
	Strand s = (Strand) si;
	int dir = (s == plusstrand)? 1 : -1;
	for (unsigned rf = 0; rf < 3; rf++){
	    int from = (s == plusstrand)? rf : n - 3 - rf;
	    orf = ORF();
	    orf.strand = s;
	    orf.complete5prime = false; // first ORF is always 5' incomplete
	    for (int pos = from; pos >= 0 && pos <= n-3; pos += dir * 3){
		if ((!orf.complete5prime && orf.start < 0) || 
		    ((startCodonProb(dna + pos, s == minusstrand) > startPthresh &&
		      ((orf.start >= 0 && startCodonProb(dna + orf.start, s == minusstrand) <= startPthresh)
		       || orf.start < 0)))){
		orf.start = pos; // let ORF start with most 5' ATG (if it exist`s)
	    }
	    if (orf.start >= 0){
		    if ((s == plusstrand && isStopcodon(dna + pos))
			|| (s == minusstrand && isRCStopcodon(dna + pos))){
			orf.end = pos;
			orf.complete3prime = true;
			if (orf.len() > longestOrf.len() && 
			    (startCodonProb(dna + orf.start, s == minusstrand) > startPthresh || !orf.complete5prime)){
			    longestOrf = orf;
			}
			orf = ORF();
			orf.complete5prime = true; // every ORF but the first is complete
			orf.strand = s;
		    } else if (pos + dir*3 >= n || pos + dir*3 < 0) {
			orf.end = pos;
			orf.complete3prime = false;
			if (orf.len() > longestOrf.len()){
			    longestOrf = orf;
			}
		    }
		}
	    }
	}
    }
    if (longestOrf.strand == minusstrand)
	swap(longestOrf.start, longestOrf.end);
    // shift thickEnd +2 as above coordinates are codon positions
    if (longestOrf.start >= 0){
	if (longestOrf.strand == plusstrand){
	    if (longestOrf.complete3prime)
		longestOrf.end += 2;
	    else 
		longestOrf.end = n-1;
	    if (!longestOrf.complete5prime)
		longestOrf.start = 0;
	} else {
	    if (longestOrf.complete5prime || startCodonProb(dna + longestOrf.end, minusstrand) > startPthresh)
		longestOrf.end += 2;
	    else 
		longestOrf.end = n-1;
	     if (!longestOrf.complete3prime)
		longestOrf.start = 0;
	}
    }
    return longestOrf;
}

/*
 * getSampledSeq
 * ****|***|***|***|***|***|
 * k bases from uniform distribution, then numCodons codons
 * k is the order of the Markov chain
 */
char *getSampledCDS(vector<Double> *emiprobs, int k, int numCodons){
    int dnaLength = k + 3*numCodons;
    char bases[5] = "acgt"; // 0->a, 1->c, 2->g, 3->t
    char *cds = new char[dnaLength+1];
    cds[dnaLength] = '\0';
    Seq2Int s2i(k+1);

    // first k bases from uniform distribution
    for (int i=0; i<k; i++) {
        cds[i] = bases[(int) (4.0 * rand() / (1.0 + RAND_MAX))];
    }
    Double sumProb, cumProb;
    double r;
    // sample one codon and one frame position at a time
    for (int i=0; i<numCodons; i++) {
	for (int f=0; f<3; f++) {// reading frame
	    sumProb = 0.0;
	    for (int base = 0; base < 4; base++) {
		cds[k+3*i+f] = bases[base];
		if (!(f==2 && GeneticCode::isStopcodon(cds + k+3*i+f-2)))
		    sumProb += emiprobs[f][s2i(cds + 3*i+f)];
	    }
	    r = rand() / (1.0 + RAND_MAX);
	    cumProb = 0.0;
	    for (int base = 0; base < 4 && cumProb <= r*sumProb; base++) {// try new base at position i, keep that base when cumulative prob. exceeds r
		cds[k+3*i+f] = bases[base];
		if (!(f==2 && GeneticCode::isStopcodon(cds + k+3*i+f-2)))
		    cumProb += emiprobs[f][s2i(cds + 3*i+f)];
	    }
	}	
    }
    return cds;
}
