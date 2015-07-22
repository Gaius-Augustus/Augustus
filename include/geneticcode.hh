/*****************************************************************************\
 * Filename : geneticcode.hh
 * Authors  : Mario Stanke
 * 
 *
 * Description: The genetic code maps codons to amino acids
 *
 *
 * Date       |   Author              |  Changes
 *------------|-----------------------|------------------------------------------
 * 15.09.2002 | Stanke                | creation of the class
\******************************************************************************/
#ifndef _GENETIC_CODE_HH
#define _GENETIC_CODE_HH

// project includes
#include "types.hh"
#include "projectio.hh"

// standard C/C++ includes
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>

using namespace std;


#define NUM_AA 20
#define NUM_TRANSTABS 24

/*
 * === am I at a possible splice site or translation start site ===
 */
#define DECLARE_ON(NAME, PATTERN, COUNT)		\
    inline bool NAME(const char* dna) {			\
	return strncmp(dna, PATTERN, COUNT) == 0;	\
    }
// 'gt'
#define DSS_SEQUENCE "gt"
#define RDSS_SEQUENCE "ac"
DECLARE_ON(onDSS,     DSS_SEQUENCE, 2)
DECLARE_ON(onRDSS,    RDSS_SEQUENCE, 2)

// 'gc'
#define ALT_DSS_SEQUENCE "gc"
#define ALT_RDSS_SEQUENCE "gc"
DECLARE_ON(onAltDSS,  ALT_DSS_SEQUENCE, 2)
DECLARE_ON(onAltRDSS, ALT_RDSS_SEQUENCE, 2)

// 'gt' or 'gc'
inline bool onGenDSS(const char* dna) {
    return 
	onDSS(dna) || (Constant::dss_gc_allowed && onAltDSS(dna));
}
inline bool onGenRDSS(const char* dna) {
    return 
	onRDSS(dna) || (Constant::dss_gc_allowed && onAltRDSS(dna));
}

// 'ag'
#define ASS_SEQUENCE "ag"
#define RASS_SEQUENCE "ct"
DECLARE_ON(onASS,     ASS_SEQUENCE, 2)
DECLARE_ON(onRASS,    RASS_SEQUENCE, 2)

// 'atg'
#define STARTCODON "atg"
#define RCSTARTCODON "cat"
DECLARE_ON(onStart,   STARTCODON, 3)
DECLARE_ON(onRStart,  RCSTARTCODON, 3)

// stop codons
DECLARE_ON(onOchre,   OCHRECODON, 3)
DECLARE_ON(onAmber,   AMBERCODON, 3)
DECLARE_ON(onOpal,    OPALCODON, 3)
DECLARE_ON(onROchre,  RCOCHRECODON, 3)
DECLARE_ON(onRAmber,  RCAMBERCODON, 3)
DECLARE_ON(onROpal,   RCOPALCODON, 3)

/*
 * wcComplement
 *
 * Watson/Crick complement of character. Preserves case.
 * Doesn't do anything if character is non-acgt
 *
 */
inline char wcComplement(char c) {
    switch (c) {
        case 'a': return 't';
        case 'c': return 'g';
        case 'g': return 'c';
        case 't': return 'a';
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default:
	    return c;
    }
}


/*
 * putReverseComplement
 *
 * copies the reverse complement of [dna, dna+len) to result
 * result may be a char* or a string iterator
 *
 * no validation is performed and no null character is appended
 */
template <class Iterator>
inline void putReverseComplement(Iterator result, const char* dna, int len) {
    const char* s = dna + len;
    Iterator t = result;
    while (--s >= dna) 
	*t++ = wcComplement(*s);
}

/*
 * reverseComplement
 * 
 * allocates a new character array and fills it with the reverse complement of 'dna'
 * eg agcngt -> acngct
 * doesn't change uppercase or lowercase property
 */
inline char* reverseComplement(const char* dna) {
    if (dna == NULL)
	return NULL;
    int len = strlen(dna);
    char* result = new char[len+1];
    putReverseComplement(result, dna, len);
    result[len] = '\0';
    return result;
}

inline void reverseComplementString(string &text) {
    int n = text.length();
    for (int i=0; i < n/2; i++) {
        char c;
        c = text[i];
        text[i] = wcComplement(text[n-i-1]);
        text[n-i-1] = wcComplement(c);
    }
    if (n%2) // n odd, must complement the middle character
	text[n/2] = wcComplement(text[n/2]);
}


inline void reverseString(string &text) {
    int i = 0;
    int n = text.length();
    while (i < (n/2)) {
        char c;
        c = (text[i]);
        text[i] = text[n-i-1];
        text[n-i-1] = c;
        i++;
    }
}



/*
 * a class for converting sequence into integer
 * replacing Base4Int
 */
class Seq2Int {
public:
    Seq2Int(int s) : size(s) {}
    int operator() (const char* s) const {
	int erg=0;
	for( int i = 0; i < size; i++ ){
	    erg <<= 2;
	    erg |= base2int(s[i]);
	}
	return erg;
    }
    int rc(const char* s) const {
	int erg=0;
	for ( int i=0; i<size; i++) 
	    erg |= base2int(wcComplement(s[i])) << (2*i);
	return erg;
    }
    int rev(const char* s) const {
	int erg=0;
	for ( int i=0; i<size; i++) 
	    erg |= base2int(s[i]) << (2*i);
	return erg;
    }
    string inv(int pn) const {
	string result(size, 'n');
	for(int i = size-1; i>=0; i--) {
	    result[i] = int2base(pn%4);
	    pn >>= 2; // == k /= 4;
	}
	return result;
    }
    string INV(int pn) const {
	string result(size, 'N');
	for(int i = size-1; i>=0; i--) {
	    result[i] = int2BASE(pn%4);
	    pn >>= 2; // == k /= 4;
	}
	return result;
    }
 	
    int read(istream& strm) const {
	char* buf = new char[size];
	// Read the character from the stream until the internal
	// representation size is reached, or a character
	// cannot be recognised.
	//------------------------------------------------------------
	for( int i = 0; i < size; i++ ){
	    char c;
	    strm.get(c);
	    c = toupper(c);
	    if( c == 'A' || c == 'C' || c == 'G' || c == 'T' )
		buf[i] = c;
	}
	int result = (*this)(buf);
	delete[] buf;
	return result;
    }

private:
    static int base2int(char c);
    static char int2base(int i);
    static char int2BASE(int i);
    int size;
};

inline int Seq2Int::base2int(char c) {
    switch (c) {
        case 'a': case 'A':
            return 0;
        case 'c': case 'C':
            return 1;
        case 'g': case 'G':
            return 2;
        case 't': case 'T':
            return 3;
        default:
	    throw InvalidNucleotideError(c);
    }
}

inline char Seq2Int::int2base(int i) {
    switch (i) {
        case 0:
            return 'a';
        case 1:
            return 'c';
        case 2:
            return 'g';
        case 3:
            return 't';
        default:
            throw ProjectError("Seq2Int::int2base: internal error: i=" + itoa(i));
    }
}

inline char Seq2Int::int2BASE(int i) {
    switch (i) {
        case 0:
            return 'A';
        case 1:
            return 'C';
        case 2:
            return 'G';
        case 3:
            return 'T';
        default:
            throw ProjectError("Seq2Int::int2base: internal error: i=" + itoa(i));
    }
}

class ORF {
public:
    ORF() { start = end = -1; complete5prime = complete3prime = 0; strand = plusstrand; }
    int len() {return abs(end - start) + 1;}
    int start;
    int end;
    Strand strand;
    bool complete5prime;
    bool complete3prime;
};

class GeneticCode {
public:
     static void init() {
	chooseTranslationTable(1);
    }
    static void printReverseGeneticMap();
    
private:
    GeneticCode() {}
    static void reverseMap();
    static const char aa_symbols_with_stop[];
    static bool start_codons[];
    static Double start_codon_probs[];
    static const char * const TranslationTables[];
    static const char * const StartCodons[];
    static Seq2Int codon;
    static int translationtable;
    static int numStartCodons;
public:
    static void chooseTranslationTable( int );
    static const char* const aa_symbols;
    static const char* const aa_names[];
    static int map[];
    static int **syncodons;
    static int *codonsOfAA;
    static int get_aa_from_symbol(char c) {
	// this one returns -1 on '*' and string::npos-1 on invalid symbols
	// if (int)string::npos is -1, we get -2 here for invalid
	return string(GeneticCode::aa_symbols_with_stop).find(c)-1;
    }
    static char translate(int n) {
	return aa_symbols_with_stop[map[n]+1];
    }
    static char translate(const char* t);
    static char revtranslate(const char* t);
    static bool isStopcodon(const char* t) {
	return translate(t)=='*';
    }
    static bool isStartcodon(const char* t, bool rc=false){
	try {
	    if (rc)
		return start_codons[codon.rc(t)];	
	    else
		return start_codons[codon(t)];
	} catch (InvalidNucleotideError){}
	return false;
    }
    static bool isStartcodon(int pn){
	return start_codons[pn];
    }
    static Double startCodonProb(const char* t, bool rc=false){
	try {
	    int pn = rc? codon.rc(t) : codon(t);
	    if (start_codons[pn])
		return start_codon_probs[pn];
	} catch (InvalidNucleotideError){}
	return 0;
    }
    static Double startCodonProb(int pn){
	return start_codon_probs[pn];
    }
    static bool isRCStopcodon(const char* t) {
	return revtranslate(t)=='*';
    }
    static bool containsInFrameStopcodon(const char*, int, int, bool, int);
    static void printStartCodons();
    static void trainStartCodonProbs(int startcounts[]);
    static void writeStart(ofstream &out); // write start codon probs to file
    static void readStart(ifstream &in); // read start codon probs from file
    static bool is_purine(int b){
	return (b==0 || b==2); // 0=a, 1=c, 2=g, 3=t, a and g are purines
    }
    static ORF longestORF(const char* dna);
};

// getSampledCDS: sample a coding region from emission probabilities
char *getSampledCDS(vector<Double> *emiprobs, int k, int numCodons);

#endif    //  _GENETIC_CODE_HH

