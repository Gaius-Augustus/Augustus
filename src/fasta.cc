/*
 * fasta.cc
 *
 * License: Artistic License, see file LICENSE.TXT or 
 *          https://opensource.org/licenses/artistic-license-1.0
 * 
 * Description: FASTA sequence input
 */

// project includes
#include "fasta.hh"

// standard C/C++ includes
#include <iostream>
#include <sstream>
#include <string>
#include <string.h>

using namespace std;

bool isFastaCommentLine(const string &line){
    size_t p = line.find_first_not_of(" \t\r\f\v");
    return p != string::npos && line[p] == ';';
}

/*
 * Skip empty lines and comment lines, i.e. lines whose first non-whitespace
 * character is ';' (comments in the original Pearson FASTA format, also used
 * e.g. for FAIR Header Reference genome (FHR) metadata).
 */
static void skipBlankAndCommentLines(istream &is){
    string line;
    while (is >> ws && is.peek() == ';')
	getline(is, line);
}

/*
 * Read the next sequence from the input stream and store the sequence,
 * its name and length in the arguments. Memory for sequence and name are allocated here.
 * Comment lines are skipped anywhere, but a '>' header line is never a comment.
 */
static void readOneFastaSeqImpl(istream &is, char* &sequence, char* &name, long &length){
    string line;
    string seq("");
    static int unnamedcount=1;
    skipBlankAndCommentLines(is);
    if (!is)
	return;
    if (is.peek() == '>') {// sequence name specified
	getline(is, line);
	// go up to the first white space, i.e. interpret only the first word as sequence identifier
	int endpos = 1;
	while (endpos < line.length() && !isspace(line[endpos]))
	    endpos++;
	endpos--;
	name = new char[endpos+1];
	strncpy(name, line.c_str()+1, endpos);
	name[endpos] = '\0';
    } else { // not correct fasta: unnamed sequence
	name = new char[14]; // at most 100000 sequences
	sprintf(name, "unnamed-%d", unnamedcount);
    }
    if (!is)
	return;
    while(is && is.peek( ) != '>'){
	if (getline(is, line) && !isFastaCommentLine(line))
	    seq.append(line);
    }
    sequence = new char[seq.length()+1];

    // now filter out any characters that are not letters
    long pos = 0;
    for (long i=0; i < seq.length(); i++)
	if (isalpha( seq[i] ))
	    sequence[pos++] = seq[i]; // tolower now postponed to after softmasking detection
    sequence[pos] = '\0';
    length = pos;
    if (length == 0){
	delete [] sequence;
	sequence = NULL;
    }
}

#ifdef ZIPINPUT
#include <boost/iostreams/filtering_stream.hpp>
using boost::iostreams::filtering_istream;

void readOneFastaSeq(filtering_istream &ifstrm, char* &sequence, char* &name, long &length){
    readOneFastaSeqImpl(ifstrm, sequence, name, length);
}
#endif

void readOneFastaSeq(std::stringstream &ifstrm, char* &sequence, char* &name, long &length){
    readOneFastaSeqImpl(ifstrm, sequence, name, length);
}

void readOneFastaSeq(std::stringstream &ifstrm, char* &sequence, char* &name, int &length){
    long llen;
    readOneFastaSeq(ifstrm, sequence, name, llen);
    length = (int) llen; // this function should be used for short sequences only, e.g. protein seqs
}

void readOneFastaSeq(ifstream &ifstrm, char* &sequence, char* &name, long &length){
    readOneFastaSeqImpl(ifstrm, sequence, name, length);
}

void readFastaHeader(ifstream &ifstrm, char* &name){
    string line;
    char   c;
    static int unnamedcount=1;
    // skip empty lines
    ifstrm >> ws;
    if (!(ifstrm))
        return;
    c = ifstrm.peek();
    if (c == '>') {// sequence name specified
        getline(ifstrm, line);
        // go up to the first white space, i.e. interpret only the first word as sequence identifier
        int endpos = 1;
        while (endpos < line.length() && !isspace(line[endpos]))
            endpos++;
        endpos--;
        name = new char[endpos+1];
        strncpy(name, line.c_str()+1, endpos);
        name[endpos] = '\0';
    } else { // not correct fasta: unnamed sequence
        name = new char[14]; // at most 100000 sequences
        sprintf(name, "unnamed-%d", unnamedcount);
    }
}

bool isFasta(ifstream &ifstrm){
    ifstrm >> ws;
    if (!(ifstrm))
        return false;
    char c = ifstrm.peek();
    if (c == '>')
        return true;
    return false;
}
