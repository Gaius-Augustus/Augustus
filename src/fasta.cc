/**********************************************************************
 * file:    fasta.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  FASTA sequence input
 * authors: Mario Stanke, mario.stanke@uni-greifswald.de
 *
 * date     |   author      |  changes
 * ---------|---------------|------------------------------------------
 * 01.09.12 | Mario Stanke  | creation of the file
 **********************************************************************/

// project includes
#include "fasta.hh"

// standard C/C++ includes
#include <iostream>
#include <sstream>
#include <string>
#include <string.h>

using namespace std;

#ifdef ZIPINPUT
#include <boost/iostreams/filtering_stream.hpp>
using boost::iostreams::filtering_istream;

/*
 * Read the next sequence from input file stream and store the sequence,
 * its name and length in the arguments. Memory for sequence and name are allocated here.
 * 
 */
void readOneFastaSeq(filtering_istream &ifstrm, char* &sequence, char* &name, int &length){
    string line;
    string seq("");
    char   c;
    static int unnamedcount=1;
    // skip empty lines
    ifstrm >> ws;
    if (!ifstrm)
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
    if (!ifstrm)
	return;
    while(ifstrm && ifstrm.peek( ) != '>'){
	if (getline(ifstrm, line))
	    seq.append(line);
    }
    sequence = new char[seq.length()+1];
     
    // now filter out any characters that are not letters
    int pos = 0;
    for (int i=0; i < seq.length(); i++) 
	if (isalpha( seq[i] ))
	    sequence[pos++] = seq[i]; // tolower now postponed to after softmasking detection
    sequence[pos] = '\0';
    length = pos;
    if (length == 0){
	delete sequence;
	sequence = NULL;
    }
}
#endif

/*
 * This is an exact copy of above.
 * For some reason a template solution does not work (because of boost class structure?).
 * However, filtering_istream is supposed to generalize basic_istream.
 */
void readOneFastaSeq(std::stringstream &ifstrm, char* &sequence, char* &name, int &length){
    string line;
    string seq("");
    char   c;
    static int unnamedcount=1;
    // skip empty lines
    ifstrm >> ws;
    if (!ifstrm)
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
    if (!ifstrm)
	return;
    while(ifstrm && ifstrm.peek( ) != '>'){
	if (getline(ifstrm, line))
	    seq.append(line);
    }
    sequence = new char[seq.length()+1];
     
    // now filter out any characters that are not letters
    int pos = 0;
    for (int i=0; i < seq.length(); i++) 
	if (isalpha( seq[i] ))
	    sequence[pos++] = seq[i]; // tolower now postponed to after softmasking detection
    sequence[pos] = '\0';
    length = pos;
    if (length == 0){
	delete sequence;
	sequence = NULL;
    }
}

void readOneFastaSeq(ifstream &ifstrm, char* &sequence, char* &name, int &length){
    string line;
    string seq("");
    readFastaHeader(ifstrm,name);
    if (!ifstrm)
        return;
    while(ifstrm && ifstrm.peek( ) != '>'){
        if (getline(ifstrm, line))
            seq.append(line);
    }
    sequence = new char[seq.length()+1];
     
    // now filter out any characters that are not letters
    int pos = 0;
    for (int i=0; i < seq.length(); i++) 
        if (isalpha( seq[i] ))
            sequence[pos++] = seq[i];// tolower now postponed to after softmasking detection
    sequence[pos] = '\0';
    length = pos;
    if (length == 0){
        delete sequence;
        sequence = NULL;
    }
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
