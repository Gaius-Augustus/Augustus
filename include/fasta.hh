/*
 * fasta.hh
 *
 * License: Artistic License, see file LICENSE.TXT or 
 *          https://opensource.org/licenses/artistic-license-1.0
 */

#ifndef _FASTA_HH
#define _FASTA_HH

#include <iostream>
#include <fstream>

#ifdef ZIPINPUT
#include <boost/iostreams/filtering_stream.hpp>
using boost::iostreams::filtering_istream;
void readOneFastaSeq(filtering_istream &ifstrm, char* &sequence, char* &name, int &length);
#endif

void readOneFastaSeq(std::stringstream &ifstrm, char* &sequence, char* &name, int &length);
void readOneFastaSeq(std::ifstream &ifstrm, char* &sequence, char* &name, int &length);
void readFastaHeader(std::ifstream &ifstrm, char* &name);
bool isFasta(std::ifstream &ifstrm);

#endif   //  _FASTA_HH
