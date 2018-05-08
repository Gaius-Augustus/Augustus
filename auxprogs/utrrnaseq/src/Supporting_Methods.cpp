/**
 * \file Supporting_Methods.cpp
 */

#include "Global.hpp"

#include <vector>
#include <string>
#include <iterator>

#include <fstream>
#include <iostream>

#include <getopt.h>

#include <boost/tokenizer.hpp>
#include <boost/assign/std/vector.hpp>

using namespace std;
using namespace boost;
using namespace boost::assign;



vector<string> tokenize(string token, string str) {
	vector<string> tokens;
	char_separator<char> sep(token.c_str()); //without c_str() not working, changes string into c-string
	tokenizer< char_separator<char> > tok(str, sep);

	for (tokenizer< char_separator<char> >::iterator i = tok.begin(); i != tok.end(); ++i) {
		tokens += *i;
	}
	return tokens;
}



void reset_getopt() {
	#ifdef __GLIBC__
		optind = 0;
	#else
		optind = 1;
	#endif
	#ifdef HAVE_OPTRESET
		optreset = 1;
	#endif
}



bool files_identical (string first_file, string second_file) {
	ifstream first_ifstream( first_file.c_str() );
	ifstream second_ifstream( second_file.c_str() );

	if (!first_ifstream) {
		ERR_STRM << "Error (in Test Mode): Could not open '" << first_file << "'!" << endl;
		cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
		abort();
	}

	if (!second_ifstream) {
		ERR_STRM << "Error (in Test Mode): Could not open '" << second_file << "'!" << endl;
		cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
		abort();
	}

	string str;
	vector <string> first;
	vector <string> second;

	while ( getline(first_ifstream, str) ) {
		first += str;
	}
	while ( getline(second_ifstream, str) ) {
		second += str;
	}

	if ( first.size() != second.size()) {
		return false;
	}
	for ( unsigned i = 0; i < first.size(); i++) {
		if ( first[i] != second[i] ) {
			return false;
		}
	}
	return true;
}
