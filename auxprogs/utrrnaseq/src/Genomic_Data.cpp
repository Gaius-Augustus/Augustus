 /*
 * \file Genomic_Data.cpp
 */

#include "Genomic_Data.hpp"
#include "Global.hpp"
#include "Supporting_Methods.hpp"
#include "Splice_Sites.hpp"
#include "Error.hpp"

#include <string>
#include <vector>
#include <utility>

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <iterator>

#include <boost/assign/std/vector.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

std::vector<Genomic_Data::Scaff_plus_Gen_Data> Genomic_Data::s_all_genes;
std::vector< std::pair <std::string, std::string> > Genomic_Data::s_sp_sites;

using namespace std;
using namespace boost;
using namespace boost::assign;



// Determine the start and end positions of runs (i.e. stretches) of lower case letters
vector< pair<unsigned, unsigned> > lower_case_runs(string str) {
	// Iterate through the string, storing runs on the fly
	vector< pair<unsigned, unsigned> > ret;
	unsigned i = 0;
	bool prev_lc = false;   // Was the previous character lower case?
	unsigned start;
	while (i < str.length()) {
		bool curr_lc = islower( char( str[i] ) );   // Is the current character lower case
		if (prev_lc) {
			if (!curr_lc) {
				unsigned end = i - 1;
				ret += make_pair(start, end);
				prev_lc = false;
			}
		}
		else {
			if (curr_lc) {
				start = i;
				prev_lc = true;
			}
		}
		++i;
	}
	if (prev_lc) {   // Is the last character of the string  lower case?
		unsigned end = i - 1;
		ret += make_pair(start, end);
	}

	return ret;
}



void Genomic_Data::initialize(string scaff_fname, string crb_fname, string intron_fname, string repeat_fname,
		vector< pair<string, string> > sp_sites, bool use_repeat_file) {

	s_all_genes.clear();
	s_sp_sites = sp_sites;

	read_scaff_file(scaff_fname);
	read_crb_file(crb_fname);
	read_intron_file(intron_fname);
	if (use_repeat_file)
		read_repeat_file(repeat_fname);

	cout << "Input Data processing finished successfully!" << endl;
	//cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
}



void Genomic_Data::read_scaff_file(string scaff_fname) {
	ifstream scaff_ifstream(scaff_fname.c_str());

	if (!scaff_ifstream) {
		ERR_STRM << "Error in scaffold file: Could not open '" << scaff_fname << "'!" << endl;
		cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
		abort();
	}

	Scaff_plus_Gen_Data seq;
	string str;
	string gene;
	string scaffold_name;

	bool first_scaff = true;
	bool corrupted_scaff = false;
	unsigned curr_scaff_row_idx = 0;

	while (getline(scaff_ifstream, str)) {
		++curr_scaff_row_idx;

		if (str.empty()) {
			continue;
		}

		if (str[0] == ';') {
			continue; //comment
		}

		int position = str.find(">");
		switch (position) {
			case 0:
				scaffold_name = str.substr(1);	//line only includes ">Name", uses string from Position 1 until the end of scaffold name

				if (first_scaff) {
					seq.name = scaffold_name;
					first_scaff = false;
				}
				else {
					if (gene.empty()) {
						ERR_STRM << "Error in scaffold file (line " << curr_scaff_row_idx - 1
							     << "): Empty scaffold found in '" << seq.name << "'!" << endl;
						cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
						seq.name = scaffold_name;
					}
					else {
						if (!corrupted_scaff) {
							// One whole scaffold finished -> saved and new scaffold started -> scaffold name saved
							seq.sequence = gene;
							s_all_genes += seq;

							gene.clear();
							seq.name = scaffold_name;
						}
						else {
							corrupted_scaff = false;	//scaffold is invalid, it will not be saved
							gene.clear();
							seq.name = scaffold_name;
						}
					}
				}
				break;

			case -1:	//> not found, scaffold in this line
				gene += str;
				break;

			default:	//> in position other than the first -> Error
				if (gene.empty()) {
					ERR_STRM << "Error in scaffold file (line " << curr_scaff_row_idx - 1
						     << "): Empty scaffold found in '" << seq.name << "'!" << endl;
					cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
				}
				else {
					if (!corrupted_scaff) {
						seq.sequence = gene;
						s_all_genes += seq;

						gene.clear();
					}
					else {
						corrupted_scaff = false;	//scaffold is invalid, it will not be saved
						gene.clear();
					}
				}

				ERR_STRM << "Error in scaffold file (line " << curr_scaff_row_idx << "): > in '"
					     << str << "' is not first character in its line!" << endl;
				cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
				corrupted_scaff = true;
				//continue;
		}
	}

	if (gene.empty()) {
		ERR_STRM << "Error in scaffold file (line " << curr_scaff_row_idx - 1
		     	 << "): Empty scaffold found in '" << seq.name << "'!" << endl;
		cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
	}
	else {
		if (!corrupted_scaff) {
			vector< pair<unsigned, unsigned> > v = lower_case_runs(gene);
			typedef pair<unsigned, unsigned> pr_type;
			BOOST_FOREACH(pr_type pr, v) {
				Repeat r = { pr.first + 1, pr.second + 1 };   // Strings are 0-based, gene coordinates are 1-based
				seq.repeats += r;
			}
			seq.sequence = gene;
			s_all_genes += seq;

			gene.clear();
			seq.name = scaffold_name;
		}
	}
	//else unnecessary, because it was the last scaffold in the file and this last gene was corrupted

	scaff_ifstream.close();
	cout << "Read in of scaffold file finished successfully!" << endl;
}


void Genomic_Data::read_crb_file(string crb_fname) {
	ifstream crb_ifstream(crb_fname.c_str());

	if (!crb_ifstream) {
		ERR_STRM << "Error (in coding region file): Could not open '" << crb_fname << "'!" << endl;
		abort();
	}

	//Fields of the GFF Format
	const unsigned COLUMN_NO = 9;
	const unsigned SEQNAME_IDX = 0;
	const unsigned FEATURE_IDX = 2;
	const unsigned STOP_CODON_IDX = 4;
	const unsigned STRAND_IDX = 6;
	const unsigned GROUP_IDX = 8;	//including transcript_id and gene_id

	string str;
	CRB curr_crb_row;
	unsigned idx_curr_crb_row = 0;

	while (getline(crb_ifstream, str)) {
		++idx_curr_crb_row;

		if (str.empty()) {
			continue;
		}

		if ( (str[0] == '#') && (str[1] == '#') ) {
			continue; //comment
		}

		vector<string> tokens = tokenize("\t", str);
		if (tokens.size() != COLUMN_NO) {
			ERR_STRM << "Error in coding region file (line " << idx_curr_crb_row << "): '"
				     << str << "' does not have the right number of columns!" << endl;
			cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
			continue;
		}

		//=======SEQ_NAME=======
		bool seqname_in_scaff_file = false;
		unsigned scaff_idx_curr_crb;
		unsigned idx_curr_scaff = 0;

		while ( ( !seqname_in_scaff_file ) && ( idx_curr_scaff < s_all_genes.size() ) ) {
			const Scaff_plus_Gen_Data& curr_scaff = s_all_genes[idx_curr_scaff];

			if (curr_scaff.name == tokens[SEQNAME_IDX] ) {
				seqname_in_scaff_file = true;
				scaff_idx_curr_crb = idx_curr_scaff;
			}
			++idx_curr_scaff;
		}

		if (!seqname_in_scaff_file) {
			ERR_STRM << "Error in coding region file (line " << idx_curr_crb_row
				     << "): '" << str << "' contains a sequence name, which is not in the scaffold file!" << endl;
			cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
			continue;
		}

		//======FEATURE========
		if (tokens[FEATURE_IDX] == "start_codon") {
			curr_crb_row.feature = "start";
		}
		else {
			if (tokens[FEATURE_IDX] == "stop_codon") {
				curr_crb_row.feature = "stop";
			}
			else {
				ERR_STRM << "Error in coding region file (line " << idx_curr_crb_row << "): In '"
					     << str << "' only start_codon or stop_codon are allowed as features!" << endl;
				cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
				continue;	//Feature is not "stop_codon" nor "start_codon"
			}
		}

		//======CODON/STRAND========
		try {
			curr_crb_row.codon_pos =  lexical_cast<unsigned>( tokens[STOP_CODON_IDX] );
		}
		catch (bad_lexical_cast const&) {
			ERR_STRM << "Error in coding region file (line " << idx_curr_crb_row << "): '"
				     << str << "' contains a codon position which is not an integer!" << endl;
			cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
			continue;
		}
		//could also use start_codon => tokens[3]

		if ( (tokens[STRAND_IDX] != "+" ) && ( tokens[STRAND_IDX] != "-" ) ) {
			ERR_STRM << "Error in coding region file (line " << idx_curr_crb_row << "): '"
				     << str << "' contains a strand other than '+' and '-'" << endl;
			cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
			continue;
		}

		curr_crb_row.strand = tokens[STRAND_IDX];
		if ( ( (curr_crb_row.strand == "+") && (curr_crb_row.feature == "start") ) || ( (curr_crb_row.strand == "-") && (curr_crb_row.feature == "stop") ) ) {
			curr_crb_row.codon_pos -= 2;
		}

		//=============GROUP===========
		curr_crb_row.grp_att = tokens[GROUP_IDX];
		s_all_genes[scaff_idx_curr_crb].crbs += curr_crb_row;
	}

	crb_ifstream.close();
	cout << "Read in of coding region file finished successfully!" << endl;
}


void Genomic_Data::read_intron_file(string intron_fname) {

	ifstream intron_ifstream(intron_fname.c_str());

	if (!intron_ifstream) {
		ERR_STRM << "Error in intron file: Could not open '" << intron_fname << "'!" << endl;
		cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
		abort();
	}

	//necessary fields of GFF Format
	const unsigned COLUMN_NO = 9;
	const unsigned SEQNAME_IDX = 0;
	const unsigned FEATURE_IDX = 2;
	const unsigned START_IDX = 3;
	const unsigned END_IDX = 4;
	const unsigned STRAND_IDX = 6;
	const unsigned MULT_IDX = 8;

	string str;
	Intron curr_intron_row;
	unsigned curr_intron_row_idx = 0;
	bool intron_accepted;

	while (getline(intron_ifstream, str)) {
		++curr_intron_row_idx;
		intron_accepted = true;

		if (str.empty()) {
			continue;
		}

		if ( (str[0] == '#') && (str[1] == '#')) {
			continue; //comment
		}

		vector<string> tokens = tokenize("\t", str);
		if (tokens.size() != COLUMN_NO) {
			ERR_STRM << "Error in intron file (line " << curr_intron_row_idx << "): '"
				     << str << "' does not have the right number of columns!" << endl;
			cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
			continue;
		}

		//=========SEQ_NAME=============
		bool seqname_in_scaff_file = false;
		unsigned scaff_idx_curr_intron;
		unsigned idx_curr_scaff = 0;

		while ( ( !seqname_in_scaff_file ) && ( idx_curr_scaff < s_all_genes.size() ) ) {
			const Scaff_plus_Gen_Data& curr_scaff = s_all_genes[idx_curr_scaff];

			if (curr_scaff.name == tokens[SEQNAME_IDX] ) {
				seqname_in_scaff_file = true;
				scaff_idx_curr_intron = idx_curr_scaff;
			}
			++idx_curr_scaff;
		}

		if (!seqname_in_scaff_file) {
			ERR_STRM << "Error in intron file (line " << curr_intron_row_idx << "): '"
				     << str << "' contains a sequence name, which is not in the scaffold file!" << endl;
			cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
			continue;
		}

		//======FEATURE (only test) =======
		if (tokens[FEATURE_IDX] != "intron") {
			ERR_STRM << "Error in intron file (line " << curr_intron_row_idx
				     << "): Feature in line '" << str << "' does not have the feature type 'intron'!" << endl;
			cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
			continue;
		}

		//=======START/END/STRAND==========
		unsigned start;
		unsigned end;
		try {
			start = lexical_cast<unsigned>(tokens[START_IDX]);
			end = lexical_cast<unsigned>(tokens[END_IDX]);
		}
		catch (bad_lexical_cast const&) {
			ERR_STRM << "Error in intron file (line " << curr_intron_row_idx
				     << "): '" << str << "' contains a coordinate which is not an integer!" << endl;
			cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
			continue;
		}

		if ( (tokens[STRAND_IDX] != "+" ) && ( tokens[STRAND_IDX] != "-" ) && ( tokens[STRAND_IDX] != "." ) ) {
			ERR_STRM << "Error in intron file (line " << curr_intron_row_idx << "): '"
				     << str << "' contains a strand other than '+', '-' and '.' ." << endl;
			cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
			continue;
		}

		pair<bool, string> intron_valid;
		bool strand_defined = tokens[STRAND_IDX] != "." ? true : false;

		if (!s_sp_sites.empty() && !strand_defined) { //if empty no splice site filtering, intron_accepted always true
			intron_valid = check_sp_sites(s_sp_sites, s_all_genes[scaff_idx_curr_intron].sequence,
							start, end, strand_defined, tokens[STRAND_IDX] );
		}
		else {
			if (s_sp_sites.empty() && !strand_defined) {
				ERR_STRM << "Error in intron file (line " << curr_intron_row_idx << "): '"
					     << str << "' contains an undefined strand '.'. In this case at least one splice site "
					     << "must be given to check for splice sites and define the strand." << endl;
				cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
				throw Error();
			}
			else {
				intron_valid.first = true;
			}
		}

		if (!intron_valid.first) {
			continue; //intron not accepted
		}

		if (tokens[STRAND_IDX] == ".")  {
			curr_intron_row.strand = intron_valid.second;
		}
		else {
			curr_intron_row.strand = tokens[STRAND_IDX];
		}

		if (curr_intron_row.strand == "+") {
			curr_intron_row.start = start;
			curr_intron_row.end = end;
		}
		else { //minus strand
			curr_intron_row.start = end;
			curr_intron_row.end = start;
		}

		//============MULT===================
		string mult_str = tokens[MULT_IDX];
		string::size_type pos_grp, pos_group;
		pos_grp = mult_str.find("grp=");
		pos_group = mult_str.find("group=");

		string mult_value_str;

		if ( (pos_grp != mult_str.npos) || (pos_group != mult_str.npos) ) {
			mult_value_str = "1";
		}
		else {
			unsigned pos = mult_str.find("mult=");
			unsigned mult_length = 0;
			unsigned mult_pos = pos + 5;
			bool mult_found = false;

			while ( (!mult_found) && ( (mult_pos + mult_length) < mult_str.length()) ) {
				char tested_char = mult_str[mult_pos + mult_length];
				if ( isdigit(tested_char) ) {
					mult_length++;
				}
				else {
					mult_found = true;
				}
			}
			mult_value_str = mult_str.substr(mult_pos, mult_length);
		}

		try {
			curr_intron_row.mult =  lexical_cast<unsigned>( mult_value_str );
		}
		catch (bad_lexical_cast const&) {
			ERR_STRM << "Error in intron file (line " << curr_intron_row_idx << "): '"
				     << str << "' contains a multiplicity which is not an integer!" << endl;
			cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
			continue;
		}

		if (intron_accepted) {
			s_all_genes[scaff_idx_curr_intron].introns += curr_intron_row;
		}
		//else intron not accepted, because of a wrong or unaccepted splice site
	}

	intron_ifstream.close();
	cout << "Read in of intron file finished successfully!" << endl;
}


void Genomic_Data::read_repeat_file(string repeat_fname) {
	ifstream repeat_ifstream( repeat_fname.c_str() );

	if (!repeat_ifstream) {
		ERR_STRM << "Error in repeat file: Could not open '" << repeat_fname << "'!" << endl;
		cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
		abort();
	}

	//necessary fields of GFF Format
	const unsigned COLUMN_NO = 9;
	const unsigned SEQNAME_IDX = 0;
	const unsigned FEATURE_IDX = 2;
	const unsigned START_IDX = 3;
	const unsigned END_IDX = 4;


	string str;
	Repeat curr_repeat_row;
	unsigned curr_repeat_row_idx = 0;

	while (getline(repeat_ifstream, str)) {
		++curr_repeat_row_idx;

		if (str.empty()) {
			continue;
		}

		if ( (str[0] == '#') && (str[1] == '#')) {
			continue; //comment
		}

		vector<string> tokens =	tokenize("\t", str);
		if (tokens.size() != COLUMN_NO) {
			ERR_STRM << "Error in repeat file (line " << curr_repeat_row_idx << "): '"
				     << str << "' does not have the right number of columns!" << endl;
			cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
			continue;
		}

		//=========SEQ_NAME=============
		bool seqname_in_scaff_file = false;
		unsigned scaff_idx_curr_repeat;
		unsigned idx_curr_scaff = 0;

		while ( ( !seqname_in_scaff_file ) && ( idx_curr_scaff < s_all_genes.size() ) ) {
			const Scaff_plus_Gen_Data& curr_scaff = s_all_genes[idx_curr_scaff];

			if (curr_scaff.name == tokens[SEQNAME_IDX] ) {
				seqname_in_scaff_file = true;
				scaff_idx_curr_repeat = idx_curr_scaff;
			}
			++idx_curr_scaff;
		}

		if (!seqname_in_scaff_file) {
			ERR_STRM << "Error in repeat file (line " << curr_repeat_row_idx << "): '"
				     << str << "' contains a sequence name, which is not in the scaffold file!" << endl;
			cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
			continue;
		}

		//======FEATURE (only test) =======
		if (tokens[FEATURE_IDX] != "nonexonpart") {
			ERR_STRM << "Error in repeat (line " << curr_repeat_row_idx << "): Feature in line '"
				     << str << "' is not 'nonexonpart'!" << endl;
			cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
			continue;
		}

		//===============START/END===================
		try {
			curr_repeat_row.start = lexical_cast<unsigned>(tokens[START_IDX]);
			curr_repeat_row.end = lexical_cast<unsigned>(tokens[END_IDX]);
		}
		catch (bad_lexical_cast const&) {
			ERR_STRM << "Error in repeat file (line " << curr_repeat_row_idx
				     << "): '" << str << "' contains a coordinate which is not an integer!" << endl;
			cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
			continue;
		}

		s_all_genes[scaff_idx_curr_repeat].repeats += curr_repeat_row;
	}

	repeat_ifstream.close();
	cout << "Read in of repeat file finished successfully!" << endl;
}



bool Genomic_Data::Intron::is_overlapping(Intron i, Intron j) {
  int i_first = 
	  (i.strand == "+") ? i.start : i.end;
  int i_second =
	  (i.strand == "+") ? i.end : i.start;
  int j_first =
	  (j.strand == "+") ? j.start : j.end;
  int j_second = 
	  (j.strand == "+") ? j.end : j.start;       													  
  bool first_opt = (j_first <= i_second) && (j_first >= i_first);
  bool second_opt = (j_second <= i_second) && (j_second >= i_first);
  bool third_opt = (j_first <= i_first) && (j_second >= i_second);
  // introns are not allowed to be immediate neighbours without spacer (this is not an overlap)
  bool fourth_opt = (i_second == (j_first - 1)) || (j_second == (i_first - 1)); 
  return first_opt || second_opt || third_opt || fourth_opt;
}



bool Genomic_Data::Repeat::operator==(const Repeat &r) const {
	return (start == r.start) && (end == r.end);
}



ostream& operator<<(ostream& os, const Genomic_Data::Repeat& r) {
	os << "[" << r.start << "," << r.end << "]";
	return os;
}
