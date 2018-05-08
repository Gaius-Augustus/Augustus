/**
 * \file Compute_UTRs.cpp
 */

#include "Genomic_Data.hpp"
#include "Supporting_Methods.hpp"
#include "Error.hpp"
#include "Compute_UTRs.hpp"
#include "UTRs.hpp"
#include "Global.hpp"

#include <map>
#include <vector>
#include <string>

#include <iostream>
#include <fstream>
#include <sstream>

#include <stdlib.h>
#include <stdio.h>
#include <iterator>

#include <boost/lexical_cast.hpp>
#include <boost/tuple/tuple.hpp>

using namespace std;
using namespace boost;



unsigned get_idx_curr_scaff(string block_scaff_name) {
	bool is_in_scaff_file = false;
	unsigned idx_former_scaff;

	for (unsigned i = 0; i < Genomic_Data::get_all_genes().size(); ++i) {
		const Genomic_Data::Scaff_plus_Gen_Data& test_if_in_scaff_file = Genomic_Data::get_all_genes()[i];

		if (test_if_in_scaff_file.name == block_scaff_name) {
				is_in_scaff_file = true;
				idx_former_scaff = i; //
				break;	//Assumption: Only one scaffold with same name
		}
	}

	if (is_in_scaff_file) {
		return idx_former_scaff;
	}

	ERR_STRM << "Error in wiggle file! Scaffold '" << block_scaff_name
			 << "' in wiggle file does not have a sequence in the fasta file! "
		 << "Further information regarding line oder row not available!" << endl;
	cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
			 //the row index or the line itself are no longer available
	throw Error();
}



string process_variable_step(const vector<string>& tokens, unsigned idx_curr_row, string curr_row) {
	const unsigned COLUMN_NO = 3;
	const unsigned SCAFF_NAME_IDX = 1;

	if (tokens.size() == COLUMN_NO) {
		ERR_STRM << "Error in wiggle file (line " << idx_curr_row << "): '"
			 << curr_row << "' includes an optional span=WindowSize!" << endl;
		cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
		throw Error();
	}

	unsigned pos = tokens[SCAFF_NAME_IDX].find("=");
	string row_scaff_name = tokens[SCAFF_NAME_IDX].substr(pos+1);

	return row_scaff_name;
}



boost::tuple<string, unsigned, unsigned> process_fixed_step(const vector<string>& tokens, unsigned idx_curr_wiggle_row, string curr_row) {
	const unsigned COLUMN_NO = 5;
	const unsigned SCAFF_NAME_IDX = 1;
	const unsigned START_IDX = 2;
	const unsigned STEP_IDX = 3;

	if (tokens.size() == COLUMN_NO) {
		ERR_STRM << "Error in wiggle file (line " << idx_curr_wiggle_row << "): '"
			 << curr_row << "' includes an optional span=WindowSize!" << endl;
		cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
		throw Error(); //have to find a way to do this without exit if possible
	}

	unsigned pos = tokens[SCAFF_NAME_IDX].find("=");
	string row_scaff_name = tokens[SCAFF_NAME_IDX].substr(pos+1);

	unsigned scaff_pos;
	unsigned step;

	pos = tokens[START_IDX].find("=");
	string start_pos_str = tokens[START_IDX].substr (pos+1);

	try {
		scaff_pos = lexical_cast<unsigned>(start_pos_str);
	}
	catch (bad_lexical_cast const&) {
		ERR_STRM << "Error in wiggle file (line " << idx_curr_wiggle_row
			 << "): Scaffold Position in '" << curr_row << "' is not an integer!" << endl;
		cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
		throw Error();
	}

	pos = tokens[STEP_IDX].find("=");
	string curr_step_str = tokens[STEP_IDX].substr (pos+1);

	try {
		step =  lexical_cast<unsigned>(curr_step_str);
	}
	catch (bad_lexical_cast const&) {
		ERR_STRM << "Error in wiggle file (line " << idx_curr_wiggle_row
			 << "): Step in '" << curr_row << "' is not an integer!" << endl;
		cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
		throw Error();
	}

	return boost::make_tuple(row_scaff_name, scaff_pos, step);
}



void compute_UTRs(string wiggle_fname) {
	ifstream wiggle_ifstream( wiggle_fname.c_str() );

	if (!wiggle_ifstream) {
	  ERR_STRM << "Error in wiggle file: Could not open '" << wiggle_fname << "'!" << endl;
	  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
		abort();
	}

	//only declaration line
	const unsigned OPTION_IDX = 0; 	//tokens[0] does not have to be an option, could also be some wiggle data (data line),
									//but if the line is a declaration line than this option is in tokens[0]

	unsigned chrom_pos;	//position for fixedStep only
	unsigned step;	//Step for fixedStep only

	string curr_row;
	map<unsigned,double> wiggle_data;	//Key: Position, Value: Coverage
	string option = "NULL";

	unsigned idx_curr_row = 0;
	unsigned idx_former_scaff;
	string block_scaff_name; //A blocks ends with a declaration line with a different scaffold name.
	bool is_first_scaffold = true;

	//extracting and analyzing every line in the wiggle file, e.g. testing if declaration or data line, and
	//then processing the information to save all the coverage data independent from the wiggle format
	while ( getline(wiggle_ifstream, curr_row) ) {
		++idx_curr_row;

		if (curr_row.empty()) {
			continue;	//empty line, I don't think that this is worth an error message
		}

		if ( curr_row[0] == '#') {
			continue; //comment, so this line will be ignored
		}

		vector<string> tokens = tokenize(" ", curr_row);

		char test_char = curr_row.at(0);
		bool is_data_line = isdigit(test_char);	//if data line, than the first element in str in an integer

		if (!is_data_line) { //not integer therefore this is not a data line
			string curr_row_option = tokens[OPTION_IDX];
			if (curr_row_option == "track") {
				continue; //error message or other notification not necessary, because this case is not a real error
			}

			if ((curr_row_option == "variableStep") || (curr_row_option == "fixedStep")) {
				string row_scaff_name;
				if (curr_row_option == "variableStep") {
					option = "variableStep";
					row_scaff_name = process_variable_step(tokens, idx_curr_row, curr_row);
				}
				else {
					if (curr_row_option == "fixedStep") {
						option = "fixedStep";
						tie(row_scaff_name, chrom_pos, step) = process_fixed_step(tokens, idx_curr_row, curr_row);
					}
				}

				if (!is_first_scaffold) {
					idx_former_scaff = get_idx_curr_scaff(block_scaff_name);
					//idx of the former scaffold, so block_scaff_name instead of row_scaff_name
				}

				bool do_computation = (row_scaff_name != block_scaff_name) && !is_first_scaffold;
				block_scaff_name = row_scaff_name;
				if (do_computation) {
					const Genomic_Data::Scaff_plus_Gen_Data& curr_scaff = Genomic_Data::get_all_genes()[idx_former_scaff];
					UTRs::compute_UTRs(curr_scaff, wiggle_data);
					do_computation = false;
					wiggle_data.clear();
				}
				is_first_scaffold = false;
			}
			else {
				ERR_STRM << "Error in wiggle file (line " << idx_curr_row
					 << "): Neither variableStep nor fixedStep, a track line or a comment!" << endl;
				cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
				throw Error();
			}
		}
		else { //not declaration line, hence it has to be a data line (normally)
			if (option == "variableStep") {

				const unsigned POSITION = 0;
				const unsigned DATA_VALUE = 1;

				unsigned chrom_pos;
				double valueData;

				try {
					chrom_pos = lexical_cast<unsigned>(tokens[POSITION]);
				}
				catch (bad_lexical_cast const&) {
					ERR_STRM << "Error in wiggle file (line " << idx_curr_row
					         << "): Position in '" << curr_row << "' is not an integer! Row will be ignored!" << endl;
					cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
					continue;
				}

				try {
					valueData = lexical_cast<double>(tokens[DATA_VALUE]);
				}
				catch (bad_lexical_cast const&) {
					ERR_STRM << "Error in wiggle file (line " << idx_curr_row
						 << "): Coverage Value in '" << curr_row << "' is not a double! Row will be ignored!" << endl;
					cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
					continue;
				}

				wiggle_data[chrom_pos] = valueData;
			}
			else {
				if (option == "fixedStep") {

					const unsigned DATA_VALUE = 0;
					try {
						double valueData = lexical_cast<double>(tokens[DATA_VALUE]);

						for (int i = 0; i < (int)step; ++i) {
							wiggle_data[chrom_pos+i] = valueData;
						}
						chrom_pos+=step;
					}
					catch (bad_lexical_cast const&) {
						ERR_STRM << "Error in wiggle file (line " << idx_curr_row
						     	 << "): Coverage Value '" << curr_row << "' is not a double!" << endl;
						cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
						throw Error(); //because otherwise coordination problemes
					}
				}
				else {
					if (option == "NULL") {
						ERR_STRM << "Error in wiggle file (line " << idx_curr_row
						    	 << "): A data line occured before a declaration line!" << endl;
						cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
						throw Error(); 	//could be done without exit, e.g. with a bool variable
					} //If unnecessary, because it could only be NULL in this case, but maybe it's easier to understand it in this way
				}
			}
		}
	}

	//for last scaffold in wiggle file
	idx_former_scaff = get_idx_curr_scaff(block_scaff_name);
	const Genomic_Data::Scaff_plus_Gen_Data& curr_scaff = Genomic_Data::get_all_genes()[idx_former_scaff];
	UTRs::compute_UTRs(curr_scaff, wiggle_data);

	wiggle_ifstream.close();
	cout << "Read in of wiggle file finished successfully!" << endl;
}

