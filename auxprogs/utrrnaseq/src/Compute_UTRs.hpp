/**
 * \file Compute_UTRs.hpp
 * Reads in and stores the coverage information for the scaffolds.
 * Computes the %UTRs for all the scaffoldswith this coverage information and the other input data.
 */

#ifndef COMPUTE_UTRS_HPP_
#define COMPUTE_UTRS_HPP_

#include "Genomic_Data.hpp"

#include <vector>
#include <string>
#include <map>

#include <boost/tuple/tuple.hpp>

/**
 * @brief Determines if the scaffold of a block in a WIG file also is in the genomic sequence file and if it is
 * finds the index number of this scaffold in the genomic sequence file (starting from 0). Else an error occurs and the user is notified. A block are all
 * entries belonging to the scaffold, independent on the WIG format.
 * @param block_scaff_name Name of the scaffold of the current block in the WIG file.
 * @return Number of the scaffold in the FASTA file.
 */
unsigned get_idx_curr_scaff(std::string block_scaff_name);

/**
 * @brief Processes definition line in a WIG file for the variableStep format. Extracts scaffold name.
 * @param tokens All substrings of the declaration line, split via tabulator.
 * @param idx_curr_row Line number of the declaration line in the WIG file. Only used for error message.
 * @param curr_row  Whole declaration line of the WIG file. Only used for error message.
 * @return Scaffold name in the declaration line.
 */
std::string process_variable_step(const std::vector<std::string>& tokens, unsigned idx_curr_row, std::string curr_row);

/**
 * @brief Processes of definition line in the WIG file for fixedStep format.
 * @param tokens All substrings of the declaration line, split via tabulator.
 * @param idx_curr_row Line number of the declaration line in the WIG file. Only used for error message.
 * @param curr_row  current Whole declaration line of the WIG file. Only used for error message.
 * @return Scaffold name, start position in scaffold and step size.
 */
boost::tuple<std::string, unsigned, unsigned>
	process_fixed_step(const std::vector<std::string>& tokens, unsigned idx_curr_row, std::string curr_row);

/**
 * @brief Processes the WIG file and determines %UTRs from WIG file.
 * @param wiggle_fname Name of the file containing the wiggle data in WIG format.
 */
void compute_UTRs(std::string wiggle_fname);

#endif /* COMPUTE_UTRS_HPP_ */
