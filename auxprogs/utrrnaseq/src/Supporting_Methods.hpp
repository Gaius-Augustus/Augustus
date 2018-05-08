/**
 * \file Supporting_Methods.hpp
 * Contains supporting methods, small methods used globally in the tool.
 */

#ifndef SUPPORTING_METHODS_HPP
#define SUPPORTING_METHODS_HPP

#include <string>
#include <vector>

/**
 * @brief splits input string into its substrings
 * @param token string trough which the input string is split
 * @param str input string, which should be split into its substrings
 * @return the substrings of the input string, which are split through expression \b token
 */
std::vector<std::string> tokenize(std::string token, std::string str);

/**
 * @brief reset getopt, so getopt could be called an additional time
 */
void reset_getopt();

/**
 * @brief Comparing two files, decides if these two files have an identical content
 * @param first_file First input file (normally its path)
 * @param second_file Second input file (normally its path)
 * @return true, if the two files are identical or false, if there are not identical
 */
bool files_identical (std::string first_file, std::string second_file);

#endif /* SUPPORTING_METHODS_HPP_ */
