/**
 * \file Splice_Sites.hpp
 * Parses of splice sites and also filters the introns and through this is able to determine the intron strand.
 */

#ifndef SPLICE_SITES_HPP_
#define SPLICE_SITES_HPP_

#include <vector>
#include <string>
#include <utility>

/**
 * @brief Filters the introns and determines the strand if undefined.
 * @param sp_sites Splice sites used for splice site filtering.
 * @param sequence Genomic sequence of the scaffold.
 * @param start Start position of intron. Strand independent. Smaller than end.
 * @param end End position of intron. Strand independent. Greater than start.
 * @param strand_defined True, if the strand was defined ( '+', '-'). False, if the strand was undefined ('.').
 * @param strand Strand of the intron as stated in the GFF file with the introns.
 * @return True, if intron is accepted; False, if intron is not accepted. And strand of the intron, only important when the
 * strand was undefined.
 */
std::pair<bool, std::string> check_sp_sites (const std::vector< std::pair <std::string, std::string> >& sp_sites, const std::string& sequence,
		unsigned start, unsigned end, bool strand_defined, std::string strand);

/**
 * @brief Converts splice site on minus strand into complementary splice site so it fits the format
 * of the splice sites on the plus strand, e.g. ...CT...AC... on minus strand, with donor site CA and
 * acceptor site TC, is converted into donor site GT respectively acceptor site AG.
 * @param sequence Genomic sequence of the scaffold.
 * @param start Start coordinate of splice site in the sequence.
 * @param end End coordinate of splice site in the sequence.
 * @returns Splice site of minus strand in plus strand format.
 */
std::string sp_site_conversion (const std::string& sequence, unsigned start, unsigned end);

/**
 * @brief Parses the splice sites from command line format.
 * @param raw_sp_sites Unparsed and raw splice sites. Original input format, all splice sites in one string.
 * @return All splice sites separated into donor and acceptor site.
 */
std::vector< std::pair <std::string, std::string> > parse_sp_sites (std::string raw_sp_sites);


#endif /* SPLICE_SITES_HPP_ */
