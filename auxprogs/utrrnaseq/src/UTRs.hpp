/*
 * \file UTRs.hpp
 */

#ifndef UTRS_HPP_
#define UTRS_HPP_

#include "Genomic_Data.hpp"

#include <vector>
#include <string>

/**
 * @brief Determines the UTRs and writes them into an output file.
 */
class UTRs {

	public:

		/**
		 * @brief All parameters that can be set by the user.
		 *
		 * All parameters that can be set by the user. For more information for some of the used parameters see
		 * external documentation: "Finding likely transcript ends based on RNA-Seq coverage".
		 */
		struct Parameters {
			std::string output_fname; ///< Name of the output file.
			unsigned read_length; ///< Length of the reads of the used RNA-Seq Data.
			unsigned window_size; ///< W Smoothing window size (see external literature).
			unsigned limit; ///< Greatest possible length of an UTR, excluding introns.
			unsigned drop_window_size; ///< V window size after UTR end (see external documentation).
			unsigned min_length; ///< Minimal length of UTR.
			unsigned min_average_cov; ///< Minimal average coverage of UTR.
			double p_win; ///< Percentage of drop window coverage in comparison to average UTR exon coverage (see external documentation).
			double p_int; ///< Percentage of intron coverage in comparison to average UTR exon coverage.
			double p_mult; ///< Percentage of multiplicity in comparison to average UTR exon coverage.
			bool zero_cov; ///< False: UTR end due to steepest decline in coverage; True: UTR end where the coverage drops to 0.
		};

		/**
		 * @brief Gets reference to parameters that can be set by the user.
		 * @returns Reference to all parameters that can be set by the user.
		 */
		static Parameters& parameters() { return s_parameters; };

		/**
		 * @brief Deletes all previous UTRs. Important for test cases.
		 */
		static void clear();

		/**
		 * @brief Determines the UTRs from the RNA-Seq coverage data \b wiggle_data for the current scaffold \b curr_scaffold.
		 * @param curr_scaffold Sequence plus genomic data of the scaffold the UTRs are determined for.
		 * @param wiggle_data Coverage information for the current scaffold \b curr_scaffold; Key: Position; Value: Coverage.
		 */
		static void compute_UTRs(const Genomic_Data::Scaff_plus_Gen_Data& curr_scaffold,
								const std::map<unsigned,double>& wiggle_data);

		/**
		 * @brief Returns the number of all computed UTRs.
		 * @return Number of UTRs.
		 */
		static unsigned size();

		/**
		 * @brief Writes all the computed UTRs into output file.
		 */
		static void output();

		struct UTR {
			std::string name; ///< Name of the scaffold the UTR belongs to.
			unsigned start; ///< Start position of the UTR.
			unsigned end; ///< End position of the UTR.
			std::string feature; ///< Either "3'-UTR" or "5'-UTR".
			std::string strand; ///< Strand of the UTR. Either '+' (plus) or '-' (minus).
			std::string group; ///< Group attributes, gene id and transcript id, extracted from GTF file with the coding region information.
			std::vector<Genomic_Data::Intron> all_introns; ///< All introns in the UTR.
		};

		/**
		 * @brief Gets all the computed UTRs.
		 * @returns All computed UTRs.
		 */
		static std::vector<UTR>& get_all_UTRs() { return s_all_UTRs; };

	private:

		/**
		 * @brief Removes less likely introns until no introns overlap, which would be biological unlikely.
		 * @param introns All introns on the same strand and belonging to the same scaffold.
		 */
		static void remove_overlapping_introns(std::vector<Genomic_Data::Intron>* introns);

		/**
		 * @brief Computes greatest genome positon for UTR end from the genome position at which the computation of the UTR starts.
		 * @param pos Position in the scaffold where the determination of the UTR starts.
		 * @param dir  Direction of the computation. Either 1 (starting from stop codon on + strand or start codon on - strand)
		 * or -1 (starting from start codon on + strand or stop codon on - strand).
		 * @param scaffold_size Length of the scaffold. Number of base pairs in the scaffold.
		 * @param read_length Length of the reads of the used RNA-Seq Data.
		 * @param coding_region_boundaries All gene boundaries on the same strand as the UTR.
		 * @returns Greatest possible position in genome. Either end of scaffold or position before the first gene after \b pos.
		 */
		static unsigned get_max_pos(unsigned pos, int dir, unsigned scaffold_size, unsigned read_length,
				                    std::vector<Genomic_Data::CRB>* coding_region_boundaries);

		/**
		 * @brief Tests if a repeat overlaps with a UTR.
		 * @param start Start position of the range of the UTR.
		 * @param end End position of the range of the UTR.
		 * @param repeats All repeats of the scaffold the UTR belongs to.
		 * @returns True: if a repeat exists in the UTR; False: if no repeat is in the UTR range.
		 */
		static bool repeat_in_UTR(unsigned start, unsigned end, std::vector<Genomic_Data::Repeat>* repeats);

		/**
		 * @brief Finds all introns in the range from start to end.
		 * @param start Start position of the range.
		 * @param end End position of the range.
		 * @param dir Direction of the computation. Either 1 (starting from stop codon on + strand or start codon on - strand)
		 * or -1 (starting from start codon on + strand or stop codon on - strand).
		 * @param all_introns All introns of the same strand as the current start/stop codon, or the UTR that is going to be computed.
		 * @returns All introns in the range from start to end.
		 */
		static std::vector<int> find_introns_in_range(unsigned start, unsigned end, int dir, std::vector<Genomic_Data::Intron>* introns);

		static std::vector<UTR> s_all_UTRs; ///< All computed UTRs.
		static Parameters s_parameters; ///< All parameters that can be set by the user.
};

#endif /* UTRs_HPP_ */
