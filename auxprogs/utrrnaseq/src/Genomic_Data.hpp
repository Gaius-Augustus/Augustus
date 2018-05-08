/**
 * \file Genomic_Data.hpp
 */

#ifndef GENOMIC_DATA_HPP_
#define GENOMIC_DATA_HPP_

#include <string>
#include <vector>
#include <utility>

/**
 * @brief Stores the genomic sequences and start and stop codon properties, plus intron and repeat positions, of the scaffolds.
 */
class Genomic_Data {

	public:

		/**
		 * @brief Extracts all necessary informations from the input files.
		 *
		 * Extracts necessary informations from the input files, with the scaffold sequences, start and stop codons,
		 * introns and repeats.
		 *
		 * @param scaff_fname Filename of the FASTA file with the scaffold sequences.
		 * @param crb_fname Filename of the GTF file with the coding region boundaries (start and stop codons).
		 * @param intron_fname Filename of the GFF file with the intron hints.
		 * @param repeat_fname Filename of the GFF file with the repeat coordinates.
		 * @param sp_sites All valid splice sites for introns.
		 * @param use_repeat_file Is the file with repeat coordinates used?
		 */
		static void initialize(std::string scaff_fname, std::string crb_fname, std::string intron_fname, std::string repeat_fname,
				std::vector< std::pair< std::string, std::string> > sp_sites, bool use_repeat_file );

		/**
		 * @brief Stores coordinates and group of a start or stop codon.
		 */
		struct CRB {
			std::string feature;	///< Either start or stop (codon).
			unsigned codon_pos; 	///< Position upstream (3' UTR) or downstream (5' UTR) from UTR start position. Last Position in the gene the UTR belongs to.
			std::string strand;	///< Either plus strand or minus strand.
			std::string grp_att;	///< List of feature attributes in the format tag=value; necessary for the output file.
		};

		/**
		 * @brief Stores coordinates and multiplicity of an intron.
		 */
		struct Intron {
			unsigned start;	///< Start position of the intron. Upstream of \b end on its strand.
			unsigned end;	///< End position of the intron. Downstream of \b start on its strand.
			std::string strand;	///< Either plus strand or minus strand.
			unsigned mult;	///< Multiplicity = copy number; absolute number of reads that span the intron.
			/**
			 * @brief Tests if two introns overlap.
			 * @param i First intron.
			 * @param j Second intron.
			 * @returns True: intron i and j overlap; False: no overlapping.
			 */
			static bool is_overlapping(Intron i, Intron j);
		};

		/**
		 * @brief Stores coordinates of a repeat.
		 */
		struct Repeat {
			unsigned start;	///< Start position of the repeat. Upstream of \b end on its strand.
			unsigned end;	///< End position of the repeat. Downstream of \b start on its strand.

			bool operator==(const Repeat& r) const;   ///< Comparison operator
		};

		/**
		 * @brief Stores the informations from the sequence file and the coding region boundaries, possible introns and repeats, which belong to the scaffold.
		 */
		struct Scaff_plus_Gen_Data {
			std::string name;	///< Name of the scaffold.
			std::string sequence;	///< Genomic sequence of the scaffold.
			std::vector<CRB> crbs;	///< Start and stop coordinates of the scaffold.
			std::vector<Intron> introns;	///< %Intron coordinates and multiplicities of the scaffold (empty if no introns exist).
			std::vector<Repeat> repeats;	///< %Repeat coordinates of the scaffold (empty of no repeats exist).
		};


		/**
		 * @brief Gets all the scaffolds and their genomic features (start and stop codons, introns and repeats).
		 * @return All scaffolds and their genomic features from the input files.
		 */
		static std::vector<Scaff_plus_Gen_Data>& get_all_genes() { return s_all_genes; }

	private:

		/**
		 * @brief Reads in the sequences and names of the scaffolds.
		 * @param scaff_fname Filename of the FASTA file with the scaffold sequences.
		 */
		static void read_scaff_file(std::string scaff_fname);

		/**
		 * @brief Reads in the start and stop codon coordinates and group attributes.
		 * @param crb_fname Filename of the GTF file with the coding region boundaries (start and stop codons).
		 */
		static void read_crb_file(std::string crb_fname);

		/**
		 * @brief Reads in the intron coordinates and multiplicity.
		 * @param intron_fname Filename of the GFF file with the intron hints.
		 */
		static void read_intron_file(std::string intron_fname);

		/**
		 * @brief Reads in the repeat coordinates.
		 * @param repeat_fname Filename of the GFF file with the repeats.
		 */
		static void read_repeat_file(std::string repeat_fname);

		//member variables
		/**
		 * @brief All information from all scaffolds, including repeats, introns, sequence and start and stop codons.
		 */
		static std::vector<Scaff_plus_Gen_Data> s_all_genes;

		/**
		 * @brief All valid splice sites for introns.
		 */
		static std::vector< std::pair< std::string, std::string> > s_sp_sites;
};



std::ostream& operator<<(std::ostream& os, const Genomic_Data::Repeat& r);   ///< Output operator

#endif /* Genomic_Data_HPP_ */
