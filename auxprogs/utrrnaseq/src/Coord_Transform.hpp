/**
 * \file Coord_Transform.hpp
 */

#ifndef COORD_TRANSFORM_H_
#define COORD_TRANSFORM_H_

#include "Genomic_Data.hpp"
#include "flex_vec.hpp"

#include <string>
#include <vector>

/**
 * @brief Class to transform the scaffold positions into positions in a relative coordinate system
 * where 0 corresponds to the start position of the determination of the UTR in the scaffold.
 * The UTR is at most, independently from the used determination method, \b limit long or ends at \b max_pos,
 * which is either the end of the scaffold or the position before the boundary of the next gene on the scaffold.
 */
class Coord_Transform {

	public:

		/**
		 * @brief Compute and save the corresponding positions of the scaffold in the relative
		 * coordinate system in a container. The positions in the container correspond to the ones in the relative coordinate
		 * system. Returns container with positions in the scaffold corresponding to the ones in the relative coordinate system, where \b comp_start
		 * is the index 0 in this container.
		 * @param comp_start Position in the scaffold where the UTR determination starts. Corresponding position in the container is 0.
		 * @param dir Direction of the computation. Either 1 (starting from stop codon on + strand or start codon on - strand) or -1
		 * (starting from start codon on + strand or stop codon on - strand).
		 * @param limit Greatest possible length of UTR, excluding introns.
		 * @param max_pos Largest end position of the UTR in the Scaffold. Either end of scaffold or beginning of next gene.
		 * @param W Smoothing window size (for more information see external documentation: "Finding likely transcript ends based
		 * on RNA-Seq coverage").
		 * @param scaffold_size Length of the scaffold. Number of base pairs in the scaffold.
		 * @param curr_introns Current set of introns of the scaffold. Only introns of the scaffold on the strand of the start or stop
		 * codon and therefore of the later computed UTR.
		 */
		Coord_Transform(unsigned comp_start, int dir, unsigned limit, unsigned max_pos, unsigned W, unsigned scaffold_size,
						std::vector<Genomic_Data::Intron>* curr_introns);

		/**
		 * @brief Transforms a coordinate in the relative coordinate system into the corresponding scaffold position.
		 * @param pos Coordinate in the relative coordinate system, which has a corresponding position in the scaffold.
		 * @return Returns the corresponding scaffold position for coordinate \b pos in the relative coordinate system.
		 */
		unsigned pos_transform(int pos);

		/**
		 * @brief Get the last coordinate in the relative coordinate system with a corresponding scaffold position.
		 * @return Returns the last coordinate in the relative coordinate system.
		 */
		int get_last_vec_pos() { return m_last_vec_pos; }

	private:

		/**
		 * @brief Scaffold positions for the corresponding positions in the relative coordinate system.
		 *
		 * This container contains the scaffold positions for the corresponding positions in the relative
		 * coordinate system, which correspond directly to the positions in this container. First index in
		 * this container is \b -W and the last index \b limit.
		 */
		Flex_Vec<unsigned> m_transform_vec;

		/**
		 * @brief Smoothing window size.
		 *
		 * Smoothing window size W (for more information see external documentation: "Finding likely transcript ends based on
		 * RNA-Seq coverage").
		 */
		unsigned m_win_size;

		/**
		 * @brief Last coordinate in the relative coordinate system.
		 *
		 * Last coordinate in the relative coordinate system with a corresponding scaffold position.
		 * This position is at most \b limit.
		 */
		int m_last_vec_pos;
};


#endif /* COORD_TRANSFORM_H_ */
