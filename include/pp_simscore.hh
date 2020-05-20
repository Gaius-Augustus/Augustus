/*
 * Name:    pp_simscore.hh
 * 
 * License: Artistic License, see file LICENSE.TXT or 
 *          https://opensource.org/licenses/artistic-license-1.0
 * 
 * Project: similarity-score algorithm for block-profile and protein sequence
 * with intron information
 * Author:  Lars Gabriel
 *
 */

#ifndef __PP_SIMSCORE_HH
#define __PP_SIMSCORE_HH

#include <iostream>
#include <fstream>
#include <limits>
#include <math.h>
#include <string>
#include "properties.hh"
#include "pp_profile.hh"

 namespace SS {
    /**
    * @brief structure representing one cell of a similarity matrix
    *
    * @details
    * Contains the highest similarity score of a position in the
    * similarity matrix, the indices of the cell from which the score
    * was computed and the corresponding type of score.
    *
    * <pre>
    * types: m: match without gaps
    *        s: gap in the protein sequence
    *        p: gap in the protein profile
    *        _: match with gap in the protein sequence in the interblock region
    * </pre>
    *
    * @author Lars Gabriel
    */
    struct Cell {
        double score;
        std::vector < std::tuple <int, int, char > > prev;
    };

    /**
    * @brief class representing one row of a SimilarityMatrix
    *
    * @details
    * Each row represents a column of a block in the protein profile.
    * Each row of the similarity matrix consists of a section of connected cells,
    * the other cells are empty. The length of the section depends on the assumption,
    * that the whole protein profile has to be included in the protein sequence.
    * The class Row contains only such cells, which are not empty.
    * These cells are stored in an array with fixed length.
    * The array can be accessed by [ ] with an index j {0, ..., length-1}.
    * The position in the similarity of the first cell of the row is stored,
    * such that the position in the similarity matrix of each cell
    * can be determined by position + j.
    *
    * @param length: number of not empty cells
    * @param position: position of the first non empty cell in the row
    *
    * @author Lars Gabriel
    */
    class Row {
        private:
        int len;
        public:
        Cell* row;
        int position;
        Row (int length, int position);
        Cell& operator[] (int n);
        int length();
    };

    /**
    * @brief class representing a similarity matrix, for the
    *   SimilarityScore algorithm
    *
    * @details
    * Example of the structur of a similarity matrix:
    *<pre>
    *
    *                |         protein sequence         |
    *
    *             –   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    *             ┬   ----XXXXXXXXX---------------------
    *             │   -----XXXXXXXXX--------------------
    *  protein-   ┴   ------XXXXXXXXX-------------------
    *  profile    –
    *             ┬   ------------XXXXXXXXXXXXXXX-------
    *             │   -------------XXXXXXXXXXXXXXX------
    *             │   --------------XXXXXXXXXXXXXXX-----
    *             ┴   ---------------XXXXXXXXXXXXXXX----
    *             –   ---------------------------------X
    *
    * </pre>
    *
    * Contains a vector of Row, each row can be accessed by [ ].
    * A new row can be added with the method addRow.
    *
    * @author Lars Gabriel
    */
    class SimilarityMatrix {
        private:
        std::vector<Row> matrix;
        public:
        Row& operator[] (int n);
        void addRow (int l, int p);
        int length ();
        bool empty ();
        ~SimilarityMatrix();
    };

    /**
    * @brief class representing a protein sequence with optional
    *informations about intron postions
    *
    * @details
    * The protein sequence has to be in fasta file format.
    * The optional intron informations must be attached to the protein sequence.
    * Format structure:
    *<pre>
    * >protein sequence header
    * XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    * XXXXXXX protein sequence XXXXXXXXXXX
    * XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    *
    * [Introns]
    * # index of the position after which an intron occures | residual nucleotides before the intron
    * 2 0
    * 5 1
    * 30 2
    * 104 1
    *</pre>
    *
    * The amino acids of the protein sequence can be accessed by []
    * with their position in the sequence.
    * With intronAt() it can be checked if there is an intron at a position
    * The number of introns in a specific frame can be determined
    * with intronsInRange(start, end)
    * For example: intronsInRange(3, 104) = 2
    *
    * @author Lars Gabriel
    */

    class ProteinSequence {
        private:
        char* sequence = NULL;
        std::map < int, int > introns;
        std::vector < int > num_prev_introns;
        int len;
        char* name;
        public:
        ProteinSequence (const char* fileName);
        ~ProteinSequence ();
        char& operator[] (int n);
        int length();
        int intronAt(int i);
        int intronsInRange (int start, int end);
    };

    /**
    * @brief Information of the representation of the profile in an alignment
    * for a readable output to stdout
    *
    * @details Contains the .

    *
    * @author Lars Gabriel
    */
    struct PrflAlignment {
        //amino acid sequence that has the highest likelihood to belong
        //to the profile
        std::string argmax_aminoacids;
        //for each column the frequency of the aligned amino acid from
        //the protein sequence in the profile
        std::vector < double > match_prob;
    };

    /**
    * @brief structure representing one element
    * of a final alignment produced by the backtracking algorithm
    *
    * @details
    * stores for every alignment the connected regions that are of the same
    * alignment type and are part of the alignment of the same block
    * List of:    - starting position of the first amino acid of the protein
    *                sequence that is included in the alignment frame
    *            - block number in which the alignment frame is located
    *            - index of the first block column that is included in
                   the alignment frame
    *            - length of the frame (number of alignment columns)
    *            - alignment type: 'm', 's'. 'p' or '-'
    *
    * @author Lars Gabriel
    */
    struct AlignmentElement {
        int protSeq_pos;
        int block_pos;
        int col_pos;
        int length;
        char type;
    };

    /**
    * @brief structure representing all information internally needed for
    * backtracking one final alignment
    *
    * @author Lars Gabriel
    */
    struct BacktrackAlignDB {
        int i;
        int j;
        int j_prev;
        int block_count;
        int col_count;
        int col_posStart;
        AlignmentElement align_element;
        std::vector <AlignmentElement > alignment;
    };

    /**
    * @brief algorithm to compare the similarity of a protein sequence
    * and a protein profile and compute an optimal alignment
    *
    * @details
    * Based on the Needleman-Wunsch algorithm, it computes with a
    * SimilarityMatrix the optimal global alignment score.
    * For every column of the protein profile and every position in the
    * protein sequence the optimal alignment score up
    * to that position will be computed.
    * Here, the columns of the SimilarityMatrix are the proteins of the
    * protein sequence and the rows of the SimilarityMatrix are the
    * columns of the protein profile.
    * At the moment, gaps in the alignment are allowed.
    *
    * @param    g: gap cost in inter-block region,
                b: gap cost for an intra-block gap,
                g_i: gap cost for a gap in the alignment of intron positions,
                iw1: intron weight for the intron score in an intra-block region,
                iw2: intron weight for the intron score in an intra-block region,
                e_i: pseudocount parameter for epsilon1,
                e_n pseudocount parameter for epsilon2
    *
    * @author Lars Gabriel
    */
    class SimilarityScore {
        private:
        //Similarity matrix:
        //one column represents one entry of the sequence
        //one row represents an entry of the blockprofile
        SimilarityMatrix S;
        //data structure for the protein sequence P
        ProteinSequence* seq;
        //data structure for the block profile B
        PP::Profile* prfl;
        //parameter for the gap costs
        double gap_cost_inter;
        double gap_cost_intra;
        double gap_cost_intron;
        //parameter for the intron weights
        double intron_weight_inter;
        double intron_weight_intra;
        //parameter for the pseudocount
        //the pseudocount is added to a relative frequency (v/w) of an
        //intron position with (v+epsi_intron)/(w+epsi_intron+epsi_noIntron)
        double epsi_intron;
        double epsi_noIntron;

        //number of protein sequences used to create the
        //intron profile of the block profile
        int intronPrfl_noProt;

        //background frequency of an introns at an intrablock intron position,
        //estimated from the intron positions of 15799 protein sequences
        const double intronIntraBlock_bfreq = 1.706e-3;
        //background frequency of the number of introns in
        //an intrablock intron position, estimated from the intron positions
        //of 15799 protein sequences in inter-block section
        //the length an inter-block section that was used for this estimation
        //is the mean of the inter-block distance interval
        const double intronAvgInterBlock_bfreq = 9.599e-3;

        //computes the intra-block intron score for a match of intron positions
        //with f residual nucleotides, if an intron of P is aligned to
        //block column s of block k
        //if f=-1: ths score of a mismatch in intron position
        //for f in {0,1,2} is computed
        double intraBlock_iscore (int k, int s, int f, int intron_frame);
        //Poisson distribution, which is used as the background distribution
        //of k intron in an inter-block section
        //lambda is estimated for I=[d_min, d_max] as
        //3 * ((d_min+d_max)/2 +1) * (intronAvgInterBlock_bfreq/intronPrfl_noProt)
        double PoiDist (int k, double lambda);
        //method for comparison, whether a score is the new best score
        //for a cell
        //and stores the information of the current best score for backtracking
        double compareScores(double new_score, double old_score, int current_i,\
            int current_j, int prev_i, int prev_j, char score_type);
        //calculates the intron score
        double intronScore (double q, double q_b);

        public:
        SimilarityScore (double g, double b, double g_i, double iw1,\
            double iw2, double e_i, double e_n);
        ~SimilarityScore ();
        //prints the inter block distance intervals of B to stdout
        void printInterBlock();
        bool SimilarityMatrix_empty();
        //prints an optimal alignment as the coordinates of connected
        //alignment sections of the same alignment type
        //Format: List of AlignmentElement
        void printAlignmentDB();
        //prints a list of a translations from the indice of a block to
        //the number of the block in the .prfl file.
        void printBlockPos ();
        //prints the average of the argmax of the block columns for
        //the complete profile
        pair<double, double> avgArgMaxProb();
        //read files:
        //sequence input file has to be in FASTA
        //intron positions are optional as '[Introns]' section
        //in the file format (see description of class ProteinSequence)
        //blockprofile input file has to be in .prfl
        //intron information are optional in the file format
        //(see README file)
        void readFiles (const char* seqFileName, const char* profileFileName);

        //algorithm to fill the similarity matrix and compute
        //the final similarity score
        void fillSimilarityMatrix ();

        //traces back the alignments corresponding to the final score
        //alignment are stored as list of AlignmentElements
        void backtrackingDB (int number_Alignments);

        //print final similarity score
        double score();

        //final alignments of the sequence and the profile
        std::vector <std::vector <AlignmentElement > > alignments;

        //datastructure for a readable alignment, which is used
        //for printing an alignment
        std::vector < pair <std::string, PrflAlignment > > alignmentsReadable;
        //prints similarity matrix and alignments with cout
        void printSimilarityMatrix ();

        //print readable alignment to stdout
        //FORMAT: - Alignment representation of P (AminoAcid, gap symbol
        //        or number of amino acids in inter-block)
        //        - Alignemnt representation of argmax of B
        //        (argmax AminoAcid for aligned block column,
        //        gap symbol or inter-block length)
        //        - frequency of amino acid of P in aligned block column of B,
        //        if alignment type is a match
        void printAlignment ();
    };
}
//namespace SS
#endif