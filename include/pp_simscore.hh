/*
 * pp_simscore.hh
 * Project: Similarity-Score for Protein-Profile and Protein-Sequence  
 * Author:  Lars Gabriel
 *
 *
 * USAGE: ../src/ ./pp_simScore protein_sequence_file.fa protein_profile_file_prfl
 */
 
 
 //TO DO: markdown, descriptions, store alignments in file, options, errors, test
 
 
#ifndef _PP_SIMSCORE_HH
#define _PP_SIMSCORE_HH
 
#include <iostream>
#include <fstream>
#include <limits>
#include <math.h>
#include "properties.hh"
#include "pp_profile.hh"

 namespace SS { 
    
    /** 
    * @brief structure representing one cell of a similarity matrix
    * 
    * @details 
    * Contains the highest similarity score of this position in the matrix and the indices of the cell from which the score was computed with corresponding type of score.
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
    * Each row represents a column of a block in the protein profile. \n
    * In each row of the similarity matrix only a connected part of cells are not empty, since the whole protein profile has to be included in the protein sequence. \n
    * The class Row contains only such cells, which are not empty. \n
    * These cells are stored in an array with fixed length. \n
    * The array can be accessed by [ ] with an index j {0, ..., length-1}. \n
    * Also the position of the first cell of the array is stored, such that the position of each cell can be determined by position + j. \n
    *
    * @param length: number of not empty cells \n
    * @param position: position of the first not empty cell in the row
    * 
    * @author Lars Gabriel
    */    
    class Row {
        private:
            
        public:
            Cell* row;       
            int position;
            Row (int length, int position);
            Cell& operator[] (int n);
            int length;   //method or const!!        
    };
    
    
    
    /**
    * @brief class representing a similarity matrix, for the SimilarityScore algorithm
    *
    * @details
    * Example of the structur of a similarity matrix: \n
    *<pre>
    *
    *                |         protein sequence         |
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
    * Contains a vector of Row, each row can be accessed by [ ]. \n
    * A new row can be added with the method addRow.
    * 
    * @author Lars Gabriel
    */
    class SimilarityMatrix {
        private:
            std::vector<Row> matrix;
        public:            
            SimilarityMatrix ();
            Row& operator[] (int n);
            void addRow (int l, int p);
            int length ();  
            bool empty ();     
            ~SimilarityMatrix(); 
    };
    
    /**
    * @brief structure representing the information of a column from the protein profile for a position of an alignment
    *
    * @details Contains the protein which is most likely at this position in the protein profile and the probability of the aligned protein of the sequence at that position. \n
    * 
    * @author Lars Gabriel
    */
    struct PrflAlignmentElement { 
        char argmax;
        double match_prob;    
    };
    
    
    /**
    * @brief structure representing all information internally needed for backtracking one final alignment
    *
    * @author Lars Gabriel
    */
    struct BacktrackAlign {//Name
        int i;
        int j;
        int j_prev;
        int block_count;
        int col_count;
        std::vector <pair <char, PrflAlignmentElement > > alignment; 
    };
    

    /**
    * @brief algorithm to compare the similarity of a protein sequence and a protein profile and compute an optimal alignment
    *
    * @details 
    * Based on the Needleman-Wunsch algorithm, it computes with a SimilarityMatrix the optimal global alignment score. \n
    * For every column of the protein profile and every position in the protein sequence the optimal alignment score up to that position will be computed. \n
    * Here are the columns of the SimilarityMatrix the proteins of the protein sequence and the rows of the SimilarityMatrix are the columns of the protein profile. \n
    * At the moment, gaps in the profile are still allowed. \n
    *
    *@param g: gap_cost
    *
    * @author Lars Gabriel
    */
    class SimilarityScore {
        private:
            //Similarity matrix:
                //one column represents one entry of the sequence
                //one row represents an entry of the blockprofile
            SimilarityMatrix S;
            double gap_cost;
            int seq_length = 0;
            char* seq = NULL;
            PP::Profile* prfl;
            
            //method for comparison, whether a score is the new best score for a cell and store the information of the current best score for backtracking 
            double compareScores(double new_score, double old_score, int current_i, int current_j, int prev_i, int prev_j, char score_type);
            //int aa2int(char c);
            
        public:
        
            //final alignments of the sequence and the profile
            std::vector <std::vector <pair <char, PrflAlignmentElement > > > alignments;
            SimilarityScore (double g);            
            ~SimilarityScore();
            
            //read files::
                //sequence input file has to be in FASTA
                //blockprofile input file has to be in .prfl 
            void readFiles (char* seqFileName, char* profileFileName);
            //algorithm to fill the similarity matrix and compute the final similarity score
            void fillSimilarityMatrix ();
            //traces back the alignments corresponding to the final score
            void backtracking ();
            double score();
            //prints similarity matrix and alignments with cout
            void printSimilarityMatrix ();
            void printAlignment ();
    };

}//namespace SS

#endif //_PP_SIMSCORE_HH
