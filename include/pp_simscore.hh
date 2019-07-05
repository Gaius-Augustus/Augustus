/*
 * Name:    pp_simscore.hh
 * Project: Similarity-Score for Protein-Profile and Protein-Sequence  
 * Author:  Lars Gabriel
 *
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
    
    //represents one cell in the similarity matrix, it contains
        //the score
        //the index of the cell from which the score was computed 
        //the type of score 
            //(match, match with gap in the interblock region, gap in sequence, gap in profile)
    struct Cell {
        double score;
        std::vector<std::tuple<int, int, char>> prev;
    };
    
    //represents a row in the similarity matrix, it contains
        //array with fixed number of cells
        //each cell can be accessed by []
        //position: number of entries before the first entry, which are empty
    class Row {
        private:
            
        public:
            Cell* row;       
            int position;
            Row (int length, int position);
            Cell& operator[] (int n);
            int length;   //method or const!!        
    };
    
    //represents a similarity matrix, it contains
        //vector of rows (can be accessed by [])
        //method, for adding a new row to the matrix
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
    
    //represents all information needed for backtracking one final alignment
    struct backtrack_align {
        int i;
        int j;
        int j_prev;
        int block_count;
        int col_count;
        std::vector<std::string > align_prfl;
        std::string align_seq;    
    };
    
   
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
            int aa2int(char c);
            
        public:
        
            //final alignments of the sequence and the profile
            std::vector <std::vector <std::string > > align_prfl;
            std::vector <std::string > align_seq;
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
