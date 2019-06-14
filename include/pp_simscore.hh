/*
 * Name:    pp_simscore.hh
 * Project: Similarity-Score for Protein-Profile and Protein-Sequence  
 * Author:  Lars Gabriel
 *
 *    g++ -o pp_simScore pp_simscore.cc types.o properties.o geneticcode.o pp_profile.o lldouble.o
 *	Usage:	./pp_simScore seq_file.fa profile.prfl
 *
 */
 
 
 
#include <iostream>
#include <fstream>
#include <limits>
#include <math.h>
#include "properties.hh"
#include "pp_profile.hh"

 namespace SS {
 
 
    struct Cell {
        double score;
        list<struct Cell *> prev;
    };

    class Column {
        private:
            int position;
            int length;
        public:
            Cell* col;

            Column (int length, int position);

            Cell& operator[] (int n);

            int getElementPosition (int n);	
            
    };


    class SimilarityMatrix {
        private:
            int length;
    
        public:
            std::vector<Column> matrix;
            SimilarityMatrix (int length);
            Column& operator[] (int n);
            void addColumn (int l, int p);
            void addColumn (Column c);
            
        
    };

    class SimilarityScore {
        private:
        
        public:
            SimilarityMatrix* S();
            std::string seq;
            PP::Profile* prfl;
            SimilarityScore(char* seqFileName, char* profileFileName);
            void readFiles(char* seqFileName, char* profileFileName);
            
        
    
    };

}//namespace SS



















