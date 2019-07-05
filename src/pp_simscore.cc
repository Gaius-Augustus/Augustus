/*
 * Name:    pp_simscore.cc
 * Project: Similarity-Score for Protein-Profile and Protein-Sequence  
 * Author:  Lars Gabriel
 *
 */


// project includes
#include "pp_simscore.hh"
#include "properties.hh"
#include "pp_profile.hh"
#include "fasta.hh"

// standard C/C++ includes
#include <iostream>
#include <fstream>
#include <limits>
#include <math.h>



using namespace SS;	

    Row::Row (int length, int position) {
        Row::length=length;
        Row::position=position;
        row = new Cell[length];
    }
    
	Cell& Row::operator[] (int n) {
        if (0 <= n && n<= length) {
            return row[n];
        } else { 
            throw out_of_range("Row out of bound. \n");
        }
     }

    SimilarityMatrix::SimilarityMatrix () {
        //SimilarityMatrix::length=length;
    }
    
    SimilarityMatrix::~SimilarityMatrix() {
        matrix.clear();
    }
    
    Row& SimilarityMatrix::operator[] (int n) {
        if (0 <= n && n< matrix.size()) {
            return matrix[n];
        } else { 
            throw out_of_range("Matrix out of bound. \n");
        }
    }
    
    void SimilarityMatrix::addRow (int l, int p) {
        Row r(l,p);
        matrix.push_back(r);    
    }
    
    
    int SimilarityMatrix::length () {
        return matrix.size();
    }
    
    bool SimilarityMatrix::empty () {
        return matrix.empty();
    }

    SimilarityScore::SimilarityScore (double g) {
        gap_cost = g;
    }
    
    SimilarityScore::~SimilarityScore() {
        delete prfl;
        for (int i = 0; i < S.length(); i++) {
            delete [] S[i].row;        
        } 
        delete [] seq;
    }
    
    void SimilarityScore::readFiles (char* seqFileName, char* profileFileName) {
    
        //read protein sequence
        ifstream ifstrm;
        char *name = NULL;
        
	    ifstrm.open(seqFileName);
	    if( !ifstrm )
	        throw string("Could not open input file \"") + seqFileName + "\"!";	        
	    //if(isFasta(ifstrm)){	    
	    readOneFastaSeq(ifstrm, seq, name, seq_length);	        
	    //}  
        
        //read protein profile
        PP::initConstants();     
        prfl = new PP::Profile(profileFileName);                
    }
    
    
    void SimilarityScore::fillSimilarityMatrix () {//only works, if profile is shorter than sequence!
        
        //row, which will be filled
        int i=0;
        //position of the current row
        int position = 0; 
        //row_length: length of a row
        //diff_pos: difference in position of current row and previous row
        //k_min: position of the first cell, which is used to compute the score of a cell
        //k_max: position of the last cell, which is used to compute the score of a cell
            //(if i represents the first position of a new block)
        int row_length, diff_pos, k_min, k_max;
        //needed for comparison of a new calculated score
        double new_score, old_score;
        //minimal space needed for the profile: length of all blocks + all min interblock distances
        int min_prfl_length = 0;
        int block_count = prfl->blockCount();
    	
	    PP::Column col; 
	    PP::DistanceType d;	
	    
	    //computes minimal profile length
	    for (int k=0; k<block_count; k++) {
	        min_prfl_length += prfl->blockSize(k) + prfl->interBlockDist(k).r.min;
	    }
	    
	    row_length = seq_length-min_prfl_length;
        S.addRow(seq_length+1,0);
        //start of the similarity score recursion 
        S[i][0].score = 0;
        //fill first row
        for (int j=1; j<seq_length+1; j++) {
            S[i][j].score = j * gap_cost;
            S[i][j].prev.push_back(make_tuple(i,j-1, 'p'));
        }
        position++;
        
        for (int t=0; t<block_count; t++) {

            PP::Position P(t,0);                         
            col = prfl->getColumn(P);
            d = prfl->interBlockDist(t);
            position += d.r.min;
            
            S.addRow(row_length, position);
            i++;                      
            
            //score for the first position in a block
            for (int j = 0; j<row_length; j++) {          
            
            diff_pos = S[i].position - S[i-1].position + j;
            k_min = ((diff_pos-d.r.max-1 < 0) ? 0 : diff_pos-d.r.max-1);
            k_max = ((row_length-1 < diff_pos) ? row_length-1 : diff_pos);            
            old_score = std::numeric_limits<int>::min();
            
            
                for (int k = k_min; k < k_max; k++) {                
                    
                    if (k < row_length && k < diff_pos-d.r.min) {
                        //match
                        new_score = S[i-1][k].score + col.L(aa2int(seq[position+j-1]));
                        old_score = compareScores(new_score, old_score, i, j, i-1, k, 'm');
                    }                    
                                   
                    else if (k < k_max) {                    
                        //gap in sequence
                        new_score = S[i-1][k].score + gap_cost * (d.r.min-k_max+k+1);
                        old_score = compareScores(new_score, old_score, i, j, i-1, k, 's');
                        
                        //match with gap
                        new_score = S[i-1][k].score + gap_cost * (d.r.min-k_max+k) + col.L(aa2int(seq[position+j]));
                        old_score = compareScores(new_score, old_score, i, j, i-1, k, '_');                        
                    }
                    else {
                        //gap in sequence
                        new_score = S[i-1][k].score + gap_cost * (d.r.min-k_max-k+1);
                        old_score = compareScores(new_score, old_score, i, j, i-1, k, 's');
                    }
                            
                }
                
                if (j>0) {
                    //gap in profile
                    new_score = S[i][j-1].score + gap_cost;
                    old_score = compareScores(new_score, old_score, i, j, i, j-1, 'p');
                }
                S[i][j].score = old_score;
            }
            
            //score for all other positions in a block
            for (int s = 1; s<prfl->blockSize(t); s++) {
                            
                PP::Position P(t,s);
                col = prfl->getColumn(P);
                position++;
                i++;
                S.addRow(row_length, position);     
                                
                for (int j = 0; j < row_length; j++) {

                    old_score = std::numeric_limits<int>::min();;
                    
                    if (j > 0) {  
                        //gap in profile  
                        new_score = S[i][j-1].score + gap_cost;
                        old_score = compareScores(new_score, old_score, i, j, i, j-1, 'p');
                    }
                    
                    if ( j < row_length-1) {
                        //gap in sequence
                        new_score = S[i-1][j+1].score +gap_cost;
                        old_score = compareScores(new_score, old_score, i, j, i-1, j+1, 's');                        
                    }
                    
                    //match
                    new_score = S[i-1][j].score + col.L(aa2int(seq[position+j-1]));
                    old_score = compareScores(new_score, old_score, i, j, i-1, j, 'm');
                    
                    S[i][j].score = old_score;                    
                }          
            }
            position++;   
        }
        
        //compute final score
        S.addRow(1,seq_length);
        old_score = std::numeric_limits<int>::min();
        d = prfl->interBlockDist(block_count);   
        for (int j=0; j<row_length; j++) {
            diff_pos = seq_length - S[i].position-j;  
            if (d.r.max < diff_pos) {
                new_score = S[i][j].score + gap_cost*(diff_pos - d.r.max);
                old_score = compareScores(new_score, old_score, i+1, 0, i, j, 'p');
            }
            else if (d.r.min <= diff_pos && diff_pos <= d.r.max) {
                new_score = S[i][j].score;
                old_score = compareScores(new_score, old_score, i+1, 0, i, j, 'm');
            }
            else if (diff_pos < d.r.min) {                
                new_score = S[i][j].score + gap_cost*(d.r.min - diff_pos);
                old_score = compareScores(new_score, old_score, i+1, 0, i, j, 's');
            }
        }
        S[i+1][0].score = old_score;
    }
    
    double SimilarityScore::score () {
        if (!S.empty()) {//!!!
            return S[S.length()][seq_length].score;
        }
        else {
            std::cerr<<"Similarity Matrix is empty. \n";
            exit(4);
        }
    }
    
    double SimilarityScore::compareScores (double new_score, double old_score, int current_i, int current_j, int prev_i, int prev_j, char score_type) {
    
        if (new_score > old_score) {
            old_score = new_score;
            S[current_i][current_j].prev.clear();
            S[current_i][current_j].prev.push_back(make_tuple(prev_i, prev_j, score_type));
        }
        else if (new_score == old_score) {
            S[current_i][current_j].prev.push_back(make_tuple(prev_i, prev_j, score_type));
        }
        return old_score;    
    }
        	
    void SimilarityScore::printSimilarityMatrix () {
        
        for (int i=0; i<S.length(); i++) {
        
            for (int j=0; j<S[i].position; j++) {
                std::cout<<"-"<<"  ";        
            }  
            
            for (int j=0; j<S[i].length; j++) {    
                for (int t=0; t<S[i][j].prev.size(); t++) {            
                    std::cout<<std::get<0>(S[i][j].prev[t])<<":"<<std::get<1>(S[i][j].prev[t])<<":"<<std::get<2>(S[i][j].prev[t])<<",";
                }
                std::cout<<" | ";
            }
            
            for (int j=0; j<seq_length+1-S[i].length-S[i].position; j++) { 
                std::cout<<"-"<<"  ";            
            }
                 
            std::cout<<endl;
        }
    }
    
    void SimilarityScore::backtracking () {
    
        std::vector <backtrack_align > bt;
        backtrack_align new_bt;
        int bt_size;
        new_bt.i = S.length()-1;
        new_bt.j = 0;        
        new_bt.block_count = prfl->blockCount()-1;
        new_bt.col_count = prfl->blockSize(new_bt.block_count)-1;
        bt.push_back(new_bt);
        
        while (!bt.empty()) { 
        
            bt_size = bt.size();
            
            for (int t = 0; t < bt_size; t++) {
            
                for (int k = S[bt[t].i][bt[t].j].prev.size() - 1; k > -1; k=k-1) {
                
                    new_bt = bt[t];
                    new_bt.j_prev = std::get<1>(S[new_bt.i][new_bt.j].prev[k]);
                    switch (std::get<2>(S[new_bt.i][new_bt.j].prev[k])) { 
                    
                        case 'm':
                            new_bt.align_seq += seq[S[new_bt.i].position + new_bt.j-1];
                            new_bt.align_seq += ' ';
                            if (new_bt.i < S.length()-1) {                        
                                new_bt.align_prfl.push_back(std::to_string(prfl[0][new_bt.block_count][new_bt.col_count][aa2int(seq[S[new_bt.i].position + new_bt.j-1])]));
                                new_bt.col_count = new_bt.col_count - 1;
                            }
                            else {
                                new_bt.align_prfl.push_back("*");
                            }
                            for (int s = S[new_bt.i].position + new_bt.j-1; s > S[new_bt.i-1].position + new_bt.j_prev; s=s-1) {
                                new_bt.align_seq += seq[s-1];
                                new_bt.align_seq += ' ';
                                new_bt.align_prfl.push_back("*");
                            }
                            
                            break;
                            
                        case 'p':
                            if (new_bt.i < S.length()-1) { 
                                new_bt.align_seq += seq[S[new_bt.i].position + new_bt.j-1];
                                new_bt.align_seq += ' ';
                                new_bt.align_prfl.push_back("_");
                            }
                            else {                            
                                for (int s = S[new_bt.i].position + new_bt.j; s > S[new_bt.i-1].position + new_bt.j_prev; s=s-1) {
                                    new_bt.align_seq += seq[s-1];
                                    new_bt.align_seq += ' ';
                                    new_bt.align_prfl.push_back("_");
                                }
                            }  
                            
                            break;
                            
                            case 's':
                                new_bt.align_seq += '_';
                                new_bt.align_seq += ' ';
                                new_bt.align_prfl.push_back(std::to_string(prfl[0][new_bt.block_count][new_bt.col_count][aa2int(seq[S[new_bt.i].position + new_bt.j-1])]));
                                new_bt.col_count = new_bt.col_count - 1;
                                for (int s = S[new_bt.i].position + new_bt.j; s > S[new_bt.i-1].position + new_bt.j_prev; s=s-1) {
                                    new_bt.align_seq += seq[s];
                                    new_bt.align_seq += ' ';
                                    new_bt.align_prfl.push_back("*");
                                }
                                            
                                break;
                                
                            case '_':
                                new_bt.align_seq += seq[S[new_bt.i].position+new_bt.j-1];
                                new_bt.align_seq += ' ';
                                new_bt.align_prfl.push_back(std::to_string(prfl[0][new_bt.block_count][new_bt.col_count][aa2int(seq[S[new_bt.i].position + new_bt.j-1])]));
                                new_bt.col_count = new_bt.col_count - 1;
                                for (int s = 0; s < S[new_bt.i].position - S[new_bt.i-1].position - 1; s++) {
                                    new_bt.align_seq += '_';
                                    new_bt.align_seq += ' ';
                                    new_bt.align_prfl.push_back("*");
                                }
                                
                                break;
                        }
                        
                        if (new_bt.col_count < 0) {
                            new_bt.block_count = (new_bt.block_count > 0) ? new_bt.block_count - 1 : 0; 
                            new_bt.col_count = prfl->blockSize(new_bt.block_count)-1;
                        }
                        new_bt.i = std::get<0>(S[new_bt.i][new_bt.j].prev[k]);            
                        new_bt.j = new_bt.j_prev; 
                        if (new_bt.i == 0 && S[new_bt.i].position+new_bt.j == 0) {
                            align_seq.push_back(new_bt.align_seq);
                            align_prfl.push_back(new_bt.align_prfl);
                            if (k == 0) {
                                bt.erase(bt.begin()+t);
                                t = t - 1;
                                bt_size = bt_size - 1;
                            }
                        }
                        else if (k > 0) {
                            bt.push_back(new_bt);
                        }
                        else {
                            bt[t] = new_bt;
                        }                            
                }
            }
      }       
    }
    
    void SimilarityScore::printAlignment () { //add number of alignments
        if (!align_prfl.empty() && !align_seq.empty()) {
            for(int i = 0; i < align_seq.size(); i++) {  
                reverse(align_seq[i].begin(), align_seq[i].end());
                reverse(align_prfl[i].begin(), align_prfl[i].end());
                std::cout<<align_seq[i]<<endl;
                std::cout<<align_prfl[i]<<endl;
                std::cout<<endl;
            }
         } 
    }   
    
    int SimilarityScore::aa2int (char c) { //use geneticcode.hh
        switch (c) {
            case '.':
                return 20;
            case 'G':
                return 0;
            case 'D':
                return 1;
            case 'E':
                return 2;
            case 'R':
                return 3;
            case 'K':
                return 4;
            case 'N':
                return 5;
            case 'Q':
                return 6;
            case 'S':
                return 7;
            case 'T':
                return 8;
            case 'A':
                return 9;
            case 'V':
                return 10;
            case 'L':
                return 11;
            case 'I':
                return 12;
            case 'F':
                return 13;
            case 'Y':
                return 14;
            case 'W':
                return 15;
            case 'H':
                return 16;
            case 'M':
                return 17;
            case 'C':
                return 18;
            case 'P':
                return 19;
        }		
        return -1;
    }
    
int main(int argc, char* argv[])
{
    if (argc != 3) 
    {
        std::cerr<<"Error (1): Wrong number of input arguments. \n";
        exit(1);
    }
    
    SimilarityScore SS(-5);
    SS.readFiles(argv[1],argv[2]);
    SS.fillSimilarityMatrix();
    SS.printSimilarityMatrix();
    std::cout<<endl;std::cout<<endl;      
    SS.backtracking();
    SS.printAlignment();

    return 0;
}
