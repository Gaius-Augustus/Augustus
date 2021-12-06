/*
 * Name:    pp_simscore.cc
 * Project: similarity-score algorithm for block-profile and protein sequence
 * with intron information
 * Author:  Lars Gabriel
 *
 */

// project includes
#include "pp_simscore.hh"
#include "properties.hh"
#include "pp_profile.hh"
#include "fasta.hh"
#include "geneticcode.hh"

// standard C/C++ includes
#include <iostream>
#include <fstream>
#include <limits>
#include <math.h>
#include <iomanip>
#include <stdio.h>
#include <cstring>
#include <getopt.h>

using namespace SS;

    //initiation of a row with its length and the position of the
    //first element in the matrix
    Row::Row (int length, int position) {
        Row::len = length;
        Row::position = position;
        row = new Cell[length];
    }

    int Row::length () {
        return len;
    }
    //access a row by the cell number relative to the row
    Cell& Row::operator[] (int n) {
        if (0 <= n && n < len)
            return row[n];
        else if (n > len)
            throw string("pp_simscore: Row out of bound (" + \
            itoa(len) + " <= " + itoa(n) + ")");
        else
            throw string("pp_simscore: Row out of bound (" + \
            itoa(n) + " < 0)");
    }

    SimilarityMatrix::~SimilarityMatrix() {
        for (unsigned i = 0; i < matrix.size(); i++) {
            delete[] matrix[i].row;
        }
    }

    //access row of matrix with column number
    Row& SimilarityMatrix::operator[] (int n) {
        if (0 <= n && n < matrix.size())
            return matrix[n];
        else if (n > matrix.size())
            throw string("pp_simscore: Similarity matrix out of bound (" \
            + itoa(matrix.size()) + " <= " + itoa(n) + ")");
        else
            throw string("pp_simscore: Similarity matrix out of bound (" \
            +  itoa(n) + " < 0)");
    }

    //add a Row to the similarity matrix
    void SimilarityMatrix::addRow (int l, int p) {
        Row r(l,p);
        matrix.push_back(r);
    }

    //returns the number of columns of the similarity matrix
    int SimilarityMatrix::length () {
        return matrix.size();
    }

    bool SimilarityMatrix::empty () {
        return matrix.empty();
    }

    //read protein sequence from a file
    //the file has to be in FASTA format and it can have an additional
    //section with [intron position]
    ProteinSequence::ProteinSequence (const char* fileName) {
        std::string fastaFile = "";
        std::string line;
        int pos, frame;
        ifstream ifstrm;
        ifstrm.open(fileName);
        if (!ifstrm)
            throw string("Could not open input file \"" + string(fileName) + "\"");
        //read sequence
        while (std::getline(ifstrm, line)) {
            if (line[0] == '[')
                break;
            if (line.find_first_not_of("\t\n\v\f\r ") != string::npos)
                fastaFile += line + '\n';
        }
        if (fastaFile != "") {
            std::stringstream s(fastaFile);
            readOneFastaSeq(s, sequence, name, len);
        }
        //read intron positions
        while (std::getline(ifstrm, line)) {
            std::istringstream iss(line);
            if (iss >> pos >> frame) {
                while (num_prev_introns.size() < pos + 1) {
                    num_prev_introns.push_back(introns.size());
                }
                introns.insert(std::pair <int, int> (pos, frame));
            }
        }
        while (num_prev_introns.size() < len) {
            num_prev_introns.push_back(introns.size());
        }
        ifstrm.close();
    }

    ProteinSequence::~ProteinSequence () {
        delete[] name;
        delete[] sequence;
    }

    //access protein sequence by index
    char& ProteinSequence::operator[] (int n) {
        if (0 <= n && n < len)
            return sequence[n];
        else if (n > len)
            throw string("pp_simscore: protein sequence out of bound (" + \
            itoa(len) + " <= " + itoa(n) + ")");
        else
            throw string("pp_simscore: protein sequence out of bound (" + \
            itoa(n) + " < 0)");
    }

    int ProteinSequence::length () {
        return len;
    }

    //returns the number of residual nucleotides f, if an intron is in the
    //protein sequence at (j,f) located
    //returns -1, if at (j,{0,1,2}) no intron is located
    int ProteinSequence::intronAt (int j) {
        std::map < int, int >::iterator it = introns.find(j);
        return (it == introns.end()) ? -1 : it->second;
    }

    //returns the number of introns with position (j,{0,1,2}), where start<=j<end
    int ProteinSequence::intronsInRange (int start, int end) {
        if (start < 0)
            start = 0;
        if (end < 0)
            end = 0;
        return num_prev_introns[end]-num_prev_introns[start];
    }

    //initiates parameter
    SimilarityScore::SimilarityScore (double g, double b, double g_i, double iw1,\
        double iw2, double e_i, double e_n) {
        gap_cost_inter = g;
        gap_cost_intra = b;
        gap_cost_intron = g_i;
        intron_weight_intra = iw1;
        intron_weight_inter = iw2;
        epsi_intron = e_i;
        epsi_noIntron = e_n;
    }

    SimilarityScore::~SimilarityScore () {
        delete prfl;
        delete seq;
    }
    bool SimilarityScore::SimilarityMatrix_empty() {
        return S.empty();
    }

    //read profile and protein sequence
    void SimilarityScore::readFiles (const char* seqFileName, \
        const char* profileFileName) {
        //read protein sequence
        seq = new ProteinSequence(seqFileName);
        //read protein profile
        PP::initConstants();
        prfl = new PP::Profile(profileFileName);
    }

    //Poisson distribution
    double SimilarityScore::PoiDist (int k, double lambda) {
        return exp(k * log(lambda) - lambda - lgamma(k+1));
    }

    //returns log10 odds score that is weighted with intron_weight
    double SimilarityScore::intronScore (double q, double q_b) {
        return intron_weight_inter * (log10(q) - log10(q_b));
    }

    //Similarity Score Algorithm
    void SimilarityScore::fillSimilarityMatrix () {
        //the current row index in the similarity matrix
        int i = 0;
        //column position of the first element of the current row
        int position = 0;
        //row_length: length of a row
        //diff_pos: difference in position of current row and previous row
        //k_max: position of the last cell, which is used to compute
        //the score of a cell,if i represents the first position of a new block
        int row_length, diff_pos, k_max;
        //new_score holds the value of a possible new highest score for the current cell
        //old_score holds the value of the current maximal score of the current cell
        double new_score, old_score;
        //minimal space needed for the profile:
        //length of all blocks + all min interblock distances
        int min_prfl_length = 0;
        //number of blocks in the profile
        int block_count = prfl->blockCount();
        //col holds the current block column that is represented by the
        //current row of the similarity matrix
        PP::Column col;
        //d holds the inter-block distance interval that is located before
        //the current block of the profile
        PP::DistanceType d;

        //minimal gaps needed for an inter-block section
        int min_gaps;
        //match score for a cell
        double match_score;
        //alignment type in {m,s,p,-}
        char type;
        //length of an inter-block in the current alignment
        int interblock_length;
        //number of residual nucleotides
        int frame;
        //intron score
        double iscore;
        //relative frequency of intron at an intron position in B
        double intronFreqPos;
        //number of introns in the protein sequence
        int noIntronsSeq;

        //computes minimal profile length
        for (int k = 0; k < block_count; k++) {
            min_prfl_length += prfl->blockSize(k);
        }
        // The algorithm only runs, if the protein sequence is at least as long
        // as the number of block columns in the protein profile.
        if (seq->length() < min_prfl_length) {
            throw string("pp_simscore: Input protein sequence too short (" \
            + itoa(seq->length()) + " < " + itoa(min_prfl_length) + ")");
        }
        //constraint for the similarity score algorithm
        //at any point of the algorithm the whole protein-profile must be
        //able to be matched
        row_length = seq->length() - min_prfl_length + 1;
        intronPrfl_noProt = prfl -> getIntronNumSeq();
        //initial case of the similarity score recursion
        S.addRow(seq->length() + 1,0);
        S[0][0].score = 0;
        //fills first row
        for (int j = 1; j < seq->length() + 1; j++) {
            S[0][j].score = S[i][j-1].score + gap_cost_inter;
            if (seq->intronAt(j-2) != -1)
                S[i][j].score += gap_cost_intron;
            S[i][j].prev.push_back(make_tuple(i, j-1, 'p'));
        }
        position++;

        //loop over all blocks of B
        for (int t = 0; t < block_count; t++) {
            PP::Position P(t,0);
            //frequencies of amino acids in the first column of the new block
            col = prfl->getColumn(P);
            //interblock distance before the block
            d = prfl->interBlockDist(t);
            //position of the first cell of current row
            //position += d.r.min;
            //cout<<position<<endl;
            S.addRow(row_length, position);
            //i iterates over all rows of S
            i++;

            //score for the first position in a block
            for (int j = 0; j < row_length; j++) {
                //index of current column in the previous row
                diff_pos = S[i].position - S[i-1].position + j;
                //calculate min and max index of the previous row
                if (t == 0)
                    k_max = diff_pos + 1;
                else
                    k_max = ((row_length < diff_pos + 1) ? row_length : diff_pos + 1);
                old_score = std::numeric_limits<int>::min();

                    for (int k = 0; k < k_max; k++) {
                        //intronscore
                        iscore = 0.0;
                        if (intron_weight_inter > 0) {
                            noIntronsSeq = seq->intronsInRange(S[i-1].position + k - 1,\
                                position + j - 1);
                            intronFreqPos = prfl->getIntronInterFreq(t, noIntronsSeq);
                            if (noIntronsSeq == 0 && intronFreqPos < 0)
                                iscore = 0.0;
                            else {
                                if (intronFreqPos < 0)
                                    intronFreqPos = 0;
                                //calculate relative frequency with pseudocount
                                intronFreqPos = (intronFreqPos + epsi_intron)/\
                                    (intronPrfl_noProt + epsi_intron + epsi_noIntron);

                                interblock_length = (diff_pos - k );
                                if (t == 0 && j == 0)
                                    interblock_length -= 1;
                                //no intron position in protein sequence
                                if (interblock_length < 1)
                                    iscore = intronScore(intronFreqPos, 1);
                                //intron match
                                else
                                    iscore = intronScore(intronFreqPos, \
                                        PoiDist(noIntronsSeq, (d.r.min+d.r.max+ 2)/2.0 * \
                                        (intronAvgInterBlock_bfreq/intronPrfl_noProt) * 3));
                            }
                        }
                        //calculate score
                        min_gaps = 0;
                        type = 'm';
                        //calculate the minimal gaps needed for this
                        //inter-block section
                        if ((d.r.min - diff_pos + k + 1) > min_gaps) {
                            min_gaps = d.r.min - diff_pos + k + 1;
                            type = '_';
                        }
                        else if ((diff_pos - k - 1 - d.r.max) > min_gaps) {
                            min_gaps = diff_pos - k - 1 - d.r.max;
                            type = '_';
                        }
                        //compute new score
                        if (k == diff_pos) {
                            //first block column of B_t as an intra-block gap
                            new_score = S[i-1][k].score + gap_cost_inter *\
                                min_gaps + gap_cost_intra;
                            type = 's';
                        }
                        else {
                            //first block column of B_t as an intra-block gap
                            match_score = col.L(GeneticCode::get_aa_from_symbol(\
                                (*seq)[position+j-1]));
                            if (gap_cost_inter > match_score) {
                                new_score = S[i-1][k].score + gap_cost_inter *\
                                    min_gaps + gap_cost_intra;
                                type = 's';
                            }
                            else {
                                //first block column of B_t matches to an amino acid of P
                                new_score = S[i-1][k].score + gap_cost_inter *\
                                    min_gaps + match_score;
                            }
                        }
                        old_score = compareScores(new_score + iscore, old_score,\
                            i, j, i-1, k, type);
                    }
                    if (j > 0) {
                    //intronscore
                        //gap in profile
                        frame = seq->intronAt(position + j - 2);
                        if (frame < 0)
                            iscore = 0.0;
                        else
                            iscore = intronScore((epsi_intron / epsi_intron + epsi_noIntron),\
                                3 * intronIntraBlock_bfreq);
                        new_score = S[i][j-1].score + gap_cost_intra + iscore;
                        old_score = compareScores(new_score, old_score, i, j, i, j-1, 'p');
                    }
                S[i][j].score = old_score;
            }

            //score for all positions in a block except the first
            for (int s = 1; s < prfl->blockSize(t); s++) {
                PP::Position P(t,s);
                col = prfl->getColumn(P);
                position++;
                i++;
                S.addRow(row_length, position);
                for (int j = 0; j < row_length; j++) {
                    old_score = std::numeric_limits<int>::min();
                    if (j > 0) {
                        //gap in B
                        if (s == prfl->blockSize(t) - 1)
                            new_score = S[i][j-1].score + gap_cost_inter;
                        else
                            new_score = S[i][j-1].score + gap_cost_intra;
                        //intron score, if intron in P at this position
                        if (seq->intronAt(position + j - 2) != -1) {
                            for (int f = 0; f < 3; f++) {
                                new_score += gap_cost_intron;
                            }
                        }
                        old_score = compareScores(new_score, old_score, i, j, i, j-1, 'p');
                    }
                    if (j < row_length - 1) {
                        //gap in sequence
                        new_score = S[i-1][j+1].score + gap_cost_intra;
                        //intron score
                        for (int f = 0; f < 3; f++) {
                            new_score += intraBlock_iscore(t, s-1, f, -1);
                        }
                        old_score = compareScores(new_score, old_score, i, j, i-1, j+1, 's');
                    }
                    //match
                    new_score = S[i-1][j].score + col.L(GeneticCode::get_aa_from_symbol(\
                        (*seq)[position+j-1]));
                    frame = seq->intronAt(position + j - 2);
                    //intron score
                    for (int f = 0; f < 3; f++) {
                        new_score += intraBlock_iscore(t, s-1, f, frame);
                    }
                    old_score = compareScores(new_score, old_score, i, j, i-1, j, 'm');
                    S[i][j].score = old_score;
                }
            }
            position++;
        }
        //compute final score
        S.addRow(1, seq->length());
        old_score = std::numeric_limits<int>::min();
        d = prfl->interBlockDist(block_count);
        diff_pos = seq->length() - S[i].position;
        for (int j = 0; j < row_length; j++) {
            //INTRONSCORE
            iscore = 0.0;
            if(intron_weight_inter > 0) {
                noIntronsSeq = seq->intronsInRange(S[i].position + j - 1, seq->length() - 1);
                intronFreqPos = prfl->getIntronInterFreq(block_count, noIntronsSeq);
                if (noIntronsSeq == 0 && intronFreqPos < 0) {
                    iscore = 0.0;
                }
                else {
                    interblock_length = diff_pos - j ;
                    if (intronFreqPos < 0) {
                        intronFreqPos = 0;
                    }
                    //calculate relative frequency with pseudocount
                    intronFreqPos = (intronFreqPos + epsi_intron)/(intronPrfl_noProt + \
                        epsi_intron + epsi_noIntron);
                    if (interblock_length < 1) {
                        if (noIntronsSeq == 0)
                            iscore = 0.0;
                        //no intron position in protein sequence
                        else
                            iscore = intronScore(intronFreqPos, 1);
                    }
                    //intron match
                    else
                        iscore = intronScore(intronFreqPos, PoiDist(noIntronsSeq, \
                            (d.r.min+d.r.max+2)/2.0 * (intronAvgInterBlock_bfreq/ \
                                intronPrfl_noProt) * 3));
                }
            }
            //compute final similarity score
            min_gaps = 0;
            type = 'm';
            if ((d.r.min - diff_pos + j) > min_gaps) {
                min_gaps = d.r.min - diff_pos + j;
                type = 's';
            }
            else if ((diff_pos - j - d.r.max) > min_gaps) {
                min_gaps = diff_pos -j - d.r.max;
                type = '_';
            }
            new_score = S[i][j].score + gap_cost_inter * min_gaps;
            old_score = compareScores(new_score + iscore, old_score, i+1, 0, i, j, type);
        }
        S[i+1][0].score = old_score;
    }

    //returns the final similarity score
    double SimilarityScore::score () {
        if (!S.empty())
            return S[S.length()-1][0].score;
        else {
            throw string("pp_simscore: Similarity matrix is empty");
        }
    }

    //computes the intra-block intron score for intron positions (k,s,f)
    //in the profile and the aligned position in P
    //intron_frame denotes the number of residual nucleotides in P,
    //after which an intron occurs or -1 if there no intron is in P located
    double SimilarityScore::intraBlock_iscore (int k, int s, int f, int intron_frame) {
        double prfl_freq;
        //relative frequency of introns in B at this intron position
        prfl_freq = prfl->getIntronIntraFreq(k, s, f);
        if (intron_frame == f)
            return intron_weight_intra * (log10((prfl_freq + epsi_intron)/ \
                (intronPrfl_noProt + epsi_intron + epsi_noIntron)) - \
                log10(intronIntraBlock_bfreq));
        else if (intron_frame != -1 && prfl_freq > 0)
            return intron_weight_intra * (log10(1-((prfl_freq + epsi_intron)/ \
                (intronPrfl_noProt + epsi_intron + epsi_noIntron))) - \
                log10(1-intronIntraBlock_bfreq));
        else if (prfl_freq > 0)
            return intron_weight_intra * log10(1 - ((prfl_freq + epsi_intron)/ \
                (intronPrfl_noProt + epsi_intron + epsi_noIntron)));
        else
            return 0;
    }

    //compares two scores of similarity matrix and processes them
    double SimilarityScore::compareScores (double new_score, double old_score, \
        int current_i, int current_j, int prev_i, int prev_j, char score_type) {
        if (new_score < old_score)
            return old_score;
        if (new_score > old_score) {
            old_score = new_score;
            S[current_i][current_j].prev.clear();
        }
        S[current_i][current_j].prev.push_back(make_tuple(prev_i, prev_j, score_type));
        return old_score;
    }

    void SimilarityScore::backtrackingDB (int number_Alignments) {
        //stores for every alignment the connected regions of the same
        //alignment type in a block of the protein profile
        //for this the start coordinates in the protein profile
        //(block number, column number) and the protein sequence,
        //the length and the type for every region is stored
        int k_max;
        if (S.empty()) {
            throw string("pp_simscore: Similarity matrix is empty");
        }

        if (number_Alignments > 0) {
            char new_type;
            std::vector <BacktrackAlignDB > bt;
            BacktrackAlignDB new_bt;
            int bt_size;
            new_bt.i = S.length() - 1;
            new_bt.j = 0;
            new_bt.block_count = prfl->blockCount() - 1;
            new_bt.col_count = prfl->blockSize(new_bt.block_count) - 1;
            new_bt.col_posStart = new_bt.col_count;
            new_bt.align_element.length = 0;

            //initial case for the recursion: start a new alignment for every
            //predecessor of the final similarity score
            k_max = (number_Alignments < S[new_bt.i][new_bt.j].prev.size()) ? \
                number_Alignments : S[new_bt.i][new_bt.j].prev.size();
            for (int k = k_max - 1; k > -1; --k) {
                new_bt.align_element.type = std::get<2>(S[S.length() - 1][0].prev[k]);
                new_bt.j_prev = std::get<1>(S[S.length() - 1][0].prev[k]);
                new_bt.i = std::get<0>(S[S.length() - 1][0].prev[k]);
                new_bt.j = new_bt.j_prev;
                bt.push_back(new_bt);
            }
            //while there are unfinished alignments
            while (!bt.empty()) {
                bt_size = bt.size();
                //for every current alignment
                for (int t = 0; t < bt_size; t++) {
                    k_max = (number_Alignments < S[bt[t].i][bt[t].j].prev.size()) ? \
                        number_Alignments : S[bt[t].i][bt[t].j].prev.size();
                    //for the number of next backtracking paths
                    for (int k = k_max - 1; k > -1; --k) {
                        new_bt = bt[t];
                        new_bt.j_prev = std::get<1>(S[new_bt.i][new_bt.j].prev[k]);

                        new_type = std::get<2>(S[new_bt.i][new_bt.j].prev[k]);
                        //skip first interblock region
                        if (new_bt.i < S.length() - 2 && \
                        (new_bt.align_element.type != 'p' || new_bt.i < S.length() - 3)) {
                            //add new Element to current alignment, if the
                            //new element type is not equal to the last
                            //element type or a block is finished
                            if (new_bt.align_element.type != new_type || \
                                (new_bt.align_element.col_pos == 0 && \
                                    new_bt.align_element.type != 'p')) {
                                new_bt.alignment.push_back(new_bt.align_element);
                                new_bt.align_element.type = new_type;
                                new_bt.col_posStart = new_bt.col_count;
                                new_bt.align_element.length = 0;
                            }
                        }
                        else
                            new_bt.align_element.type = new_type;
                        //update alignment element parameters
                        new_bt.align_element.protSeq_pos = S[new_bt.i].position + new_bt.j - 1;
                        new_bt.align_element.block_pos = new_bt.block_count;
                        new_bt.align_element.col_pos = new_bt.col_count;
                        new_bt.align_element.length++;

                        if (new_type == 's' || new_type == '_' || \
                        (new_type == 'm' && new_bt.i < S.length() - 1))
                            new_bt.col_count = new_bt.col_count - 1;
                        //update current position in the similarity matrix and
                        //corresponding position in the protein profile
                        if (new_bt.col_count < 0) {
                            if (new_bt.block_count > 0) {
                                new_bt.block_count = new_bt.block_count - 1;
                                new_bt.col_count = prfl->blockSize(new_bt.block_count) - 1;
                            }
                            else {
                                new_bt.block_count = 0;
                                new_bt.col_count = 0;
                            }
                        }
                            new_bt.i = std::get<0>(S[new_bt.i][new_bt.j].prev[k]);
                            new_bt.j = new_bt.j_prev;
                        //an alignment is finished
                        if (new_bt.i == 0 && S[new_bt.i].position + new_bt.j == 0) {
                            new_bt.alignment.push_back(new_bt.align_element);
                            std::reverse(std::begin(new_bt.alignment), std::end(new_bt.alignment));
                            alignments.push_back(new_bt.alignment);

                            if (k == 0) {
                                bt.erase(bt.begin() + t);
                                t = t - 1;
                                bt_size = bt_size - 1;
                            }
                        }
                        //create new alignments, if there is more than one predecessor
                        //in the similarity matrix
                        else if (k > 0)
                            bt.push_back(new_bt);
                        //update current alignment
                        else
                            bt[t] = new_bt;
                    }
                }
            }
        }
    }

    //print Alignment in list of alignment elements
    void SimilarityScore::printAlignmentDB () {
        for (int a = 0; a < alignments.size(); a++) {
            for (int i = 0; i < alignments[0].size(); i++) {
                cout << a << '\t' << i << '\t' << alignments[a][i].protSeq_pos \
                << '\t' << alignments[a][i].block_pos << '\t' << \
                alignments[a][i].col_pos << '\t' << alignments[a][i].length << \
                '\t' << alignments[a][i].type << endl;
            }
        }
    }

    //translation from block indices of data structure to block position in file
    void SimilarityScore::printBlockPos () {
        for (int i = 0; i < prfl->blockCount(); i++) {
                cout << i << '\t' << prfl[0][i].blockNumbInFile << endl;
        }
    }

    //prints the similarity matrix
    void SimilarityScore::printSimilarityMatrix () {
        if (S.empty()) {
            throw string("pp_simscore: Similarity matrix is empty");
        }
        std::cout.precision(2);
        std::cout << std::fixed;
        for (int i = 0; i < S.length(); i++) {
            for (int j = 0; j < S[i].position; j++) {
                std::cout << std::left << setw(7) << "-";
            }
            for (int j = 0; j < S[i].length(); j++) {
                //print the score for every cell
                std::cout << std::left << setw(7) << S[i][j].score;
                //print the coordinates of the predecessor cells
                for (int t=0; t<S[i][j].prev.size(); t++) {
                    std::cout << std::get<0>(S[i][j].prev[t]) << ":" << \
                    std::get<1>(S[i][j].prev[t]) << ":" << \
                    std::get<2>(S[i][j].prev[t]) << ",";
                }
                std::cout << " | ";
            }
            for (int j = 0; j < seq->length() + 1 - S[i].length() - S[i].position; j++) {
                std::cout << std::left << setw(7) << "-";
            }
            std::cout << endl;
        }
    }

    //computes an readable alignment representation to stdout
    void SimilarityScore::printAlignment () {
        if (alignments.empty()) {
            throw string("pp_simscore: No alignment was determined");
        }
        int number_of_decimals = 4;
        std::vector <int> prob_index;
        std::string decimals;
        std::string seq_align;
        PrflAlignment prfl_align;
        int prot_index;
        int pos = 0;
        int row_numb;
        int char_numb = 60;

        //compute alignment
        for (int i = 0; i < alignments.size(); i++) {
            seq_align = "";
            prfl_align.argmax_aminoacids = "";
            prfl_align.match_prob.clear();
            for (int j = 0; j < alignments[i].size(); j++) {
                if (j == 0) {
                    seq_align += to_string( alignments[i][j].protSeq_pos );
                    prfl_align.argmax_aminoacids += to_string( alignments[i][j].protSeq_pos);
                }
                else if (alignments[i][j].block_pos != alignments[i][j-1].block_pos) {
                    seq_align += to_string( alignments[i][j].protSeq_pos - pos);
                    prfl_align.argmax_aminoacids += to_string( alignments[i][j].protSeq_pos - pos);
                }
                pos = alignments[i][j].protSeq_pos + alignments[i][j].length;
                for (int k = 0; k < alignments[i][j].length; k++) {
                    if (alignments[i][j].type =='s') {
                        seq_align += "_";
                        prfl_align.argmax_aminoacids += prfl[0][alignments[i][j].block_pos][alignments[i][j].col_pos + k].argmax();
                    }
                    else if (alignments[i][j].type =='p') {
                        seq_align += (*seq)[alignments[i][j].protSeq_pos + k];
                        prfl_align.argmax_aminoacids += "_";
                    }
                    else {
                        seq_align += (*seq)[alignments[i][j].protSeq_pos + k];
                        prfl_align.argmax_aminoacids += prfl[0][alignments[i][j].block_pos][alignments[i][j].col_pos + k].argmax();
                        prfl_align.match_prob.push_back(prfl[0][alignments[i][j].block_pos][alignments[i][j].col_pos+k][GeneticCode::get_aa_from_symbol((*seq)[alignments[i][j].protSeq_pos + k])]);
                    }
                }
                if (j == alignments[i].size() - 1) {
                    seq_align += to_string(seq->length() - pos);
                    prfl_align.argmax_aminoacids += to_string( seq->length() - pos);
                }
            }
            alignmentsReadable.push_back(std::make_pair(seq_align, prfl_align));
        }

        //print alignment
        //output format see description in header file
        for (int i = 0; i < alignmentsReadable.size(); i++) {
            prob_index.clear();
            for (int i = 0; i < number_of_decimals + 1; i++){
                prob_index.push_back(0);
            }
            std::cout << "ALIGNMENT "<< i << ":" << endl;
            row_numb = alignmentsReadable[i].first.length() / char_numb + 1;
            for (int k = 0; k < row_numb; k++) {
                std::cout << std::left << std::setw(25) << "protein sequence:";
                std::cout << alignmentsReadable[i].first.substr(k * char_numb, char_numb);
                std::cout << endl;
                std::cout << std::left << std::setw(25) << "protein profile argmax:";
                std::cout << alignmentsReadable[i].second.argmax_aminoacids.substr(k * char_numb, char_numb);
                std::cout << endl;

                for (int d = 1; d < number_of_decimals + 1; d++) {
                    if  (d == 2)
                        std::cout << std::left << std::setw(25) << \
                        "protein profile prob.";
                    else if (d == 3)
                        std::cout << std::left << std::setw(25) << \
                        "of prot. seq. character:";
                    else
                        std::cout << std::left << std::setw(25) << "";

                    prot_index = k * char_numb;
                    while (prot_index < (k+1) * char_numb && prob_index[d-1] < \
                            alignmentsReadable[i].second.match_prob.size()) {
                        if (std::isdigit(alignmentsReadable[i].first[prot_index]) || \
                            alignmentsReadable[i].first[prot_index] == '_' || \
                            alignmentsReadable[i].first[prot_index] == 's' || \
                            alignmentsReadable[i].second.argmax_aminoacids[prot_index] == '_')
                            std::cout << ' ';
                        else {
                            std::cout << to_string(alignmentsReadable[i].second.match_prob[prob_index[d-1]])[d];
                            prob_index[d-1]++;
                        }
                        prot_index++;
                    }
                    std::cout << endl;
                }
                std::cout << endl;
                std::cout << endl;
            }
            std::cout << endl;
            std::cout << endl;
        }
    }

    pair<double, double> SimilarityScore::avgArgMaxProb () {
        int prfl_length = 0;
        double consents = 0;
        for (int k = 0; k < prfl->blockCount();k++){
            for (int s = 0; s < prfl->blockSize(k); s++) {
                consents += prfl[0][k][s].max_val();
                prfl_length += 1;
            }
        }
        return make_pair(consents, prfl_length);
    }

    void SimilarityScore::printInterBlock() {
        PP::DistanceType d;
        for (int k = 0; k < prfl->blockCount() + 1; k++){
            d = prfl->interBlockDist(k);
            cout << d.r.min << '\t' << d.r.max << endl;
        }
    }



void printHelp() {

    std::cout <<
" Algorithm for calculating the similarity score and the optimal alignments between a protein profile and a protein sequence: \n\
The algorithm can optional take intron positions into account.\n\n\
USAGE: ./pp_simScore [OPTIONS] --fasta protein_sequence_file --prfl protein_profile_file\n\n\
 [OPTIONS]:\n\
    --gap_inter:                gap costs for an alignment column that is a gap in an inter-block section\n\
                                default setting: -5\n\n \
    --gap_intra:                gap costs for an alignment column that is a gap in a block \n\
                                default setting: -50\n\n\
    --gap_intron:               gap costs for a gap in intron positions\n\
                                default setting: -5\n\n\
    --epsilon_intron:           pseudocount parameter epsilon1, the pseudocount is added to a relative intron frequency v/w with (v+epsilon1)/(w+epsilon1+epsilon2)\n\
                                default setting: 0.0000001\n\n\
    --epsilon_noIntron:         pseudocount parameter epsilon2, the pseudocount is added to a relative intron frequency v/w with (v+epsilon1)/(w+epsilon1+epsilon2)\n\
                                default setting: 0.1\n\n\
    --intron_weight_intra:      value that is added to an intron score for a match of intron positions in a block\n\
                                default setting: 5\n\n\
    --intron_weight_inter:      value that is added to an intron score for a match of intron positions in an inter-block\n\
                                default setting: 5\n\n\
    --alignment:                number of optimal alignments that are computed\n\
                                default setting: 1\n\n\
    --help:                     print USAGE\n\n\
    --out:                      denotes the output format, the following output options, between \" \", are implemented:\n\n\
                                       \"score\" : output is the similarity score\n\
                                       \"matrix\" : output are similarity matrix and similarity score\n\
                                       \"alignment\" : output are the computed alignments to the console as\n\n\
                               Alignment representation of P as symbols of {AminoAcid, gap symbol or number of amino acids in inter-block}\n\
                               Alignment representation of argmax of B as symbols of {argmax AminoAcid for aligned block column, gap symbol or inter-block length}\n\
                               Frequency of amino acid of P in aligned block column of B, if alignment type is a match\n\n\
                                        \"matrix+alignment\": output are similarity matrix, \n\
                                                            similarity score and the computed \n\
                                                            alignments in the format described above\n\
                                        \"db\" : output are the computed alignment \n\
                                               as list of alignment frames, \n\
                                               an element of the list consists of:\n\n\
                               - starting position of the first amino acid of the protein sequence that is included in the alignment frame\n\
                               - block number in which the alignment frame is located\n\
                               - index of the first block column that is included in the alignment frame\n\
                               - length of the frame (number of alignment columns)\n\
                               - alignment type: 'm', 's'. 'p' or '-'\n\n\
                                        \"bp\" : output is a list of translations from the index of a block \n\
                                               to the number of the block in the .prfl file\n\
                                        \"consents\" : output is the average of the argmax \n\
                                                     of the block columns for the complete profile\n\
                                        \"interblock\" : output is a list of all inter-block distance intervals\n\
                                default setting: \"score\" \n";
}

int main(int argc, char* argv[])
{
    std::string fasta, prfl, db;
    std::string out_option = "score";
    double gapCost_inter = -5;
    double gapCost_intra = -50;
    double gapCost_intron = -5;
    double intron_weight_inter = 5;
    double intron_weight_intra = 5;
    double epsilon_intron = 0.0000001;
    double epsilon_noIntron = 0.1;
    int alignmentNo = 1;
    const char* const short_opts = "g:b:r:e:n:o:f:p:i:t:a:h";
    pair <int, int> count_miss;

    const option long_opts[] = {
        {"gap_inter", required_argument, nullptr, 'g'},
        {"gap_intra", required_argument, nullptr, 'b'},
        {"gap_intron", required_argument, nullptr, 'r'},
        {"epsilon_intron", required_argument, nullptr, 'e'},
        {"epsilon_noIntron", required_argument, nullptr, 'n'},
        {"out", required_argument, nullptr, 'o'},
        {"fasta", required_argument, nullptr, 'f'},
        {"prfl", required_argument, nullptr, 'p'},
        {"intron_weight_intra", required_argument, nullptr, 'i'},
        {"intron_weight_inter", required_argument, nullptr, 't'},
        {"alignment", required_argument, nullptr, 'a'},
        {"help", no_argument, nullptr, 'h'},
    };

    while (true)
    {
        const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);
        if (-1 == opt)
            break;
        switch (opt) {
            case 'g':
                gapCost_inter = std::stod(optarg);
                break;
            case 'b':
                gapCost_intra = std::stod(optarg);
                break;
            case 'r':
                gapCost_intron = std::stod(optarg);
                break;
            case 'e':
                epsilon_intron = std::stod(optarg);
                break;
            case 'n':
                epsilon_noIntron = std::stod(optarg);
                break;
            case 'o':
                out_option = optarg;
                break;
            case 'f':
                fasta = optarg;
                break;
            case 'p':
                prfl = optarg;
                break;
            case 'a':
                alignmentNo = std::stoi(optarg);
                break;
            case 'i':
                intron_weight_intra = std::stod(optarg);
                break;
            case 't':
                intron_weight_inter = std::stod(optarg);
                break;
            case 'h':
            case '?':
            default:
                printHelp();
                break;
        }
    }

    if (!fasta.empty() && !prfl.empty()) {
        SimilarityScore SS(gapCost_inter, gapCost_intra, gapCost_intron, \
            intron_weight_intra, intron_weight_inter, epsilon_intron, epsilon_noIntron);
        SS.readFiles(fasta.c_str(), prfl.c_str());
        if (out_option == "consents") {
            pair < double, double > output = SS.avgArgMaxProb();
            cout << output.first << '\t' << output.second;
        }
        else if (out_option == "interblock")
            SS.printInterBlock();
        else if (out_option == "bp")
            SS.printBlockPos();
        else {
            SS.fillSimilarityMatrix();
            if (out_option == "score")
                std::cout << SS.score() << endl;
            else {
                SS.backtrackingDB(alignmentNo);
                if (out_option == "matrix") {
                    SS.printSimilarityMatrix();
                    std::cout << endl << endl << "Final Score:  " << SS.score() << endl;
                }
                else if (out_option == "alignment") {
                    std::cout << endl << endl << "Final Score:  " << SS.score() << endl;
                    SS.printAlignment();
                }
                else if (out_option == "matrix+alignment") {
                    SS.printSimilarityMatrix();
                    std::cout << endl << endl << "Final Score:  " << SS.score() << endl;
                    SS.printAlignment();
                }
                else if (out_option == "db")
                    SS.printAlignmentDB();
            }
        }
    }
    else {
        std::cout << "invalid input arguments\n\n";
        printHelp();
    }
    return 0;
}
