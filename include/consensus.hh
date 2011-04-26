/**********************************************************************
 * file:    consensus.hh
 *
 * descr.:  header file for consensus.cc which finds the consensus pattern
 * authors: Sukadeb Acharya, sukadeb@gmail.com
 *
 **********************************************************************/

#ifndef _CONSENSUS_HH
#define _CONSENSUS_HH

#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

#include "geneticcode.hh"

using namespace std;

class consensus;

//stores the input sequences in this structure
typedef struct sequence{
	string seq;
	int length;
} sequence;


//stores the final output in this structure
typedef struct histogram_data{
  string consensus_pattern;
  vector<string> neighbours;
  vector<int> freq_neighbours;
  vector<float> mean_neighbours;
  int f_value;
  float mean_value;
  vector<int> position_of_occurence;
} histogram_data;


class consensus{
  /* pattern size is the input pattern size
   * max_allowed_mismatch is the number of mismatches to consider while calculating the neighbours
   * number_seq stores the number of sequences
   * starting and ending store the starting and ending point in the sequence to consider for our window
   * no_iterations specifies the number of consensus patterns we need to store for the output
   * delta is used to calculate whether a certain string is relevant or not
   * p_value is the threshold used to identify the significant strings
   * powers is the array used to store the powers 4 for caculating the character to int conversion and vice versa
   * tpm is the transition probabiliry matrix
   * var is the vector of structure sequence containing all the input sequences
   * t_values is a flag array used to store the actuall frequencies of the patterns
   * f_values stores the actual frequencies of the patterns
   * back_probs store the background probabilities
   * significant_strings and relevant_strings store the significant and relevant strings respectively
   */
  int pattern_size,max_allowed_mismatch,number_seq,starting,ending,no_iterations;
  float delta,p_value;
  vector<double> powers;
  float tpm[4][4];
  vector<sequence> var;
  vector<int> t_values,f_values;
  vector<float> r_values,mean_values;
  float back_probs[4];
  vector<int> significant_strings,relevant_strings;
  string analyse_pattern;

  /* consensus_data stores the consensus patterns as their corresponding integers
   * final_list is a vector of struccture histogram_data
   * max_string_length is the maximum length of the input string sequences
   * max_freq is the maximum frequency of occurence at a particular position and is used to scale the histogram
   */
  vector<int> consensus_data;
  vector<histogram_data> final_list;
  int max_string_length;
  int max_freq;

public:

  consensus(int starting1,int ending1);
  /* sets the parameters pattern size, p_value to calculate the significant strings, delta to calculate the relevant 
   * strings which is the ratio of the actual frequency to the observed frequency, the maximum mismatch allowed 
   * to find the neighbours and the number of iterations n we want to perform to get the consensus patterns
   */
  void set_values(int pattern_size1,float p_value1, float delta1, int max_allowed_mismatch1,int n);

  /* takes the file name to store the sequences
   */
  void set_file_name(string filename);
  /* adds a particular string to consider for caculation
   */
  void add_string(string fileline);
  /* takes as input the pattern that we want to analyse and also the required parameters and does all the necessary 
   * calculation only for the given pattern
   */
  void analyse(string pattern1,float p_value1,float delta1, int max_allowed_mismatch1);
  /* calculates the cumulative value using normal distribution
   */
  float calculate_cumulative_value(float z);
  /* calculates the poisson threshold for a given frequency occurrence i and threshold value p 
   */
  double poisThresh(int i, float p);
  /* calculates the poisson inverse between values a and b with threshold value p , number of iterations to perform n and 
   * frequency of occurence i
   */
  double poisInv(double a,double b, float p,int n,int i);
  /* starts the calculations
   */
  void start();
  /* plots the histogram
   */
  void plot_histogram();
  /* returns the consensus pattern and other reltaed information
   */
  vector<histogram_data> get_consensus();
  /* prints the consensus pattern and other related information
   */
  void print_consensus();
  /* returns the factorial of m
   */
  int factorial( int m);
  /* returns n choose r
   */
  int nCr(int n, int r);
  /* converts a character from A C G T to coresponding number value 0 1 2 3 respectively
   */
  int char2num(char m);
  /* converts the numbers to corresponding character
   */
  char num2char(int a);
  /* recursive function used by find_neighbours function to find k mismatch neighbours which uses the newly formed number
   * the actual number for which we want the neighbours, number of mismatches k, length of pattern, object Seq2Int for
   * performing the calculations and the vector where the new neighbours is to be stored
   */
  void k_mismatch(int new_m, int old_m, int k, int length,Seq2Int s2i,vector<int> &neighbours_index);
  /* returns a vector containing k mismatch neighbours of pattern with interger value m and length as specified
   */
  vector<int> find_neighbours(int m,int k, int length);
};

#endif
