#ifndef _EXON_SEG_HH
#define _EXON_SEG_HH

#include<iostream>
#include<fstream>
#include<math.h>
#include<getopt.h>
#include<vector>
#include<map>
#include<stdlib.h>
#include<limits.h>
using namespace std;


/* Minus_infinity is a very large number to use as a sign of impossible transitions
 * in the transition probability matrix
 */
const long int minus_infinity=LONG_MIN;
/* maxcov stores the maximum coverage depth that we are considering in storing the distributions
 * numstates store the number of states 
 */
const int maxcov =10000;
#define NUMSTATES 4
/* EXONP denotes positive exon EXONM denotes the negative exon INTRONJ and INTRONI denote the 
 * J and I lables of introns Together these four are teh states that we are considering
 * in our HMM
 */
const int EXONP =0;
const int EXONM =1;
const int INTRONI =2;
const int INTRONJ =3;
/* STRANDP STRANDM and STRANDB denote the positive negative and unknown strand respectively
 */
const int STRANDP =0;
const int STRANDM =1;
const int STRANDB =2;
/* pott_gamma denotes the gamma value that we will be using in calulating the
 * Pott's functionals
 */
const int pott_gamma =1500;
/* moving_window stores the window size that we will be using while finding the 
 * splice sites using the template matching technique
 */
const int moving_window=70;
/* The convergence limit that we use for calculating the labda values using the iterative technique
 */
const double convergence_limit=0.01;
/* L stores the threshold betweeen introns and exons which we use when we try to estimate the
 * the distribution using the train_function which calculates the distribution using the exons
 * predicted from a gff file
 */
const int L=4;

/* lengths stores the length for the states 
 */
const int lengths[]={1,1,38,1};

/* fragment is a structure which stores small streches of the hints generated in the first part and also
 * the necessary information associated with them
 */
typedef struct fragment{
  int state,start,end;
  double avg_cov;
  /* cov_exp stores the coverage depth for that stretch
   * in all the experiments
   */
  vector<double> cov_exp;
} fragment;


/* It is a small class to store the information extracted from the wig file 
 */
class dataset{
  
public:
  /* define the map used to store the data
   * The top most vector stores different strands 
   * the next one stores the different tracks
   * The next vector contains the score at each position in the track
   */
  map<string, vector< vector< vector<int> > > > database;
  /* define the iterator to loop through the map
   */
  map<string, vector< vector< vector<int> > > >::iterator it;
  int max_depth;
  int max_len;
  int no_of_tracks;
  vector<string> name_of_tracks;
  double avg_depth;
} ;

class exon_segmentation{
/* Variable lambda1 stores the parameter for the distribution of intron
   * and lambda2 stores the parameter for the distribution of exons
   */
  double lambda1;
  double lambda2;
  /* stores the probability distribution for each experiment
   */
  vector<vector<double> > prob_dist_exon;
  vector<vector<double> > prob_dist_intron;


public:
  //constructor
  exon_segmentation();
  /* function to calculate the geometric distribution parameters (lambdas)
   * which will be used later to calculate the distributions. 
   */
  void calculate_lambda(dataset &coverage_info, int c,float s);
  /* calculates the geometric distribution values using the lambda values and stores them
   */
  void prob_density(int max,int no_of_tracks);
  /* this function calculates the probability distribution values using an input gff file for training
   */
  void train_function(dataset &coverage_info,string inputfile,string key);
  /* functions which return the distribution value for a particular coverage depth
   * by accessing the corresponding index in the distribution table
   */
  double prob_value_exon(int num,int track);
  double prob_value_intron(int num,int track);
  /* Find the splice sites
   * INCOMPLETE
   */
  void find_splice(vector< vector< vector<int> > > &input_set,dataset &coverage_info);
  /* reads the file and stores it in the the object of the class dataset
   */
  void read_file(dataset &coverage_info,string filename);
  /* it takes as input pointer to a 2-D vecctor, 
   * the input dataset for a chromosome and
   * stores the emission probability values in it
   */
  void calculate_emissions( vector< vector<double> > &emission_probs,int max_len,vector< vector< vector<int> > > &input_set,int chunksize,int iteration);
  /* takes as input a row and column and calculates the value to be filled 
   * in that element in the gamma matrix
   * there is also a pointer spred to store the predecessor state in it
   */
  void fill_viterbi(vector< vector<double> > &gamma, int numstates,int j, int s,double emission_val,double a[][NUMSTATES],int* spred);
  /* this function implements the viterbi algorithm
   */
  vector<int> viterbi(vector< vector< vector<int> > > &input_set,int max_score,int max_len,int chunksize,int iteration);
  /* function to convert the state sequence into segments and also to store the avg coverage depth
   */
  vector< fragment > segment( vector<int> state_seq,vector< vector< vector<int> > > &input_set,int include_intron,int no_of_tracks );
  /* This function takes input the averrage coverage depth and transforms them
   * into some usable form in the pott's functional
   */
  double pott_convert(double d);
  /* calculate m[a,b]=(psum[b]-psum[a-1])/(psuml[b]-psuml[a-1])
   * to be used for calculating pott's functionals
   */
  double calculate_m(int a,int b,vector<vector<int> > &psuml,vector<vector<double> > &psum,int exp);
  /* Calculate S(a,b) to be used in the pott's functional
   */
  double calculate_s(int a,int b,vector<vector<int> > &psuml,vector<vector<double> > &psum,vector<vector<double> > &psqrsum,vector<vector<double> > &pcubesum,vector<fragment> &segments,int exp);
  /* the function which takes as input the vector of fragments and clusters the exons together to form genes
   */
  vector<int> cluster_exons(vector<fragment> segments,int no_of_tracks);

  /* Takes the dataset, the name of output file and the threshold c and the ratio of complement 
   * of the paramters(s) as input and makes the gff file
   */
  void make_gff(dataset &coverage_info,string outputfile,int c,float s,string train_file,int flag_r,string bed_file,string augustus_hints,int chunksize);
  /* takes the dataset and vector of fragments as input and
   * stores the distribution of number of positions with the 
   * coverage depths for both introns and exons 
   * in separate files for plotting
   */
  void make_graph(dataset &coverage_info,vector<fragment > &segments);
};

#endif
