/*****************************************************************************\
 * Filename : train_logReg_param.hh
 * Authors  : Lizzy Gerischer
 *
 * Description: implementation of logistic regression for training OE and EC 
 *              feature parameters using GSL 
 *              
 *
 * Date       |   Author              |  Changes
 *------------|-----------------------|------------------------------------------
 * 09.12.2016 | Lizzy Gerischer       | creation of the file
\******************************************************************************/
#ifndef _TRAINLOGREG_HH
#define _TRAINLOGREG_HH

// project includes

#include "types.hh"
#include "properties.hh"

// standard C/C++ includes
#include <iostream>
#include <vector>
#include <unordered_map>

#include <gsl/gsl_multimin.h>

using namespace std;

class train_data{

public:
  train_data(unordered_map<string, pair<int, vector<double> > > *samples, unordered_map<string,int> *ref_class, int numSp);
  vector<pair<int, vector<double> > > *exon_samples;
  vector<pair<int, vector<double> > > *intron_samples;
  int num_features;
  double prob_true_false_exons;
  double prob_true_false_introns;
  vector<int> exon_feature_idx;
  vector<int> intron_feature_idx;
  void set_mean_std(vector<pair<int, vector<double> > > *samples, vector<int> feature_idx);

};

double activation_f(const gsl_vector *theta, vector<double> *f);
double get_feature(vector<double> *f, int i);
double cross_entropy_error_f(const gsl_vector *theta, void *features);
void cross_entropy_error_df(const gsl_vector *theta, void *features, gsl_vector *df);
void cross_entropy_error_fdf(const gsl_vector *params, void *features, double *f, gsl_vector *df);
double rob_cross_entropy_error_f(const gsl_vector *theta, void *features);
void rob_cross_entropy_error_df(const gsl_vector *theta, void *features, gsl_vector *df);
void rob_cross_entropy_error_fdf(const gsl_vector *params, void *features, double *f, gsl_vector *df);
void optimize_parameters(train_data *data);
void reference_from_file(unordered_map<string,int> *ref_class);
void trainFeature_from_file();
void train_OEscore_params(int numSpecies);


#endif
