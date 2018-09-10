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

struct stdValues {
  stdValues(): mean(0), std(0), max(numeric_limits<double>::min()), min(numeric_limits<double>::max()), numSamples(0) {}
  double mean;
  double std;
  double max;
  double min;
  int numSamples;
};

class TrainData{

public:
  TrainData(unordered_map<string,int> *ref_class, int numSp, string type);

  bool isExon;
  vector<vector<double> > samples;
  unordered_map<string,int> *ref_class;
  vector<vector<bool> > traitExists;
  int num_features;
  double prob_true_false_samples;
  void setFeatures(vector<double> *f, Traits *t);
  void calculateBinQuantils();
  vector<stdValues> std_values; // identified by feature idx, containing mean and standard deviation  
  vector<bool> to_be_std;
  void standardizeFeatures();
  void setNumberOfFeatures();
  void printParams();
};


class LogReg{
public:

  static double activation_f(const gsl_vector *theta, vector<double> *f);
  //static double get_feature(vector<double> *f, int i);
  //static double get_binning_feature(int i, double s);
  static double cross_entropy_error_f(const gsl_vector *theta, void *features);
  static void cross_entropy_error_df(const gsl_vector *theta, void *features, gsl_vector *df);
  static void cross_entropy_error_fdf(const gsl_vector *params, void *features, double *f, gsl_vector *df);
  static double rob_cross_entropy_error_f(const gsl_vector *theta, void *features);
  static void rob_cross_entropy_error_df(const gsl_vector *theta, void *features, gsl_vector *df);
  static void rob_cross_entropy_error_fdf(const gsl_vector *params, void *features, double *f, gsl_vector *df);
  static void optimize_parameters(TrainData *data);
  static void reference_from_file(unordered_map<string,int> *ref_class);
  static void trainFeature_from_file();
  static void train_OEscore_params(int numSpecies);
  static void sgd(vector<vector<double> > *samples, gsl_vector *theta, int n);
  static void gsl_minimizer(vector<vector<double> > *samples, gsl_vector *theta, int n);
  static void updateWeights(const gsl_vector *theta, TrainData *data);
  static void backtransformWeigths(gsl_vector *theta, TrainData *data);

};
#endif
