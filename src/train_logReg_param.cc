 /*****************************************************************************\
 * Filename : train_logReg_param.cc
 * licence  : Artistic Licence, see file LICENCE.TXT or 
 *            http://www.opensource.org/licenses/artistic-license.php
 * Authors  : Mario Stanke, mario.stanke@uni-greifswald.de
 *
 *
 * Date       |   Author              |  Changes
 *------------|-----------------------|------------------------------------------
 * 09.12.2016 | Lizzy Gerischer       | creation of the file
 *
 ******************************************************************************/

// project includes

// standard C/C++ includes
#include <iostream>
#include <fstream>
#include <iomanip>
#include <regex>
#include <random>

#include "train_logReg_param.hh"

int num_species;
double prob_true_false;


void LogReg::train_OEscore_params(int numSpecies){

  try{
    Constant::alpha = Properties::getdoubleProperty("/CompPred/alpha");
  }catch(...){}
    
  unordered_map<string,int> ref_class;
  reference_from_file(&ref_class);

  TrainData exon_data(&ref_class, numSpecies, "exon");
  optimize_parameters(&exon_data);
  TrainData intron_data(&ref_class, numSpecies, "intron");
  optimize_parameters(&intron_data);

  // output file for trained log_reg parameter
  ofstream paramfile;
  if(Properties::hasProperty("param_outfile"))
    paramfile.open(Properties::getProperty("param_outfile"));
  else if(Constant::cgpParsFile != "")
    paramfile.open(Constant::cgpParsFile);
  else
    paramfile.open("log_reg_train.pars");

  for(int i = 0; i < Constant::lr_features.size(); i++){
    Constant::lr_features[i]->print(&paramfile);
  }

  paramfile.close();
}

void TrainData::printParams(){
    // test if exon score parameters are read in properly
    cout << "features used for CGP predition and training" << endl << "id\tdescription\tstd?\tnumBins\tbinning_boundaries\tweights" << endl;
    for(int i=0; i<Constant::lr_features.size(); i++){
      Constant::lr_features[i]->print();
    }
}



TrainData::TrainData(unordered_map<string,int> *rc, int numSp, string type){

  //  number of features
  
  if(type == "exon")
    isExon = true;
  else 
    isExon = false;

  ref_class = rc;

  calculateBinQuantils();
  setNumberOfFeatures();
    
  printParams();

  int numRef = 0;
  num_species = numSp;
  
  std_values.resize(num_features, stdValues());

  // store all features in huge data matrix and standardize
  for(auto samp = Constant::logReg_samples.begin(); samp != Constant::logReg_samples.end(); samp++){
    
    if((samp->first.at(0) == 'C' && isExon) || (samp->first.at(0) != 'C' && !isExon) ){
      
      vector<double> features;
     
      auto got = ref_class->find(samp->first);
      if(got == ref_class->end()){
	features.push_back(0);
      }else{
	features.push_back(1);
	numRef++;
      }
      traitExists.push_back(vector<bool>(1,false));
      setFeatures(&features, &samp->second);
      samples.push_back(features);
      
    }
  }

  /*   
    cout << "samples:" << endl;
  for(int i=0; i<num_features; i++)
    cout << i << "\t";
  cout << endl;
  for(int i=0; i<samples.size(); i++){
    cout << type <<"_sample\t";
    for(int j=0; j<samples[i].size(); j++){
      cout << samples[i][j] << "\t";
    }
    cout << endl;
  }
  */

  standardizeFeatures();

  cout << "# " << numRef << " reference " << type << " samples in training set." << endl;
  cout << "# " << samples.size() << " " << type << " candidates are used for training." << endl;

  prob_true_false_samples = (double) numRef / samples.size();
  prob_true_false = 1;

  /*
  cout << "samples standardized:" << endl;
    for(int i=0; i<num_features; i++)
    cout << i << "\t";
  cout << endl;
  for(int i=0; i<samples.size(); i++){
    cout << type <<"_samp_std\t";
    for(int j=0; j<samples[i].size(); j++){
      cout << samples[i][j] << "\t";
    }
    cout << endl;
  }
  */
}

bool sort_samples(pair<double, double> p1, pair<double, double> p2){
  return(p1.first < p2.first);
}


void TrainData::calculateBinQuantils(){

  // calculate bin quantils

  for(int i=0; i < Constant::lr_features.size(); i++){
    if( (isExon && Constant::lr_features[i]->id < 20) || (!isExon && Constant::lr_features[i]->id >= 20) ){
      int numBins = Constant::lr_features[i]->num_bins;
      if(numBins > 1){
	int num_samples = 0;
	int numRef = 0;
	vector<pair<double,double> > samples;
	LRbinnedFeatureGroup *lrb = static_cast<LRbinnedFeatureGroup*>(Constant::lr_features[i]);
	for(auto lr_it = Constant::logReg_samples.begin(); lr_it != Constant::logReg_samples.end(); lr_it++){	
	  if( (lr_it->first[0] == 'C' && isExon) || (lr_it->first[0] != 'C' && !isExon) ){
	    bool hasFeature = true;
	    double x = Properties::getFeature(lrb->id, &lr_it->second, &hasFeature);
	    if(hasFeature){ 
	      double widx;
	      auto got = ref_class->find(lr_it->first);
	      if(got == ref_class->end())
		widx = 0;
	      else{
                widx = 1;
		numRef++;
	      }
	      samples.push_back(pair<double, double>(x,widx));	      
	      num_samples++;
	    }
	  }
	}
	lrb->bb.clear();
	sort(samples.begin(), samples.end(), sort_samples);
	double w_true = (double) 1/numRef/2;
	double w_false = (double) 1/(num_samples-numRef)/2;
	cout << "numSamples = " << num_samples << "\tnumRef = " << numRef << "\tw_true = " << w_true << "\tw_false = " << w_false << endl;
	if(samples[0].second == 1)
	  samples[0].second = w_true;
	else
	  samples[0].second = w_false;

	for(int j=1; j<samples.size(); j++){
     	  if(samples[j].second == 1)
	    samples[j].second = samples[j-1].second + w_true;
	  else
	    samples[j].second = samples[j-1].second + w_false;
	}
	double qdiff = (double)1/numBins;
	double qmax = 1 - qdiff;
	double qidx = 0;
	for(double q = qdiff; q <= qmax; q += qdiff){
	  cout << "qmax = " << qmax << "; qidx =" << qidx << "; q = " << q << endl;
	  while(qidx < num_samples-1 && q > samples[qidx].second)
	    qidx++;
	  
	  double b = samples[qidx].first;
	  if(lrb->bb.empty() || lrb->bb.back() != b){ // skip quantile if it falls onto the previous quantile value
	    lrb->bb.push_back(b);
	    cout << "featureGroup " << lrb->id << " with bin quantil " << q << ": " << lrb->bb.back() << endl;
	  }else
	    cout << "skiped quantil " << q << endl;
	}
	if(lrb->bb.size() != lrb->weights.size()){ // resize weights vector if neccessary (quantile has been skiped)
	  lrb->weights.resize(lrb->bb.size(), 0);
	  Constant::lr_features[i]->num_features = lrb->weights.size();
	  Constant::lr_features[i]->num_bins = lrb->weights.size()+1;
	}
      }
    }
  }
}


void TrainData::setFeatures(vector<double> *f, Traits *t){

  for(int i = 0; i < Constant::lr_features.size(); i++){
    LRfeatureGroup *lr_feature = Constant::lr_features[i];
    if(lr_feature->id == 0 || lr_feature->id == 20) // skip intercept
      continue;
    if( (lr_feature->id < 20 && isExon) || (lr_feature->id >= 20 && !isExon) ){
      bool hasTrait = true;
      double x = Properties::getFeature(lr_feature->id, t, &hasTrait);
      if(lr_feature->num_bins == 1){
	f->push_back(x);
	traitExists.back().push_back(hasTrait);
      }else{
	traitExists.back().insert(traitExists.back().end(), lr_feature->num_features, hasTrait);
	if(!hasTrait){ // if sample does not have the required traits insert zeros and process next feature group
	  f->insert(f->end(), lr_feature->num_features, 0);
	  continue;
	}
	LRbinnedFeatureGroup* lr_bf = static_cast<LRbinnedFeatureGroup*>(lr_feature);
	vector<double> bf = Properties::getBinnedFeature(lr_bf, x);
	f->insert(f->end(), bf.begin(), bf.end());
      }	
    }
  }
  //cout << std_values.size() << ";" << f->size() << ";" << traitExists.back().size() << endl;
  if(std_values.size() != 0 && (std_values.size() != f->size() || f->size() != traitExists.back().size()) )
    throw ProjectError("Error in TrainData::setFeatures: feature vector sizes differ");

  for(int i = 1; i < f->size(); i++){
    double fe = (*f)[i];
    //if(fe == 0)
    //traitExists.back()[i] = false;
    if(traitExists.back()[i] && to_be_std[i]){
      std_values[i].mean += fe;
      std_values[i].std += pow(fe, 2);
      std_values[i].numSamples++;
      if(std_values[i].max < fe)
	std_values[i].max = fe;
      if(std_values[i].min > fe)
	std_values[i].min = fe;
    }
  }
}



/*
 * get mean and std of sample features for and standardize
 */
void TrainData::standardizeFeatures(){

  for(int i=0; i<std_values.size(); i++){
    if(to_be_std[i]){
      std_values[i].mean /= std_values[i].numSamples;
      std_values[i].std = std_values[i].std / std_values[i].numSamples - pow(std_values[i].mean, 2);
      if(std_values[i].std < -0.001)
	throw ProjectError("Error in TrainData::standardizeFeatures(): attempt to build the squareroot of a negative number in calculation of standard deviation.");
      else if(std_values[i].std <= 0){
      	cerr << "Warning: In TrainData::standardizeFeatures(): standard deviation is zero. Cannot standardize this feature. Continue without standardization." << endl;
	std_values[i].mean = 0;
	std_values[i].std = 1;
      }
    }else{
      std_values[i].mean = 0;
      std_values[i].std = 1;
    }
  }

  for(int i=0; i<std_values.size(); i++)
    cout << "std_value " << i << ": mean = " << std_values[i].mean << "\tstd = " << std_values[i].std << "\tmax = " << std_values[i].max <<  "\tmin = " << std_values[i].min << "\tnumSamples = " << std_values[i].numSamples << endl;

  for(int i=0; i<samples.size(); i++){
    for(int j=0; j<samples[i].size(); j++){
      if(traitExists[i][j])
	samples[i][j] = (samples[i][j] - std_values[j].mean) / std_values[j].std;
    }
  } 
}

void TrainData::setNumberOfFeatures(){

  num_features = 0;
  to_be_std.clear();
  for(int i=0; i<Constant::lr_features.size(); i++){
    if( (isExon && Constant::lr_features[i]->id < 20) || (!isExon && Constant::lr_features[i]->id >= 20) ){
      num_features += Constant::lr_features[i]->num_features;
      if(Constant::lr_features[i]->if_std)
	to_be_std.insert(to_be_std.end(), Constant::lr_features[i]->num_features, true);
      else
	to_be_std.insert(to_be_std.end(), Constant::lr_features[i]->num_features, false);
    }
  }
}

void LogReg::backtransformWeigths(gsl_vector *theta, TrainData *data){

  int n = theta->size;
  double theta_sum = 0;
  for(int i=1; i<n; i++){
    theta_sum += gsl_vector_get(theta, i) * data->std_values[i].mean / data->std_values[i].std;
    gsl_vector_set(theta, i, gsl_vector_get(theta, i) / data->std_values[i].std);
  }
  gsl_vector_set(theta, 0, gsl_vector_get(theta,0) - theta_sum);  
}
  

void LogReg::updateWeights(const gsl_vector *theta, TrainData *data){

  int idx = 0;
  for(int i = 0; i < Constant::lr_features.size(); i++){
    if( (data->isExon && Constant::lr_features[i]->id < 20) || (!data->isExon && Constant::lr_features[i]->id >= 20) ){
      LRfeatureGroup *lr_feature = Constant::lr_features[i];
      if(lr_feature->num_bins == 1){
	lr_feature->weight = gsl_vector_get(theta, idx);
	idx++;
      }else{
	LRbinnedFeatureGroup *lr_bf = static_cast<LRbinnedFeatureGroup*>(lr_feature);
	for(int j = 0; j < lr_bf->weights.size(); j++){
	  lr_bf->weights[j] = gsl_vector_get(theta, idx);
	  idx++;
	}
      }
    }
  }
}
 
double LogReg::activation_f(const gsl_vector *theta, vector<double> *s){

  double theta_x = gsl_vector_get(theta,0);
  int n = theta->size;
  
  for(int i=1; i<n; i++){
    theta_x += (*s)[i] * gsl_vector_get(theta,i);
  }
  return theta_x;
}

/*
double LogReg::get_feature(vector<double> *f, int i){

  double x;
  switch(i){
  case 0: x = 1;                                          // intercept
    break;
  case 1: x = ( ((*f)[3] > 0) ? 0 : 1 );                  // not having omega
    break;
  case 2: x = ( 1 - (*f)[9] );                            // is not an OE 
    break;
  case 3: x = log((*f)[0]);                               // log exon length
    break;
  case 4: x = (*f)[1];                                    // posterior probability  
    break;
  case 5: x = (*f)[2];                                    // mean base probability 
    break;
  case 6: x = (*f)[3];                                    // omega
    break;
  case 7: x = (*f)[4];                                    // variance of omega
    break;
  case 8: x = (*f)[5];                                    // conservation
    break;
  case 9: x = (*f)[7];                                    // containment
    break;
  case 10: x = (*f)[6];                                   // diversity
    break;
  case 11: x = (*f)[8];                                   // number of exons involved in OE 
    break;
  case 12: x = ( ((*f)[1] <= 0) ? 1 : 0 );                // is not sampled 
    break;
  case 13: x = (*f)[5] * (*f)[6];                         // conservation * diversity
    break; 
  case 14: x = (*f)[3] * (*f)[6];                         // omega * diversity 
    break;
    //  case 15: x = (*f)[8]/num_species;                       // number of species involved in this ortho exon divided by clade size
    //break;
  case 15: x = ((*f)[10]>0 && (*f)[11]>0 && (*f)[12]>0 && (*f)[13]>0) ? ( (log((*f)[11]) - log((*f)[10])) * exp(Constant::lambda * (log((*f)[11]) - log((*f)[10]))) + (log((*f)[12]) - log((*f)[13])) * exp(Constant::lambda * (log((*f)[12]) - log((*f)[13]))) )/( exp(Constant::lambda * (log((*f)[11]) - log((*f)[10]))) + exp(Constant::lambda * (log((*f)[12]) - log((*f)[13]))) ) : 0;
    break;
    // case 16: x = ((*f)[14] > 0 && (*f)[14] < 1)? -log(1/(*f)[14] - 1) : ( ((*f)[14] == 0)? -50 : 50);
    //break;
  default: 
    if(Constant::binQuantiles.size() > 0 && i < 15 + Constant::binQuantiles.size()){
      return get_binning_feature(i-16, (*f)[14]);
    }else
      throw ProjectError("feature number " + itoa(i) + " does not exist!");
    break;
  }
  return x;
}


double LogReg::get_binning_feature(int i, double s){
  
  cout << "calling get_binning_feature with params i=" << i<< "; s=" << s << endl;

  int n = Constant::binQuantiles.size();
  if(n == 0){
    cout << "binQuantiles vector is empty" << endl;
    return 0;
  }
  double curr_bin = n-1;
  int j = n-1;
  while(s < Constant::binQuantiles[j]){
    j--;
  }
  curr_bin = j+1;
  cout << "return value: ";
  if( (i==0 && curr_bin==0) || (i==n-1 && curr_bin==n) ){
    cout << "1" << endl;
    return 1;
  }else if(curr_bin==i+1 && i >= 0 && i < n-1){
    cout << ( (Constant::binQuantiles[i+1] - s) / (Constant::binQuantiles[i+1] - Constant::binQuantiles[i]) ) << endl;
    return ( (Constant::binQuantiles[i+1] - s) / (Constant::binQuantiles[i+1] - Constant::binQuantiles[i]) );
  }else if(curr_bin==i && i > 0 && i <= n-1){
    cout << ( (s - Constant::binQuantiles[i-1]) / (Constant::binQuantiles[i] - Constant::binQuantiles[i-1]) ) << endl;
    return ( (s - Constant::binQuantiles[i-1]) / (Constant::binQuantiles[i] - Constant::binQuantiles[i-1]) );
  }else{
    cout << "0" << endl; 
    return 0;
  }
  
}

*/

double LogReg::cross_entropy_error_f(const gsl_vector *theta, void *samples){

  double cross_entropy_error = 0;
  double curr_CEE = 0;
  auto *samp  = (vector<vector<double> > *)samples;
  int n = theta->size;
  int m = samp->size();
  //  cout << "########### extreme candidates ################" << endl;
  for(int i=0; i<m; i++){
    int y = (*samp)[i][0];
    double theta_x = activation_f(theta, &(*samp)[i]);
    double exp_theta_x = exp(-theta_x);
    curr_CEE = -y * log(1 + exp_theta_x) + (1 - y) * (-theta_x - log(1 + exp_theta_x));
    cross_entropy_error += curr_CEE;

    //cout.precision(10); 
    //cout << "cross entropy error of sample " << i << " = " << curr_CEE << endl;
 
    // test output of "extreme" candidates: labels contradict prediction
    // double sigma = 1/(1+exp(-theta_x));
     //    if(abs(sigma - (*x)[i].first) > 0.8){
     //     cout << "feature_" << i << "\t" << (*x)[i].first << "\t" << theta_x << "\t" << sigma << "\t" << curr_CEE << endl;
     //}
  }
  // regularization
  double reg = 0;
  //cout << "theta in CEE calculation: (";
  //cout << gsl_vector_get(theta, 0) << ",";
  for(int i = 1; i < n; i++){
    reg += pow(gsl_vector_get(theta, i), 2);
    //cout << gsl_vector_get(theta, i) << ",";
  }
  //cout << ")" << endl;
  reg *= Constant::alpha;

  //cout << "##############################################" << endl;

  cross_entropy_error = -cross_entropy_error/m + reg/(2*m);
  //cout << "CEE+reg: " <<cross_entropy_error<< endl;
  return (cross_entropy_error);
}

// The gradient of f, df = (df/dx_1, df/dx_2, ..., df/dx_n). 
void LogReg::cross_entropy_error_df(const gsl_vector *theta, void *samples, gsl_vector *df){
  
  vector<vector<double> > *samp = (vector<vector<double> > *)samples;
  int m = samp->size();
  int n = theta->size;
 
  vector<double> grad(n,0);
  for(int i=0; i<m; i++){
    double sigma = 1/(1 + exp(-activation_f(theta, &(*samp)[i])));
    double s_y = sigma - (*samp)[i][0];
    if(s_y == 0){
      cout << "gradient vanishes in sample " << i << ". Try to solve numerical problem : ";
      sigma = ( 1 - exp(-activation_f(theta, &(*samp)[i])) ) / ( 1 - exp(-activation_f(theta, &(*samp)[i]) * 2) ); 
      s_y = sigma - (*samp)[i][0];
      if(s_y != 0)
	cout << "success" << endl;
      else
	cout << "failure" << endl;
    }
    grad[0] += s_y;
    for(int j=1; j<n; j++){
      grad[j] += (*samp)[i][j]  * s_y;
    }
  }
  //cout << "(theta_i|gradient_i) = (" << gsl_vector_get(theta, 0) << "|" << grad[0] << "\t";
  gsl_vector_set(df, 0, grad[0]/m);  
  for(int j=1; j<n; j++){ // with regularization
    gsl_vector_set(df, j, (grad[j] +  Constant::alpha * gsl_vector_get(theta, j)) / m );
    //cout << gsl_vector_get(theta, j) << "|" << grad[j] + 2 * Constant::alpha * gsl_vector_get(theta, j) << "\t";
  }
  //cout << ")" << endl; 
}

// Compute both f and df together. 
void LogReg::cross_entropy_error_fdf(const gsl_vector *params, void *features, double *f, gsl_vector *df){

  *f = cross_entropy_error_f(params, features);
  cross_entropy_error_df(params, features, df);
}


// label-noise robust version 

double LogReg::rob_cross_entropy_error_f(const gsl_vector *theta, void *samples){

  double cross_entropy_error = 0;
  double curr_CEE = 0;
  auto *samp  = (vector<vector<double> > *)samples;
  double epsilon = Constant::label_flip_prob; 
  double gamma = prob_true_false;
  int m = samp->size();
  for(int i=0; i<m; i++){
    double sigma = 1/(1 + exp(-1 * activation_f(theta, &(*samp)[i])));
    int y = (*samp)[i][0];
    curr_CEE = y * log( (1 - sigma)*epsilon*gamma + sigma*(1 - epsilon) ) + (1 - y) * log( (1 - sigma)*(1 - epsilon*gamma) + sigma*epsilon );
    //cout << "feature_" << i << "\t" << (*x)[i].first << "\t" << sigma << "\t" << curr_CEE << endl;
    cross_entropy_error += curr_CEE;
  }
  return (-cross_entropy_error);
}

void LogReg::rob_cross_entropy_error_df(const gsl_vector *theta, void *samples, gsl_vector *df){

  auto *samp = (vector<vector<double> > *)samples;
  int m = samp->size();
  int n = theta->size;
  double epsilon = Constant::label_flip_prob;
  double gamma = prob_true_false;
  vector<double> grad(n,0);
  for(int i=0; i<m; i++){
    double sigma = 1/(1 + exp(-1 * activation_f(theta, &(*samp)[i])));
    int y = (*samp)[i][0];
    double g = - ( y * (1 - epsilon - epsilon*gamma) / (epsilon*gamma + sigma * (1 - epsilon - epsilon*gamma)) + (1 - y)*(epsilon + epsilon*gamma - 1) / (1 - epsilon*gamma + sigma*(epsilon + epsilon*gamma - 1)) ) * sigma * (1 - sigma);
    grad[0] += g;
    for(int j=1; j<n; j++){
      grad[j] += (*samp)[i][j] * g;
      //cout << "gradient of exon 0 ... " << i << " for feature " << j << ": " << grad[j] << endl;
    }
   
  }
  for(int j=0; j<n; j++){
    gsl_vector_set(df, j, grad[j]);
  }
}

// Compute both f and df together. 
void LogReg::rob_cross_entropy_error_fdf(const gsl_vector *params, void *features, double *f, gsl_vector *df){

  *f = rob_cross_entropy_error_f(params, features);
  rob_cross_entropy_error_df(params, features, df);
}


void LogReg::optimize_parameters(TrainData *data){

  int n = data->num_features;

  if(data->samples.size() <= 1 || n == 0){ // no samples or no features for the sample type
    string type;
    if(data->isExon)
      type = "exon";
    else
      type = "intron";
    cerr << "Warning: Could not train " << type << " score: no " << type << " samples or relevant features found!" << endl;
    return;
  }

  gsl_vector *theta;
  theta = gsl_vector_alloc(n);
  gsl_vector_set_all(theta, 0);
  //  gsl_vector_set(theta, 0, -2.375175516);

  bool useTFprob = 0;
    
  if(Constant::rLogReg){
    cout << "# using lable-noise robust logistic regression with label flip probability " << Constant::label_flip_prob << endl;
    try{
      useTFprob = Properties::getBoolProperty("useTFprob");
    } catch(...){
      useTFprob = 0;
    }
    if(useTFprob)
      prob_true_false = data->prob_true_false_samples;
  }

  if(data->isExon)
    cout << "# exon parameter training" << endl;
  else
    cout << "# intron parameter training" << endl;

  // random allocation
  default_random_engine generator(time(NULL));
  uniform_real_distribution<double> distribution(-1,1);
  cout << "random theta: ";
  for(int i=0; i<n; i++){
    double r = distribution(generator);
    //double r = ( (double) rand() / RAND_MAX );
    gsl_vector_set(theta, i, r);
    cout << r << "\t";
  }
  cout << endl;


  bool use_sgd;
  try{
    use_sgd = Properties::getBoolProperty("use_sgd");
  } catch(...){
    use_sgd = 0;
  }
  
  if(std::isnan(cross_entropy_error_f(theta, &data->samples)))
    throw ProjectError("Error in optimizeParameters(): cannot calculate cross entropy error. NaN's produced.");
  
  if(use_sgd)
    sgd(&data->samples, theta, n);
  else
    gsl_minimizer(&data->samples, theta, n);

  cout << "cross entropy test error after training: " << cross_entropy_error_f(theta, &data->samples) << endl;
  
  backtransformWeigths(theta, data);
  updateWeights(theta, data);
  data->printParams();  

  gsl_vector_free(theta);

  
  
}

void LogReg::gsl_minimizer(vector<vector<double> > *samples, gsl_vector *theta, int n){

  int iter = 0;
  int status;

  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;


  try{
    Constant::rLogReg = Properties::getBoolProperty("rLogReg");
    Constant::label_flip_prob = Properties::getdoubleProperty("label_flip_prob");
  }catch(...){}

  gsl_multimin_function_fdf CE_error_f;
  CE_error_f.n = n;
  if(Constant::rLogReg){
    CE_error_f.f = &rob_cross_entropy_error_f;
    CE_error_f.df = &rob_cross_entropy_error_df;
    CE_error_f.fdf = &rob_cross_entropy_error_fdf;
  }else{
    CE_error_f.f = &cross_entropy_error_f;
    CE_error_f.df = &cross_entropy_error_df;
    CE_error_f.fdf = &cross_entropy_error_fdf;  
  }
  
  CE_error_f.params = (void *)samples;


  T = gsl_multimin_fdfminimizer_vector_bfgs2;
  //T = gsl_multimin_fdfminimizer_conjugate_fr;
  s = gsl_multimin_fdfminimizer_alloc (T, n);

  try{
    Constant::GD_stepsize = Properties::getdoubleProperty("GD_stepsize");
  }catch(...){}
  gsl_multimin_fdfminimizer_set (s, &CE_error_f, theta, Constant::GD_stepsize, 0.1);
  
  do{
    iter++;
    status = gsl_multimin_fdfminimizer_iterate (s);
    
    if (status){
      cout << "gsl error_code: " << status << endl; 
      break;
    }
    status = gsl_multimin_test_gradient (s->gradient, 1e-6);
   
    if (status == GSL_SUCCESS)
      cout << "Minimum found at:\n";
    
    /*============================================================
    double epsilon = 1e-4;
    gsl_vector *approx_g = gsl_vector_alloc(n);
    gsl_vector *calc_g = gsl_vector_alloc(n);

    for(int i=0; i<n; i++){
      gsl_vector_set(s->x, i, gsl_vector_get(s->x, i) + epsilon);
      double L1 = cross_entropy_error_f(s->x, samples);
      gsl_vector_set(s->x, i, gsl_vector_get(s->x, i) - 2 * epsilon);
      double L2 = cross_entropy_error_f(s->x, samples);
      gsl_vector_set(s->x, i, gsl_vector_get(s->x, i) + epsilon);
      
      double res = (L1 - L2) / (2*epsilon);
      gsl_vector_set(approx_g, i, res );
    }
    cross_entropy_error_df(s->x, samples, calc_g);

    for(int i=0; i<n; i++){
      cout << "approx grad theta_" << i << " = " << gsl_vector_get(approx_g, i) << "\t calculated grad = " << gsl_vector_get(calc_g, i) << endl; 
    }
    //============================================================*/

    cout << setw(10) << iter;
    for(int i=0; i<n; i++){   
      cout << setw(20) << gsl_vector_get(s->x, i) << "|" << gsl_vector_get(s->gradient, i);
    }
    cout << setw(10) << s->f << endl;
    

  }while (status == GSL_CONTINUE && iter < 1000);

  for(int i=0; i<n; i++)
    gsl_vector_set(theta, i, gsl_vector_get(s->x, i));


  gsl_multimin_fdfminimizer_free (s);  
  
}


void LogReg::sgd(vector<vector<double> > *samples, gsl_vector *theta, int n){
  
  
  double learning_rate;
  try{
    learning_rate = Properties::getdoubleProperty("learning_rate");
  }catch(...){
    learning_rate = 0.01;
  }
  double start_lr = learning_rate;

  int maxEpoch;
  try{
    maxEpoch = Properties::getIntProperty("max_sgd_epoch");
  }catch(...){
    maxEpoch = 100;
  }
  int iter = 0;
  int num_samples = samples->size();
  double curr_lowest_error = numeric_limits<double>::max();
  gsl_vector *best_theta = gsl_vector_alloc(n);

  cout << "num_samples: " << num_samples << endl;
  cout << "maxEpoch: " << maxEpoch << endl;
  
  while(iter/num_samples < maxEpoch){
    random_shuffle ( samples->begin(), samples->end() );
        for(int s = 0; s < num_samples; s++){
      double sigma = 1/(1 + exp(-activation_f(theta, &(*samples)[s])));
      double s_y = sigma - (*samples)[s][0];
      for(int i=0; i<n; i++){
	double grad;
	if(i==0)
	  grad = s_y;
	else
	  grad = s_y * (*samples)[s][i] + Constant::alpha/num_samples * gsl_vector_get(theta, i);
	double new_theta_i =  gsl_vector_get(theta, i) - learning_rate * grad;
	//cout << "theta[" << i << "] = " << gsl_vector_get(theta, i) << " has Grad " << grad << " and learning rate " << learning_rate << endl;
	gsl_vector_set(theta, i, new_theta_i);
      }
      iter++;
       /*
      if(iter%1000==0){
	double cee = cross_entropy_error_f(theta, samples);
	cout << setw(15) << "iter" << iter;
	for(int i=0; i<n; i++)
	  cout << setw(15) << gsl_vector_get(theta, i);
	cout << setw(15) << cee << setw(15) << learning_rate << endl;
      }
      */
    }
    learning_rate = start_lr * 1/(iter/num_samples);
    //learning_rate *= 1/(iter/num_samples); 

    double cee = cross_entropy_error_f(theta, samples);
    if(cee < curr_lowest_error){
      gsl_vector_memcpy(best_theta, theta);
      curr_lowest_error = cee;
    }
    gsl_vector *grad = gsl_vector_alloc(n);
    cross_entropy_error_df(theta, samples, grad);
    cout << setw(15) << iter/num_samples;
    for(int i=0; i<n; i++)
      cout << setw(20) << gsl_vector_get(theta, i) << "|" << gsl_vector_get(grad, i);
    cout << setw(15) << cee << endl; 
    gsl_vector_free(grad);
  }
  gsl_vector_memcpy(theta, best_theta);
  
  cout << setw(15) << iter/num_samples;
  for(int i=0; i<n; i++)
    cout << setw(15) << gsl_vector_get(theta, i);
  cout << setw(15) << curr_lowest_error << endl;

  gsl_vector_free(best_theta);
}

void LogReg::trainFeature_from_file(){

  ifstream trainfile;
  string buffer;
  string filename = Properties::getProperty("trainFeatureFile");
  trainfile.open(filename, ifstream::in);
  if(!trainfile){
    string errmsg = "Could not open the train feature file " + filename + ".";
    throw PropertiesError(errmsg);
  }
  while(getline(trainfile,buffer)){
    if(buffer[0] == '#')
      continue;
    stringstream linestream(buffer);
    string chr;
    string type;
    int startPos;
    int endPos;
    char strand;
    char frame;
    string attributes;
    linestream >> chr >> buffer >> type >> startPos >> endPos >> buffer >> strand >> frame >> attributes;

    if(type == "exon")
      type = "CDS";
    if(type == "intron" || type == "CDS"){
      stringstream str;
      str << type << "\t" << chr << "\t" << startPos<< "\t" << endPos << "\t" << strand << "\t" << frame;
      string key = str.str();
      
      auto got = Constant::logReg_samples.find(key);
      if ( got == Constant::logReg_samples.end() ){
	OEtraits t(1); 
	t.setLength(endPos - startPos + 1);   // exon length	
	pair<string, OEtraits> entry(key, t);
	got = Constant::logReg_samples.insert(entry).first;
      }

      stringstream s(attributes);
      string ap;
      vector<string> attr;
      while(std::getline(s, ap, ';')){
	attr.push_back(ap);
      }
      for(int i=0; i<attr.size(); i++){
	vector<string> attr_pair;
	stringstream ss(attr[i]);
	while(std::getline(ss, ap, '=')){
	  attr_pair.push_back(ap);
	}

	if (attr_pair[0] == "postProb")
	  got->second.setPostProb(stod(attr_pair[1]));
	else if (attr_pair[0] == "avgBaseProb" || attr_pair[0] == "mbp")
	  got->second.setMBP(stod(attr_pair[1]));
	else if (attr_pair[0] == "Eomega" || attr_pair[0] == "PMomega")
	  got->second.setOmega(stod(attr_pair[1]));
	else if (attr_pair[0] == "VarOmega" || attr_pair[0] == "varOmega")
	  got->second.setVarOmega(stod(attr_pair[1]));
	else if (attr_pair[0] == "cons")
          got->second.setConsScore(stod(attr_pair[1]));
	else if (attr_pair[0] == "div")
          got->second.setDiversity(stod(attr_pair[1]));
	else if (attr_pair[0] == "containment" || attr_pair[0] == "contain")
          got->second.setContainment(stoi(attr_pair[1]));
	else if (attr_pair[0] == "n" || attr_pair[0] == "numSpecies")
          got->second.setNumExons(stoi(attr_pair[1]));
	else if (attr_pair[0] == "leftBoundaryExtOmega")
	  got->second.setLeftExtOmega(stod(attr_pair[1]));
	else if (attr_pair[0] == "leftBoundaryIntOmega")
          got->second.setLeftIntOmega(stod(attr_pair[1]));
	else if (attr_pair[0] == "rightBoundaryIntOmega")
	  got->second.setRightIntOmega(stod(attr_pair[1]));
	else if (attr_pair[0] == "rightBoundaryExtOmega")
	  got->second.setRightExtOmega(stod(attr_pair[1]));
	else if (attr_pair[0] == "codingPostProb")
	  got->second.setCodingPostProb(stod(attr_pair[1]));
      }
      if(got->second.getNumExons() == 1)
	got->second.setIsOE(false);
      else
	got->second.setIsOE(true);
      //      got->second.print();
    }
  }
}


void LogReg::reference_from_file(unordered_map<string,int> *ref_class){

  ifstream reffile;
  string buffer;
  bool found_other_feature = 0;
  reffile.open(Constant::referenceFile, ifstream::in);
  if(!reffile){
    string errmsg = "Could not open the reference exon file " + Constant::referenceFile + ".";
    throw PropertiesError(errmsg);
  }
  while(getline(reffile,buffer)){
    if(buffer[0] == '#')  
      continue;
    stringstream linestream(buffer);
    string chr;
    string type;
    int startPos;
    int endPos;
    char strand;
    char frame;
    linestream >> chr >> buffer >> type >> startPos >> endPos >> buffer >> strand >> frame >> buffer;
    /*size_t pos = buffer.find("class=");
    int ref_class = buffer[pos+6] - '0';
    if(ref_class < 0){
      ref_class = buffer[pos+7] - '0';
      ref_class*= -1;
    }
    */
   
    if(type != "intron" && type != "CDS"){
      found_other_feature = 1;
    }else{
      stringstream str;
      str << type << "\t" << chr << "\t" << startPos<< "\t" << endPos << "\t" << strand << "\t" << frame;
      string key = str.str();
      //      cout << "ref_feature: " << key << endl;
      pair<string, int> k(key,1);
      ref_class->insert(k);
    }
  }
  if(found_other_feature)
    cerr << "Warning: reference feature(s) of type other than CDS or intron found! Make sure that relevant features have the correct type definition, e.g. intron or CDS. Note, that the program is case sensitive." << endl;

}

