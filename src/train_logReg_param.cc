/*
 * train_logReg_param.cc
 *
 * License: Artistic License, see file LICENSE.TXT or 
 *          https://opensource.org/licenses/artistic-license-1.0
 */

// project includes

// standard C/C++ includes
#include <iostream>
#include <fstream>
#include <iomanip>
#include <regex>

#include "train_logReg_param.hh"

int num_species;
double prob_true_false;
vector<double> mean;
vector<double> se;


void train_OEscore_params(int numSpecies){

  unordered_map<string,int> ref_class;
  reference_from_file(&ref_class);
  train_data data(&Constant::logReg_feature, &ref_class, numSpecies);
  optimize_parameters(&data);

}


train_data::train_data(unordered_map<string, pair<int, vector<double> > > *s, unordered_map<string,int> *ref_class, int numSp){
  
  int numExonRef = 0;
  int numIntronRef = 0;
  num_species = numSp;
  exon_samples = new vector<pair<int, vector<double> > >;
  intron_samples = new vector<pair<int, vector<double> > >;
  while(!s->empty()){
    unordered_map<string, pair<int, vector<double> > >::iterator samp = s->begin();
    unordered_map<string, int>::iterator got = ref_class->find(samp->first);
    if(got == ref_class->end())
      samp->second.first = 0;
    else{
      samp->second.first = 1;
      if(samp->first.at(0) == 'C')
	numExonRef++;
      else
	numIntronRef++;
    }
    if(samp->first.at(0) == 'C')
      exon_samples->push_back(samp->second);
    else
      intron_samples->push_back(samp->second);
    s->erase(samp);
  }
  cout << "# " << numExonRef << " reference exon samples and " << numIntronRef << " reference intron samples in training set." << endl;
  cout << "# " << exon_samples->size() << " exon and " << intron_samples->size() << " intron candidates are used for training." << endl;

  int num_samples = exon_samples->size();
  for(int i=0; i<num_samples-1; i++){
    try{
      if( (*exon_samples)[i].second.size() != (*exon_samples)[i+1].second.size() ){
	throw ProjectError("Train data vectors for logistic regression have inconsistent length!");
      }
    }catch(...){
      throw ProjectError("Train data vectors for logistic regression have inconsistent length!");
    }
  }

  /*
   *  number of features
   */

  try{
    num_features = Properties::getIntProperty("/CompPred/num_features");
  }catch(...){
    num_features = 15;
  }
  mean.resize(num_features);
  se.resize(num_features);

  // feaures that are included in standardization
  exon_feature_idx   = {0,0,0,1,1,1,1,1,1,1,1,1,0,1,1,1};
  intron_feature_idx = {0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0};

  prob_true_false_exons = (double) numExonRef / exon_samples->size();
  prob_true_false_introns = (double) numIntronRef / intron_samples->size();
  prob_true_false = 1;

  /*      
    cout << "exon samples:" << endl << "class\t";
  for(int i=0; i<num_features; i++)
    cout << i << "\t";
  cout << endl;
  for(int i=0; i<exon_samples->size(); i++){
    cout << (*exon_samples)[i].first << "\t";
    for(int j=0; j<(*exon_samples)[i].second.size(); j++){
      cout << (*exon_samples)[i].second[j] << "\t";
    }
    cout << endl;
  }

  cout << "intron samples:" << endl << "class\t" << endl;
  for(int i=0; i<num_features; i++)
    cout << i << "\t";
  cout << endl;
  for(int i=0; i<intron_samples->size(); i++){
    cout << (*intron_samples)[i].first << "\t";
    for(int j=0; j<(*intron_samples)[i].second.size(); j++){
      cout << (*intron_samples)[i].second[j] << "\t";
    }
    cout << endl;
  }
  
  */

}


/*
 * get mean and std of sample features for standardization
 */
void train_data::set_mean_std(vector<pair<int, vector<double> > > *samples, vector<int> feature_idx){

  fill(mean.begin(), mean.end(), 0);
  fill(se.begin(), se.end(), 0);

  int m = samples->size();
  for(int i=0; i<m; i++){
    for(int j=0; j<num_features; j++){
      if(feature_idx[j]){
	double f = get_feature(&(*samples)[i].second, j);
	mean[j] += f;
	se[j] += pow(f, 2);
      }else{
	se[j] = 1;
      }
    }
  }
  for(int j=0; j<num_features; j++){
    if(feature_idx[j]){
      mean[j] /= m;
      se[j] = sqrt( se[j]/m - pow(mean[j], 2) );
    }
  }
  
  cout << "mean and std of features" << endl;
  for(int i=0; i<num_features; i++){
    cout << i << "\t" << mean[i] << "\t" << se[i] << endl;
  } 
   
}

double activation_f(const gsl_vector *theta, vector<double> *f){

  double theta_x = gsl_vector_get(theta, 0);
  int n = theta->size;

  for(int i=1; i<n; i++){
    theta_x += gsl_vector_get(theta, i) * ( (get_feature(f, i) - mean[i]) / se[i] ) ;
    //cout << "f" << i << ":" << setw(10) << gsl_vector_get(theta, i) << ":" << get_feature(f, i) << "(" << ( (get_feature(f, i) - mean[i]) / se[i] ) << ")\t"; 
  }
  //cout << endl;
  return theta_x;
}

double get_feature(vector<double> *f, int i){

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
  default: throw ProjectError("feature number " + itoa(i) + " does not exist!");
    break;
  }
  return x;
}


double cross_entropy_error_f(const gsl_vector *theta, void *features){

  double cross_entropy_error = 0;
  double curr_CEE = 0;
  vector<pair<int, vector<double> > > *x  = (vector<pair<int, vector<double> > > *)features;

  int m = x->size();
  //  cout << "########### extreme candidates ################" << endl;
  for(int i=0; i<m; i++){
    double theta_x = activation_f(theta, &(*x)[i].second);
    curr_CEE = -1 * (*x)[i].first * log(1 + exp(-theta_x)) + (1-(*x)[i].first) * (-theta_x - log(1 + exp(-theta_x)));
    cross_entropy_error += curr_CEE;
    //cout.precision(10); 
    //   cout << "cross entropy error of sample " << i << " = " << ( -1 * (*x)[i].first * log(1 + exp(-theta_x)) + (1-(*x)[i].first) * (-theta_x - log(1 + exp(-theta_x)))) << endl;
 
    // test output of "extreme" candidates: labels contradict prediction
    // double sigma = 1/(1+exp(-theta_x));
     //    if(abs(sigma - (*x)[i].first) > 0.8){
     //     cout << "feature_" << i << "\t" << (*x)[i].first << "\t" << theta_x << "\t" << sigma << "\t" << curr_CEE << endl;
     //}
  }
  //cout << "##############################################" << endl;
  //  cout << "CEE: " << (-1 * cross_entropy_error) << endl;  
  return (-1 * cross_entropy_error);
}

// The gradient of f, df = (df/dx_1, df/dx_2, ..., df/dx_n). 
void cross_entropy_error_df(const gsl_vector *theta, void *features, gsl_vector *df){
  
  vector<pair<int, vector<double> > > *x = (vector<pair<int, vector<double> > > *)features;
  int m = x->size();
  int n = theta->size;

  vector<double> grad(n,0);
  for(int i=0; i<m; i++){
    double sigma = 1/(1 + exp(-1 * activation_f(theta, &(*x)[i].second)));
    for(int j=0; j<n; j++){
      grad[j] += ( (get_feature(&(*x)[i].second, j) - mean[j]) / se[j] ) * (sigma - (*x)[i].first);
      //cout << "gradient of exon 0 ... " << i << " for feature " << j << ": " << grad[j] << endl;
    }
  }
  for(int j=0; j<n; j++){
    gsl_vector_set(df, j, grad[j]);
  }
}

// Compute both f and df together. 
void cross_entropy_error_fdf(const gsl_vector *params, void *features, double *f, gsl_vector *df){

  *f = cross_entropy_error_f(params, features);
  cross_entropy_error_df(params, features, df);
}


// label-noise robust version 

double rob_cross_entropy_error_f(const gsl_vector *theta, void *features){

  double cross_entropy_error = 0;
  double curr_CEE = 0;
  vector<pair<int, vector<double> > > *x  = (vector<pair<int, vector<double> > > *)features;
  double epsilon = Constant::label_flip_prob; 
  double gamma = prob_true_false;
  int m = x->size();
  for(int i=0; i<m; i++){
    double sigma = 1/(1 + exp(-1 * activation_f(theta, &(*x)[i].second)));
    curr_CEE = (*x)[i].first * log( (1 - sigma)*epsilon*gamma + sigma*(1 - epsilon) ) + (1-(*x)[i].first) * log( (1 - sigma)*(1 - epsilon*gamma) + sigma*epsilon );
    //cout << "feature_" << i << "\t" << (*x)[i].first << "\t" << sigma << "\t" << curr_CEE << endl;
    cross_entropy_error += curr_CEE;
  }
  return (-1 * cross_entropy_error);
}

void rob_cross_entropy_error_df(const gsl_vector *theta, void *features, gsl_vector *df){

  vector<pair<int, vector<double> > > *x = (vector<pair<int, vector<double> > > *)features;
  int m = x->size();
  int n = theta->size;
  double epsilon = Constant::label_flip_prob;
  double gamma = prob_true_false;
  vector<double> grad(n,0);
  for(int i=0; i<m; i++){
    double sigma = 1/(1 + exp(-1 * activation_f(theta, &(*x)[i].second)));
    double g = -1 * ( (*x)[i].first * (1 - epsilon - epsilon*gamma) / (epsilon*gamma + sigma * (1 - epsilon - epsilon*gamma)) + (1 - (*x)[i].first)*(epsilon + epsilon*gamma - 1) / (1 - epsilon*gamma + sigma*(epsilon + epsilon*gamma - 1)) ) * sigma * (1 - sigma);
    for(int j=0; j<n; j++){
      grad[j] += ( (get_feature(&(*x)[i].second, j) - mean[j]) / se[j] ) * g;
      //cout << "gradient of exon 0 ... " << i << " for feature " << j << ": " << grad[j] << endl;
    }
   
  }
  for(int j=0; j<n; j++){
    gsl_vector_set(df, j, grad[j]);
  }
}

// Compute both f and df together. 
void rob_cross_entropy_error_fdf(const gsl_vector *params, void *features, double *f, gsl_vector *df){

  *f = rob_cross_entropy_error_f(params, features);
  rob_cross_entropy_error_df(params, features, df);
}


void optimize_parameters(train_data *data){

  // output file for trained log_reg parameter
  ofstream paramfile;
  if(Properties::hasProperty("param_outfile"))
    paramfile.open(Properties::getProperty("param_outfile"));
  else
    paramfile.open(Constant::configPath + "/cgp/log_reg_parameters_trained.cfg");

  int n = data->num_features;

  gsl_vector *theta;
  theta = gsl_vector_alloc(n);
  gsl_vector_set_all(theta, 0);
  
  bool useTFprob = 0;
    
  if(Constant::rLogReg){
    cout << "# using lable-noise robust logistic regression with label flip probability " << Constant::label_flip_prob << endl;
    try{
      useTFprob = Properties::getBoolProperty("useTFprob");
    } catch(...){
      useTFprob = 0;
    }
    if(useTFprob)
      prob_true_false = data->prob_true_false_exons;
  }

  cout << "# exon parameter training" << endl;
 
  // random allocation
  srand(time(NULL));
  cout << "random theta: ";
  for(int i =0; i<n; i++){
    double r = ( (double) rand() / RAND_MAX );
    gsl_vector_set(theta, i, r);
    cout << r << "\t";
  }
  cout << endl;

  data->set_mean_std(data->exon_samples, data->exon_feature_idx);

  bool use_sgd;
  try{
    use_sgd = Properties::getBoolProperty("use_sgd");
  } catch(...){
    use_sgd = 0;
  }
  
  if(use_sgd)
    sgd(data->exon_samples, theta, n);
  else
    gsl_minimizer(data->exon_samples, theta, n);

  // output back transformed coefficients
  double theta_sum = 0;
  for(int i=1; i<n; i++){
    theta_sum += gsl_vector_get(theta, i) * mean[i] / se[i];
  }
  paramfile << "/CompPred/exon_score0\t" << ( gsl_vector_get(theta, 0) - theta_sum ) << endl; 
  for(int i=1; i<n; i++)
    paramfile << "/CompPred/exon_score" << i << "\t" << ( gsl_vector_get(theta, i) / se[i] ) << endl;
  

  cout << "# intron parameter training" << endl;

  // random allocation                                                                                                                
  srand(time(NULL));
  cout << "random theta: ";
  for(int i =0; i<n; i++){
    double r = ( (double) rand() / RAND_MAX );
    gsl_vector_set(theta, i, r);
    cout << r << "\t";
  }
  cout << endl;

  data->set_mean_std(data->intron_samples, data->intron_feature_idx);

  if(Constant::rLogReg && useTFprob)
    prob_true_false = data->prob_true_false_introns;

  if(use_sgd)
    sgd(data->intron_samples, theta, n);
  else
    gsl_minimizer(data->intron_samples, theta, n);
  
  // output back transformed coefficients
  theta_sum = 
      (gsl_vector_get(theta, 4) * mean[4] / se[4])
    + (gsl_vector_get(theta, 5) * mean[5] / se[5]) 
    + (gsl_vector_get(theta, 3) * mean[3] / se[3]);
  

  paramfile << "/CompPred/intron_score" << 0 << "\t" << ( gsl_vector_get (theta, 0) - theta_sum ) << endl;
  paramfile << "/CompPred/intron_score" << 1 << "\t" << ( gsl_vector_get (theta, 4) / se[4] ) << endl;
  paramfile << "/CompPred/intron_score" << 2 << "\t" << ( gsl_vector_get (theta, 5) / se[5] ) << endl;
  paramfile << "/CompPred/intron_score" << 3 << "\t" << ( gsl_vector_get (theta, 3) / se[3] ) << endl;
  
  paramfile.close();

  
}

void gsl_minimizer(vector<pair<int, vector<double> > > *samples, gsl_vector *theta, int n){

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

  cout << "cross entropy test error: " << cross_entropy_error_f(theta, samples) << endl;

  T = gsl_multimin_fdfminimizer_vector_bfgs2;
  s = gsl_multimin_fdfminimizer_alloc (T, n);

  try{
    Constant::GD_stepsize = Properties::getdoubleProperty("GD_stepsize");
  }catch(...){}
  gsl_multimin_fdfminimizer_set (s, &CE_error_f, theta, Constant::GD_stepsize, 1e-4);
  
  do{
    iter++;
    status = gsl_multimin_fdfminimizer_iterate (s);
    
    if (status){
      cout << "gsl error_code: " << status << endl; 
      break;
    }
    status = gsl_multimin_test_gradient (s->gradient, 1e-3);
   
    if (status == GSL_SUCCESS)
      cout << "Minimum found at:\n";
    
    cout << setw(10) << iter;
    for(int i=0; i<n; i++)
      cout << setw(10) << gsl_vector_get(s->x, i);
    cout << setw(10) << s->f << endl;
    

  }while (status == GSL_CONTINUE && iter < 10000);

  for(int i=0; i<n; i++)
    gsl_vector_set(theta, i, gsl_vector_get(s->x, i));

  gsl_multimin_fdfminimizer_free (s);  
  
}


void sgd(vector<pair<int, vector<double> > > *samples, gsl_vector *theta, int n){

  random_shuffle ( samples->begin(), samples->end() );

  double learning_rate;
  try{
    learning_rate = Properties::getdoubleProperty("learning_rate");
  }catch(...){
    learning_rate = 0.001;
  }
  int maxInter;
  try{
    maxInter = Properties::getIntProperty("max_sgd_inter");
  }catch(...){
    maxInter = 10;
  }
  int iter = 0;
  int num_samples = samples->size();
  double curr_lowest_error = numeric_limits<double>::max();
  gsl_vector *best_theta = gsl_vector_alloc(n);

  if(iter < maxInter*num_samples){
    for(vector<pair<int, vector<double> > >::iterator it = samples->begin(); it != samples->end(); it++){
      //cout << "++++++++++ exon with class " <<  it->first << endl;
      for(int i=0; i<n; i++){
	double sigma = 1/(1 + exp(-1 * activation_f(theta, &it->second)));
	//cout << "x_" << i << " = " << ( (get_feature(&it->second, i) - mean[i]) / se[i] ) << "; pred = " << sigma << endl;
	double grad = ( it->first - sigma ) * ( (get_feature(&it->second, i) - mean[i]) / se[i] );
	//      cout << "gradient = " << grad << endl;
	double new_theta_i =  gsl_vector_get(theta, i) + learning_rate * grad;
	gsl_vector_set(theta, i, new_theta_i);
	//      cout << "theta[" << i << "] = " << gsl_vector_get(theta, i) << endl;
      }
      
      if(iter % 1000 == 0){
	double cee = cross_entropy_error_f(theta, samples);
	if(cee < curr_lowest_error){
	  gsl_vector_memcpy(best_theta, theta);
	  curr_lowest_error = cee;
	}
	cout << setw(10) << iter;
	for(int i=0; i<n; i++)
	  cout << setw(10) << gsl_vector_get(theta, i);
	cout << setw(10) << cee << endl;
      }
      iter++;
    }
  }
  gsl_vector_memcpy(theta, best_theta);
  gsl_vector_free(best_theta);
}

void trainFeature_from_file(){

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
            
            unordered_map<string, pair<int, vector<double> > >::iterator got = Constant::logReg_feature.find(key);
            if ( got == Constant::logReg_feature.end() ){
                vector<double> feature(14,0);
                feature[0] = endPos - startPos + 1;   // exon length	
                pair<int, vector<double> > p;
                p = make_pair(-1, feature);
                pair<string, pair<int, vector<double> > > entry;
                entry = make_pair(key, p);
                got = Constant::logReg_feature.insert(entry).first;
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

                if (attr_pair[0] == "postProb"){
                    try {
                        got->second.second[1] = stod(attr_pair[1]);
                    }
                    catch (...) {
                    }
                }     
                else if (attr_pair[0] == "avgBaseProb" || attr_pair[0] == "mbp"){
                    try {
                        got->second.second[2] = stod(attr_pair[1]);
                    }
                    catch (...) {
                    }
                }
                else if (attr_pair[0] == "Eomega" || attr_pair[0] == "PMomega"){
                    try {
                        got->second.second[3] = stod(attr_pair[1]);
                    }
                    catch (...) {
                    }
                }
                else if (attr_pair[0] == "VarOmega" || attr_pair[0] == "varOmega"){
                    try {
                        got->second.second[4] = stod(attr_pair[1]);
                    }
                    catch (...) {                        
                    }
                }
                else if (attr_pair[0] == "cons"){
                    try {
                        got->second.second[5] = stod(attr_pair[1]);
                    }
                    catch (...) {
                    }
                }
                else if (attr_pair[0] == "div"){
                    try {
                        got->second.second[6] = stod(attr_pair[1]);
                    }
                    catch (...) {
                    }
                }
                else if (attr_pair[0] == "containment" || attr_pair[0] == "contain"){
                    try {
                        got->second.second[7] = stod(attr_pair[1]);
                    }
                    catch (...) {
                    }
                }
                else if (attr_pair[0] == "n" || attr_pair[0] == "numSpecies")
                    got->second.second[8] = stod(attr_pair[1]);
                else if (attr_pair[0] == "oescore")
                    got->second.second[9] = 1;
                else if (attr_pair[0] == "leftBoundaryExtOmega"){
                    try {
                        got->second.second[10] = stod(attr_pair[1]);
                    }
                    catch (...) {
                    }
                }
                else if (attr_pair[0] == "leftBoundaryIntOmega"){
                    try {
                        got->second.second[11] = stod(attr_pair[1]);
                    }
                    catch (...) {
                    }
                }
                else if (attr_pair[0] == "rightBoundaryIntOmega"){
                    try {
                        got->second.second[12] = stod(attr_pair[1]);
                    }
                    catch (...) {
                    }
                }
                else if (attr_pair[0] == "rightBoundaryExtOmega"){
                    try {
                        got->second.second[13] = stod(attr_pair[1]);
                    }
                    catch (...) {
                    }
                }
            }
        }        
    }
}


void reference_from_file(unordered_map<string,int> *ref_class){

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

