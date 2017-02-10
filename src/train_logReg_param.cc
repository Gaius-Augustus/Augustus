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
\******************************************************************************/

// project includes

// standard C/C++ includes
#include <iostream>
#include <fstream>

#include "train_logReg_param.hh"


train_data::train_data(unordered_map<string, pair<int, vector<double> > > *s, unordered_map<string,int> *ref_class){
  
  exon_samples = new vector<pair<int, vector<double> > >;
  intron_samples = new vector<pair<int, vector<double> > >;
  while(!s->empty()){
    unordered_map<string, pair<int, vector<double> > >::iterator samp = s->begin();
    unordered_map<string, int>::iterator got = ref_class->find(samp->first);
    if(got == ref_class->end())
      samp->second.first = 0;
    else
      samp->second.first = 1;
    if(samp->first.at(0) == 'C')
      exon_samples->push_back(samp->second);
    else
      intron_samples->push_back(samp->second);
    s->erase(samp);
  }

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

  num_features = 14;

  /*

  cout << "exon samples:" << endl << "class\t" << endl;
  for(int i=0; i<10;i++)
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
  for(int i=0; i<10;i++)
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

double activation_f(const gsl_vector *theta, vector<double> *f){

  double theta_x = 0;
  int n = theta->size;
  for(int i=0; i<n; i++){
    theta_x += gsl_vector_get(theta, i) * get_feature(f, i);
  }
  return theta_x;
}

double get_feature(vector<double> *f, int i){

  switch(i){
  case 0: return 1;                                          // intercept
    break;
  case 1: return log((*f)[0]);                               // log exon length
    break;
  case 2: return (*f)[1];                                    // posterior probability  
    break;
  case 3: return (*f)[2];                                    // mean base probability 
    break;
  case 4: return (*f)[3]; // ( ((*f)[3] > 0) ? log((*f)[3]) : (*f)[3] ); // omega
    break;
  case 5: return (*f)[4];                                    // variance of omega
    break;
  case 6: return (*f)[5];                                    // conservation
    break;
  case 7: return (*f)[6];                                    // diversity
    break;
  case 8: return log((*f)[7]+1);                             // containment
    break;
  case 9: return (*f)[8];                                    // number of exons involved in OE 
    break;
  case 10: return (*f)[5] * (*f)[6];                         // conservation * diversity
    break; 
  case 11: return (*f)[3] * (*f)[6];  // ( ((*f)[3] > 0) ? log((*f)[3]) : (*f)[3] ) * (*f)[6];   // omega * diversity 
    break;
  case 12: return (1 - (*f)[9]);                             // is not an OE
    break;
  case 13: return ( ((*f)[1] <= 0) ? 1 : 0 );                // is not sampled
    break;
  default: throw ProjectError("feature number " + itoa(i) + " does not exist!");
    break;
  }
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
    //    cout << "cross entropy error of sample " << i << " = " << ( -1 * (*x)[i].first * log(1 + exp(-theta_x)) + (1-(*x)[i].first) * (-theta_x - log(1 + exp(-theta_x)))) << endl;
 
    // test output of "extreme" candidates: labels contradict prediction
    /* double sigma = 1/(1+exp(-theta_x));
    if(abs(sigma - (*x)[i].first) > 0.8){
      cout << "feature_" << i << "\t" << (*x)[i].first << "\t" << sigma << "\t" << curr_CEE << endl;
    }
    */
  }
  //cout << "##############################################" << endl;
  return (-1 * cross_entropy_error);
}

/* The gradient of f, df = (df/dx_1, df/dx_2, ..., df/dx_n). */
void cross_entropy_error_df(const gsl_vector *theta, void *features, gsl_vector *df){
  
  vector<pair<int, vector<double> > > *x = (vector<pair<int, vector<double> > > *)features;
  int m = x->size();
  int n = theta->size;

  vector<double> grad(n,0);
  for(int i=0; i<m; i++){
    double sigma = 1/(1 + exp(-1 * activation_f(theta, &(*x)[i].second)));
    for(int j=0; j<n; j++){
      grad[j] += get_feature(&(*x)[i].second, j) * (sigma - (*x)[i].first);
    }
  }
  for(int j=0; j<n; j++){
    gsl_vector_set(df, j, grad[j]);
  }
}

/* Compute both f and df together. */
void cross_entropy_error_fdf(const gsl_vector *params, void *features, double *f, gsl_vector *df){

  *f = cross_entropy_error_f(params, features);
  cross_entropy_error_df(params, features, df);
}

void optimize_parameters(train_data *data){

  // output file for trained log_reg parameter
  ofstream paramfile;
  if(Properties::hasProperty("param_outfile"))
    paramfile.open(Properties::getProperty("param_outfile"));
  else
    paramfile.open(Constant::configPath + "/cgp/log_reg_parameters_trained.cfg");
  
  int iter = 0;
  int status;
  int n = data->num_features;  
  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;

  gsl_vector *theta;
  theta = gsl_vector_alloc(n);
  gsl_vector_set_all(theta, 1);

  gsl_multimin_function_fdf CE_error_f;
  CE_error_f.n = n;
  CE_error_f.f = &cross_entropy_error_f;
  CE_error_f.df = &cross_entropy_error_df;
  CE_error_f.fdf = &cross_entropy_error_fdf;

  // exon parameter optimization
  CE_error_f.params = (void *)data->exon_samples;

  //cout << "cross entropy test error: " << cross_entropy_error_f(theta, data->exon_samples) << endl;

  T = gsl_multimin_fdfminimizer_conjugate_fr;
  s = gsl_multimin_fdfminimizer_alloc (T, n);

  try{
    Constant::GD_stepsize = Properties::getdoubleProperty("GD_stepsize");
  }catch(...){}
  //  cout << "# initial step size = " << Constant::GD_stepsize << endl;
  gsl_multimin_fdfminimizer_set (s, &CE_error_f, theta, Constant::GD_stepsize, 1e-4);
 
  do{
    iter++;
    status = gsl_multimin_fdfminimizer_iterate (s);
    
    if (status){
      //  cout << "gsl error_code: " << status << endl; 
      break;
    }
    
    status = gsl_multimin_test_gradient (s->gradient, 1e-3);
    
    /*
    if (status == GSL_SUCCESS)
      printf ("Minimum found at:\n");
    
    printf("%5d ", iter);
    for(int i=0; i<n; i++)
      printf("%.5f ", gsl_vector_get (s->x, i));
    printf("%10.5f\n", s->f);
    */

  }while (status == GSL_CONTINUE && iter < 10000);

  for(int i=0; i<n; i++)
    paramfile << "/CompPred/exon_score" << i << "\t" << gsl_vector_get (s->x, i) << endl;
  
  gsl_multimin_fdfminimizer_free (s);

  // intron parameter optimization
  CE_error_f.params = (void *)data->intron_samples;

  //cout << "cross entropy test error: " << cross_entropy_error_f(theta, data->intron_samples) << endl;

  T = gsl_multimin_fdfminimizer_conjugate_fr;
  s = gsl_multimin_fdfminimizer_alloc (T, n);

  gsl_multimin_fdfminimizer_set (s, &CE_error_f, theta, Constant::GD_stepsize, 1e-4);

  iter = 0;
  do{
    iter++;
    status = gsl_multimin_fdfminimizer_iterate (s);

    if (status){
      //      cout << "gsl error_code: " << status << endl;
      break;
    }

    status = gsl_multimin_test_gradient (s->gradient, 1e-3);

    /*
    if (status == GSL_SUCCESS)
      printf ("Minimum found at:\n");

    printf("%5d ", iter);
    for(int i=0; i<n; i++)
      printf("%.5f ", gsl_vector_get (s->x, i));
    printf("%10.5f\n", s->f);
    */

  }while (status == GSL_CONTINUE && iter < 10000);


  for(int i=0; i<4; i++)
    paramfile << "/CompPred/intron_score" << i << "\t" << gsl_vector_get (s->x, i) << endl;
  paramfile.close();

  gsl_multimin_fdfminimizer_free (s);  
  
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
    int frame;
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

