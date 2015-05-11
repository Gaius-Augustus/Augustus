#include<iostream>
#include<fstream>
#include<math.h>
#include<getopt.h>
#include<vector>
#include<map>
#include<stdlib.h>
#include "exon_seg.hh"

using namespace std;

  exon_segmentation::exon_segmentation(){
    /* Initialise lambda1 and lambda2
     */
    lambda1=0.2;
    lambda2=0.3;
  }
  /* function to calculate the geometric distribution parameters (lambdas)
   * which will be used later to calculate the distributions. 
   */
  void exon_segmentation::calculate_lambda(dataset &coverage_info, int c,float s){
    cout<<"Estimating the parameters of the geometric distribution."<<endl;
    /*int i,strand,track,l,sum=0;
    double prev_lambda1,prev_lambda2,sum1=0.0,sum2=0.0,sum11=0.0,sum22=0.0;
    double log_lambda1,log_lambda2,lob_comp_lambda1,log_comp_lambda2,log_r;
    vector<long double> p1,p2;
    vector<double> count;*/
    ///* Make the count vector
    //*/
    /*count.resize(maxcov+1,0);
    cout<<"Started making count vectors"<<endl;
    for(coverage_info.it=coverage_info.database.begin();coverage_info.it!=coverage_info.database.end();coverage_info.it++)
      for(strand=0;strand<coverage_info.it->second.size();strand++)// loop over strands +,-,b
	for(track=0;track<coverage_info.it->second[strand].size();track++)// loop over tracks
	  for(l=0;l<coverage_info.it->second[strand][track].size();l++)// loop over position on track
	    if(coverage_info.it->second[strand][track][l]<maxcov+1){
	      count[coverage_info.it->second[strand][track][l]]++;
	      sum++;
	      }*/
    /* Normalize the count vector to a distribution
     */
    //for(i=0;i<count.size();i++)
    //count[i]=count[i]/sum;
    /* Find the threshold value between introns and exons using the 
     * information about the % of exons in the sequence
     */
    /*sum1=0.0;
    c=-1;
    while(sum1<(1-r)){
      c++;
      sum1+=count[c];      
    }
    cout<<" The threshold is :"<<c<<'\t'<<count[c]<<'\t'<<r<<endl;
    for(i=0;i<5;i++)
    cout<<count[i]<<endl;*/
    //exit(1);

    /* Here we calculate the values of lambda1 and lambda2 by solving the two equations
     * (1-lambda1)/(1-lambda2)= s  (predefined value by the user ) and
     * lambda1(1-lambda1)^c=lambda2(1-lambda2)^c
     */
    lambda1=(s-1)/(pow(s,c+1)-1);
    lambda2=(pow(s,c+1)-pow(s,c))/(pow(s,c+1)-1);



    /* Initialise the values of lambda1 and lambda2 in the order of
     * 1/(mean coverage depth)
     */
    /*lambda2=(1/coverage_info.avg_depth)*1.1;
    lambda1=(1/coverage_info.avg_depth)*0.9;

     Initialise the p1 and p2 vectors
     */
    //p1.resize(maxcov+1,0.0);
    //p2.resize(maxcov+1,0.0);
    /* Choose the initial values of prev_lambda away from the values of
     * lambdas so that it does not satisfy the convergence condition
     * in the first iteration itself.
     */
    /*prev_lambda1=lambda1+2*convergence_limit;
      prev_lambda2=lambda2-2*convergence_limit;
    cout<<"The initial guesses  :"<<lambda1<<'\t'<<lambda2<<endl;*/
    /* Start the iteration untill the values converge
     */
    /*while(fabs(prev_lambda1-lambda1)>convergence_limit || fabs(prev_lambda2-lambda2)>convergence_limit){
      sum1=0.0;
      sum2=0.0;
      for(i=0;i<maxcov+1;i++){*/
	/* It is of the form a/(a+b) )=1/(1+b/a)
	 * b=(1-r)*lambda2*(1-lambda2)^i
	 * a=r*lambda1*(1-lambda1)^i
	 */
	/*p1[i]=1/(1+((1-r)/r)*(lambda2/lambda1)*pow(((1-lambda2)/(1-lambda1)),i));
	p2[i]=1-p1[i];
	//cout<<i<<'\t'<<r<<'\t'<<lambda1<<'\t'<<lambda2<<endl;
	//if(i>-1 && i<15)
	//cout<<i<<'\t'<<p1[i]<<'\t'<<p2[i]<<endl;
	//}*/      
      /* Start re-estimating the values of lambda
       */
    /*sum1=0.0;
      sum2=0.0;
      sum11=0.0;
      sum22=0.0;
      for(i=0;i<maxcov+1;i++){
	sum1+=i*p1[i]*count[i];
	sum11+=p1[i]*count[i];
	sum2+=i*p2[i]*count[i];
	sum22+=p2[i]*count[i];*/
	//if(i>-1 && i<25)
	//cout<<i<<"iterating sums  :"<<'\t'<<sum1<<'\t'<<sum2<<endl;
    //}
      /* Store the current value to prev_lambda
       */
      /*prev_lambda1=lambda1;
	prev_lambda2=lambda2;*/
      /* Update the values of lambda
       */
      //cout<<"sum1  "<<sum1<<'\t'<<"sum2  "<<sum2<<endl;
      /*lambda1=1/(1+sum1/sum11);
      lambda2=1/(1+sum2/sum22);
      cout<<"The iteratedd values :"<<lambda1<<'\t'<<lambda2<<endl;
    }*/

  //cout<<"The estimated values of lambda's are\nlambda1 :"<<lambda1<<"\nlambda2 :"<<lambda2<<endl;
    cout<<"lambda1 "<<lambda1<<"  lambda2 "<<lambda2<<endl;
    //exit(1);
  }
  /* calculates the geometric distribution values using the lambda values and stores them
   */
  void exon_segmentation::prob_density(int max,int no_of_tracks){
    int i,size=maxcov,j;
    vector<double> temp;
    /* calculate the log(lambda) values and log(1-lambda) values
     * beforehand to use them in the next step
     */
    double log_lambda1=log(lambda1);
    double log_comp_lambda1=log(1-lambda1);
    double log_lambda2=log(1-lambda2);
    double log_comp_lambda2=log(1-lambda2);
    /* If the maximum value of caverage observed(max) is less than the maximum coverage depth 
     * we have used for calculations then we resize the distribution vectors to the value max+1
     * (We add one to include the value of distribution point 0 also)
     */
    if(max<maxcov)
      size=max;
    temp.resize(size+1,0.0);
    prob_dist_exon.resize(no_of_tracks,temp);
    prob_dist_intron.resize(no_of_tracks,temp);
    for(j=0;j<no_of_tracks;j++){
      for(i=0;i<size+1;i++){
	prob_dist_exon[j][i]=log_lambda1+i*log_comp_lambda1;
	prob_dist_intron[j][i]=log_lambda2+i*log_comp_lambda2;
      }
    }
  }
  /* this function calculates the probability distribution values using an input gff file for training
   */
  void exon_segmentation::train_function(dataset &coverage_info,string inputfile,string key){
    int find1,i,start=0,end=0,seg_ptr,track,strand,l,depth;
    double sum_exons=0,sum_introns=0,sum_q_exon_num=0,sum_q_exon_den=0,sum_q_intron_num=0,sum_q_intron_den=0;
    /* r stores the fraction of postitions with coverage depths > L
     * q is a parameter such that L+1/q= "mean of coverages that are >L "
     * log_comp_q stores the value of log(1-q) for introns and exons
     */
    double q_intron,q_exon,log_comp_q_intron,log_comp_q_exon,r_exon=0,r_intron=0;
    ifstream infile(inputfile.c_str());
    string fileline,feature,strand1,line;
    /* stores the number of positions with a parrticular coverage depth
     */
    vector<int> count_introns,count_exons;
    /* stores the normalized value of the count vector
     */
    vector<double> f_introns,f_exons,temp1;
    /*structure where the data inormation from the training gff file is stored
     */
    vector<fragment> segments;
    /* a temporary fragment is defined to store the data initially and then it is pushed to the vector
     */
    fragment temp;
    /* initialise the vectors
     */
    count_introns.resize(maxcov+1,0);
    count_exons.resize(maxcov+1,0);
    f_exons.resize(maxcov+1,0.0);
    f_introns.resize(maxcov+1,0.0);
    temp1.resize(maxcov+1);
    prob_dist_exon.resize(coverage_info.no_of_tracks,temp1);
    prob_dist_intron.resize(coverage_info.no_of_tracks,temp1);

    /* We read the input gff file and store the information about introns and exons in a vector os fragments
     */
    getline(infile,line);
    while(!infile.eof() ){
      fileline=line;
      /* initially ignore all the lines containing other information
       */
      if(fileline[0]=='b'||fileline[0]=='#'||fileline[0]=='t'){
	getline(infile,line);
	continue;
      }
      /* analyse the line containing the information on feature and its position in the sequence
       */
      else{      
	/* we assume that the entries are separated by tabs in the gff file
	 */
	find1=fileline.find('\t');      
	/* loops through the number of tabs possible in a line
	 */
	for(i=0;i<8;i++){
	  /* extract the relevant information from the corrsponding i
	   */
	  if(i==2)
	    feature=fileline.substr(0,find1);
	  if(i==3)
	    start=atoi(fileline.substr(0,find1).c_str());
	  if(i==4)
	    end=atoi(fileline.substr(0,find1).c_str());
	  if(i==6){
	    strand1=fileline[0];
	  }
	  /* after a particular entry is processed remove that from the string fileline
	   */
	  fileline=fileline.substr(find1+1);
	  find1=fileline.find('\t');	
	}
	/* store the state in the fragment according to the feature found in the file
	 */
	if( (feature=="exon" || feature=="CDS") && strand1=="+"){
	  temp.state=EXONP;
	}
	else if( (feature=="exon" || feature=="CDS") && strand1=="-"){
	  temp.state=EXONM;
	}
	else{
	  temp.state=INTRONI;
	}
	temp.start=start;
	temp.end=end;
	/* push this fragment to the vector of fragments
	 */
	segments.push_back(temp);
      }	   
      getline(infile,line);
    }
    /* Make the probability distribution vector for each track
     */
    for(track=0;track<coverage_info.no_of_tracks;track++){
      for(strand=0;strand<coverage_info.database[key].size();strand++){
	/* ignore if the strand does not contain any experiments
	 */
	if(coverage_info.database[key][strand].size()==0)
	  continue;
	for(l=0;l<coverage_info.database[key][strand][track].size();l++){
	  if(l+1>segments[seg_ptr].end)
	    seg_ptr++;
	  depth=coverage_info.database[key][strand][track][l];
	    /* check for compatiability while storing the counts and ignore coverage depths greater than
	     * the maximum predefined coverage depth allowed for consideration
	     */
	  if( (segments[seg_ptr].state==EXONP && depth <maxcov+1 && strand!=STRANDM)|| (segments[seg_ptr].state==EXONM && depth <maxcov+1 && strand!=STRANDP) ){
	    count_exons[depth]++;
	    sum_exons++;
	    if(depth>L)
	      r_exon++;
	  }
	  else
	    if(depth<maxcov+1){
	      count_introns[depth]++;
	      sum_introns++;
	      if(depth>L)
		r_intron++;
	    }
	} 
	seg_ptr=0;
      }
      /* normalize the r value for exons and introns
       */
      r_exon=r_exon/sum_exons;
      r_intron=r_intron/sum_introns;
      /* calculate the relative frequencies in the loop
       */
      for(i=0;i<maxcov+1;i++){
	f_exons[i]=count_exons[i]/sum_exons;
	f_introns[i]=count_introns[i]/sum_introns;
	/* for i> L we have to calculate the numerator as the sum over count[i]*i 
	 * and the denominator as the sum over count[i] from L to maximum coverage 
	 * depth allowed ffor both inrons and exons
	 * then we take the ratio of denomminator to numerator to calculate the
	 * corresponding q values
	 */
	if(i>L){
	  sum_q_exon_num+=count_exons[i]*i;
	  sum_q_intron_num+=count_introns[i]*i;
	  sum_q_exon_den+=count_exons[i];
	  sum_q_intron_den+=count_introns[i];
	}
      }
      q_exon=sum_q_exon_den/sum_q_exon_num;
      q_intron=sum_q_intron_den/sum_q_intron_num;
      /* initialise the probability distribution vectors
       */
      log_comp_q_intron=log(1-q_intron);
      log_comp_q_exon=log(1-q_exon);
      for(i=0;i<maxcov+1;i++){
	/* for depths <= L we directly use the relative frequency values
	 */
	if(i<L+1){
	  prob_dist_exon[track][i]=log(f_exons[i]);
	  prob_dist_intron[track][i]=log(f_introns[i]);
	}
	/* For depth of L+1 we calculate the log value once
	 */
	else if(i==L+1){
	  prob_dist_exon[track][i]=log(r_exon)+log(q_exon);
	  prob_dist_intron[track][i]=log(r_intron)+log(q_intron);
	}
	/* we use the recursive definition to calculate the log values for depths > L+1
	 */
	else{
	  prob_dist_exon[track][i]=prob_dist_exon[track][i-1]+log_comp_q_exon;
	  prob_dist_intron[track][i]=prob_dist_intron[track][i-1]+log_comp_q_intron;
	}
      }
      r_exon=0;
      r_intron=0;
      sum_exons=0;
      sum_introns=0;
      sum_q_exon_num=0;
      sum_q_intron_num=0;
      sum_q_exon_den=0;
      sum_q_intron_den=0;
    }
    /* part to store the distribution data into a file for analysis
     */
    /*ofstream outfile("distribution data.dat",ios::out);
    outfile<<"r_exon:"<<r_exon<<" r_intron:"<<r_intron<<" q_exon:"<<q_exon<<" q_intron:"<<q_intron<<endl;
    outfile<<"log of the probability distribution table"<<endl;
    outfile<<"exon"<<'\t'<<"intron"<<endl;
    for(i=0;i<maxcov+1;i++)
      outfile<<prob_dist_exon[i]<<'\t'<<prob_dist_intron[i]<<endl;
      outfile.close();*/
  }
  
  /* functions which return the distribution value for a particular coverage depth
   * by accessing the corresponding index in the distribution table
   */
  double exon_segmentation::prob_value_exon(int num,int track){
    if(num>=prob_dist_exon[track].size())
      return prob_dist_exon[track][prob_dist_exon.size()-1];
    else
      return prob_dist_exon[track][num];
  }
  double exon_segmentation::prob_value_intron(int num,int track){
    if(num>=prob_dist_intron[track].size())
      return prob_dist_intron[track][prob_dist_intron.size()-1];
    else
      return prob_dist_intron[track][num];
  }
  /* Find the splice sites
   * INCOMPLETE
   */
  void exon_segmentation::find_splice(vector< vector< vector<int> > > &input_set,dataset &coverage_info){
    cout<<"Started finding splice sites."<<endl;
    int i,j,strand,track,l,template_size=90;
    double mean_template=0.0,f_mean=0,sigma_template=0,sigma_f=0;
    vector<double> smooth_data,corr_value;
    /* It will store the template which we have to find
     */
    double pattern_to_find[]={0.242857,0.257143,0.285714,0.314286,0.342857,0.371429,0.414286,0.457143,0.5,0.542857,0.585714,0.628571,0.671429,0.714286,0.785714,0.857143,0.928571,1.07143,1.21429,1.35714,1.5,1.65714,1.84286,2.07143,2.31429,2.57143,2.9,3.27143,3.71429,4.18571,4.74286,5.28571,5.9,6.57143,7.32857,8.17143,9.08571,10.0143,10.9714,11.8714,12.8,13.7143,14.6429,15.6143,16.6,17.6143,18.6429,19.7,20.7286,21.8571,22.9429,24.0429,25.1571,26.2143,27.3286,28.4286,29.6,30.7571,31.8429,32.9429,34.0429,35.1429,36.2429,37.3857,38.4143,39.4429,40.4143,41.4,42.2714,43.1429,44.0143,44.9286,45.7429,46.5571,47.3714,48.1714,48.9429,49.7571,50.5571,51.3571,52,52.6429,53.2714,53.9143,54.5571,55.1857,55.6714,56.1429,56.5,56.8571,};    
    ofstream splice_place("splice_site_positions.txt",ios::out);    
    /* Testing the method on first track of the +ve strand only
     */
    strand=STRANDP;
    track=0;
    /* Initialize the vector to store the smooth data
     */
    smooth_data.resize(input_set[0][0].size(),0.0);
    /* For the begining the code is written just to wrok on the
     * first track of the + strand only
     */
    for(l=0;l<smooth_data.size();l++){// loop over position on track
      /* Dont change the values for the begining values and end values
       * which cannot be covered by the window
       */
      if(l<moving_window/2){
	smooth_data[l]=input_set[strand][track][l];
	//s_test<<l+1<<'\t'<<smooth_data[l]<<'\n';
      }
      else if(l>=input_set[strand][track].size()-moving_window/2-1){
	smooth_data[l]=input_set[strand][track][l];
	//s_test<<l+1<<'\t'<<smooth_data[l]<<'\n';
      }
      /* Calculate the moving average value for the first element in window
       */
      else if(l==moving_window/2){
	for(i=l-moving_window/2;i<l+moving_window/2+1;i++)
	  smooth_data[l]=smooth_data[l]+input_set[strand][track][i]/moving_window;
	//s_test<<l+1<<'\t'<<smooth_data[l]<<'\n';
      }
      /* Use the previous average value to calculate the average value for the rest of the indices
       */
      else{
	smooth_data[l]=(smooth_data[l-1]*moving_window-input_set[strand][track][l-moving_window/2-1]+input_set[strand][track][l+moving_window/2])/moving_window;
	//s_test<<l+1<<'\t'<<smooth_data[l]<<'\n';
      }
    }  

    /* Make the template pattern that we have to search for
     */
    mean_template=0;
    //pattern_to_find.resize(template_size,0.0);
    /* Selecting a pattern which has steep rise
     */
    for(j=0;j<template_size;j++){
      mean_template+=pattern_to_find[j];
    }
    mean_template=mean_template/template_size;  
    /* Find the standard deviation of the template
     */
    for(i=0;i<template_size;i++){
      sigma_template+=(pattern_to_find[i]-mean_template)*(pattern_to_find[i]-mean_template);
    }
    sigma_template=sqrt(sigma_template);
    /* Testing the template matching using Normalized Cross Correlation
     */
    corr_value.resize(smooth_data.size(),0.0);
    f_mean=0;
    for(l=0;l<smooth_data.size()-template_size;l++){
      
      /* Find the standard deviation and mean 
       */
      sigma_f=0;
      f_mean=0;
      for(i=l;i<l+template_size;i++)
	f_mean+=smooth_data[i];
      f_mean=f_mean/template_size;
      for(i=l;i<l+template_size;i++)
	sigma_f+=(smooth_data[i]-f_mean)*(smooth_data[i]-f_mean);
      sigma_f=sqrt(sigma_f);
      corr_value[l]=0;
      for(i=0;i<template_size;i++)
	corr_value[l]+=(smooth_data[l+i]-f_mean)*(pattern_to_find[i]-mean_template);
      corr_value[l]=corr_value[l]/(sigma_f*sigma_template*(template_size-1) );
      if(corr_value[l]>0.01123)
	splice_place<<l<<'\t'<<corr_value[l]<<endl;
    }
    splice_place.close();
    cout<<"Finished finding splice sites."<<endl;
  }

  /* reads the file and stores it in the the object of the class dataset
   */
  void exon_segmentation::read_file(dataset &coverage_info,string filename){
    vector<int> temp_vector;//temporary vectors to initialize higher dimension vectors
    vector< vector<int> > temp_vector1;
    size_t found,found1;
    string fileline,temp_string,chrom,sign;
    int i,j,pos=0,val,flag=2,max_depth=0,max_len=0,no_of_tracks=0;
    int sum=0,total_size=0,start=0,step=0,span=1,count=0,strand=0;//0 for plus, 1 for - and 2 for both
    
    ifstream infile (filename.c_str()); 
    getline(infile,fileline);  
    while(!infile.eof()){
      /* if the line begins with f or v then store the necessary information from there
       */
      if(fileline[0]=='f'||fileline[0]=='v'){
	count=0;//count stores the number of steps been taken in case of fixedstep type
	/* store the name of the chromosome from the file
	 */
	found=fileline.find("chrom");
	if(found==string::npos){
	  cout<<"Incorrect Wiggle Format"<<endl;
	  exit(0);
	}
	else{
	  temp_string=fileline.substr(found);
	  found=temp_string.find(' ');
	  chrom=temp_string.substr(6,found-6);	
	  if(coverage_info.database[chrom].size()==0){
	    temp_vector1.clear();
	    cout<<"Found new chromosome: "<< chrom<<endl;
	    /* alocate space for the strands for this chromosome
	     */
	    for(i=0;i<3;i++)
	      coverage_info.database[chrom].push_back(temp_vector1);
	  }   
	  /* otherwise directly make space in the corresponding strand 
	   */
	  coverage_info.database[chrom][strand].push_back(temp_vector);
	}
	/* store the value of span from the string fileline else store the default value
	 */
	found=fileline.find("span");
	if(found!=string::npos){
	  temp_string=fileline.substr(found);
	  found=temp_string.find(' ');
	  span=atoi(temp_string.substr(5,found-5).c_str());
	}
	else
	  span=1;
	/* store the value of step from the string fileline else store the default value
	 */
	found=fileline.find("step");
	if(found!=string::npos){
	  temp_string=fileline.substr(found);
	  found=temp_string.find(' ');
	  step=atoi(temp_string.substr(5,found-5).c_str());
	}
	else
	  step=0;
	/* store the value of start from the string fileline else store the default value
	 */
	found=fileline.find("start");
	if(found!=string::npos){
	  temp_string=fileline.substr(found);
	  //check for separation by both space and tab
	  found1=temp_string.find('\t');
	  found=temp_string.find(' ');
	  if(found>0 && found>found1)
	    found=found1;
	  start=atoi(temp_string.substr(6,found-6).c_str());
	}
	else
	  start=0;

	/* if span=0 then set the flag for the simplest case
	 */
	if(span==0)
	  flag =0;
	else if (step==0)
	  /* set the flag for the case of variable start with span
	   */
	  flag=1;
	else
	  /* set the flag for fixed step with span
	   */
	  flag=2;	
	cout<<"start="<<start<<'\t'<<"step="<<step<<'\t'<<"span="<<span<<'\t'<<"flag="<<flag<<endl;
      }
      /* ignore the lines that are not useful
       */
      else if(fileline[0]=='#'||fileline[0]=='b'){
	getline(infile,fileline);
	continue;
      }
      /* we use the track line and determine whether the strand is +(0),-(1) or both(2) and store it in strand
       */
      else if(fileline[0]=='t'){
	/* store the maximum genome length by comparing with the previous track length
	 */
	if(pos>max_len)
	  max_len=pos;      
	found=fileline.find('(');
	if(found!=string::npos){
	  if(fileline.substr(found+1,1)=="+")
	    strand=STRANDP;	 
	  else if( fileline.substr(found+1,1)=="-" )
	    strand=STRANDM;
	  else
	    strand=STRANDB;	  
	}
	/* assuming it to be unknown strand otherwise
	 */
	else
	  strand=STRANDB;
      }
      /* Now we use the line where score are written and store them in proper position
       */
      else{
	/* the integers can be separated by either tab or space so check for both
	 */
	found1=fileline.find('\t');
	found=fileline.find(' ');
	if(found>0 && found>found1)
	  found=found1;
	/* in case of other than fixedstep type store the position and value
	 */
	if(flag!=2){
	  pos=atoi(fileline.substr(0,found+1).c_str());
	  val=atoi(fileline.substr(found).c_str());
	}
	/* in case of fixedstep we only have the values to store
	 */
	else
	  val=atoi(fileline.c_str());   
	if(val>max_depth)
	  max_depth=val;
	/* store the corresponding vector index in i, to which values should be pushed into
	 */
	i=coverage_info.database[chrom][strand].size()-1;
	/* calculate the position if the track is of fixedstep type
	 */
	if(flag==2){
	  pos=start+count*step;
	  count++;
	}
	/* fill in 0 to the values not present in the data set
	 */
	for(j=coverage_info.database[chrom][strand][i].size();j<pos-1;j++)
	  coverage_info.database[chrom][strand][i].push_back(0);
	/* fil in the obtained score to all the positions specified
	 * We also add the coverage depth value to sum which will be 
	 * used to calculate the overall mean of the coverage depths
	 */
	for(j=pos;j<pos+span;j++){
	  coverage_info.database[chrom][strand][i].push_back(val);
	  sum+=val;
	}
      }
      getline(infile,fileline);
    }
    /* store the maximum genome length by comparing with the previous length for the last track
     */
    if(pos>max_len)
      max_len=pos;
    coverage_info.max_len=max_len;
    coverage_info.max_depth=max_depth;
    cout<<"maximum length is :"<<max_len<<endl;
    /* resize all the vectors to maximum size
     */
    for(coverage_info.it=coverage_info.database.begin();coverage_info.it!=coverage_info.database.end();++coverage_info.it)
      for(i=0;i<coverage_info.it->second.size();i++){//loop over all the strands
	if(no_of_tracks < coverage_info.it->second[i].size())
	  no_of_tracks=coverage_info.it->second[i].size();
	for(j=0;j<coverage_info.it->second[i].size();j++) {//loops over all the tracks	
	  /* Add the size of the track to total_track to calculate the average
	   */
	  total_size+=coverage_info.it->second[i][j].size();
	  /* Resize the vector
	   */
	  coverage_info.it->second[i][j].resize(max_len,0);
	}
      }
    coverage_info.avg_depth=sum/total_size;
    /* Store the number of tracks which we take to be
     * maximum over the number of tracks present in each strand
     */
    coverage_info.no_of_tracks=no_of_tracks;
    cout<<"NUmber of tracks "<<no_of_tracks<<endl;
  }


  /* it takes as input pointer to a 2-D vecctor, 
   * the input dataset for a chromosome and
   * stores the emission probability values in it
   */
void exon_segmentation::calculate_emissions( vector< vector<double> > &emission_probs,int max_len,vector< vector< vector<int> > > &input_set,int chunksize,int iteration){
  int i,j,m,index,strand,track,limit;
    /* g_matrix is used to store the probaility values temporarily so as to
     * reduce the time complexity in states having length > 1
     */
    vector<vector<double> > g_matrix;
    /*define temporary vectors to initialise higher dimension vectors
     */
    vector<double> temp,temp1;
    double sum,sum1;
    temp1.resize(chunksize,0.0);
    for(i=0;i<NUMSTATES;i++){
      if(lengths[i]>1)
	/*makes the g_matrix only for states having length > 1
	 */
	g_matrix.push_back(temp1);
      else
	g_matrix.push_back(temp);
    }
    /* Calculate the size of the current chunk and store it in limit so as to 
     * avoid index overflow in tha dataset
     */
    if((iteration+1)*chunksize >max_len)
      limit=max_len-iteration*chunksize;
    else
      limit=chunksize;
    for(i=0;i<NUMSTATES;i++){
      for(j=0;j<limit;j++){
	/* ignore the case with negative indices
	 */
	if((j-lengths[i]+1)<0)
	  continue;            
	sum=0;
	sum1=0;
	/* for the states with lengths 1 we directly use this formula without affecting the time complexity
	 * we also calculate this for the first instance for storing in the g_matrix for states with length >1
	 */
	if( lengths[i]==1 || j==lengths[i]-1){
	  /* loops for all possible lengths
	   */
	  for(m=1;m<=lengths[i];m++){
	    /* loops through all possible strands
	     */
	    for(strand=0;strand<3;strand++){
	      /* loops through all the tracks in each strand
	       */
	      for(track=0;track<input_set[strand].size();track++){
		/* store the index(depth) for which probability value is required 
		 * here we also add the offset determined by the chunksize and the number of iteration
		 */
		index=input_set[strand][track][j+1-m+iteration*chunksize];
		/* if the strand amd the label are compatiable then use the
		 * probabiliity distribution for exons
		 */
		if( (i==EXONP && (strand!=STRANDM))||(i==EXONM && (strand!=STRANDP)) ){
		  sum=sum+prob_value_exon(index,track);
		}
		else{
		  sum=sum+prob_value_intron(index,track);
		}
	      }
	    }
	    /* for states with length >1 store it in g_matrix
	     */
	    if(lengths[i]>1)
	      g_matrix[i][j+1-m]=sum;
	    sum1=sum1+sum;
	    sum=0;
	  }
	  emission_probs[i][j]=sum1;
	}
	/* for the states with length > 1 we use the previously calculated g values to find the emission probability
	 */
	else{
	  sum=0;//initialise the value of sum
	  /* first we calculate the next g value
	   */
	  for(strand=0;strand<3;strand++){
	    for(track=0;track<input_set[strand].size();track++){
	      index=input_set[strand][track][j+iteration*chunksize];
	      /* check compatiability 
	       */
	      if( (i==EXONP && (strand!=STRANDM))||(i==EXONM && (strand!=STRANDP)) )
		sum=sum+prob_value_exon(index,track);
	      else
		sum=sum+prob_value_intron(index,track);	    	    
	    }
	    g_matrix[i][j]=sum;//stores the value of g
	    /* calculates emission value by adding the head value and substracting 
	     * the tail value to the previous entry in g_matrix
	     */
	    emission_probs[i][j]=emission_probs[i][j-1]+g_matrix[i][j]-g_matrix[i][j-lengths[i]];
	  }     
	}
      } 
    }   
  }

  /* takes as input a row and column and calculates the value to be filled 
   * in that element in the gamma matrix
   * there is also a pointer spred to store the predecessor state in it
   */
  void exon_segmentation::fill_viterbi(vector< vector<double> > &gamma, int numstates,int j, int s,double emission_val,double a[][NUMSTATES],int* spred){
    int k;
    double sum=0,max=minus_infinity;
    /* for negative value of j-len(s) we make it as minus infinity
     */
    if( (j-lengths[s]) <0 )
      sum=minus_infinity;
    /* else we take the maximum through all possible states k
     */
    else{
      max=minus_infinity;//initialise max
      for(k=0;k<numstates;k++){
	/* ignore the transitions that are not possible
	 */
	if(a[k][s]==minus_infinity)
	  continue;
	sum=gamma[j-lengths[s] ][k] + a[k][s] + emission_val;
	if(sum>max){
	  max=sum;
	  *spred=k;//store the predecessor
	}
      }
    }
    gamma[j][s]=max;
  }
  
  /* this function implements the viterbi algorithm
   */
vector<int> exon_segmentation::viterbi(vector< vector< vector<int> > > &input_set,int max_score,int max_len,int chunksize,int iteration){

  int i,j,s,spred=0,limit;
    double sum1,max;
    /* stores the state sequence and the reverse state sequence
     */
    vector<int> state_seq,state_seq1;
    vector<double> temp;
    /* define the vectors to store the emission probabilities
     * and the dynamic programming table
     */
    vector< vector<double> > emission_probs,gamma;//vectors to store the emission probabilities and dynamic matrix
    /* define the transition probability matrix for the states
     */
    double a[NUMSTATES][NUMSTATES]={ {0,minus_infinity,-10,minus_infinity},
				  {minus_infinity,0,-10,minus_infinity},
				  {minus_infinity,minus_infinity,minus_infinity,0},
				  {-10,-10,minus_infinity,0}
    };
    /* Calculate the value of limit which stores the maximum size of the current chunk
     * to avoid index overflow
     */
    if((iteration+1)*chunksize >max_len)
      limit=max_len-iteration*chunksize;
    else
      limit=chunksize;
    
    /* initialise the dynamic programming matrix
     */
    temp.clear();//clears the temporary vectors
    for(i=0;i<NUMSTATES;i++){
      temp.push_back(0.0);
    }
    gamma.resize(limit+1,temp);
    
    /* initialise the emission probability matrix
     */
    temp.clear();
    temp.resize(limit,0.0);
    emission_probs.resize(NUMSTATES,temp);
    calculate_emissions(emission_probs, max_len,input_set,chunksize,iteration);
    /*for(i=0;i<emission_probs[0].size();i++){
      cout<<i<<'\t';
      for(j=0;j<emission_probs.size();j++)
	cout<<emission_probs[j][i]<<'\t';
      cout<<endl;
      }
    cout<<"tesint "<<prob_dist_exon[0][20]<<'\t'<<prob_dist_intron[0][20]<<endl;  
    cout<<"tesing "<<prob_dist_exon[0][0]<<'\t'<<prob_dist_intron[0][0]<<endl;*/
    
    /* calculate the dynamic programing matrix
     */
    state_seq.clear();
    for(j=1;j<limit+1;j++)
      for(s=0;s<NUMSTATES;s++){
	sum1=emission_probs[s][j-1];
	fill_viterbi(gamma,NUMSTATES,j,s,sum1,a,&spred);
      }

    /* bactracking to get the sequence
     */
    j=limit;
    max=minus_infinity;
    for(i=0;i<NUMSTATES;i++)
      if(gamma[j][i]>max){
	max=gamma[j][i];
	s=i;
      }
    while(j>0){ 
      state_seq.push_back(s);
      sum1=emission_probs[s][j-1];
      fill_viterbi(gamma,NUMSTATES,j,s,sum1,a,&spred); 
      j=j-lengths[s];
      s=spred;      
    }
    /* reverse the state sequence obtained after backtracking
     */
    state_seq1.clear();
    for(i=state_seq.size()-1;i>=0;i--)
      state_seq1.push_back(state_seq[i]);
    return state_seq1;      
  }



  /* function to convert the state sequence into segments and also to store the avg coverage depth
   */
  vector< fragment > exon_segmentation::segment( vector<int> state_seq,vector< vector< vector<int> > > &input_set,int include_intron ,int no_of_tracks){
    int start,i,end,j,k;
    double sum=0;
    fragment temp;
    vector< fragment > segmentation;
    vector<double> temp1;
    start=1;
    end=lengths[state_seq[0]];
    cout<<"Started segmentation and calculating average coverage depths for the fragments."<<endl;
    for(i=1;i<state_seq.size();i++){
      /* if the next state is same as the previous state or they both are introns then we continue
       */
      if( (state_seq[i]==state_seq[i-1]) || (state_seq[i]==3&&state_seq[i-1]==2)||(state_seq[i]==2&&state_seq[i-1]==3)){
	end=end+lengths[state_seq[i]];
	continue;
      }
      else{
	/* once we reach a different state we store the stretch of information
	 */
	/* Initialise the vector to store the average coverage depth
	 * for each experiment
	 */
	temp1.resize(no_of_tracks,0.0);
	temp.state=state_seq[i-1];
	temp.start=start;
	temp.end=end;
	if(temp.state==EXONP){
	  for(k=0;k<input_set[STRANDP].size();k++){//loops over the tracks in + strand for exon+ state
	    for(j=start;j<end+1;j++){//loops in the positions from start to end in the track
	      sum=sum+input_set[STRANDP][k][j-1];//sum over to get the total average
	      temp1[k]+=input_set[STRANDP][k][j-1];//sum to get average for particular track
	    }
	    
	  }
	  for(k=0;k<input_set[STRANDB].size();k++){//loops over the tracks in . strand for exon+ state
	    for(j=start;j<end+1;j++){//loops in the positions from start to end in the track
	      sum=sum+input_set[STRANDB][k][j-1];
	      temp1[k]+=input_set[STRANDB][k][j-1];//sum to get average for particular track
	    }	    
	  }
	  /* calculate the average coverage depth for this fragment
	   */
	  temp.avg_cov=sum/(end-start+1);	  
	  for(k=0;k<no_of_tracks;k++)
	    temp1[k]=temp1[k]/(end-start+1);
	  temp.cov_exp=temp1;
	  sum=0;
	}
	if(temp.state==EXONM){
	  for(k=0;k<input_set[STRANDM].size();k++){//loops over the tracks in - strand for exon- state
	    for(j=start;j<end+1;j++){//loops in the positions from start to end in the track
	      sum=sum+input_set[STRANDM][k][j-1];
	      temp1[k]+=input_set[STRANDM][k][j-1];//sum to get average for particular track
	    }
	  }
	  for(k=0;k<input_set[STRANDB].size();k++){//loops over the tracks in . strand for exon- state
	    for(j=start;j<end+1;j++){//loops in the positions from start to end in the track
	      sum=sum+input_set[STRANDB][k][j-1];
	      temp1[k]+=input_set[STRANDB][k][j-1];//sum to get average for particular track
	    }
	  }
	  temp.avg_cov=sum/(end-start+1);
	  for(k=0;k<no_of_tracks;k++)
	    temp1[k]=temp1[k]/(end-start+1);
	  temp.cov_exp=temp1;
	  sum=0;
	}
	if(temp.state==EXONP || temp.state==EXONM ||( (temp.state==INTRONJ||temp.state==INTRONI) &&include_intron==1)){
	  segmentation.push_back(temp);
	}
	start=end+1;
	end=end+lengths[state_seq[i]];
      }
    } 
    /* Pushing the last state
     */ 
    if(state_seq[i-1]==EXONP || state_seq[i-1]==EXONM ||( (state_seq[i-1]==INTRONJ||state_seq[i-1]==INTRONI) &&include_intron==1)){
      temp.state=state_seq[i-1];
      temp.start=start;
      temp.end=end;
      temp.cov_exp=temp1;
      if(state_seq[i-1]==INTRONJ||state_seq[i-1]==INTRONI)
	temp.avg_cov=0;
      segmentation.push_back(temp);
    }    
    return segmentation;		     
  }

  /* calculate m[a,b]=(psum[b]-psum[a-1])/(psuml[b]-psuml[a-1])
   * to be used for calculating pott's functionals
   */
  double exon_segmentation::calculate_m(int a,int b,vector<vector<int> > &psum1,vector<vector<double> > &psum,int exp){
    return (psum[exp][b]-psum[exp][a-1])/(psum1[exp][b]-psum1[exp][a-1]);
  }
  /* Calculate S(a,b) to be used in the pott's functional
   */
  double exon_segmentation::calculate_s(int a,int b,vector<vector<int> > &psum1,vector<vector<double> > &psum,vector<vector<double> > &psqrsum,vector<vector<double> > &pcubesum,vector<fragment> &segments,int exp){
    /* we calculate s[a,b]= summation ( l[i]*c[i]*(c[i]-m[a,b])^2=l[i]*c[i]^3+l[i]*c[i]*m[a,b]^2-2*l[i]*c[i]^2*m[a,b]
     * =(psqrsum[b]-psqrsum[a-1])+(m[a,b])^2*(psuml[b]-psuml[a-1])-2*m[a,b]*(psum[b]-psum[a-1])
     * We use a-1 instead of a because we have to include the index a itself
     */
    return ( (pcubesum[exp][b]-pcubesum[exp][a-1])+pow(calculate_m(a,b,psum1,psum,exp),2)*(psum[exp][b]-psum[exp][a-1])-2*calculate_m(a,b,psum1,psum,exp)*(psqrsum[exp][b]-psqrsum[exp][a-1]) );

  }
/* Function to convert the average coverage depths into some usable form
 * in teh pott's functional
 */
double exon_segmentation::pott_convert(double d){
  if(d<1)
    return d;
  else 
    return log(d);
}
  /* the function which takes as input the vector of fragments and clusters the exons together to form genes
   */
  vector<int> exon_segmentation::cluster_exons(vector<fragment> segments,int no_of_tracks){
    cout<<"Started clustering the exons."<<endl;
    int i,k,n,argmin,m,exp;
    /* define the vectors to store the partial sums and partials square sums
     */
    double min,temp,sum;
    vector<vector<double> > x,psum,psqrsum,pcubesum;
    vector<vector<double> > trans_avg_covs;
    vector<vector<int> > psum1;
    vector<int> p,temp1;
    vector<double> h,temp2;
    vector<int> J,J_temp;        
    /* store the number of exons in the variable n;
     */
    n=segments.size();

    temp2.resize(n,0.0);
    temp1.resize(n+1,0);

    /* Convert the average scores to the log values
     */
    trans_avg_covs.resize(no_of_tracks,temp2);
    for(k=0;k<no_of_tracks;k++)
      for(i=0;i<n;i++)
	trans_avg_covs[k][i]=pott_convert(segments[i].cov_exp[k]);
    /* initnialize the partial sum vectors
     * psum stores the summation of l_i*c_i
     * psqrsum stores the summation l_i*c_i^2
     * pcubesum stores the summation l_i*c_i^3
     * psuml stores the summattion l_i
     */
    temp2.resize(n+1,0.0);
    psum.resize(no_of_tracks,temp2);
    psqrsum.resize(no_of_tracks,temp2);
    pcubesum.resize(no_of_tracks,temp2);
    psum1.resize(no_of_tracks,temp1);

    /* initialize the vectors required for the dynamic programing matrix
     */
    h.resize(n+1,0.0);
    x.resize(no_of_tracks,temp2);
    p.resize(n+1,0);
    /* calculate the first element of the partial sum vector
     */
    for(k=0;k<no_of_tracks;k++){
      psum1[k][1]=(segments[0].end-segments[0].start+1);
      psum[k][1]=(segments[0].end-segments[0].start+1)*trans_avg_covs[k][0];
      psqrsum[k][1]=(segments[0].end-segments[0].start+1)*trans_avg_covs[k][0]*trans_avg_covs[k][0];
      pcubesum[k][1]=(segments[0].end-segments[0].start+1)*trans_avg_covs[k][0]*trans_avg_covs[k][0]*trans_avg_covs[k][0];
    }
    for(k=0;k<no_of_tracks;k++){
      for(i=2;i<=n;i++){
	/* Calculate partial sums by using the previous entry in the vector
	 */
	psum1[k][i]=psum1[k][i-1] + (segments[i-1].end-segments[i-1].start+1);
	psum[k][i]=psum[k][i-1]+(segments[i-1].end-segments[i-1].start+1)*trans_avg_covs[k][i-1];
	psqrsum[k][i]=psqrsum[k][i-1]+(segments[i-1].end-segments[i-1].start+1)*pow(trans_avg_covs[k][i-1],2);
	pcubesum[k][i]=pcubesum[k][i-1]+(segments[i-1].end-segments[i-1].start+1)*pow(trans_avg_covs[k][i-1],3);
      }
    }
    k=0;
    /* Intialiazation for the first element of the dynamic programing table
     */
    /* Change the value of pott_gamma to no_of_tracks*pott_gamma
     * to account for multiple number of experiments in the file
     */
    h[0]=-1*pott_gamma*no_of_tracks;
    /* Start the dynamic programing algorithm
     */
    for(k=1;k<=n;k++){
      /* Initialize min value to the first entry
       */
      sum=0;
      /* sum the calue of s over all the experimetns
       */
      for(exp=0;exp<no_of_tracks;exp++)
	sum+=calculate_s(1,k,psum1,psum,psqrsum,pcubesum,segments,exp);
      min=h[0]+sum;
      argmin=0;
      for(i=0;i<k;i++){	
	sum=0;
	for(exp=0;exp<no_of_tracks;exp++)
	  sum=sum+calculate_s(i+1,k,psum1,psum,psqrsum,pcubesum,segments,exp);
	temp=h[i] + sum;
	if(temp<min){
	  min=temp;
	  argmin=i;
	}
      }
      p[k]=argmin;
      h[k]=pott_gamma*no_of_tracks + min;
    }

    /* Code to do the backtracking
     */
    k=n;
    J.clear();
    while(k>0){
      for(exp=0;exp<no_of_tracks;exp++)
	x[exp][k]=calculate_m(p[k]+1,k,psum1,psum,exp);
      k=p[k];
      if(k>0)
	J.push_back(k);
    }
    /* Reverse the J vector and store it in J_temp vector
     */ 
    J_temp.resize(J.size(),0);
    k=J.size();
    for(i=0,m=k-1;i<k;i++,m--)
      J_temp[i]=J[m];
    cout<<"Finshed clustering."<<endl;
    cout<<"J vector size :"<<J_temp.size()<<endl;
    return J_temp;
  }

  /* Takes the dataset, the name of output file and the threshold c and the ratio of complement 
   * of the paramters(s) as input and makes the gff file
   * train_file specifies the gff file to be used for training the distributions
   * flag_r specifies whether the function will be iterated by training thhe distribtuion from the
   * output produced
   */
void exon_segmentation::make_gff(dataset &coverage_info,string outputfile,int c,float s,string train_file,int flag_r,string bed_file,string augustus_hints,int chunksize){
    int i,J_pointer=0,include_intron=0;
    vector<int> state_seq,temp_state_seq,J;
    vector< vector< vector<int> > > input_set;
    vector<fragment> segmentation,splice_sites;
    /* Assign the default values to the file names
     */
    string file_augustus="augustus_hints.gff",file_bed="output.bed";
    ofstream outfile(outputfile.c_str(),ios::out);
    /* Check for the names of the output bed file and augustus hints file
     * specified.Else give the default values
     */
    if(augustus_hints!=" ")
      file_augustus=augustus_hints;
    if(bed_file!=" ")
      file_bed=bed_file;       
    ofstream bed_format(file_bed.c_str(),ios::out);
    ofstream augustus_hint(file_augustus.c_str(),ios::out);
    /* Call the function to find the splice sites
     */
    //find_splice(coverage_info.database["chr2L"],coverage_info);

    for(coverage_info.it=coverage_info.database.begin();coverage_info.it!=coverage_info.database.end();++coverage_info.it){
      cout<<"chromosome considering for calculation: "<<coverage_info.it->first<<endl;

      /* Calculate the probaility distribution for introns and exons
       */
      if(train_file==" "){//no file for training is specified
	/* call the function to calculate the parameters for the two geometric distributions
	 */
	calculate_lambda(coverage_info,c,s);
	/* Call the function to make the geometric distribution table using the lambdas
	 */
	prob_density(coverage_info.max_depth,coverage_info.no_of_tracks);
      }
      else{//if the file for training is specified	
	/* Call the function to make the geometric distributtion vector after training 
	 * from the specified gff file and the given data set and also it takes the 
	 * key to the chromosome as input
	 */    
	train_function(coverage_info,train_file,coverage_info.it->first);
      }



      input_set.clear();
      state_seq.clear();
      segmentation.clear();
      /* store the information in input_set for precossing
       */
      input_set=coverage_info.it->second;
      /* call the viterbi function to get the state sequence for each chunk
       * We call the viterbi function for each chunk present in the given data
       */
      for(i=0;i<coverage_info.max_len/chunksize+1;i++){
	temp_state_seq.clear();
	temp_state_seq=viterbi(input_set,coverage_info.max_depth,coverage_info.max_len,chunksize,i);
	/* Now append the temporary state sequence to the state sequence vector
	 */
	state_seq.reserve( state_seq.size() + temp_state_seq.size());
	state_seq.insert( state_seq.end(), temp_state_seq.begin(), temp_state_seq.end());
      }

      
      if(flag_r==1){//if the flag to do the train the distribution from output is specified
	/* give the state sequence as input to the segment function
	 * and get the vector of fragments as output including the intron as it is required for
	 * training the distribution
	 */
	include_intron=1;
	segmentation=segment(state_seq,input_set,include_intron,coverage_info.no_of_tracks);
	include_intron=0;
	//for(i=0;i<segmentation.size();i++)
	//cout<<segmentation[i].start<<'\t'<<segmentation[i].end<<'\t'<<segmentation[i].state<<'\t'<<segmentation[i].avg_cov<<endl;
	/* Make the temporary output file
	 */
	ofstream outfile1("temp.gff",ios::out);
	outfile1<<"browser position "<<coverage_info.it->first<<":1-100000\n";
	outfile1<<"track name=seg description=\"Segmentation first round\" visibility=2\n";
	for(i=0;i<segmentation.size();i++){
	  if(segmentation[i].state==EXONP)
	    outfile1<<coverage_info.it->first<<"\tsource\t"<<"exon\t"<<segmentation[i].start<<'\t'<<segmentation[i].end<<'\t'<<"1\t+\t0\tgrp"<<i<<"-"<<segmentation[i].avg_cov<<"\n";
	  else if(segmentation[i].state==EXONM)
	    outfile1<<coverage_info.it->first<<"\tsource\t"<<"exon\t"<<segmentation[i].start<<'\t'<<segmentation[i].end<<'\t'<<"1\t-\t0\tgrp"<<i<<"-"<<segmentation[i].avg_cov<<"\n";
	  else
	    outfile1<<coverage_info.it->first<<"\tsource\t"<<"intron\t"<<segmentation[i].start<<'\t'<<segmentation[i].end<<'\t'<<"1\t-\t0\tgrp"<<i<<"-"<<segmentation[i].avg_cov<<"\n";
	  
	}
	outfile1.close();	
	/* Call function to train the datasets using the output file
	 */
	train_function(coverage_info,"temp.gff",coverage_info.it->first);
	input_set.clear();
	state_seq.clear();
	segmentation.clear();
	/* store the information in input_set for precossing
	 */
	input_set=coverage_info.it->second;
	/* call the viterbi function to get the state sequence
	 */
	/* call the viterbi function to get the state sequence for each chunk
	 * We call the viterbi function for each chunk present in the given data
	 */
	for(i=0;i<coverage_info.max_len/chunksize+1;i++){
	  temp_state_seq.clear();
	  temp_state_seq=viterbi(input_set,coverage_info.max_depth,coverage_info.max_len,chunksize,i);
	  /* Now append the temporary state sequence to the state sequence vector
	   */
	  state_seq.reserve( state_seq.size() + temp_state_seq.size());
	  state_seq.insert( state_seq.end(), temp_state_seq.begin(), temp_state_seq.end());
	}	
      }           
      /* give the state sequence as input to the segment function
       * and get the vector of fragments as output
       */
      /* We only want the exons in the segment vectors as we have to cluster them together
       */
      //include_intron=1;
      segmentation=segment(state_seq,input_set,include_intron,coverage_info.no_of_tracks);
      //for(i=0;i<segmentation.size();i++)
      //cout<<segmentation[i].start<<'\t'<<segmentation[i].end<<'\t'<<segmentation[i].state<<'\t'<<segmentation[i].avg_cov<<endl;
      /* call the function to cluster the exons
       */
      J=cluster_exons(segmentation,coverage_info.no_of_tracks);
      /*make the output gff file from the vector of fragments
       */
      outfile<<"browser position "<<coverage_info.it->first<<":1-100000\n";
      outfile<<"track name=seg description=\"Segmentation first round\" visibility=2\n";
      augustus_hint<<"browser position "<<coverage_info.it->first<<":1-100000\n";
      augustus_hint<<"track name=seg description=\"Segmentation first round\" visibility=2\n";
      bed_format<<"browser position "<<coverage_info.it->first<<":1-100000\n";
      bed_format<<"track name=\"Exon positions\" description=\"Clustering First Round\" visibility=2 itemRgb=\"On\" useScore=1\n";
      for(i=0;i<segmentation.size();i++){
	/* If we encounter a position to jump then we increase the value of J_pointer
	 * Also we store the relevant information for approximate gene boundaries
	 * in the augustus hint file
	 */	
	if(J[J_pointer]==i+1){
	  J_pointer++;
	  augustus_hint<<coverage_info.it->first<<"\twig2hints\t"<<"diffgene\t"<<(segmentation[i].end+segmentation[i].start)/2<<'\t'<<(segmentation[i+1].end+segmentation[i+1].start)/2<<'\t'<<"-"<<"\t+\t.\tsrc=W;left="<<segmentation[i].start<<".."<<segmentation[i].end<<";right="<<segmentation[i+1].start<<".."<<segmentation[i+1].end<<"\n";
	}
	if(segmentation[i].state==EXONP){
	  outfile<<coverage_info.it->first<<"\tsource\t"<<"exon\t"<<segmentation[i].start<<'\t'<<segmentation[i].end<<'\t'<<"1\t+\t0\tgrp"<<i<<"-"<<J_pointer<<"-"<<segmentation[i].avg_cov<<"\n";
	  augustus_hint<<coverage_info.it->first<<"\twig2hints\t"<<"exonpart\t"<<segmentation[i].start<<'\t'<<segmentation[i].end<<'\t'<<segmentation[i].avg_cov<<"\t+\t.\tsrc=W"<<"\n";
	  bed_format<<coverage_info.it->first<<"\t"<<segmentation[i].start<<'\t'<<segmentation[i].end<<'\t'<<"exon-"<<segmentation[i].avg_cov<<"\t"<<segmentation[i].avg_cov<<"\t+\t"<<segmentation[i].start<<'\t'<<segmentation[i].end<<"\t0,255,0\n";
	}
	else if(segmentation[i].state==EXONM){
	  outfile<<coverage_info.it->first<<"\tsource\t"<<"exon\t"<<segmentation[i].start<<'\t'<<segmentation[i].end<<'\t'<<"1\t-\t0\tgrp"<<i<<"-"<<J_pointer<<"-"<<segmentation[i].avg_cov<<"\n";
	  augustus_hint<<coverage_info.it->first<<"\twig2hints\t"<<"exonpart\t"<<segmentation[i].start<<'\t'<<segmentation[i].end<<'\t'<<segmentation[i].avg_cov<<"\t-\t.\tsrc=W"<<"\n";
	  bed_format<<coverage_info.it->first<<"\t"<<segmentation[i].start<<'\t'<<segmentation[i].end<<'\t'<<"exon-"<<segmentation[i].avg_cov<<"\t"<<segmentation[i].avg_cov<<"\t-\t"<<segmentation[i].start<<'\t'<<segmentation[i].end<<"\t0,0,255\n";
	}
	else
	  /* only print the exons in the gff file for better visualization
	   */
	  continue;
      //outfile<<"chr2L\tsource\t"<<"intron\t"<<segmentation[i].start<<'\t'<<segmentation[i].end<<'\t'<<"1\t.\t0\tgrp"<<i<<"-"<<segmentation[i].avg_cov<<"\n";
      }    
    /* Add the gene jumps in a separate track to the bed file
     */
      bed_format<<"browser position "<<coverage_info.it->first<<":1-100000\n";
      bed_format<<"track name=\"Gene Jumps\" description=\"Clustering First Round\" visibility=2 itemRgb=\"On\" useScore=1\n";
      for(i=0;i<J.size();i++){
	bed_format<<coverage_info.it->first<<"\t"<<(segmentation[J[i]-1].end+segmentation[J[i]-1].start)/2<<'\t'<<(segmentation[J[i]].end+segmentation[J[i]].start)/2<<'\t'<<"GENE_JUMP\t"<<"500"<<"\t+\t"<<(segmentation[J[i]-1].end+segmentation[J[i]-1].start)/2<<"\t"<<(segmentation[J[i]].end+segmentation[J[i]].start)/2<<"\t255,0,0\n";
      }      
    }
    outfile.close();  
    augustus_hint.close();
    bed_format.close();
    cout<<"Finished writing to "<<outputfile<<endl;
    
    //cout<<"Making graph data file"<<endl;
    //make_graph(coverage_info,segmentation);
  }
  /* takes the dataset and vector of fragments as input and
   * stores the distribution of number of positions with the 
   * coverage depths for both introns and exons 
   * in separate files for plotting
   */
  void exon_segmentation::make_graph(dataset &coverage_info,vector<fragment > &segments){
    int i,l,seg_ptr=0;
    int strand,track;
    vector<int>  intron_data,temp;
    vector< vector<int> > exon_data;
    temp.resize(2,0);
    exon_data.resize(coverage_info.max_depth+1,temp);
    intron_data.resize(coverage_info.max_depth+1,0);
    
    for(coverage_info.it=coverage_info.database.begin();coverage_info.it!=coverage_info.database.end();++coverage_info.it)
	for(strand=0;strand<coverage_info.it->second.size();strand++)// loop over strands +,-,b
	  for(track=0;track<coverage_info.it->second[strand].size();track++){// loop over tracks
	    for(l=0;l<coverage_info.it->second[strand][track].size();l++){// loop over position on track
	      if(l+1>segments[seg_ptr].end)
		seg_ptr++;
	      if(segments[seg_ptr].state==EXONP && strand!=STRANDM)
		exon_data[coverage_info.it->second[strand][track][l]][0]++;
	      else if(segments[seg_ptr].state==EXONM && strand!=STRANDP)
		exon_data[coverage_info.it->second[strand][track][l]][1]++;
	      else
		intron_data[coverage_info.it->second[strand][track][l]]++;
	    } 
	    seg_ptr=0;
	  }  
    ofstream outfile1("exon_graph.dat",ios::out);
    ofstream outfile2("intron_graph.dat",ios::out);
    outfile1<<"#depths"<<'\t'<<"freq for +\tfreq for-\n";
    outfile2<<"#depths"<<'\t'<<"frequency for introns\n";
    for(i=0;i<exon_data.size();i++){
      outfile1<<i<<'\t'<<exon_data[i][0]<<'\t'<<exon_data[i][1]<<'\n';
      outfile2<<i<<'\t'<<intron_data[i]<<'\n';
    }
    outfile1.close();
    outfile2.close();
  }

