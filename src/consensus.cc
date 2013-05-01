/**********************************************************************
 * file:    consensus.cc
 *
 * descr.:  finds the consensus pattern
 * authors: Sukadeb Acharya, sukadeb@gmail.com
 *
 **********************************************************************/

#include "consensus.hh"

const int max_line_len=100;
const int nRep=20;

using namespace std;

consensus::consensus(int starting1,int ending1){
  //gives default value to the paramters
  pattern_size=6;
  p_value=.01;
  delta=2;
  max_allowed_mismatch=1;
  no_iterations=10;
  starting=starting1;
  ending=ending1;
  analyse_pattern="none";
  
  //clears all the vectors used 
  var.clear();
  powers.clear();
  t_values.clear();
  f_values.clear();
  r_values.clear();
  mean_values.clear();
  significant_strings.clear();
  relevant_strings.clear();
  final_list.clear();
  consensus_data.clear();
  
  max_string_length=0;
  int i;
  // powers stores all the powers of 4 upto the pattern length for calculation purpose
  for(i=pattern_size-1;i>=0;i--) 
    powers.push_back(pow((double)4,(double)i));

  //initialises the mean value, significance value, relevant value vectors
  for(i=0;i<pow((double)4,(double)pattern_size);i++){
    t_values.push_back(0);
    f_values.push_back(0);
    r_values.push_back(0);
    mean_values.push_back(0);
  }
}

void consensus::set_values(int pattern_size1,float p_value1, float delta1, int max_allowed_mismatch1,int n){
  
  int i;
  pattern_size=pattern_size1;
  p_value=p_value1;
  delta=delta1;
  max_allowed_mismatch=max_allowed_mismatch1;
  no_iterations=n;
  powers.clear();
  for(i=pattern_size-1;i>=0;i--) 
    powers.push_back(pow((double)4,(double)i));

  for(i=0;i<pow((double)4,(double)pattern_size);i++){
    t_values.push_back(0);
    f_values.push_back(0);
    r_values.push_back(0);
    mean_values.push_back(0);
  }
}
void consensus::set_file_name(string filename){
  //reading the file and storing the strings
  string fileline,fileline1;
  int len1,temp_starting,temp_ending;
  sequence temp;
  ifstream infile (filename.c_str());
  //ignores the first line
  getline(infile,fileline1);
  //reads the second line 
  getline(infile,fileline1);
  while(!infile.eof()){
    if(fileline1[0]!='>'){
      fileline.append(fileline1);
    }
    else{
      //modifying code for negative inputs
      if(starting<0){
	temp_starting=fileline.length()+starting;
	temp_ending=fileline.length()+ending;
      }
      else{
	temp_starting=starting;
	temp_ending=ending;
      }
      char *buf=new char[fileline.length()];//stores the required length of string
      if ((starting>=0 && fileline.length()>temp_starting) || (starting<0 && temp_ending>0)){//if the length of string is less than start then ignore it
	if( (starting>=0 && (ending==0||ending>fileline.length()))|| (starting<0 && temp_starting<0)){//if end is not defined or it is shorter then take the available string
	  if(starting >=0){
	    len1=fileline.length()-starting+1;
	  }
	  else{
	    len1=fileline.length()+ending+1;
	  }
	}
	else{
	  if(starting>=0){
	    len1=ending-starting+1;
	  }
	  else{
	    len1=ending-starting+1;
	  }
	}	  
	temp.seq=buf;//string is copied to the structure
	temp.seq.erase(len1-1);//extra length is erased
	temp.length=len1;
	var.push_back(temp);
	if(len1>max_string_length)
	  max_string_length=len1;
	fileline="";
      }
    }
    
    getline(infile,fileline1);	
  }
  //to add the last string in the file
  if(starting<0){
    temp_starting=fileline.length()+starting;
    temp_ending=fileline.length()+ending;
  }
  else{
    temp_starting=starting;
    temp_ending=ending;
  }
  char *buf=new char[fileline.length()];//stores the required length of string
  if ((starting>=0 && fileline.length()>temp_starting) || (starting<0 && temp_ending>0)){//if the length of string is less than start then ignore it
    if( (starting>=0 && (ending==0||ending>fileline.length()))|| (starting<0 && temp_starting<0)){//if end is not defined or it is shorter then take the available string
      if(starting >=0){
	len1=fileline.length()-starting+1;
      }
      else{
	len1=fileline.length()+ending+1;
      }
    }
    else{
      if(starting>=0){
	len1=ending-starting+1;
      }
      else{
	len1=ending-starting+1;
      }
    }
    temp.seq=buf;
    temp.seq.erase(len1-1);
    temp.length=len1;
    var.push_back(temp);
    if(len1>max_string_length)
      max_string_length=len1;
    fileline="";
  } 
  infile.close();
  number_seq=var.size();
  cout << " finished reading the file " << endl;
}


//function add some extra strings
void consensus::add_string(string fileline){
  int len1,temp_starting,temp_ending;
  sequence temp;
  if(starting<0){
    temp_starting=fileline.length()+starting;
    temp_ending=fileline.length()+ending;
  }
  else{
    temp_starting=starting;
    temp_ending=ending;
  }
  char *buf=new char[fileline.length()];//stores the required length of string
  if ((starting>=0 && fileline.length()>temp_starting) || (starting<0 && temp_ending>0)){//if the length of string is less than start then ignore it
    if( (starting>=0 && (ending==0||ending>fileline.length()))|| (starting<0 && temp_starting<0)){//if end is not defined or it is shorter then take the available string
      if(starting >=0){
	len1=fileline.length()-starting+1;
      }
      else{
	len1=fileline.length()+ending+1;
      }
    }
    else{
      if(starting>=0){
	len1=ending-starting+1;
      }
      else{
	len1=ending-starting+1;
      }
    }
    temp.seq=buf;
    temp.seq.erase(len1-1);
    temp.length=len1;
    var.push_back(temp);
    if(len1>max_string_length)
      max_string_length=len1;
    number_seq++;
  }
}

float consensus::calculate_cumulative_value(float z){
  int k;
  float normal_table[38]={0,.0398,.0793,.1179,.1554,.1915,.2258,.2580,.2881,.3159,.3413,.3643,.3849,.4032,.4192,.4332,.4452,.4554,.4641,.4713,.4772,.4821,.4861,.4893,.4918,.4938,.4953,.4965,.4974,.4981,.4987,.4990,.4993,.4995,.4997,.4998,.4998,.4999};
  if(z>=0 && z< 3.8){
      k=(int)z*10;
      return (.5 - (normal_table[k]+ (z-.1*k)*(normal_table[k+1]-normal_table[k])/.1));
    }
    else if(z>=3.8){
      return 0;
	
    }
    else{
      k=(int)(-1*(z*10));
      if(k<37){
	return(0.5 + normal_table[k]+ (-1*z-.1*k)*(normal_table[k+1]-normal_table[k])/.1);
      } 
      else 
	return 1;
    }
}

double consensus::poisThresh(int i, float p){
  double a=0,b=i;
  return poisInv(a,b,p,nRep,i);
}

double consensus::poisInv(double a,double b,float p,int n,int i){
  double sum=0,lambda=(a+b)/2;
  int j;
  if(n<=0)
    return (a+b)/2;
  for(j=0;j<=i;j++)
    sum=sum+pow(2.73,-1*lambda)*pow(lambda,(double)j)/factorial(j);
  if(sum>1-p)
    return poisInv(a,(a+b)/2,p,n-1,i);
  else
    return poisInv((a+b)/2,b,p,n-1,i);
}


void consensus::start(){
  int i,j,k,m,sum,sum1,number;
  float L[30];//stores the poison thresholds for calculating the significance values
  for(i=0;i<30;i++)
    L[i]=poisThresh(i+1,p_value);//calls the function to calculate the threshold values
  static Seq2Int s2i(pattern_size);
  for ( i=0;i<4;i++)
    for(j=0;j<4;j++)
      tpm[i][j]=0;
  int row=0,column=0;
  sum=0;
  //formulating the transition probability matrix
  for(int i=0;i<number_seq;i++)
    {
      sum+=var[i].seq.length();
      for(int j=1;j<var[i].seq.length();j++)
	{
	  row=char2num(var[i].seq[j-1]);
	  column=char2num(var[i].seq[j]);
	  if(row>-1 && column > -1)
	    tpm[row][column]++;
	}
    }
  //calculating the background probabilites
  float k1=0;
  for(i=0;i<4;i++){
    k1=0;
    for(j=0;j<4;j++)
      k1+=tpm[j][i];
    back_probs[i]=k1/sum;
  }

  //normalizing the matrix
  for(i=0;i<4;i++){
    sum=0;
    for(j=0;j<4;j++)
      sum+=(int)tpm[i][j];
    for(j=0;j<4;j++)
      tpm[i][j]=tpm[i][j]/sum;	
  }
  //calculating the f values
  for(i=0;i<number_seq;i++){ 
    if(var[i].seq.length()<pattern_size) 
      continue;
    for(k=0;k<pow((double)4,(double)pattern_size);k++)
      t_values[k]=0;
    for( j=0;j<var[i].seq.length()- pattern_size + 1;j++){ 
      try{	  
	number = s2i(var[i].seq.c_str()+j);
	if(t_values[number]==0)
	  {
	    f_values[number]++;
	    t_values[number]=1;
	  } 
      }
      catch(exception& e ){
	//ignoring the exception
      }
    }
  }
  cout << " calculated the frequencies " << endl;
  //calculating the mean value and storing the significance patterns
  sum=0;
  //sum stores the total number of places where a pattern can begin
  for(j=0;j<number_seq;j++)
    sum+=var[j].seq.length()-pattern_size+1;
  string a;
  for(i=0;i<pow((double)4,(double)pattern_size);i++){
    a=s2i.inv(i);
    mean_values[i]=back_probs[char2num(a[0])];
      
    for(j=1;j<pattern_size;j++){
      row=char2num(a[j-1]);
      column=char2num(a[j]);
      mean_values[i]=mean_values[i]*tpm[row][column];
    }
    mean_values[i]=mean_values[i]*sum;//finished calculating the mean value for the pattern i

    //starting the significant test
    if(f_values[i]>30){//use the normal distribution table
      float z=calculate_cumulative_value((f_values[i] - mean_values[i])/sqrt(mean_values[i]));    
      if(z<p_value)
	significant_strings.push_back(i);
    }
    else{//use the poisson distribution
      if(L[f_values[i]-1]>mean_values[i])
	 significant_strings.push_back(i);
    }
  }
  cout << " calculated the significant values" << endl;
  //calculating the relevance values and storing the relevant ones
  for(i=0;i<pow((double)4,(double)pattern_size);i++)
    r_values[i]=f_values[i]/mean_values[i];

  for(i=0;i<significant_strings.size();i++)
    if( r_values[significant_strings[i]] > delta)
      relevant_strings.push_back(significant_strings[i]);
  sum=0;

  int max=0,max_index=0,counter=0;
  vector<int> prev_neighbours,new_neighbours;
  vector<int> neighbours_index;    
  //greedy method to calculate the cumulative neighbourhood significance values
  cout << "Calculating the greedy way " << endl;
  while(relevant_strings.size()>0){    
    counter++;
    if(counter>no_iterations)
      break;
    max=0;
    for(i=0;i<relevant_strings.size();i++){
      sum=0;//initialize the sum value
      new_neighbours.clear();
      
      //adds the cumulative f values for all the k neighbours
      for(j=0;j<=max_allowed_mismatch;j++){
	sum1=0;

	neighbours_index.clear();//clears the previous neighbour values
	//find all the neighbours with j mismatches
	neighbours_index=find_neighbours(relevant_strings[i],j,pattern_size);
	//store the neighbours for adding the values and erasing from the lists
	for(m=0;m<neighbours_index.size();m++)
	  for(k=0;k<significant_strings.size();k++){
	    if(neighbours_index[m]==significant_strings[k]){
	      sum1+=f_values[neighbours_index[m]];
	      new_neighbours.push_back(k);
	      break;
	    }
	  }
	sum=sum+sum1/(int)(pow((double)3,(double)j)*nCr(pattern_size,j));
      }

      if(sum>=max){
	max=sum;
	max_index=i;
	prev_neighbours.clear();
	prev_neighbours=new_neighbours;
	new_neighbours.clear();
      }
    }
    //storing details of the consensus pattern
    consensus_data.push_back(relevant_strings[max_index]);
    histogram_data temp1;
    temp1.consensus_pattern=s2i.inv(relevant_strings[max_index]);
    temp1.f_value=f_values[relevant_strings[max_index]];
    temp1.mean_value=mean_values[relevant_strings[max_index]];

    for(i=1;i<prev_neighbours.size();i++){
      temp1.neighbours.push_back(s2i.inv(significant_strings[prev_neighbours[i]]));
      temp1.freq_neighbours.push_back(f_values[significant_strings[prev_neighbours[i]]]);
      temp1.mean_neighbours.push_back(mean_values[significant_strings[prev_neighbours[i]]]);
    }

    int number=0,number1=0,temp_number=0,flag=0;
    a=s2i.inv(relevant_strings[max_index]);
    cout << counter <<'\t'<< " consensus pattern   " << a<< endl; 
    number=relevant_strings[max_index];

    //storing the left shift neighbour for deletion
    number1=(number-char2num(a[0])*((int)powers[0]))*4;
    for(i=0;i<4;i++){
      temp_number=number1+i;
      if( temp_number!=number)
	for(j=0;j<significant_strings.size();j++)
	  if(significant_strings[j]==temp_number){
	    for(k=0;k<prev_neighbours.size();k++)
	      if(prev_neighbours[k]==j){
		flag=1;
		break;
	      }
	    if(flag==0)
	      prev_neighbours.push_back(j);
	    flag=0;
	    break;
	  }
    }
    
    //storing the right shift neghbour for deletion
    number1=number/4;
    flag=0;
    for(i=0;i<4;i++){
      temp_number=number1+i*((int)powers[0]);
      if( temp_number!=number)
	for(j=0;j<significant_strings.size();j++)
	  if(significant_strings[j]==temp_number){
	    for(k=0;k<prev_neighbours.size();k++)
	      if(prev_neighbours[k]==j){
		flag=1;
		break;
	      }
	    if(flag==0)
	      prev_neighbours.push_back(j);
	    flag=0;
	    break;
	  }
    }

    //deleting the neighbours from both lists
    for(i=0;i<prev_neighbours.size();i++){
      k=significant_strings[prev_neighbours[i]];
      //erasing neighbours from the significant list
      
      significant_strings.erase(significant_strings.begin()+ prev_neighbours[i]);
      //erasing the neighours from relevant list
      for(j=0;j<relevant_strings.size();j++)
	if(relevant_strings[j]==k ){
	  relevant_strings.erase(relevant_strings.begin() + j);
	  break;
	}
      
      //adjusting the indices to accomodate for the deleted element
      for(j=i+1;j<prev_neighbours.size();j++)
	if(prev_neighbours[j]>prev_neighbours[i])
	  prev_neighbours[j]=prev_neighbours[j]-1;
    }
    
    prev_neighbours.clear();
    final_list.push_back(temp1);
  }
  cout<<" calculated the first "<< no_iterations << " consensus patterns."  << endl; 

  //storing the histogram data in a structure
  max_freq=0;
  for(k=0;k<final_list.size();k++){
    cout << " storing histogram for sequence " << final_list[k].consensus_pattern<<endl;
    for(i=0;i<max_string_length;i++)
      final_list[k].position_of_occurence.push_back(0);
    for(i=0;i<number_seq;i++){      
      if(var[i].seq.length()<pattern_size) 
	continue;
      //for(k=0;k<pow((double)4,(double)pattern_size);k++)
      //t_values[k]=0;
      for( j=0;j<var[i].seq.length()- pattern_size +1;j++){ 
	try{
	  if(consensus_data[k] == s2i(var[i].seq.c_str()+j)){
	    final_list[k].position_of_occurence[j]++;
	    if(final_list[k].position_of_occurence[j]>max_freq)
	      max_freq=final_list[k].position_of_occurence[j];
	  }
	}
	catch(exception& e ){
	  //ignoring the exception
	}
      }
    }
  }
}

void consensus::analyse(string pattern1,float p_value1,float delta1,int max_allowed_mismatch1){
  int i,j,k,sum,number,number1,flag;
  histogram_data data1;
  analyse_pattern=pattern1;
  pattern_size=analyse_pattern.length();
  p_value=p_value1;
  delta=delta1;
  max_allowed_mismatch=max_allowed_mismatch1;
  powers.clear();
  for(i=pattern_size-1;i>=0;i--) 
    powers.push_back(pow((double)4,(double)i));
  static Seq2Int s2i(pattern_size);
  number=s2i(analyse_pattern.c_str());
  vector<int> neighbours_list;
  neighbours_list=find_neighbours(number,max_allowed_mismatch,pattern_size);

  for(i=0;i<pow((double)4,(double)pattern_size);i++){
    t_values.push_back(0);
    f_values.push_back(0);
    r_values.push_back(0);
    mean_values.push_back(0);
  }

  float L[30];//stores the poison thresholds for calculating the significance values
  for(i=0;i<30;i++)
    L[i]=poisThresh(i+1,p_value);//calls the function to calculate the threshold values  
  for ( i=0;i<4;i++)
    for(j=0;j<4;j++)
      tpm[i][j]=0;
  int row=0,column=0;
  sum=0;
  //formulating the transition probability matrix
  for(int i=0;i<number_seq;i++){
    sum+=var[i].seq.length();
    for(int j=1;j<var[i].seq.length();j++){
      row=char2num(var[i].seq[j-1]);
      column=char2num(var[i].seq[j]);
      if(row>-1 && column > -1)
	tpm[row][column]++;
    }
  }
  //calculating the background probabilites
  float k1=0;
  for(i=0;i<4;i++){
    k1=0;
    for(j=0;j<4;j++)
      k1+=tpm[j][i];
    back_probs[i]=k1/sum;
  }

  //normalizing the matrix
  for(i=0;i<4;i++){
    sum=0;
    for(j=0;j<4;j++)
      sum+=(int)tpm[i][j];
    for(j=0;j<4;j++)
      tpm[i][j]=tpm[i][j]/sum;	
  }

  //calculating the f values
  for(i=0;i<number_seq;i++){ 
    if(var[i].seq.length()<pattern_size) 
      continue;
    for(k=0;k<pow((double)4,(double)pattern_size);k++)
      t_values[k]=0;
    for( j=0;j<var[i].seq.length()- pattern_size + 1;j++){ 
      try{	  
	number1 = s2i(var[i].seq.c_str()+j);
	if(t_values[number]==0)
	  {
	    f_values[number1]++;
	    t_values[number1]=1;
	  } 
      }
      catch(exception& e ){
	//ignoring the exception
      }
    }
  }

  cout << " calculated the frequenecy values" << endl;
  data1.consensus_pattern=analyse_pattern;
  data1.f_value=f_values[number];

  sum=0;
  //sum stores the total number of places where a pattern can begin
  for(j=0;j<number_seq;j++)
    sum+=var[j].seq.length()-pattern_size+1;
  flag=0;
  string a;
  neighbours_list.push_back(number);//to calculate the expectance of the given pattern
  for(i=0;i<neighbours_list.size();i++){
    //test whether the neighbour is significant or not
    number1=neighbours_list[i];
    a=s2i.inv(number1);
    mean_values[number1]=back_probs[char2num(a[0])];
    for(j=1;j<pattern_size;j++){
      row=char2num(a[j-1]);
      column=char2num(a[j]);
      mean_values[number1]=mean_values[number1]*tpm[row][column];
    }
    mean_values[number1]=mean_values[number1]*sum;//finished calculating the mean value for the neighbour pattern

    if(f_values[number1]>30){//use the normal distribution table
      float z=calculate_cumulative_value((f_values[number1] - mean_values[number1])/sqrt(mean_values[number1]));    
      if(z<p_value)
	flag=1;
    }
    else{//use the poisson distribution
      if(L[f_values[number1]-1]>mean_values[number1])
	flag=1;
    }
    //store in the neighbour list if they are significant
    if(flag==1&& i!=(neighbours_list.size()-1)){
      data1.neighbours.push_back(s2i.inv(number1));
      data1.freq_neighbours.push_back(f_values[number1]);
      data1.mean_neighbours.push_back(mean_values[number1]);
      flag=0;
    }
  }  
  data1.mean_value=mean_values[number];
  final_list.push_back(data1);
  //storing the histogram
  max_freq=0;
  for(k=0;k<final_list.size();k++){
    cout << " storing histogram for sequence " << final_list[k].consensus_pattern<<endl;
    for(i=0;i<max_string_length;i++)
      final_list[k].position_of_occurence.push_back(0);
    for(i=0;i<number_seq;i++){      
      if(var[i].seq.length()<pattern_size) 
	continue;
      for( j=0;j<var[i].seq.length()- pattern_size +1;j++){ 
	try{
	  if(number == s2i(var[i].seq.c_str()+j)){
	    final_list[k].position_of_occurence[j]++;
	    if(final_list[k].position_of_occurence[j]>max_freq)
	      max_freq=final_list[k].position_of_occurence[j];
	  }
	}
	catch(exception& e ){
	  //ignoring the exception
	}
      }
    }
  }
}

void consensus::plot_histogram(){
  int i,j,k,length;
  for(i=0;i<final_list.size();i++){
    cout<<" The pattern is " << final_list[i].consensus_pattern<< endl;
    cout<<max_line_len<<'\t'<<max_freq<<endl;
    cout<< " The scale is "<< (float)max_line_len/max_freq <<" * = one occurence"<< endl;
    if(starting==0)
      cout <<" and it starts from begining of sequence." <<endl;
    else
      cout <<" and it begins from position "<< starting <<"."<<endl;
    for(j=0;j<max_string_length;j++){      
      cout <<j +1<< '\t';
      if(final_list[i].position_of_occurence[j]==0){
	cout<<endl;
	continue;
      }
      length=(final_list[i].position_of_occurence[j]*max_line_len)/max_freq;
      if(length==0)
	length=1;
      for(k=0;k<length;k++)
	cout <<"*";
      cout<<endl;
    }
  }
}

vector<histogram_data> consensus::get_consensus(){    
  return final_list;
}

void consensus::print_consensus(){
  int i,j;
  cout <<" Printing the consensus patterns, f values, expected values and neighbours" << endl;
  for(i=0;i<final_list.size();i++){
    cout<<final_list[i].consensus_pattern<<'\t'<<final_list[i].f_value<< '\t'<< final_list[i].mean_value<< endl;
    cout<< "The mismatch neighbours are:" << endl;
    for(j=0;j<final_list[i].neighbours.size();j++)
      cout<<final_list[i].neighbours[j]<<'\t'<<final_list[i].freq_neighbours[j]<<'\t'<<final_list[i].mean_neighbours[j]<<endl;
    cout <<endl<<endl;
  }
}

int consensus::factorial( int m){
  if(m<=1)
    return 1;
  else
    return m*factorial(m-1);
}

int consensus::nCr(int n, int r){
  return factorial(n)/( factorial(n-r)*factorial(r) ) ;
}

int consensus::char2num(char m)
{
  switch(m){
  case 'a' : case 'A':
    return 0;
    break;
  case 'c' : case 'C':
    return 1;
    break;
  case 'g' : case 'G':
    return 2;
    break;
  case 't' : case 'T':
    return 3;
    break;
  default:
    return -1;
    break;
  }
  return -1;
}

char consensus::num2char(int a){
  switch(a){
  case 0:
    return 'A';
  case 1:
    return 'C';
  case 2:
    return 'G';
  case 3 :
    return 'T';
  default:
    return 0;
  }
}

//k mismatch function starts here
void consensus::k_mismatch(int new_m, int old_m, int k, int length,Seq2Int s2i,vector<int> &neighbours_index){
  int i,temp_length,temp_k;
  //return if number of mismatches to do is more than the available length
  if(k> length){
    return;
  }
  //push if k_mismatches are done
  if(k==0){
    neighbours_index.push_back(new_m);
    return;
  }

  //without changing the last element
  temp_length=length-1;
  k_mismatch(new_m,old_m,k,temp_length,s2i,neighbours_index);
  //change the last element
  temp_k=k-1;
  for(i=0;i<4;i++){
    int current_int=char2num(s2i.inv(new_m)[length-1]);
    if(current_int==i)
      continue;
    k_mismatch((int)(new_m+(i-current_int)*powers[length-1]),old_m,temp_k,temp_length,s2i,neighbours_index);
  }
    
}

vector<int> consensus::find_neighbours(int m,int k, int length){
  vector<int> neighbours_index;
  static Seq2Int s2i(length);
  neighbours_index.clear();
  k_mismatch(m,m,k,length,s2i,neighbours_index);
  return neighbours_index;
}
//mismatch function ends here

























