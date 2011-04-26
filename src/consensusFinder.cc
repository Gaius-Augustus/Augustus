#include "consensus.hh"

#include <getopt.h>
#include <cstdlib>

using namespace std;

void print_consensus(vector<histogram_data> final_list){
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

void usage(){
  cout<<"Usage:\n ./consensusFinder input fasta file [options] \n -l:pattern length\n -p:p value for finding the significant strings \n -d:delta to find the relevant strings\n -s:starting position\n -e:ending position\n -m:no. of mismatches to consider while finding neighbours\n -c:string to consider for analysis\n -n:to specify the number of consensus patterns to store\n -t:to print the histogram\n -h:help\n"; 

}
int main(int argc, char** argv){
  vector<histogram_data> final_list;
  string filename;
  int pattern_size=6,no_iterations=10;
  float p_value=.01;
  float delta=2;
  int max_allowed_mismatch=1,start=0,end=0,plot_hist=0;
  string pattern="none";
  int c,index,longIndex;
  static const struct option longOpts[]={
    {"pattern_size",required_argument,NULL,'l'},
    {"p_value",required_argument,NULL,'p'},
    {"start",required_argument,NULL,'s'},
    {"end",required_argument,NULL,'e'},
    {"delta",required_argument,NULL,'d'},
    {"mismatch",required_argument,NULL,'m'},
    {"no_consensus",required_argument,NULL,'n'},
    {"print_hist",no_argument,NULL,'t'},
    {"string_to_analyse",required_argument,NULL,'c'},
    {"help",no_argument,NULL,'h'},
    { NULL, no_argument, NULL, 0 }
    };
  opterr = 0;  
  while ((c = getopt_long (argc, argv, ":thl:p:d:m:s:n:e:c:",longOpts, &longIndex)) != -1){
    switch (c)
      {
      case 'l':
	pattern_size = atoi(optarg);
	break;
      case 'p':
	p_value = atof(optarg);
	break;
      case 's':
	start = atoi(optarg);
	break;
      case 'e':
	end = atoi(optarg);
	break;
      case 'd':
	delta = atof(optarg);
	break;
      case 'm':
	max_allowed_mismatch = atoi(optarg);
	break;
      case 'n':
	no_iterations = atoi(optarg);
	break;
      case 't':
	plot_hist=1;
	break;
      case 'c':
	pattern=optarg;
	break;
      case 'h':    
	usage();
	exit(0);
	break;
      case 0:     /* long option without a short arg */
	cout<<" Wrong option"<<endl;
	usage();
	break;

      default:
	cout<<" Wrong option"<<endl;
	usage();
	exit(0);
      } 
  }  
  if((argc-optind)==0){
    usage();
    exit(0);
  }
  for (index = optind; index < argc; index++)
    filename=argv[index];

  consensus consensus1(start,end);
  consensus1.set_file_name(filename);

  if(pattern=="none"){
    consensus1.set_values(pattern_size,p_value,delta,max_allowed_mismatch,no_iterations);    
    consensus1.start();
  }
  else
    consensus1.analyse(pattern,p_value,delta,max_allowed_mismatch);
  final_list=consensus1.get_consensus();
  if(plot_hist==1)
    consensus1.plot_histogram();
  print_consensus(final_list);
  return 0;
}


