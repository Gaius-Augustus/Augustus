#include <iostream>
#include <string>
#include <vector>
#include <getopt.h>
#include <cstdlib>
#include "exon_seg.hh"

using namespace std;


void usage(){
  cout<<"Program to make exon segmentation and grouping the exons together to form the hints for possible position of genes to be used by augustus. \nUsage:\n curve2hints --in=file --outexons=out.gff --c=int --s=float\noptions:\n--in=s\t\tinput file is single file containing tracks in wig format\n--outexons=s\tfile for the output segmentaion\n-h\t\thelp\n--c\t\tspecify the threshold to distinguish between introns and exons in terms of coverage depths (default 4)\n--s\t\tspecify the ratio(1-lambda1)/(1-lambda2) It should be >1 (default 1.12)\n--z\t\tspecify the chunk size (default value 50000)\n--r\t\tTo enable re-estimation of distribution from the output produced\n--train_file\t\tSpecify the gff file to be used for training the distributions.It should contain the position of both introns and exons\n--bed_file\t\tSpecify the output bedfile which contains the position of gene jumps(default output.bed)\n--augustus_hints\t\tspecify the output augustus hint file(default augustus_hints.gff\n"; 
}

/* define the main file 
 */
int main(int argc, char** argv){
  int c_thresh=4,chunksize=50000;
  string inputfile=" ",outputfile=" ",train_file=" ",augustus_hints=" ",bed_file=" ";
  int c,longIndex,flag_r=0;
  float s=1.12;

  dataset coverage_info;
  static const struct option longOpts[]={
    {"in",required_argument,NULL,'i'},
    {"outexons",required_argument,NULL,'o'},
    {"train_file",required_argument,NULL,'t'},
    {"bed_file",required_argument,NULL,'b'},
    {"augustus_hints",required_argument,NULL,'a'},
    {"help",no_argument,NULL,'h'},
    {"re-estimate",no_argument,NULL,'r'},
    {"cthreshold",required_argument,NULL,'c'},
    {"sratio",required_argument,NULL,'s'},
    {"zchunksize",required_argument,NULL,'z'},
    { NULL, no_argument, NULL, 0 }
  };
  opterr = 0;
  while ((c = getopt_long (argc, argv, ":rhi:o:c:s:z",longOpts, &longIndex)) != -1){
    switch (c)
      {
      case 'i':
	inputfile = optarg;
	break;
      case 'o':
	outputfile =optarg;
	break;
      case 't':
	train_file =optarg;
	break;
      case 'b':
	bed_file =optarg;
	break;
      case 'a':
	augustus_hints =optarg;
	break;
      case 'c':
	c_thresh=atoi(optarg);
	break;
      case 'z':
	chunksize=atoi(optarg);
	break;
      case 's':
	s=atof(optarg);
	break;
      case 'r':    	
	flag_r=1;
	break;
      case 'h':    
	usage();
	exit(0);
	break;
      case 0:     /* long option without a short arg */
	cout<<" Wrong option"<<endl;
	usage();
	exit(0);
	break;
      default:
	cout<<" Wrong option"<<endl;
	usage();
	exit(0);
      } 
  }  
  if (  ((argc-optind)==0 && inputfile==" ")||  ( (argc-optind)==0 && outputfile==" ")) {
    usage();
    exit(0);
  }
  exon_segmentation seg1;
  seg1.read_file(coverage_info,inputfile.c_str());
  seg1.make_gff(coverage_info,outputfile,c_thresh,s,train_file,flag_r,bed_file,augustus_hints,chunksize);
}
