/*
 * Name:    pp_simscore.cc
 * Project: Similarity-Score for Protein-Profile and Protein-Sequence  
 * Author:  Lars Gabriel
 *
 *    g++ -o pp_simScore pp_simscore.cc types.o properties.o geneticcode.o pp_profile.o lldouble.o
 *	Usage:	./pp_simScore seq_file.fa profile.prfl
 *
 */


#include <iostream>
#include <fstream>
#include <limits>
#include <math.h>
#include "properties.hh"
#include "pp_profile.hh"

using namespace PP;
	


//std::string subs_file = "../config/profile/blosum62.qij";
//int subs_matrix[20][20];
/*
void read_subs_matrix(std::string file_name){
	std::fstream subs_file;
	subs_file.open(file_name);
	std::string line;
	
	while(std::getline(subs_file, line))
	{
		if(line[0] != '#')
		{
			std::istringstream iss(line);
			
			
		
		}
		
	}

}
	*/
	
	

int aa2int(char c)
{
	switch (c)
	{
		case '.':
			return 20;
		case 'G':
			return 0;
		case 'D':
			return 1;
		case 'E':
			return 2;
		case 'R':
			return 3;
		case 'K':
			return 4;
		case 'N':
			return 5;
		case 'Q':
			return 6;
		case 'S':
			return 7;
		case 'T':
			return 8;
		case 'A':
			return 9;
		case 'V':
			return 10;
		case 'L':
			return 11;
		case 'I':
			return 12;
		case 'F':
			return 13;
		case 'Y':
			return 14;
		case 'W':
			return 15;
		case 'H':
			return 16;
		case 'M':
			return 17;
		case 'C':
			return 18;
		case 'P':
			return 19;
	}		
}	

int main(int argc, char* argv[])
{	
	if (argc != 3) 
	{
		std::cerr<<"Error (1): Wrong number of input arguments. \n";
		exit(1);
	}
	double match_cost = 0.9;
	double gap_cost = 0.1;

	//read protein sequence
	std::fstream input_seq;
	

	input_seq.open(argv[1]);
	if(!input_seq)
	{
		std::cerr<<"Error (2): Sequence file not found. \n";
		exit(2);	
	}	
	std::string seq, line;
	while(std::getline(input_seq, line))
	{	
		seq += line;
	}
	
	//read protein profile
	initConstants();
	Profile prfl(argv[2]);
	
	//std::map<int, vector<vector<Double>>> similarity_matrix;
	
	int seq_length=seq.length();
	int prfl_length=0;	
	
	for(int i=0; i<prfl.blockCount();i++)
	{
		prfl_length+=prfl.blockSize(i);	
	}
	
	
	double similarity_matrix [seq_length+1][prfl_length+1];
	similarity_matrix[0][0]=1;
	
	for(int i=1; i<seq_length+1; i++)
	{
		similarity_matrix[i][0]=similarity_matrix[i-1][0]*gap_cost;
	}
	
	for(int j=1; j<prfl_length+1; j++)
	{
		similarity_matrix[0][j]=similarity_matrix[0][j-1]*gap_cost;
	}

	int j=0;
	Column col;
	DistanceType d;
	double help, max=-1;
	for(int t=0; t<prfl.blockCount(); t++)
	{
		Position P(t,0);
		col=prfl.getColumn(P);
		d=prfl.interBlockDist(t);
		for(int i=1; i<seq_length+1; i++)
		{
			//Rekursion 2	
			for(int k=0;k<d.r.max+2;k++)
			{	
				if(k>d.r.min && k<=i)
				{

					help=similarity_matrix[i-k][j]*col[aa2int(seq[i])]*match_cost;
					max = (help>max) ? help : max;		
								
				}
				if(k<=d.r.min && k<=i)
				{
					help=similarity_matrix[i-k][j]*pow(gap_cost,(d.r.min-k+1));
					max = (max<help) ? help:max;	
					if(k>0)
					{
						help=similarity_matrix[i-k][j]*pow(gap_cost,(d.r.min-k))*col[aa2int(seq[i])]*match_cost;
						max = (max<help) ? help:max;
					}				
				}

				similarity_matrix[i][j+1]=max;
				
			}	
		}		
		
		for(int s=2; s<prfl.blockSize(t)+1;s++)
		{
			Position P(t,s);
			col=prfl.getColumn(P);
			for(int i=1; i<seq_length+1; i++)
			{
				//Rekursion 1				
				similarity_matrix[i][j+s]= std::max({similarity_matrix[i-1][j+s-1]*col[aa2int(seq[i])]*match_cost, similarity_matrix[i-1][j+s]*gap_cost, similarity_matrix[i][j+s-1]*gap_cost});
			}	
			
		}
		j+=prfl.blockSize(t);
		max=-1;
	}
	
	
	//print similarity matrix
	/*for(int i=0;i<seq_length+1;i++)
	{
		for(int j=0;j<prfl_length+1;j++)
		{
			std::cout<<similarity_matrix[i][j];
			std::cout<<"  ";
		}
		std::cout<<endl;
	}*/
	//calculate final score
	
	//
	std::cout<<similarity_matrix[seq_length][prfl_length]<<endl;
	return 0;
}

