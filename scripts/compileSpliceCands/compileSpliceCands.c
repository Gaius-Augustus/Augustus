/*******************************************************************************/
/* compileSpliceCands: 			                                       */
/*                                                                             */
/* Autor: Ralph Krimmel                                                        */
/* email: thirsty@milk-and-cookies.net                                         */
/*******************************************************************************/




#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include "list.h"


typedef int bool;
#define TRUE   (1)
#define FALSE  (0)

#define FILEBUFFER 1048576
#define LINEBUFFER 256

/*variable for debug output*/
static int debug=0;

 /* chomp a string*/
void chomp(const char *s)
{
	char *p;
	p = strchr(s, '\n'); if (p != NULL) *p = '\0';
}

/*prints out an error message when allocating memory fails and exits*/
void printOomError()
{
	printf("Could not allocate enough memory, exiting..\n");
	fflush(stdout);
	exit(1);
}

/*returns the tokens of a string split by the separator in an array of char **/
int strtoken(char ***target, char *in, char sep) {
	int ret = 0;
	int i = 0;
	int si = 0;
	char **out = *target;
	while (in[i] != '\0') 
	{
		if ((in[i] != sep)) 
		{
   			out[ret][si++] = in[i];
 		} 
		else 
		{
   			out[ret][si] = '\0';
   			si = 0;
   			ret++; 
		}
 		i++;
	}
	out[ret][si] = '\0';
	if (i > 0) 
	{
		ret++;
	}
	

 	return ret;
}

/*checks if the input is still sorted. not used at the moment*/
int stillSorted(char ***run,char *chkstr,int strcount)
{
	int i;
	char **target = *run;
	for (i=0; i < strcount; i++) 
	{
		if (strcmp(target[i], chkstr) == 0)
		{
		        fprintf(stderr, "Found %s!. This should never happen on sorted input data.\nPlease sort the input.\nYou can do this by using the GNU sort tool.\n", *target);
			return FALSE;
		}
        }
	return TRUE;
}

/*checks two strands which can be +,- and . if they are compatible.
 * compatible strands are '++', '--', '+.', '.+', '-.', '.-' and '..'  */
int isCompatibleStrand(char strand1, char strand2)
{
	if((strand1==strand2) || (strand1=='.') || (strand2 == '.'))
	{
		return TRUE;
	}
	return FALSE;
}

/*parses a double out of a string formatted like this: grp0-0-22.1373 
 * The return value would be 22.1373 then */
double getAverageCoverage(char *s)
{
	int i=0;
	char *pch;

	pch = strtok (s,"-");
	while (pch != NULL)
	{	
		i++;
		pch = strtok (NULL, "-");
		if(i==2)
		{	
			return atof(pch);
		}
		
 	}
	return 0.0;
}

/*returns the reverse complement of a char*/
char rc(char in)
{
	switch(in)
	{
		case 'a':
		  return 't';
		case 't':
		  return 'a';
		case 'c':
		  return 'g';
		case 'g':
		  return 'c';
		default:
		  return in;
	}
}

/*searches on a genome in two ranges for splice sites and prints the resulting introns in gff format
 * if there are unknown strands, it searches in both directions. This function should just be called with compatibleStrand arguments
 * (see isCompatibleStrand function) */
void findPossibleSpliceSites(char *chromosomeName,char *genomePart, long int startPosLeft, long int endPosLeft, long int startPosRight, long int endPosRight, char strandLeft, char strandRight)
{
	long int i,j;
	char resultingStrand;
	if(debug) printf("startPosLeft: %ld: | endPosLeft: %ld | startPosRight: %ld | endPosRight: %ld |strandleft: %c | strandRight: %c\n",startPosLeft,endPosLeft,startPosRight,endPosRight,strandLeft,strandRight);
		
	resultingStrand = (strandLeft==strandRight) ? strandLeft : '.';
	if((strandLeft=='+') || strandRight=='+' || (( strandLeft=='.') && (strandRight=='.' )))
	{
		for(i=startPosLeft;i<endPosLeft-1;i++)
			{
			if((genomePart[i]=='c') && (genomePart[i+1] =='t') )
			{
				for (j=startPosRight;j<endPosRight-1;j++)
				{
					if((genomePart[j] == 'a') && (genomePart[j+1] == 'g'))
					{
						printf("%s\t%s\t%s\t%ld\t%ld\t%i\t%c\t%i\t%s\n",chromosomeName,"source","intron",i,j,99,resultingStrand,1,"");
					}
				}	
			}
		}
	}
	if((strandLeft=='-') || (strandRight=='-') || ((strandLeft=='.') && (strandRight=='.' )))
	{
		for(i=endPosLeft;i>startPosLeft+1;i--)
		{
			if(((genomePart[i])=='g') && (genomePart[i-1] == 'a'))
			{
				for (j=endPosRight;j>startPosRight+1;j--)
				{
					if((genomePart[j]=='t') && (genomePart[j-1]=='c'))
					{
						printf("%s\t%s\t%s\t%ld\t%ld\t%i\t%c\t%i\t%s\n",chromosomeName,"source","intron",i,j,99,resultingStrand,1,"");
					}
				}
			}
		}
	}
}


/*searches a needed genome given by chromosomename and returns is as string*/
char *getChromosome(char *fn,char* chromosomename)
{
	/*the file handle that hold the genome file*/
	FILE *genome_handle=NULL;
	/*the buffer that holds a line while reading the genome*/
        char line[LINEBUFFER];
	/*the return string*/
        char *s=NULL;
	/*the length of the actual read string*/
        int line_length=0;
	/*when the fasta description of the next chromosome is found done is set to true*/
        int done=FALSE;
	/*true when the wanted chromosome is found. from this point on, the file is read into ram*/
	int startreading=FALSE;
	/*the string we are searching for in the genome (means: ">chromosomename")*/
	char *checkstring=NULL;
	/* the actual position we want to write the next char of the genome to*/ 
	long int position=0;
	int i=0;
	if(debug) 
	{
		fflush(stdout);
		printf("Trying to get %s from file\n",chromosomename);	
	}
	/*assemble the checkstring*/
	checkstring=calloc(strlen(chromosomename)+2,1);
	if(checkstring==NULL) printOomError();

	checkstring[0]='>';
	strncat(checkstring,chromosomename,strlen(chromosomename));
	
	/*open the genome file*/
	genome_handle = fopen ( fn, "r" );
        if (genome_handle == NULL)
        {
                printf("Could not open genome file\n");
        }
        
	/*get enough memory to hold the chromosome*/
	s = malloc(3500000000);
	if(s == NULL) printOomError();
	
	/*read in the chromosome*/
	while ((fgets ( line, sizeof line, genome_handle ) != NULL) && !(done))
	{
		chomp(line);
		if(strcmp(line,checkstring)==0)	
		{
			startreading=TRUE;
			if(debug) 
			{
				printf("Found %s, reading it into ram\n",chromosomename);
				fflush(stdout);
			}
			continue;
		}	
		if(startreading)
		{

			line_length=strlen(line);
			if(line[0]=='>')
			{
				done=TRUE;
			}
			else
			{
				for(i=0;i<line_length;i++)
				{
					s[i+position]=tolower(line[i]);
				}
				position+=line_length;
			}
		}
	}
	free(checkstring);
	fclose(genome_handle);
	if(debug) printf("Successfully loaded chromosome %s into memory\n",chromosomename);
	fflush(stdout);
	return s;
}

/*parses a line and returns a filled Data struct */
Data parseLine(char *ln)
{
	Data d = malloc(sizeof(struct LineData));

	char **tokens=NULL;
	int i;	
	tokens = malloc(9*sizeof(char*));
	if (tokens ==NULL) 
	{
		printf("Couln't allocate memory, exiting\n");
		exit(1);
	}
			for (i = 0; i < 9; i++) 
	{
		tokens[i] = malloc(LINEBUFFER);
		if(tokens[i] == NULL)
		{	
			printf("Couln't allocate memory, exiting\n");
			exit(1);
		}
	}

	strtoken(&tokens,ln,'\t');
	strncpy(d->name,tokens[0],LINEBUFFER);
	d->startPos=atoi(tokens[3]);
	d->endPos=atoi(tokens[4]);
	d->strand=tokens[6][0];
	d->averageCoverage=getAverageCoverage(tokens[8]);
	
	for(i=0;i<9;i++) 
	{
		free(tokens[i]);
		tokens[i]=NULL;
	}
	free(tokens); 
	tokens=NULL;
return d;
}

/*parses the linked list of Data structs */
int parseList(List *L, char *actualChromosome,long int chromosomeLength,char* oldname,int maxSpliceSiteDiff,float threshold,int maxIntronLength)
{
	List linkedList = *L;
	Position P;

	Data left;
	Data right;
	int isCandidate;
	signed int leftBorder;
	long int rightBorder;
	

	P=Header(linkedList);
	if(debug) printf("parsing list\n"),fflush(stdout);
	while  (!IsLast(P,linkedList))
	{
		P=Advance(P);
		if(P == First(linkedList))
		{
			left=Retrieve(P);	
		}
		else
		{
			right=Retrieve(P);
			isCandidate = ((right->averageCoverage>=left->averageCoverage*threshold) || 
				       (left->averageCoverage*((1/threshold))<right->averageCoverage));
			if(isCandidate && (maxIntronLength > (right->endPos-left->startPos)))
			{	
				if(isCompatibleStrand(left->strand,right->strand))
				{
					leftBorder=left->startPos-maxSpliceSiteDiff;
					
					rightBorder=right->endPos+maxSpliceSiteDiff;
					if ((leftBorder>0) && (rightBorder<chromosomeLength))
					{
						findPossibleSpliceSites(left->name,actualChromosome,left->startPos-maxSpliceSiteDiff,left->endPos+maxSpliceSiteDiff,right->startPos-maxSpliceSiteDiff,right->startPos+maxSpliceSiteDiff,left->strand,right->strand); 
					}
					else
					{
						if((left->startPos>0) && (right->endPos<chromosomeLength))
						{
							findPossibleSpliceSites(left->name,actualChromosome,left->startPos,left->endPos,right->startPos,right->endPos,left->strand,right->strand);
						}
					}

				}
			}
		}
	}
	return 0;
}


void printList(const List L) 
{
	Position P = Header(L);
	if ( IsEmpty(L) ) 
	{
		printf("[Empty list]\n");
	} 
	else 
	{
		do 
		{
      			P = Advance(P);
			printf("name:%s| avgc:%f | start: %ld | end: %ld\n", Retrieve(P)->name,Retrieve(P)->averageCoverage,Retrieve(P)->startPos,Retrieve(P)->endPos);
    		} while (!IsLast(P, L));
	}
	fflush(stdout);
}



int main (int argc, char **argv)
{
	char *gff_filename=NULL;
	char *genome_filename=NULL;
	/*char *trackname=NULL; */
	
	FILE *gff_handle=NULL; 

	char usage[]="\ncompileSpliceCands: Find introns parsing a set of exon candidates in gff format (output of curve2hints)  USAGE: compileSpliceCands -f <splice candidates filename> -g <genome filename>\n\nOption:\t\tArgument:\tDescription:\n\t-f\t<filename>\tThe potential splice sites in gff format\n\t-g\t<filename>\tThe genome file in fasta format)\n\t-c\t<integer>\tDefines the number of potential splice site every potential splice site itself is compared with\n";
	char usage1[]="\t-t\t<float>\t\thas to be between 0 and 1 and  defines how much the average coverage may differ\n\t-m\t<integer>\tdefines how many bases arround the splice site should be checked.\n\t-d\t<no argument>\tenables debugging output\n\t-i\t<integer>\tThe maximum length an Intron can be. Default 500000(human genome).\n";
	
	char buffer[FILEBUFFER];
	
	char line[LINEBUFFER];

	/* the name of the last chromosome */	
        char oldname[LINEBUFFER];
	
	int maxSpliceSiteDiff=10;
	
	float threshold=0.6;

	int checkWidth=3;

	int i=0;
	
	int c=0;

	char *actualChromosome=NULL;

	int listCount=0;

	int sameChromosome=TRUE;
	
	long int chromosomeLength;

	int maxIntronLength=500000;

	List L=NULL;
	
	Position P=NULL; 

	Data tmpData;

	opterr = 0;
	
	if ( argc < 5  ) 
	{
		printf ("%s%s",usage,usage1);
		exit(2);
	}

	/*Parse command line arguments*/
	while ((c = getopt (argc, argv, "f:g:c:t:m:i:d")) != -1)
	        switch (c)
		{
		case 'f':
			gff_filename = optarg;
			break;
		case 'g':
			genome_filename = optarg;
			break;
		case 'c':
			checkWidth = atoi(optarg);
			break;
		case 'm':
			maxSpliceSiteDiff = atoi(optarg);
			break;
		case 'd':
			debug=1;
			break;
		case 'i':
			maxIntronLength = atoi(optarg);
			break;
		case 't':
			threshold = atof(optarg);
			break;
		case '?':
			if (optopt == 'c')
		  		fprintf (stderr, "Option -%c requires an argument.\n", optopt);
			else if (isprint (optopt)) 
			{
		  		fprintf (stderr, "Unknown option `-%c'.\n", optopt);
				fprintf (stderr, "%s%s",usage,usage1);
		  	}
		  	else
		  	{
				fprintf (stderr,"Unknown option character `\\x%x'.\n",optopt);
				fprintf (stderr, "%s%s",usage,usage1);
				return 1;
		  	}
		default:
			fprintf(stderr, "Dont know what to do, exiting\n");
			abort ();
           	}
	
	
	/* invalid command line arguments */
	if(checkWidth == 0) 
	{
		checkWidth=4;
	}
	if(maxSpliceSiteDiff == 0)
	{
		maxSpliceSiteDiff=10;
	}
	if(threshold == 0.0)
	{
		threshold = 0.6;
	}
	if(maxIntronLength == 0)
	{

	}	maxIntronLength = 500000;

	/*cause we start counting from zero ;) */
	checkWidth++;
	gff_handle = fopen ( gff_filename, "r" );
	if ( gff_handle != NULL )
	{
		/*create an empty linked list and set the position variable to the header*/	
		 L = MakeEmpty(NULL);
		 P = Header(L);
		
		/*we want buffered input */
		setbuf(gff_handle, buffer);

		i=0;
		/*initial fill of the linked list*/
		while((fgets ( line, sizeof line, gff_handle ) != NULL)) 
		{
			chomp(line);
			tmpData=parseLine(line);
			if(i==0)
			{
				actualChromosome=getChromosome(genome_filename,tmpData->name);
				chromosomeLength=strlen(actualChromosome);
				strncpy(oldname,tmpData->name,LINEBUFFER);
				i++;
			}
			sameChromosome=((strcmp(tmpData->name,oldname)==0));
			if(debug) 
			{
				printf("tmpData->name=%s tempData->averageCovergage:%f| oldname:%s\n| listcount: %i | checkWidth: %i  | sameChromosome: %i\n",tmpData->name,tmpData->averageCoverage,oldname,listCount,checkWidth,sameChromosome);
				printList(L);
				printf("\n");
				fflush(stdout);
			}
			if((listCount<checkWidth && sameChromosome))
			{
				Insert(tmpData,L,P);
				P=Advance(P);
				listCount++;
			}
			else
			{
				if((listCount==checkWidth && sameChromosome))
				{
					parseList(&L,actualChromosome,chromosomeLength,oldname,maxSpliceSiteDiff,threshold,maxIntronLength); 
					Delete(Retrieve(First(L)),L);
					Insert(tmpData,L,P);
					P=Advance(P);
				}

				else
				{
					parseList(&L,actualChromosome,chromosomeLength,oldname,maxSpliceSiteDiff,threshold,maxIntronLength); 
					free(actualChromosome);
					actualChromosome=getChromosome(genome_filename,tmpData->name);
					chromosomeLength=strlen(actualChromosome);
					strcpy(oldname,tmpData->name);
					L=MakeEmpty(L);
					P=Header(L);
					Insert(tmpData,L,P);
					P=Advance(P);
					listCount=1;
				}
			}
		}

		if(debug) printList(L);	
		parseList(&L,actualChromosome,chromosomeLength,oldname,maxSpliceSiteDiff,threshold,maxIntronLength); 
		printf("\n");
	 	fclose ( gff_handle );
		DeleteList(L); 
	}
	else
	{
		perror ( gff_filename ); /* why didn't the file open? */
		printf("%s%s",usage,usage1);
	}
	free(actualChromosome);
	return 0;
}
