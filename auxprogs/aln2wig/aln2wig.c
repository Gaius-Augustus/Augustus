/*******************************************************************************/
/* aln2wig: converts psl and shrimp files to wig format                                   */
/* see http://genome.ucsc.edu/goldenPath/help/customTrack.html#PSL             */
/* and http://genome.ucsc.edu/goldenPath/help/wiggle.html for a description    */
/* of the formats                                                              */
/*                                                                             */
/* Autor: Ralph Krimmel                                                        */
/* email: krimmel@math.uni-goettingen.de                                       */
/*******************************************************************************/




#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>




typedef int bool;
typedef char** stringArray;
#define TRUE   (1)
#define FALSE  (0)
#define FILEBUFFER 1048576
#define SHRIMP_CHROM_LENGTH 250000000
#define START_TOKEN_SIZE 128
#define SIZE_TOKEN_SIZE 128
#define TOKEN_SIZE 512


 /* chomp */
void chomp(const char *s)
{
	char *p;
	while (NULL != s && NULL != (p = strrchr(s, '\n')))
	{
		*p = '\0';
	}
}


/*returns the tokens of a string split by the separator in an array of char **/
int strtoken(char ***target, char *in, char sep) {
	int ret = 0;
	int i = 0;
	int si = 0;
	char **out = *target;

 	while (in[i] != '\0') 
	{
		if (in[i] != sep) 
		{
   			out[ret][si++] = in[i];
 		} 
		else 
		{
   			out[ret][si] = '\0';
   			si = 0;
   			ret++; }
 		i++;
	}
	out[ret][si] = '\0';
	if (i > 0) 
	{
		ret++;
	}
 	return ret;
}

/*creates and writes the wigfile*/
void createWig(char *name, int *coveragefield, long int fieldlength, bool span)
{
	int i;
	int oldindex;
	int oldcoverage;
	if (span)
	{
		i=0;
		oldindex=0;
		oldcoverage=coveragefield[0];
			
		do
		{
			if(oldcoverage==coveragefield[i])
			{
				i++;
			}
			else
			{	if (coveragefield[i] != 0) 
				{
					printf("variableStep chrom=%s span=%i\n",name,i-oldindex);
					printf("%i %i\n",i+1,coveragefield[i]);
				}
				oldcoverage=coveragefield[i];
				oldindex=i;
				i++;
			}

		} while(i<fieldlength);
	}
	else
	{
		printf("variableStep chrom=%s\n",name);
		for (i=0;i< fieldlength; i++)
		{
		        if(coveragefield[i]!=0)
                        {
			        printf("%i %i\n",i+1,coveragefield[i]);
			}
		}
	}
}

/*checks if the input is still sorted*/
int stillSorted(char ***run,char *chkstr,int strcount)
{
	int i;
	char **target = *run;
	for (i=0; i < strcount; i++) 
	{
		if (strcmp(target[i], chkstr) == 0)
		{
		        fprintf(stderr, "Found %s!. This should never happen on sorted input data.\nPlease sort the input.\nYou can do this by using the GNU sort tool: sort -k 14,14 <filename>", *target);
			return FALSE;
		}
        }
	return TRUE;
}

int detectTokenCount(char *input, char tok)
{
        int i;
        int count=0;
        for (i=0;i < strlen(input); i++)
        {
                if (input[i] == tok)
                {
                        count ++;
                }
        }
        return count + 1;
}



int main (int argc, char **argv)
{
	int c;
	char *filename=NULL;
	char *trackname=NULL;
	/*Our filehandle*/
	FILE *handle; 
	
	char usage[]="USAGE: aln2wig -f <filename>\nInput file can be a psl file or a shrimp file\nOutput goes to STDOUT\n\nOptions:\n\t-f\t<filename>\n\t-s\tUse span notation\n\t-t\tName of the track\n";
	
	char buffer [FILEBUFFER];
	
	/* line buffer of 1024 bytes */
	char line [ 1024 ];

	/* we need information about the length of our target sequence from the first line*/
 	char firstline [TOKEN_SIZE]; 
	
	/*the list of completed chromosomes. we need this to check if the input is sorted or not*/
	char **completedChromosomes=NULL;

	int completedChromosomeCount=0;

	/* the name of the actual chromosome */	
        char name[TOKEN_SIZE];

	/*we write the wigfile on every chromosome change and free the memory after it. To determine the change, we need to save to old name. */	
	char oldname[TOKEN_SIZE];

	/*The tokens of the psl file we split by \t*/
	char **tokens;
	
	/* we save the chromosomeLength to create a coverageArray big enough to hold all positions. */
	long int chromosomeLength; 

	/* the coverage array */
	int *coverage;
	
	bool span=FALSE;

	int i,j;
	
	
	/*****************/
	/* PSL variables */
	/*****************/
	
	/*number of tokens in the file. We need this to determine if input is in psl or shrimp format*/
	int tokenCounts;

	/*The tokens of the psl file we split by ",". Namely,this is the block size*/
	char **sizeTokens;

	/* The tokens of the psl file we split by ",". This is the absolute starting Positions of the block. */
	char **startTokens; 

	/*the string we split into the  startTokens*/
	char *blockStarts; 
	
	/*the string we split into the sizeTokens*/
	char *blockSizes; 

	/* the number of blocks we have in each line */
	int blockCounts;
	

	/*********+++********/
	/* shrimp variables */
	/************+++*****/

	long int oldChromosomeLength;

	int contigStart=0;
	
	int contigEnd=0;

	opterr = 0;
	
	if ( argc == 1 ) 
	{
		printf ("%s",usage);
		exit(2);
	}
		
	while ((c = getopt (argc, argv, "t:sf:")) != -1)
	        switch (c)
		{
		case 'f':
			filename = optarg;
			break;
		case 't':
			trackname = optarg;
			break;
		case 's':
			span=TRUE;
			break;
		case '?':
			if (optopt == 'c')
		  		fprintf (stderr, "Option -%c requires an argument.\n", optopt);
			else if (isprint (optopt)) 
			{
		  		fprintf (stderr, "Unknown option `-%c'.\n", optopt);
				fprintf (stderr, "%s",usage);
		  	}
		  	else
		  	{
				fprintf (stderr,"Unknown option character `\\x%x'.\n",optopt);
				fprintf (stderr, "%s",usage);
				return 1;
		  	}
		default:
			abort ();
           	}

	handle = fopen ( filename, "r" );

	if ( handle != NULL )
	{
		/*print track name*/
		completedChromosomes = malloc(sizeof(char*));
		printf("track name=%s type=wiggle_0\n",trackname==NULL? filename : trackname);
		

		/*we want buffered input */
		setbuf(handle, buffer);

		/*read in the first line*/
		if (fgets (firstline, sizeof firstline, handle) == NULL){
		    perror ( filename );
		    printf ("Could not read first line.");
		    exit(2);
		}
		chomp(firstline);
		
		tokenCounts = detectTokenCount(firstline,'\t');
		if (tokenCounts == 21) 
		{	
			/*psl detected*/
			tokens = malloc(21*sizeof(char*));
			for (i = 0; i < 21; i++) {
				tokens[i] = malloc(TOKEN_SIZE);
				if (tokens[i] == NULL) exit(2);
			}
			
			strtoken(&tokens,firstline,'\t');
			chromosomeLength=atoi(tokens[14]);

			/* The 14th token is the name of the actual chromosome */
			strncpy(oldname,tokens[13],TOKEN_SIZE-1);
			
			coverage =  (int *) calloc(chromosomeLength,sizeof(int));
			if (coverage==NULL) exit (2);
			
			/*set the file pointer back to start because we also need to parse the first line too */
			fseek(handle, 0, SEEK_SET);
	
			while (fgets ( line, sizeof line, handle ) != NULL )
			{
				chomp(line);
				strtoken(&tokens,line,'\t');
				strncpy(name,tokens[13],TOKEN_SIZE -1);
				/*if this happens the chromosome changes*/
				if (strcmp(name,oldname)!=0) 
				{
					/*create the wigfile for the parsed chromosome*/	
					createWig(oldname,coverage,chromosomeLength,span);
					completedChromosomeCount++;
					
					/*remember the parsed chromosome to check if the input is sorted*/
					completedChromosomes = realloc(completedChromosomes,completedChromosomeCount*sizeof(char*));
					completedChromosomes[completedChromosomeCount-1] = malloc(TOKEN_SIZE*sizeof(char));
					if (completedChromosomes[completedChromosomeCount-1] == NULL )
					{
						printf("Couldnt allocate Memory, exiting\n");
						exit(1);
					}
					/*add the old chromosome to the list of completed chromosomes */
					strncpy(completedChromosomes[completedChromosomeCount-1],oldname,TOKEN_SIZE-1);
					
					free(coverage);
					/*The 15th token is the target sequence size*/
					chromosomeLength=atoi(tokens[14]);
					coverage = (int *) calloc(chromosomeLength,sizeof(int));
					if (coverage==NULL) exit (2);
					strncpy(oldname,name,TOKEN_SIZE-1);
					
					if (!stillSorted(&completedChromosomes,name,completedChromosomeCount))
					{
						exit (1);
					}
				}
			
				blockCounts = atoi(tokens[17])+1;
				blockSizes =  tokens[18];
				blockStarts = tokens[20];
				startTokens = malloc(blockCounts*sizeof(char*));
				if( startTokens == NULL) exit(2);
				for (i = 0; i < blockCounts; i++) 
				{
					startTokens[i] = malloc(START_TOKEN_SIZE);
					if (startTokens[i] == NULL) exit(2);
				}
				
				sizeTokens = malloc(blockCounts*sizeof(char*));
				if (sizeTokens == NULL) exit(2);
				for (i = 0; i < blockCounts; i++) 
				{
					sizeTokens[i] = malloc(SIZE_TOKEN_SIZE);
					if (sizeTokens[i] == NULL) exit(2);
				}
			
				strtoken((char***) &startTokens,blockStarts,',');
				strtoken((char***) &sizeTokens,blockSizes,',');
			
			
				for (i = 0; i < blockCounts-1; i++)
				{

					for (j = atoi(startTokens[i]); j < atoi(startTokens[i])+atoi(sizeTokens[i]); j++)
					{
				
						coverage[j]++;

					}
				}
			
				/*Cleaning up memory needed for parsing this line */
				
				/*freeing startTokens array */
				for (i=0; i< blockCounts; i++)
				{
					free(startTokens[i]);
					startTokens[i]=NULL; 
					
				}
				free(startTokens);
				startTokens=NULL;
				
				/*freeing sizeTokens array*/
				for (i=0; i< blockCounts; i++)
				{
					free(sizeTokens[i]);
					sizeTokens[i]=NULL; 
				}
				free(sizeTokens);
				sizeTokens=NULL;
			
			}
			createWig(oldname,coverage,chromosomeLength,span);	
			/*All done, clean up*/
			for(i=0;i<21;i++) 
			{
				free(tokens[i]);
				tokens[i]=NULL;
			}
			free(tokens);
			free(coverage);
		
			for(i=0;i<completedChromosomeCount;i++) 
			{
				free(completedChromosomes[i]);
				completedChromosomes[i]=NULL;
			}
			free(completedChromosomes);
		}
		else
		{
			/*Shrimp detected */
			tokens = malloc(10*sizeof(char*));
			oldChromosomeLength=SHRIMP_CHROM_LENGTH;
			coverage =  (int *) calloc(oldChromosomeLength,sizeof(int));
			for (i = 0; i < 11; i++) {
				tokens[i] = malloc(TOKEN_SIZE);
			}
			/*Get the name of the first chromosome*/	
			strtoken(&tokens,firstline,'\t');
			strncpy(oldname,tokens[1],TOKEN_SIZE -1);
			
			/*set the file pointer back to start because we also need to parse the first line too */
			fseek(handle, 0, SEEK_SET);
			
			while (fgets ( line, sizeof line, handle ) != NULL )
			{
				strtoken(&tokens,line,'\t');
				strncpy(name,tokens[1],TOKEN_SIZE-1);
				chromosomeLength=atoi(tokens[4]);
				if (chromosomeLength > oldChromosomeLength) 
				{
					oldChromosomeLength += (chromosomeLength-oldChromosomeLength+(SHRIMP_CHROM_LENGTH/5));
					coverage = realloc(coverage,oldChromosomeLength*sizeof(int));
				}
				if (strcmp(name,oldname)!=0) 
				{
					/*create the wigfile for the parsed chromosome*/	
					createWig(oldname,coverage,oldChromosomeLength,span);
					completedChromosomeCount++;
					/*remember the parsed chromosome to check if the input is sorted*/
					completedChromosomes = realloc(completedChromosomes,completedChromosomeCount*sizeof(char*));
					completedChromosomes[completedChromosomeCount-1] = malloc(TOKEN_SIZE*sizeof(char));
					if (completedChromosomes[completedChromosomeCount-1] == NULL )
					{
						printf("Couldnt allocate Memory, exiting\n");
						exit(1);
					}
					strncpy(completedChromosomes[completedChromosomeCount-1],oldname,TOKEN_SIZE-1);
				
					/*free(coverage); *
					chromosomeLength=atoi(tokens[4]);
					coverage = (int *) calloc(oldChromosomeLength,sizeof(int)); */
					memset ( coverage, 0x00, oldChromosomeLength); 
					strncpy(oldname,name,TOKEN_SIZE-1);
					if (coverage==NULL) exit (2);
					if (!stillSorted(&completedChromosomes,name,completedChromosomeCount))
					{
						exit (1);
					}
				}
				contigStart=atoi(tokens[3]);
				contigEnd=atoi(tokens[4]);
				for (i=contigStart;i<contigEnd+1;i++)
				{
					coverage[i]++;
				}
			}
			createWig(oldname,coverage,oldChromosomeLength,span);
			/*All done, clean up*/
			for(i=0;i<10;i++) 
			{
				free(tokens[i]);
				tokens[i]=NULL;
			}
			free(tokens);
			free(coverage);
		
			for(i=0;i<completedChromosomeCount;i++) 
			{
				free(completedChromosomes[i]);
				completedChromosomes[i]=NULL;
			}
			free(completedChromosomes);
		}

	 	fclose ( handle );
	}
	else
	{
		perror ( filename ); /* why didn't the file open? */
		printf("%s",usage);
	}
	return 0;
}
