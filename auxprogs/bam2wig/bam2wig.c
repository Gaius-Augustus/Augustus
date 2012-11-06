/* 
   Creates a Wiggle file with coverage information coming from a BAM file 	
 																			
   NOTE: 
   Depending on the version of the compiler,  the call to "-lcurses" might have 
   to be replaced to "-lncurses"
 					
   Tonatiuh Pena-Centeno														
   Created: 12-June-2012 													
   Last modified:   6-November-2012												
*/


#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include "sam.h"  
#include "bam.h"

// Auxiliary data structure
typedef struct {     
	bamFile fp;      // the file handler
	bam_iter_t iter; // NULL if a region not specified
	int min_mapQ;    // mapQ (for filtering purposes but not used in this app)
} aux_t;


// Reads a BAM alignment from a BAM file.
static int read_bam(void *data, bam1_t *b) 
{// If necessary, include UNMAPPED filtering in here to avoid pileup (use min_mapQ)

	// data is a pointer to the auxiliary structure
	aux_t *aux = (aux_t*)data; 

	// Compute coverage according to the specified region, if one has been provided (i.e. bam_iter_read)
	// or compute coverage of the complete alignment otherwise (bam_read1)
	int ret = aux->iter? bam_iter_read(aux->fp, aux->iter, b) : bam_read1(aux->fp, b);

	return ret;
}

extern bam_index_t *bam_index_core(bamFile fp);

void usage()
{
  printf("\n");
  fprintf(stderr, "Usage: bam2wig [-r region] [-t trackname] <in.bam> \n");
  printf("------------------------------------------------------------\n");
  printf(" -r   Allows to specify a target region, e.g. 'chr3L:10-250'\n");	
  printf("      This option can only be used if an index file exists\n");
  printf("      See: samtools index \n");	
  printf(" -t   A string might be provided as track name\n");
  printf("\n");	
  printf("NOTE:");
  printf("File needs to be sorted by Reference ID (i.e. target name)\n");
  printf("Use 'samtools sort <in.bam>' to such effect.\n"); 
  printf("\n");
  exit(1);
}

int main(int argc, char *argv[])
{

	char *filename=NULL;
	char *trackname=NULL;

	int n, tid, beg, end, pos, *n_plp;
	const bam_pileup1_t **plp;
	char *reg = 0; // specified region
	bam_header_t *h = 0; // BAM header of the 1st input
	aux_t **data;
	bam_mplp_t mplp;
	/* bam_index_t *idx;  */

	// Parsing the command line
	while ((n = getopt(argc, argv, "r:t:")) >= 0) 
	  {
		switch (n) 
			{
			  case 'r': reg = strdup(optarg); break;   // parsing a region requires a BAM header
			  case 't': trackname = strdup(optarg); break;    // mapping quality threshold
			  default :
				usage();
			}
	  }

	if (optind == argc || (argc-optind) != 1)
	  {
		usage();
	  }


	// Initialising auxiliary data structures
	data = calloc(1, sizeof(void*)); // data[0] is array for just one BAM file
	// set the default region. left-shift "end" by appending 30 zeros (i.e. end=1073741824) 
	beg = 0; end = 1<<30; tid = -1;  

	// Opening BAM file
	char *oldTargetName = "", *newTargetName;
	filename = argv[optind];
	data[0] = calloc(1, sizeof(aux_t));
	data[0]->fp = bam_open(filename, "r"); 			// file handler of BAM
	data[0]->min_mapQ = 0;                    		// mapQ is not used by this app
	// Reading BAM header
	bam_header_t *htmp = 0;							 
	htmp = bam_header_read(data[0]->fp);         	


	// parsing region
	if (reg) 
		{ 
		  bam_parse_region(htmp, reg, &tid, &beg, &end); 
		}

	if (tid >= 0) 
	  { // if a region is specified and parsed successfully
		bam_index_t *idx = bam_index_load(argv[optind]);  // load the index
		data[0]->iter = bam_iter_query(idx, tid, beg, end); // set the iterator
		bam_index_destroy(idx); // the index is not needed any more; phase out of the memory
	  }


	// The set of multi-pileup functions are used to obtain the coverage information 
	mplp = bam_mplp_init(1, read_bam, (void**)data); // initialization
	n_plp = calloc(1, sizeof(int)); // n_plp[0] contains the number of covering reads 
	plp = calloc(1, sizeof(void*)); // plp[0] points to the array of covering reads (internal in mplp)

	// Print default trackname or use the one specified by the user
	printf("track name=%s type=wiggle_0\n", trackname==NULL? filename : trackname);


	while (bam_mplp_auto(mplp, &tid, &pos, n_plp, plp) > 0)
	  { // come to the next covered position

		// If requested region is of range, skip
		if (pos < beg || pos >= end) 
		 	continue;

		// Print the reference name as the header of each track
		newTargetName = htmp->target_name[tid];
		if (strcmp(oldTargetName, newTargetName))
		  { 
			printf("variableStep chrom=%s\n", newTargetName); 
		  }

		// Verifying whether the array of covering reads corresponds to "del" or "refskip"
		int j, m = 0;
		for (j = 0; j < n_plp[0]; ++j)
		  {
			const bam_pileup1_t *p = plp[0] + j; 	// DON'T modify plp[][] unless you really know
			if (p->is_del || p->is_refskip) 
			  {// having dels or refskips at tid:pos
				++m;
			  } 	
		  }

		// Estimates the depth, discounting "m" from the pile if: dels/refskips were found
		int coverage = n_plp[0] - m;

		// Prints position and coverage
		if (coverage > 0) 
		  {
			printf("%d %d\n", pos+1, coverage);
		  }

		// Update reference name
		oldTargetName = newTargetName;

	  } // end while


	// Freeing used memory and closing handles
	free(n_plp); 
	free(plp);
	bam_mplp_destroy(mplp);
	bam_header_destroy(h);
	bam_close(data[0]->fp);

	// Iterator is used only when a region was provided
	if (data[0]->iter) 
	  { 
		bam_iter_destroy(data[0]->iter); 
	  }
	free(data[0]); 
	free(data); 
	free(reg);

	return 0;
}
