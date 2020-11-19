/* 
    Creates a Wiggle file with coverage information taken from a BAM file 	
										
*/


#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include "bgzf.h"
#include "sam.h"

// Auxiliary data structure
typedef struct {
	htsFile *hts_fp; // the file handler
	BGZF *bgzf_fp;   // the file handler
	hts_itr_t *iter; // NULL if a region not specified
	int min_mapQ;    // mapQ (for filtering purposes but not used in this app)
} aux_t;


// Reads a BAM alignment from a BAM file.
static int read_bam(void *data, bam1_t *b) 
{// If necessary, include UNMAPPED filtering in here to avoid pileup (use min_mapQ)

	// data is a pointer to the auxiliary structure
	aux_t *aux = (aux_t*)data; 

	// Compute coverage according to the specified region, if one has been provided (i.e. bam_iter_read)
	// or compute coverage of the complete alignment otherwise (bam_read1)
	int ret = aux->iter? hts_itr_next(aux->bgzf_fp, aux->iter, b, aux->hts_fp) : bam_read1(aux->bgzf_fp, b);

	return ret;
}

void usage()
{
  printf("\n");
  fprintf(stderr, "Usage: bam2wig [-r region] [-t trackname] <in.bam> \n");
  printf("-----------------------------------------------------------------\n");
  printf(" -r   Allows user to specify a target region, e.g. 'chr3L:10-250'\n");
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

	int n, tid, *n_plp;
	hts_pos_t beg, end, pos;
	const bam_pileup1_t **plp;
	char *reg = 0; // specified region
	aux_t **data;
	bam_mplp_t mplp;

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


	// Initializing auxiliary data structures
	data = calloc(1, sizeof(void*)); // data[0] is array for just one BAM file
	// set the default region to the maximum value of hts_pos_t
	beg = 0; end = HTS_POS_MAX; tid = -1;

	// Opening BAM file
	char *oldTargetName = "", *newTargetName;
	filename = argv[optind];
	data[0] = calloc(1, sizeof(aux_t));
	data[0]->hts_fp = hts_open(filename, "r");              // file handler of BAM
	if (data[0]->hts_fp == NULL) {
		fprintf(stderr, "Failed to open file \"%s\" : No such file or directory or not a bam file.\n", filename);
		exit(1);
	}
	if (data[0]->hts_fp->format.format != bam) {
		fprintf(stderr, "File \"%s\" is not in bam file format.\n", filename);
		exit(1);
	}
	data[0]->bgzf_fp = data[0]->hts_fp->fp.bgzf;            // file handler of BAM
	data[0]->min_mapQ = 0;                    		// mapQ is not used by this app
	// Reading BAM header
	sam_hdr_t *htmp = bam_hdr_read(data[0]->bgzf_fp);

	// parsing region
	if (reg) 
		{ 
			hts_parse_region(reg, &tid, &beg, &end, (hts_name2id_f)bam_name2id, htmp, 0);
		}

	if (tid >= 0) 
	  { // if a region is specified and parsed successfully
		hts_idx_t *idx = hts_idx_load(argv[optind], HTS_FMT_BAI);  // load the index

		if (idx == NULL)
		  {
		    fprintf(stderr, "Missing indexed BAM file!\n");
		    fprintf(stderr, "See: samtools index\n");
		    fprintf(stderr, "Do `samtools index <in.bam>` and an index file with extension \".bai\" will be generated.\n");
		    exit(1);
		  }

		data[0]->iter = sam_itr_queryi(idx, tid, beg, end); // set the iterator
		hts_idx_destroy(idx); // the index is not needed any more; phase out of the memory
	  }


	// The set of multi-pileup functions are used to obtain the coverage information 
	mplp = bam_mplp_init(1, read_bam, (void**)data); // initialization
	n_plp = calloc(1, sizeof(int)); // n_plp[0] contains the number of covering reads 
	plp = calloc(1, sizeof(void*)); // plp[0] points to the array of covering reads (internal in mplp)

	// Print default trackname or use the one specified by the user
	printf("track name=%s type=wiggle_0\n", trackname==NULL? filename : trackname);

	int exitCode = 0; // less than zero if e.g. bam file is sorted
	while ((exitCode = bam_mplp64_auto(mplp, &tid, &pos, n_plp, plp)) > 0)
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
			printf("%ld %d\n", pos+1, coverage);
		  }

		// Update reference name
		oldTargetName = newTargetName;

	  } // end while


	// Freeing used memory and closing handles
	free(n_plp); 
	free(plp);
	bam_mplp_destroy(mplp);
	sam_hdr_destroy(htmp);
	hts_close(data[0]->hts_fp);

	// Iterator is used only when a region was provided
	if (data[0]->iter) 
	  { 
		hts_itr_destroy(data[0]->iter); 
	  }
	free(data[0]); 
	free(data); 
	free(reg);

	return exitCode == 0 ? 0 : 1;
}
