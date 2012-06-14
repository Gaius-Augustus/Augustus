/* This program demonstrates how to generate pileup from multiple BAMs
 * simutaneously, to achieve random access and to use the BED interface.
 * To compile this program separately, you may:
 *
 *   gcc -g -O2 -Wall -o bam2depth -D_MAIN_BAM2DEPTH bam2depth.c -L. -lbam -lz
 */

// Some other notes
/* 																			*/
/* From BAM.H, the structure BAM_PILEUP1_T is defined as follows: 			*/
/* 																			*/
/* 		typedef struct { 													*/
/*   		bam1_t *b; 														*/
/*   		int32_t qpos; 													*/
/*   		int indel, level; 												*/
/*   		uint32_t is_del:1, is_head:1, is_tail:1, is_refskip:1, aux:28; 	*/
/* 		} bam_pileup1_t; 													*/
/* 						   													*/
/* Created: 12-June-2012 */
/* Last modified: 13-June-2012 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include "bam.h"

typedef struct {     // auxiliary data structure
	bamFile fp;      // the file handler
	bam_iter_t iter; // NULL if a region not specified
	int min_mapQ;    // mapQ filter
} aux_t;

void *bed_read(const char *fn); // read a BED or position list file
void bed_destroy(void *_h);     // destroy the BED data structure
int bed_overlap(const void *_h, const char *chr, int beg, int end); // test if chr:beg-end overlaps

// This function reads a BAM alignment from one BAM file.
static int read_bam(void *data, bam1_t *b) // read level filters better go here to avoid pileup
{
	aux_t *aux = (aux_t*)data; // data in fact is a pointer to an auxiliary structure
	int ret = aux->iter? bam_iter_read(aux->fp, aux->iter, b) : bam_read1(aux->fp, b);
	if ((int)b->core.qual < aux->min_mapQ) b->core.flag |= BAM_FUNMAP;
	return ret;
}

int main(int argc, char *argv[])
{

	char *filename=NULL;
	char *trackname=NULL;

	int n, tid, beg, end, pos, *n_plp, baseQ = 0, mapQ = 0;
	const bam_pileup1_t **plp;
	char *reg = 0; // specified region
	void *bed = 0; // BED data structure
	bam_header_t *h = 0; // BAM header of the 1st input
	aux_t **data;
	bam_mplp_t mplp;

	// parse the command line
	while ((n = getopt(argc, argv, "r:b:q:Q:")) >= 0) {
		switch (n) {
			case 'r': reg = strdup(optarg); break;   // parsing a region requires a BAM header
			case 'b': bed = bed_read(optarg); break; // BED or position list file can be parsed now
			case 'q': baseQ = atoi(optarg); break;   // base quality threshold
			case 'Q': mapQ = atoi(optarg); break;    // mapping quality threshold
		}
	}
	if (optind == argc) {
		fprintf(stderr, "Usage: bam2wig [-r reg] [-q baseQthres] [-Q mapQthres] [-b in.bed] <in1.bam> [...]\n");
		return 1;
	}

	// initialize the auxiliary data structures
	
	if ((argc - optind)!=1) 
	  { // number of input BAMs should be equal to 1
		printf("Just one BAM file as input is admitted!!!\n"); 
		exit(1);
	  }

	data = calloc(1, sizeof(void*)); // data[i] for the i-th input; data is an array
	beg = 0; end = 1<<30; tid = -1;  // set the default region

	// Open BAM file, including header
	char *oldName, *newName;
	bam_header_t *htmp = 0;							// keep the header of the 1st BAM */
	data[0] = calloc(1, sizeof(aux_t));
	filename = argv[optind];
	data[0]->fp = bam_open(filename, "r"); 		// open BAM
	data[0]->min_mapQ = mapQ;                    	// set the mapQ filter
	htmp = bam_header_read(data[0]->fp);         	// read the BAM header

	if (reg)
	  { // parse the region
		bam_parse_region(htmp, reg, &tid, &beg, &end);
	  }

	if (tid >= 0) 
	  { // if a region is specified and parsed successfully
		bam_index_t *idx = bam_index_load(argv[optind]);  // load the index
		data[0]->iter = bam_iter_query(idx, tid, beg, end); // set the iterator
		bam_index_destroy(idx); // the index is not needed any more; phase out of the memory
	  }

	// the core multi-pileup loop
	mplp = bam_mplp_init(1, read_bam, (void**)data); // initialization
	n_plp = calloc(1, sizeof(int)); // n_plp[i] is the number of covering reads from the i-th BAM
	plp = calloc(1, sizeof(void*)); // plp[i] points to the array of covering reads (internal in mplp)

	// Print track name
	printf("track name=%s type=wiggle_0\n", trackname==NULL? filename : trackname);

	oldName = "";

	while (bam_mplp_auto(mplp, &tid, &pos, n_plp, plp) > 0)
	  { // come to the next covered position

		if (pos < beg || pos >= end) continue; // out of range; skip
		if (bed && bed_overlap(bed, htmp->target_name[tid], pos, pos + 1) == 0) continue; // not in BED; skip

		// prints reference name (I will not require it all the time)
		newName = htmp->target_name[tid];
		/* fputs(newName, stdout);  */
		if (strcmp(oldName, newName))
		  {
			printf("variableStep chrom=%s\n", newName);
		  }

		// Sweeping through all BAM files; just one in my case
		// base level filters have to go here
		int m = 0;
		int j;
		for (j = 0; j < n_plp[0]; ++j)
		  {
			const bam_pileup1_t *p = plp[0] + j; 	// DON'T modify plp[][] unless you really know
			if (p->is_del || p->is_refskip) {++m;} 	// having dels or refskips at tid:pos
			/* 	else if (bam1_qual(p->b)[p->qpos] < baseQ) {++m;} // low base quality */
		  }

		// Prints out the depth, discounting "m" from the pile if: dels/refskips or low quality
		// alignments were present
		int coverage = n_plp[0] - m;

		// prints position and coverage
		if (coverage > 0) printf("%d %d\n", pos+1, coverage);

		// Updating reference name
		oldName = newName;
	  } // end while


	// Freeing used memory and closing handles
	free(n_plp); 
	free(plp);
	bam_mplp_destroy(mplp);

	bam_header_destroy(h);
	bam_close(data[0]->fp);
	if (data[0]->iter) 
	  { 
		bam_iter_destroy(data[0]->iter); 
	  }
	free(data[0]); 
	free(data); 
	free(reg);
	if (bed) 
	  { 
		bed_destroy(bed); 
	  }
	return 0;
}
