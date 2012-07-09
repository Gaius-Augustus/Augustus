/* Checks whether a file has been sorted by target name	 					*/
/* 														 					*/
/* Tonatiuh Pena-Centeno 								 					*/
/* Created: 11-April-2012 								 					*/
/* Last modified:  9-July-2012                           					*/

#include <stdio.h>
#include <stdlib.h> 
#include "sam.h"  
#include "bam.h"

extern bam_index_t *bam_index_core(bamFile fp);

int main(int argc, char *argv[])  
{  
    if (argc == 1) {
        fprintf(stderr, "Usage: checkSortedness <in.bam> \n");
        return 1;
    }

  	const char *fn; 
	fn = argv[1];
	const char *_fnidx;
	_fnidx = 0;
	bamFile fp;
	bam_index_t *idx; 
	if ((fp = bam_open(fn, "r")) == 0) {
		fprintf(stderr, "[checkTargetSortedness] fail to open the BAM file.\n");
		return -1;
	}
	idx = bam_index_core(fp);
	bam_close(fp);
	if(idx == 0) {
		fprintf(stderr, "[checkTargetSortedness] file is not sorted by target name and position.\n");
		return -1; 
	} 
	printf("[checkTargetSortedness]: BAM file correctly sorted.\n");

	return 0;

}  
