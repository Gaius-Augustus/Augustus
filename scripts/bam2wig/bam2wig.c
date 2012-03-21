/**
 * Author:
 *	Pierre Lindenbaum PhD
 * Contact:
 *	plindenbaum@yahoo.fr
 * WWW:
 *	http://plindenbaum.blogspot.com
 *	http://samtools.sourceforge.net/
 *	http://samtools.sourceforge.net/sam-exam.shtml
 * Reference:
 *	http://genome.ucsc.edu/goldenPath/help/wiggle.html
 * Motivation:
 *	creates a WIGGLE file from a BAM file.
 * Usage:
 *	bam2wig <bam-file> (<region>)
 */
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <errno.h>
#include "sam.h"

static const int NUM_ZERO_ACCEPTED_DEFAULT=10;

typedef struct parameter_t
	{
	/** output stream */
	FILE* out;
	/** previous chromosome seen */
	int prev_tid;
	/** previous genomic position seen */
	int prev_pos;
	/** user's start position */
	int beg;
	/** user's end position */
	int end;
	/** number of depth=0 seen */
	int count_zero;
	/** max number of depth=0 allowed */
	int pref_zero;
	/** input BAM */
	samfile_t *in;
} Param,*ParamPtr;

// callback for bam_fetch()
static int fetch_func(const bam1_t *b, void *data)
{
	bam_plbuf_t *buf = (bam_plbuf_t*)data;
	bam_plbuf_push(b, buf);
	return 0;
}
// callback for bam_plbuf_init()
static int  scan_all_genome_func(uint32_t tid, uint32_t pos, int depth, const bam_pileup1_t *pl, void *data)
	{
	ParamPtr param = (ParamPtr)data;
	if ((int)pos >= param->beg && (int)pos < param->end)
		{
		if(depth==0) /* no coverage */
			{
			param->count_zero++;
			}
		else
			{
			
			if(param->prev_tid!=tid || /* not the same chromosome */
			   param->prev_pos+1+param->count_zero!=(int)pos || /* not the expected index  */
			   param->count_zero > param->pref_zero /* too many depth=0 */
			   )
				{
				param->count_zero=0;/* reset count depth=0 */
				/* print WIGGLE header . First base of a WIG is 1*/
				fprintf(param->out,"fixedStep chrom=%s start=%d step=1 span=1\n", param->in->header->target_name[tid],pos+1);
				}
			while(param->count_zero >0)
				{
				fputs("0\n",param->out);
				param->count_zero--;
				}
			fprintf(param->out,"%d\n",depth);
			param->prev_pos=(int)pos;
			}
		param->prev_tid=(int)tid;
		}
	
	return 0;
	}

static void usage()
	{
	fprintf(stdout, "Author: Pierre Lindenbaum PHD. 2011.\n");
	fprintf(stdout, "Last compilation:%s %s\n",__DATE__,__TIME__);
	fprintf(stdout, "Usage: bam2wig (options) <aln.bam> [chr:start-end]\n");
	fprintf(stdout, "Options:\n");
	fprintf(stdout, " -z <int> number of depth=0 accepted before starting a new WIG file (default:%d).:\n",NUM_ZERO_ACCEPTED_DEFAULT);
	fprintf(stdout, " -o <filename-out> save as... (default:stdout).\n");
	fprintf(stdout, " -t print a ucsc custom track header.\n");
	}
	
int main(int argc, char *argv[])
	{
	int optind=1;
	int header=0;
	char* fileout=NULL;
	Param parameter;
	
	parameter.out=stdout;
	parameter.prev_tid=-1;
	parameter.prev_pos=-1;
	parameter.beg = 0;
	parameter.end = INT_MAX;
	parameter.pref_zero=NUM_ZERO_ACCEPTED_DEFAULT;
	parameter.count_zero=0;
	parameter.in=NULL;
	
	while(optind < argc)
		{
		if(strcmp(argv[optind],"-h")==0)
		        {
		        usage();
		        return EXIT_FAILURE;
		        }
		else if(strcmp(argv[optind],"-o")==0 && optind+1<argc)
			{
			fileout=argv[++optind];
			}
		else if(strcmp(argv[optind],"-t")==0)
			{
			header=1;
			}
		else if(strcmp(argv[optind],"-z")==0 && optind+1<argc)
		        {
		      	parameter.pref_zero=atoi(argv[++optind]);
		        if(parameter.pref_zero<0)
		        	{
		        	parameter.pref_zero=0;
		        	}
		        }
		else if(strcmp(argv[optind],"--")==0)
		        {
		        ++optind;
		        break;
		        }
		else if(argv[optind][0]=='-')
		        {
		        fprintf(stderr,"%s: unknown option '%s'\n",argv[0],argv[optind]);
		        exit(EXIT_FAILURE);
		        }
		else
		        {
		        break;
		        }
		++optind;
		}
        
        if(optind==argc)
		{
		usage();
		return EXIT_FAILURE;
		}
	parameter.in = samopen(argv[optind], "rb", 0);
	if (parameter.in == 0)
		{
		fprintf(stderr, "Cannot open BAM file \"%s\".\n", argv[optind]);
		return EXIT_FAILURE;
		}
	
	if(fileout!=NULL)
		{
		errno=0;
		parameter.out=fopen(fileout,"w");
		if(parameter.out==NULL)
			{
			fprintf(stderr, "Cannot open \"%s\" %s.\n",fileout,strerror(errno));
			return EXIT_FAILURE;
			}
		}
	if(header!=0)
		{
		fputs( "track name=\"__TRACK_NAME__\" description=\"__TRACK_DESC__\" type=\"wiggle_0\"\n",parameter.out);
		}
	if (optind+1 == argc)
		{
		sampileup(parameter.in, -1, scan_all_genome_func, &parameter);
		}	
	else  if (optind+2 == argc)
	        {
		int ref;
		bam_index_t *idx;
		bam_plbuf_t *buf;
		idx = bam_index_load(argv[optind]); // load BAM index
		if (idx == 0)
			{
			fprintf(stderr, "BAM indexed file is not available for \"%s\".\n",argv[optind]);
			return EXIT_FAILURE;
			}
		bam_parse_region(parameter.in->header, argv[optind+1], &ref,
		                 &parameter.beg, &parameter.end); // parse the region
		if (ref < 0)
			{
			fprintf(stderr, "Invalid region %s\n", argv[optind+1]);
			return EXIT_FAILURE;
			}
		buf = bam_plbuf_init( scan_all_genome_func, &parameter); // initialize pileup
		bam_fetch(parameter.in->x.bam, idx, ref, parameter.beg, parameter.end, buf, fetch_func);
		bam_plbuf_push(0, buf); // finalize pileup
		bam_index_destroy(idx);
		bam_plbuf_destroy(buf);
		}
	else
		{
		fprintf(stderr, "illegal number of arguments.\n");
		return EXIT_FAILURE;
		}
	samclose(parameter.in);
	if(fileout!=NULL)
		{
		fflush(parameter.out);
		fclose(parameter.out);
		}
	return EXIT_SUCCESS;
	}
