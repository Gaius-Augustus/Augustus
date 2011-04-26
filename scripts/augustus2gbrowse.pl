#!/usr/bin/perl
# convert AUGUSTUS output to Gbrowse format GFF file
# Mario Stanke

@input = <STDIN>;


foreach(@input){
    next unless (/\tAUGUSTUS\t/);
    if (/\tgene\t/){
	s/\t(\S+)$/\tGene \1/;
    }
	
    s/\ttranscript\t/\tmRNA\t/;
    if (/\tmRNA\t/){
	s/\t(\S+)\.(\w+)$/\tmRNA \1.\2 ; Gene \1/;
    } else {
	s/transcript_id "(.*)"; gene_id "(.*)";/mRNA \1/;
    }
    print;
}
