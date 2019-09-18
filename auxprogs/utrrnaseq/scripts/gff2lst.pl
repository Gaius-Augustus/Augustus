#!/usr/bin/env perl

# Author: Katharina J. Hoff, December 18th 2012
#
# This script converts a gff output from Maria Hartmann's tool for identification of UTRs from RNA-seq data to a list for evaluation the identified UTRs with the overlapStat.pl script.
#
# Input format example:
# chr19   makeUTR 3'-UTR  281211  281388  .       -       .       gene_id "NM_177526"; transcript_id "NM_177526";
# chr19   makeUTR 3'-UTR  617249  617274  .       -       .       gene_id "NM_005035"; transcript_id "NM_005035";
# chr19   makeUTR 3'-UTR  647600  648128  .       -       .       gene_id "NM_194460"; transcript_id "NM_194460";
# chr19   makeUTR 3'-UTR  810694  811986  .       -       .       gene_id "NM_024888"; transcript_id "NM_024888";
# chr19   makeUTR intron  811987  812152  .       -       .       gene_id "NM_024888"; transcript_id "NM_024888";
# chr19   makeUTR 3'-UTR  812153  812570  .       -       .       gene_id "NM_024888"; transcript_id "NM_024888";
# chr19   makeUTR 3'-UTR  896663  897437  .       -       .       gene_id "NM_138774"; transcript_id "NM_138774";
# chr19   makeUTR 3'-UTR  1009799 1010349 .       -       .       gene_id "NM_001033026"; transcript_id "NM_001033026";
# chr19   makeUTR 5'-UTR  1020995 1021047 .       -       .       gene_id "NM_001033026"; transcript_id "NM_001033026";
#
# Output format example:
# chr19_-_3_NM_177526     281211  281388
# chr19_-_3_NM_005035     617249  617274
# chr19_-_3_NM_194460     647600  648128
# chr19_-_3_NM_024888     810694  811986
# chr19_-_3_NM_024888     812153  812570
# chr19_-_3_NM_138774     896663  897437
# chr19_-_3_NM_001033026  1009799 1010349
# chr19_-_3_NM_001033026  1020995 1021047

if(scalar(@ARGV) != 1){
   die("\nScript that converts gff output from Maria's UTR identification routine to a lst output for comparison with OverlapStat.pl\n\nUsage: gff2lst.pl input.gff > out.lst\n\nFormat examples are given in script source code.")
}

my @t1;
my @t2;

open(GFF, "<", $ARGV[0]) or die("Could not open gff file $ARGV[0]!\n");

while(<GFF>){
    @t1 = split(/\t/);
    @t2 = split(/"/, $t1[8]);
    if($_=!m/intron/){
	print "$t1[0]_$t1[6]_";
	if($t1[2]=~m/5/){print "5_";}elsif($t1[2]=~m/3/){print "3_";}else{print STDERR "Could not determine end! $t1[2]\n"}
	print "$t2[1]\t$t1[3]\t$t1[4]\n";
    }
}
close(GFF) or die("Could not close gff file $ARGV[0]!\n");
