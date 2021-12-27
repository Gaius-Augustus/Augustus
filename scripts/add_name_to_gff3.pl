#!/usr/bin/env perl
#
# add name of feature to a gff3 file that produced by gtf2gff.pl from Augustus/scripts
#
# Katharina J. Hoff
# November 22nd 2019

use strict;
use warnings;
use Getopt::Long;


my $help = 0;
my $outfile;
my $filter;

GetOptions(
    'out=s'=>\$outfile,
    'filter!' => \$filter,
    'help!'=>\$help);

exec("perldoc $0") if ($help || !defined($outfile));

open(OUT, ">", $outfile) or die("Could not open file $outfile for writing!\n");
while(<STDIN>){
	my @t = split(/\t/);
	if( ($t[2] =~ m/gene/ ) or ($t[2] =~ m/mRNA/) ){ # EVM might only require Name for gene and mRNA
	    $t[8] =~ m/ID=([^;]+);/;
	    my $fname = $1;
	    $t[8] =~ s/\n//;
            print OUT $t[0]."\t". $t[1]."\t".$t[2]."\t".$t[3]."\t".$t[4]."\t".$t[5]."\t".$t[6]."\t".$t[7]."\t".$t[8]."Name=$fname\n";
	}else{
	    if(defined($filter)){
		if( ($t[2] =~ m/CDS/) or ($t[2] =~ m/exon/)){
	      	    $_ =~ s/;$//;
		    print OUT $_;
		}
	    }else{
		$_ =~ s/;$//;
		print OUT $_;
	    }
	}
}
close(OUT) or die("Could not close file $outfile!\n");


__END__

=pod

=head1 NAME

add_name_to_gff3.pl   -   add name of feature of a gff3 file that produced by gtf2gff.pl from Augustus/scripts 
                          to ensure compatibilty with e.g. EvidenceModeler  

=head1 SYNOPSIS

add_name_to_gff3.pl <in.gff3 --out=out.gff3
    
=head1 OPTIONS

  --out=file          gff3 output file
  --filter            print only the features gene, mRNA, CDS and exon

=head1 DESCRIPTION
    
    example input:

    chr1 AUGUSTUS  gene                     12656   14013   0.04    +   .   ID=g50;
    chr1 AUGUSTUS  mRNA                     12656   14013   0.04    +   .   ID=g50.t1;Parent=g50;
    chr1 AUGUSTUS  transcription_start_site 12656   12656   .       +   .   ID=g50.t1.tss1;Parent=g50.t1;
    chr1 AUGUSTUS  five_prime_utr           12656   12867   0.2     +   .   ID=g50.t1.5UTR1;Parent=g50.t1;
    chr1 AUGUSTUS  exon                     12656   12993   .       +   .   ID=g50.t1.exon1;Parent=g50.t1;
    chr1 AUGUSTUS  start_codon              12868   12870   .       +   0   ID=g50.t1.start1;Parent=g50.t1;
    chr1 AUGUSTUS  CDS                      12868   12993   0.8     +   0   ID=g50.t1.CDS1;Parent=g50.t1;
    chr1 AUGUSTUS  intron                   12994   13248   1       +   .   ID=g50.t1.intron1;Parent=g50.t1;
    chr1 AUGUSTUS  CDS                      13249   13479   1       +   0   ID=g50.t1.CDS2;Parent=g50.t1;
    chr1 AUGUSTUS  exon                     13249   14013   .       +   .   ID=g50.t1.exon2;Parent=g50.t1;
    chr1 AUGUSTUS  stop_codon               13477   13479   .       +   0   ID=g50.t1.stop1;Parent=g50.t1;
    chr1 AUGUSTUS  three_prime_utr          13480   14013   0.17    +   .   ID=g50.t1.3UTR1;Parent=g50.t1;
    chr1 AUGUSTUS  transcription_end_site   14013   14013   .       +   .   ID=g50.t1.tts1;Parent=g50.t1;

    example output (with --filter):

    chr1    AUGUSTUS        gene    12656   14013   0.04    +       .       ID=g50;Name=g50
    chr1    AUGUSTUS        mRNA    12656   14013   0.04    +       .       ID=g50.t1;Parent=g50;Name=g50.t1
    chr1    AUGUSTUS        exon    12656   12993   .       +       .       ID=g50.t1.exon1;Parent=g50.t1
    chr1    AUGUSTUS        CDS     12868   12993   0.8     +       0       ID=g50.t1.CDS1;Parent=g50.t1
    chr1    AUGUSTUS        CDS     13249   13479   1       +       0       ID=g50.t1.CDS2;Parent=g50.t1
    chr1    AUGUSTUS        exon    13249   14013   .       +       .       ID=g50.t1.exon2;Parent=g50.t1

=cut
