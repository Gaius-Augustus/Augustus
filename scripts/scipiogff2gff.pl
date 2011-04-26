#!/usr/bin/perl
#
# convert a gff from scipio to a gff file in a format that
# can then be converted to genbank, e.g. with gff2gbSmallDNA.pl
# This is useful for creating a training set from a set of protein sequences.
#
# Mario Stanke, 9.06.2009

use strict;
use Getopt::Long;

my $usage .= "$0 -- convert a gff from scipio to a gff file in a format that\n";
$usage .= "can then be converted to genbank, e.g. with gff2gbSmallDNA.pl\n\n";
$usage .= "Usage: $0 --in=scipio.gff --out=gene.gff\n";
$usage .= "\n";
$usage .= "Typical usage in a pipeline:\n";
$usage .= "   scipio.pl genome.fa genes.aa > scipio.yaml\n";
$usage .= "   cat scipio.yaml | yaml2gff.pl > scipio.gff\n";
$usage .= "   $0 --in=scipio.gff --out=gene.gff\n";
$usage .= "\n";

if ($#ARGV < 1) {
    die "$usage";
}

my ($infile, $outfile);

GetOptions( 'in=s' => \$infile, 'out=s' => \$outfile);

open(INFILE, "<$infile") || die "Couldn't open $infile.\n";
open(OUTFILE, ">$outfile") || die "Could not open $outfile\n";
while(<INFILE>){
    if (/\tScipio\t/){
	s/\tprotein_match\t/\tCDS\t/;
	s/\tID=([^;]+);.*/\ttranscript_id "\1"/; 
	print OUTFILE;
    } 
}
close(OUTFILE);
close(INFILE);
