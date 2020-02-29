#!/usr/bin/env perl

# Authors: Katharina Hoff
# Contact: katharina.hoff@uni-greifswald.de
# This script converts the output of startAlign.pl to gtf format

use strict;
use warnings;

my $usage = "Usage\n\ngth2gtf.pl gth.gff3 gth.gtf\n\n";

if(scalar @ARGV != 2){
    print $usage;
    exit(1);
}

my $v = 4;

my $align = $ARGV[0];
my $out = $ARGV[1];

print STDOUT "\# " . (localtime) . ": Converting GenomeThreader file $align "
    . "to gtf format\n" if ($v > 2);

open( GTH,    "<", $align ) or die(
    "ERROR in file " . __FILE__ ." at line ". __LINE__
    . "\nCould not open file $align!\n");
open( GTHGTF, ">", $out )   or die(
    "ERROR in file " . __FILE__ ." at line ". __LINE__
    . "\nCould not open file $out!\n");

my $geneId;

# GTH may output alternative transcripts; we don't want to have any
# alternatives in training gene set, only print the first of any occuring
# alternatives
my %seen;
while (<GTH>) {
    chomp;
    my @gtfLine = split(/\t/);
    if (m/\tgene\t/) {
        my @idCol = split( /=/, $gtfLine[8] );
        $geneId = $idCol[1];
    }
    elsif (m/\tCDS\t/) {
        my @gtfLineLastCol      = split( /;/, $gtfLine[8] );
        my @gtfLineLastColField = split( /=/, $gtfLineLastCol[1] );
        if (not( defined( $seen{ "$gtfLine[0]" . "_" . $geneId . "_" } ) )
            )
        {
            $seen{ "$gtfLine[0]" . "_" . $geneId . "_" }
                = "$gtfLine[0]" . "_"
                . $geneId . "_"
                . $gtfLineLastColField[1];
        }
        if ( $seen{ "$gtfLine[0]" . "_" . $geneId . "_" } eq "$gtfLine[0]"
            . "_"
            . $geneId . "_"
            . $gtfLineLastColField[1] )
        {
            print GTHGTF "$gtfLine[0]\t$gtfLine[1]\t$gtfLine[2]\t"
                        . "$gtfLine[3]\t$gtfLine[4]\t$gtfLine[5]\t"
                        . "$gtfLine[6]\t$gtfLine[7]\tgene_id \""
                        . "$gtfLine[0]_g_" .$geneId . "_"
                        . $gtfLineLastColField[1] . "\"; transcript_id "
                        . "\"$gtfLine[0]_t" . "_" . $geneId . "_"
                        . $gtfLineLastColField[1] . "\";\n";
            print GTHGTF "$gtfLine[0]\t$gtfLine[1]\texon\t$gtfLine[3]\t"
                        . "$gtfLine[4]\t$gtfLine[5]\t$gtfLine[6]\t"
                        . "$gtfLine[7]\tgene_id \"$gtfLine[0]_g" . "_"
                        . $geneId . "_"
                        . $gtfLineLastColField[1] . "\"; transcript_id \""
                        . "$gtfLine[0]_t" . "_" . $geneId . "_"
                        . $gtfLineLastColField[1] . "\";\n";
        }
    }
}
close(GTHGTF) or die("ERROR in file " . __FILE__ ." at line ". __LINE__
    . "\nCould not close file $out!\n");
close(GTH)    or die("ERROR in file " . __FILE__ ." at line ". __LINE__
    . "\nCould not close file $align!\n");
