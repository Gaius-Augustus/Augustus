#!/usr/bin/env perl

# Authors: Katharina Hoff & Sophie Kersting
# Contact: katharina.hoff@uni-greifswald.de
# This script computes the average length of genes in a gtf file

use strict;
use warnings;
use List::Util qw[min max];
use POSIX qw/floor/;

my $usage = "Usage\n\ncomputeFlankingRegion.pl input.gtf\n\n";

my $totalgenelength=0;
my $numberofgenes=0;
my $lastgenename=0;

if(scalar @ARGV<1){
    print $usage;
    exit(1);
}

my $printedWarn = 0;

my $filename = $ARGV[0];

my %gene;

open(GEN, "<", $filename) or die("Could not open file ".$filename."\n");
while(<GEN>){
	if($_ =~ m/transcript_id/){
    	if($_=~m/^(.*)\t.*\tCDS\t(\d*)\t(\d*)\t.*\t(\+|-)\t(0|1|2)\t.*transcript_id "(.*)"/){
    		if(not(defined($gene{$1.$6}{'start'}))) {
    			$gene{$1.$6}{'start'} = min($2, $3);
    		}elsif($gene{$1.$6}{'start'} > min($2, $3)){
    			$gene{$1.$6}{'start'} = min($2, $3);
    		}
    		if(not(defined($gene{$1.$6}{'stop'}))) {
    			$gene{$1.$6}{'stop'} = max($2, $3);
    		}elsif($gene{$1.$6}{'stop'} < max($2, $3)){
    			$gene{$1.$6}{'stop'} = max($2, $3);
    		}
		}
	}elsif($_=~m/\tCDS\t/){
		if($printedWarn==0){
			print STDERR "Warning: transcript_id was not detected in the last column. Will assume that the entire last column serves as transcript identifier!\n";
			$printedWarn = 1;
		}
		if($_=~m/^(.*)\t.*\tCDS\t(\d*)\t(\d*)\t.*\t(\+|-)\t(0|1|2|\.)\t(.*)/){
    		if(not(defined($gene{$1.$6}{'start'}))) {
    			$gene{$1.$6}{'start'} = min($2, $3);
    		}elsif($gene{$1.$6}{'start'} > min($2, $3)){
    			$gene{$1.$6}{'start'} = min($2, $3);
    		}
    		if(not(defined($gene{$1.$6}{'stop'}))) {
    			$gene{$1.$6}{'stop'} = max($2, $3);
    		}elsif($gene{$1.$6}{'stop'} < max($2, $3)){
    			$gene{$1.$6}{'stop'} = max($2, $3);
    		}
		}
	}
}
close(GEN) or die("Could not close file ".$filename."\n");

foreach(keys %gene){
	$numberofgenes++;
	$totalgenelength += $gene{$_}{'stop'} - $gene{$_}{'start'} + 1;
}

my $average_length = $totalgenelength/$numberofgenes;
print "\nTotal length gene length (including introns): ".$totalgenelength.". Number of genes: ".$numberofgenes.". Average Length: ".$average_length."\n\n";
# calculate the flanking DNA value
my $flanking_DNA = min(floor($average_length/2), 10000);
print "The flanking_DNA value is: ".$flanking_DNA." (the Minimum of 10 000 and ".floor($average_length/2).")\n";
