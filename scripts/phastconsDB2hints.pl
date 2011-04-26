#!/usr/bin/perl
#
# convert phastcons track to a file with hints to AUGUSTUS
# 
# input file format:
# chromosome  begin  end  name  score
# e.g.
# chr1    21102591        21102633        lod=35  321
#
# Each phastcons interval gives rise to 5 CDSpart hints of different regins. 
# If the score is above a threshold and the CDSpart interval is not a the ends
# then the grade is 2, otherwise the grade is 1. This can be used in the 
# extrinsic.cfg file to consider the internal parts more reliable than the regions at the 
# margins of the conserved region.
#
# Mario Stanke, mstanke@gwdg.de

use strict;
use Getopt::Long;
my ($usage, $priority, $grp);
my ($chrom, $name, $chromStart, $blockCount, $blockSizes, $chromStarts, $score, $len, $grade);
my $source = 'X';
my $cutends = 5;
my $marginA = 10;
my $marginB = 10;
my $minscore = 260;
my $thresh = 400;
my $gradeoffset = 1; # to use together with exoniphy (grade 0)

#
# raw region:   --------------------------------------------------------------------------------------
#               cutends | marginA | marginB |             rest            | marginB | marginA| cutends    
#                       |  CDSpart| CDSpart |           CDSpart           | CDSpart | CDSpart|
# if raw region is shorter than 2(cutends+marginA+marginB), scale down to all 1/5
#
#default values
$priority = 4;

$usage .= "$0 -- convert phastcons conserved regions to hints\n";
$usage .= "Usage: $0 [--priority=p --minscore=a --source=s --cutends=n]\n";
$usage .= "reads from stdin writes to stdout.\n";
$usage .= "takes input in tab separated format:\n";
$usage .= "chrom\tstart\tend\tname\tscore\n";

GetOptions('priority:i'=>\$priority,
	   'source:s' =>\$source,
	   'minscore:i' =>\$minscore,
	   'cutends:i' =>\$cutends);

my ($from, $to, $newfrom, $newto);
while (<>) {
    chomp;
    my @f = split /\t/, $_;
    my @p;
    next unless (@f > 4);

    $chrom = $f[0];
    $name = $f[3];
    $score = $f[4];
    $from = $f[1];
    $to = $f[2];
    $len = $to - $from;

    next unless ($score >= $minscore);
    next unless ($len-2*$cutends>0);
    
    if ($len > 2*($cutends+$marginA+$marginB)){
	@p = ($from+$cutends, $from+$cutends+$marginA, $from+$cutends+$marginA+$marginB, $to - $marginB-$marginA-$cutends, $to-$marginA-$cutends, $to-$cutends);
    } else {
	my $fifth = int(($len-2*$cutends)/5);
	@p = ($from+$cutends, $from+$cutends+$fifth, $from+$cutends+2*$fifth, $to-$cutends-2*$fifth, $to-$cutends-$fifth, $to-$cutends);
    }
    for (my $i=0; $i<@p-1; $i++){
	$newfrom = $p[$i];
	$newto = $p[$i+1]-1;
	if ($i != 0 && $i != 4 && $score > $thresh) {
	    $grade = 2;
	} else {
	    $grade = 1;
	}
	print "$chrom\tphast\tCDSpart\t$newfrom\t$newto\t$grade\t\.\t.\tsrc=$source;pri=$priority\n";
    }
}
    
    
