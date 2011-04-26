#!/usr/bin/perl
#
#  convert the retroGenes track to a file with hints to AUGUSTUS
#

use strict;
use Getopt::Long;
my ($usage, $priority, $grp);
my ($chrom, $name, $chromStart, $blockCount, $blockSizes, $chromStarts, $score);
my $source = 'R';
my $max_small_gap_len = 30;
my $minscore=650;

#default values
$priority = 4;

$usage .= "$0 -- convert retroposed genes output to hints\n";
$usage .= "Usage: $0 [--priority=p --minscore=a --maxSmallGapLen=m --source=s]\n";
$usage .= "reads from stdin writes to stdout.\n";
$usage .= "takes input in tab separated format:\n";
$usage .= "name\tchrom\tstrand\tcdsStart\tcdsEnd\n";
$usage .= "You get this format on the UCSC system by typing e.g.\n";
$usage .= "echo \"select name, chrom, strand, cdsStart, cdsEnd from exoniphy order by chrom, cdsStart;\" | hgsql hg18 | grep -v cdsStart\n";

GetOptions('priority:i'=>\$priority,
	   'source:s' =>\$source,
	   'minscore:i' =>\$minscore,
           'maxSmallGapLen:i'=>\$max_small_gap_len);

my ($from, $to, $newfrom, $newto);
while (<>) {
    chomp;
    my @f = split /\t/, $_;
    next unless (@f > 6);

    $chrom = $f[0];
    $name = $f[1];
    $score = $f[2];
    $chromStart = $f[3];
    $blockCount = $f[4];
    $blockSizes = $f[5];
    $chromStarts = $f[6];
    
    next unless ($score >= $minscore);
    
    my @sizes = split /,/, $blockSizes;
    my @starts = split /,/, $chromStarts;
    $from = $to = -1000;
    for (my $i=0; $i<$blockCount; $i++) {
	$newfrom = $chromStart + $starts[$i] + 1;
	$newto = $chromStart + $starts[$i] + 1 + $sizes[$i]-1;
	if ($to >= 0 && $to < $newfrom - $max_small_gap_len - 1) {
	    print "$chrom\tretro\tnonexonpart\t$from\t$to\t$score\t\.\t.\tgrp=$name:$chromStart;src=$source;pri=$priority\n";
	    $from = $newfrom;
	} elsif ($from<0) {
	    $from = $newfrom;
	}
	$to = $newto;
    }
    print "$chrom\tretro\tnonexonpart\t$from\t$to\t$score\t\.\t.\tgrp=$name:$chromStart;src=$source;pri=$priority\n";
}

