#!/usr/bin/perl
#
#  convert the exoniphy exon predictions to a file with hints to AUGUSTUS
#

use strict;
use Getopt::Long;
my ($usage, $cutoff, $priority, $grp, $seqname, $strand, $begin, $end);

#default values
$cutoff = 9; 
$priority = 4;

$usage .= "$0 -- convert exoniphy output to hints file\n";
$usage .= "Usage: $0 [--priority=p --cutoff=c]\n";
$usage .= "reads from stdin writes to stdout.\n";
$usage .= "takes input in tab separated format:\n";
$usage .= "name\tchrom\tstrand\tcdsStart\tcdsEnd\n";
$usage .= "You get this format on the UCSC system by typing e.g.\n";
$usage .= "echo \"select name, chrom, strand, cdsStart, cdsEnd from exoniphy order by chrom, cdsStart;\" | hgsql hg18 | grep -v cdsStart\n";

GetOptions('priority:i'=>\$priority,
           'cutoff:i'=>\$cutoff);

while (<>) {
    chomp;
    my @f = split /\t/, $_;
    next unless (@f > 4);
    $grp     = $f[0];
    $seqname = $f[1];
    $strand  = $f[2];
    $begin   = $f[3]+1;
    $end     = $f[4];

    print "$seqname\texoniphy\tCDS\t$begin\t$end\t0\t$strand\t.\tgrp=$grp;pri=$priority;src=X\n";
    $begin += $cutoff;
    $end   -= $cutoff;
    if ($begin > $end) {
	$begin = $end = int(($begin+$end)/2);
    }
    print "$seqname\texoniphy\tCDSpart\t$begin\t$end\t0\t$strand\t.\tgrp=$grp;pri=$priority;src=X\n";
}
