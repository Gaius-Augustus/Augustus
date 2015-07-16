#!/usr/bin/perl
#
# summarize multiple identical hints to one with mult=n
#
# Mario Stanke, 4.1.2010

use strict;
use Getopt::Long;

my $usage = "$0 -- summarize multiple identical hints to one with mult=n\n";
$usage .= "\n";
$usage .= "Usage: $0 <in.psl >joined.psl\n";
$usage .= "  PREREQUISITE: input GFF file must be sorted so that hints that should be summarized are below each other\n";
$usage .= "  e.g. do a cat hints.gff | sort -n -k 4,4 | sort -s -n -k 5,5 | sort -s -k 3,3 | sort -s -k 1,1 | join_mult_hints.gff\n";


my $help=0;
GetOptions('help!'=>\$help);
if ($help) {
    print "$usage";
    exit(0);       
}

my @f;
my @lf;
my ($lm,$m);

while (<>) {
    @f = split(/\t/,$_);
    if (!(@lf)){
	@lf = @f;	
    } elsif (!(($f[0] eq $lf[0]) && ($f[2] eq $lf[2]) && ($f[3] == $lf[3]) && ($f[4] == $lf[4])  && ($f[6] eq $lf[6]) && ($f[7] eq $lf[7]))){
	print join("\t",@lf);
	@lf = @f;
    } else {
	# update lf by adding f to it
	$lf[8] =~ s/gro?u?p=[^;]*;//;
	if ($lf[8] =~ /mult=(\d+);/){
	    $lm = $1;
	    $lf[8] =~ s/mult=\d+;//;
	} else {
	    $lm = 1;
	}
	if ($f[8] =~ /mult=(\d+);/){
	    $m = $1;
	} else {
	    $m = 1;
	}
	$lf[8] = "mult=" . ($lm+$m) . ";" . $lf[8];
    }
}
print join("\t",@lf) if (@lf);
