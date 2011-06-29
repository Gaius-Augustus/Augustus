#!/usr/bin/perl
#
# choose an interval with good coverage
#
# Mario Stanke, March 2010
#
# Sebastian Adler, November 25 2010
#
use strict;
use Getopt::Long;

my $usage = "$0 -- choose an interval from each target with good coverage\n";
$usage .= "\n";
$usage .= "Usage: cat in.wig | [perl] $0 > out.bed\n";
$usage .= "Options:\n";
$usage .= "  --mincov=n     minimal coverage for each position of the interval (default: 1)\n";
$usage .= "  --minrelcov=f  minimal relative coverage (0 <= f <= 1, default: 0)\n";
$usage .= "  --maxgap=n     coverage gaps of maximal this length each are allowed (default: 0)\n";

my $help = 0;
my $mincov = 1;    # as given on the command line
my $mincoverage;   # adjusted
my $minrelcov = 0; # per default no constraint
my $maxgap = 0;

GetOptions(
    'mincov:n'=>\$mincov,
    'minrelcov:f'=>\$minrelcov,	#optional floating point value
    'maxgap:n'=>\$maxgap,
    'help!'=>\$help);

if ($help) {
    print "$usage";
    exit(0);
}

my $tname;
my ($pos, $cov);
my @ctglines = (); # all lines of the input file belonging to one contig

while (<>) {
    if (/variableStep chrom=(\S+)/){
	if (defined($tname)){
	    processCtg(\@ctglines);
	}
	$tname = $1;
	@ctglines = ();
	next;
    }
    push @ctglines, $_;
}

if (defined($tname)) { 
    processCtg(\@ctglines);
}

sub processCtg {
    my $ctglinesref = shift;
    my $r = getBestRange($ctglinesref, $mincov);
    #print "ohne Filter: $tname\t" . $r->{"a"} . "\t" . $r->{"b"} . "\t" . $r->{"score"} . "\n"; 
    if ($r->{"score"}>=0){
	my $relcov = $r->{"score"} / ($r->{"b"} - $r->{"a"} + 1);
	$mincoverage = $mincov;
	if ($mincoverage < $relcov * $minrelcov){
	    #print "increasing mincov from $mincov to " . ($relcov * $minrelcov) . " in $tname\n";
	    $mincoverage = $relcov * $minrelcov;
	    $r = getBestRange(\@ctglines, $mincoverage);
	}
	print "$tname\t" . $r->{"a"} . "\t" . $r->{"b"} . "\t" . $r->{"score"} . "\n"; # print best interval of previous target
    }
}

sub getBestRange {
    my $lines = shift;
    my $mc = shift;
    my $ma;  # left and right boundaries of the 
    my $mb;  # best interval seen so far
    my $ms;  # score of the best interval seen so far
    my $a;   # boundaries of the current/last interval
    my $b;
    my $s;
   
    $ma = $mb = $ms = $a = $b = -1;
    $s = 0;
    
    foreach my $line (@$lines){
	if ( $line =~ /(\d+)\s(\d+)/){
	    $pos = $1;
	    $cov = $2;
	    # update intervals and scores
	    next if ($cov < $mc);
	    if ($pos <= $b + 1 + $maxgap) {
		$b = $pos;
		$s += $cov;	# problem, if $s still undefined?
	    } else {
		$a = $b = $pos;
		$s = $cov; # fixed bug here on June 8th 2010: was s=pos
	    }
	    if ($ms < 0 || $ms < $s) {# current interval is best interval
		#print "new best: $a, $b, $s\n";
		($ma, $mb, $ms) = ($a, $b, $s);
	    }
	}
    }
    my %h= ("a" => $ma, "b" => $mb, "score" => $ms);
    return \%h;
}
