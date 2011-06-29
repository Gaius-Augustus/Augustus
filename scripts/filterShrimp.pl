#!/usr/bin/perl
#
# filter a file from the Shrimp alignment program
#
#
#
use strict;
use Getopt::Long;

my $usage = "$0 -- filter a SHRIMP format file\n";
$usage .= "\n";
$usage .= "Usage: $0 <in.psl >out.f.psl\n";
$usage .= "  PREREQUISITE: input file must be sorted by query name\n";
$usage .= "  \n";
$usage .= "  options:\n";
$usage .= "  --minScore=n       minimal percentage of identity (default 300)\n";
$usage .= "  --uniq             take only best match and only, when second best is much worse (default false)\n";
$usage .= "  --uniqthresh=f     the best alignment is considered uniquen when the second best has a score at most this much fraction of the best (default 0.9)\n";
$usage .= "  --best             take (all) best alignment(s) if they pass the minScore filter.\n";
$usage .= "  --commongenefile=s file name in which to write cases where one read maps to several different genes\n";
$usage .= "  --verbose          output debugging info (default false)\n";


my $minScore = 300;;
my $uniqthresh = 0.9; # a match is considered unique if the second best match has less than this percentage of the best
my $uniq = 0;
my $best = 0;
my $commongenefile;
my $verbose = 0;
my $help = 0;
my $cmdline = join(" ", @ARGV);

GetOptions(
    'help!'=>\$help,
    'minScore:i'=>\$minScore,
    'uniqthresh:f'=>\$uniqthresh,
    'uniq!'=>\$uniq,
    'best!'=>\$best,
    'commongenefile:s'=>\$commongenefile,
    'verbose!'=>\$verbose);

if ($help) {
    print "$usage";
    exit(0);
}
if ($uniq && $best){
    print "--uniq and --best cannot be used at the same time\n";
    exit(0);
}

my ($readname, $contigname, $strand, $contigstart, $contigend, $readstart, $readend, $readlength, $score, $editstring);
my $oldreadname = "";
my @f;
my ($outMinScore, $outUnique, $outBest) = (0,0,0); # number of reasons for filtering (nested, this order)
my @qali = (); # array of array references: lines for each query (pair)
open (COMMON, ">$commongenefile") or die ("Could not open $commongenefile for writing.") if (defined($commongenefile));

while (<>) {
    s/^#.*//;
    next unless /\S/;

    @f = split /\t/, $_;
    if (@f < 10) { warn "Not SHRIMP format"; next }
    # ($readname, $contigname, $strand, $contigstart, $contigend, $readstart, $readend, $readlength, $score, $editstring) = @f;
    $readname = $f[0];
    $contigname = $f[1];
    $score    = $f[8];
    
    processQuery() if ($oldreadname ne $readname && $oldreadname ne "");

    # filter for minimum score
    if ($score < $minScore){
	$outMinScore++;
	next;
    }

    push @qali, [$_, $score, $contigname];
    $oldreadname = $readname;
}
processQuery() if ($readname ne "");

close COMMON if (defined($commongenefile));


print STDERR "\n        filtered:\n";
print STDERR "----------------:\n";
print STDERR "mininum score   : $outMinScore\n";
print STDERR "uniqueness      : $outUnique\n";
print STDERR "best            : $outBest\n";
print STDERR "command line: $cmdline\n";

sub processQuery(){
    # print "processing " . scalar(@qali) . " alignments\n";
    
    if (($uniq || $best) && @qali>1){
	my %rscores;
	foreach my $ali (@qali){
	    $rscores{$ali} = scoreAli($ali); # store scores in hash so later sorting is faster
	}
	@qali = sort {-($rscores{$a} <=> $rscores{$b})} @qali;
	if ($uniq){
	    # let pass only best alignment, and only if second is significantly worse
	    my $ratio = $rscores{$qali[1]}/$rscores{$qali[0]};
	    if ($verbose){
		print "best two alignments\n";
		print "" . $qali[0]->[0] ."\nscore=$rscores{$qali[0]}\n";
		print "" . $qali[1]->[0] ."\nscore=$rscores{$qali[1]}\n";
		print "ratio = $ratio\n";
	    }
	    if ($ratio < $uniqthresh){
		print $qali[0]->[0];
	    } else {
		$outUnique++;
	    }
	} else { # best filtering, take all alignments that have the same score as the best
	    my $bestscore = $rscores{$qali[0]};
	    my @bestTnames = ();
	    while ($rscores{$qali[0]} == $bestscore){
		print $qali[0]->[0];
		push @bestTnames, $qali[0]->[2];
		shift @qali;
	    }
	    $outBest += @qali;
	    if (@bestTnames > 1){
		my %genenames = ();
		foreach my $Tname (@bestTnames) { $Tname =~ s/\.t\d+//; $genenames{$Tname}=1; }
		print COMMON $oldreadname . "\t" . join(" ", keys %genenames) . "\n" if (%genenames>1 && defined($commongenefile));
	    }
	}
    } else {
	foreach my $ali (@qali){
	    print $ali->[0];
	}
    }
    @qali = ();
}

# for comparing quality of two alignments
sub scoreAli{
    my $ali = shift;
    return $ali->[1];
}
