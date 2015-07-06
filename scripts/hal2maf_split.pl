#!/usr/bin/perl
# export a hal alignment to overlapping maf alignment chunks
# Stefanie Koenig 2.06.2015

use strict;
use warnings;

use Getopt::Long; # for parameter specification on the command line

my $usage = <<'ENDUSAGE';
hal2maf_split.pl	...
SYNOPSIS

hal2maf_split.pl --halfile aln.hal --refGenome genome

OPTIONS

    --help             			output this help message
    --keepDupes                         keep paralogs
    --keepAnncestors                    write ancestral sequences
    --refSequence                       name of reference sequence within reference genome (default: whole genome)
    --chunksize N                       size of the aligment chunk, e.g. sequence segment in the reference genome that covers the alignment chunk
    --overlap N                         length of overlap between the alignment chunks
    --cpus N                            number of cpus
    --hal_exec_dir                      path to hal executables. If not specified
                                        it must be in \$PATH environment variable.
    --no_split_list                     list of genomic intervals, in which genes are assumed.
                                        If such a list is given, the splitting of the alignmnet is outside of these regions. Format:

                                        chr2 120567671 120601255
                                        chr2 120604238 120609520
                                        chr5 65261850  65335670
                                        chr5 56530780  865308994
                                        ...

DESCRIPTION
      
  Example:
    hal2maf_split.pl --halfile flies.hal --refGenome dmel

ENDUSAGE

sub cannot_overlap;
sub find_splitting_point;

my ($halfile, $refGenome, $keepDupes, $keepAncestors, $refSequence, $no_split_list, $help); # options
my $h2m_param = "";
my $hal_exec_dir;
my $chunksize = 1500000;
my $overlap = 500000;
my $padding = 10000;
my $cpus = 1;

GetOptions('halfile:s'=>\$halfile,
           'refGenome:s'=>\$refGenome,
           'refSequence:s'=>\$refSequence,
	   'keepDupes!'=>\$keepDupes,
	   'keepAncestors!'=>\$keepAncestors,
	   'chunksize:i' =>\$chunksize,
	   'overlap:i' =>\$overlap,
	   'cpus:i' =>\$cpus,
	   'no_split_list:s' =>\$no_split_list,
	   'hal_exec_dir:s' =>\$hal_exec_dir,
           'help!'=>\$help);

if ($help){
    	print $usage;
    	exit(1);
}

if(!defined($halfile)){
    print "hal alignment missing.\n$usage";
    exit(1);
}

if(!defined($refGenome)){
    print "reference genome missing.\n$usage";
    exit(1);
}

if(!defined($refGenome)){
    print "reference genome missing.\n$usage";
    exit(1);
}
if(!defined($keepDupes)){
    $h2m_param.=" --noDupes";
}
if(!defined($keepAncestors)){
    $h2m_param.=" --noAncestors";
}

# check whether this perl module for paralell execution is installed
my $got_ForkManager = 0;
eval { require Parallel::ForkManager };
unless ($@) {
    $got_ForkManager = 1;
    Parallel::ForkManager->import();
}

if (!$got_ForkManager && $cpus > 1){
    print STDERR "The perl module Parallel::ForkManager is required to run hal2maf_split.pl in parallel.\n";
    print STDERR "Install this module first. On Ubuntu linux install with\n";
    print STDERR "sudo apt-get install libparallel-forkmanager-perl\n";
    print STDERR "Will now run sequentially (--cpus=1)...\n";
    $cpus = 1;
}

if(defined($hal_exec_dir)){
    if($hal_exec_dir =~/[^\/]$/){
	$hal_exec_dir .= '/';
    }
}
else{
    $hal_exec_dir = "";
}

my $halStats = $hal_exec_dir . "halStats";
my $hal2maf = $hal_exec_dir . "hal2maf";

# check whether halStats and hal2maf are properly installed
if (qx(which "$halStats") !~ /halStats$/){
 die ("$halStats is not executable. Please add the directory which contains the executable halStats to the PATH environment variable or specify the path with --hal_exec_dir.");
}

if (qx(which "$hal2maf") !~ /hal2maf$/){
 die ("$hal2maf is not executable. Please add the directory which contains the executable hal2maf to the PATH environment variable or specify the path with --hal_exec_dir.");
}

my @seqlist = qx($halStats --chromSizes $refGenome $halfile);
my $r_c = $?; # return code
if($r_c != 0){
    print "terminated after an error in halStats.\n";
    exit(1);
}

my %seqHash = ();
foreach(@seqlist){
    chomp;
    my($seqname,$len)=split("\t",$_);
    next if(defined($refSequence) && $refSequence ne $seqname);
    $seqHash{$seqname}=$len;
}

my @intervals = ();
if(defined($no_split_list)){
    open(IN, "<$no_split_list") or die ("Could not open $no_split_list");
    while(<IN>){
	chomp;
	my($chr,$start,$end)=split("\t",$_);
	push @intervals,[$chr, $start, $end];
    }
    close(IN);
}

# join overlapping genomic intervals
# sort intervals by 1. chromosome, 2. start
@intervals = sort {$a->[0] cmp $b->[0] || $a->[1] <=> $b->[1]} @intervals;

my @joined=(); # array of joined genomic intervals
my ($chr, $start, $end);

foreach (@intervals){
    if(defined($chr) && $_->[0] eq $chr && $_->[1] - $end <= 0 ) { # overlap between the last and the current interval
        if($end < $_->[2]){
            $end = $_->[2];
        }
    } else {
        if (defined($chr)){
            push @joined,[$chr, $start, $end];
        }
        ($chr, $start, $end) = ($_->[0], $_->[1], $_->[2]);
    }
}
push @joined,[$chr, $start, $end];

#print "no split intervals\n";
#foreach (@joined){
#    print "$_->[0]\t$_->[1]\t$_->[2]\n";
#}

# calculate start and end positions of alignment chunks
my @aln_chunks = ();

foreach my $seq (sort {$a cmp $b} keys %seqHash ){
    my $seqlen = $seqHash{$seq};
    my $start = 0;
    my $pred = 0;
    while($start + $chunksize <= $seqlen){
	my $end = $start + $chunksize - 1;
	while (@joined && cannot_overlap($joined[0], [$seq, $start, $end])){
	    $pred=$joined[0]->[2];
	    shift @joined;
	}
	
	foreach(@joined) {
	    last if($_->[0] gt $seq || ($_->[0] eq $seq && $_->[1] > $start));
	    if($_->[0] eq $seq && $_->[1] <= $start && $_->[2] >= $start){
		#print "overlapping start\n";
		#print "$seq:$start-$end\t$_->[0]:$_->[1]-$_->[2]\tpred=$pred\n";
		find_splitting_point($start, [$pred,$_->[1]]);
		last;
	    }
	    $pred=$_->[2];
	}
	foreach(@joined) {
	    last if($_->[0] gt $seq || ($_->[0] eq $seq && $_->[1] > $end));
	    if($_->[0] eq $seq && $_->[1] <= $end && $_->[2] >= $end){
		#print "overlapping end\n";
		#print "$seq:$start-$end\t$_->[0]:$_->[1]-$_->[2]\tpred=$pred\n";
		find_splitting_point($end, [$pred,$_->[1]]);
		last;
	    }
	    $pred=$_->[2];
	}
	push @aln_chunks,[$seq, $start, $end];
	$start = $end + 1 - $overlap;
    }
    my $end = $seqlen - 1;
    push @aln_chunks,[$seq, $start, $end];
}

# export alignment chunks in parallel to maf format
if($cpus > 1){
    my $pm = new Parallel::ForkManager($cpus);
    foreach(@aln_chunks){
	# this part is done in parallel by the child process
	my $pid = $pm->start and next;
	my ($seq, $start, $end) = ($_->[0], $_->[1], $_->[2]);
	my $chunksize = $end - $start + 1;
	my $cmd = "$hal2maf --refGenome $refGenome $h2m_param --refSequence $seq --start $start --length $chunksize $halfile $seq.$start-$end.maf";
	print "$cmd\n";
	qx($cmd);
	my $r_c = $?; # return code
	if($r_c != 0){
	    print "terminated after an error in hal2maf.\n";
	    exit(1);
	}
	$pm->finish; # terminate the child process
    }
    $pm->wait_all_children;
}
else{ # export alignment chunks sequentially to maf format
    foreach(@aln_chunks){
	my ($seq, $start, $end) = ($_->[0], $_->[1], $_->[2]);
	my $chunksize = $end - $start + 1;
	my $cmd = "$hal2maf --refGenome $refGenome $h2m_param --refSequence $seq --start $start --length $chunksize $halfile $seq.$start-$end.maf";
	print "$cmd\n";
	qx($cmd);
	my $r_c = $?;
	if($r_c != 0){
	    print "terminated after an error in hal2maf.\n";
	    exit(1);
	}
    }
}

sub cannot_overlap { 
    my $a = shift;
    my $b = shift;
    if ($a->[0] lt $b->[0] or ( $a->[0] eq $b->[0] and $a->[2] < $b->[1])){
        return 1;
    }
    return 0;
}

sub find_splitting_point {
    my $pos = shift;
    my $left = shift; # intergenic region left of pos
    my $left_start=$left->[0];
    my $left_end=$left->[1];
    my $split_point = $left_end - $padding;
    if($split_point < 0){
	$split_point = 0;
    }
    if($left_end > $split_point){
	$split_point = ($left_end + $left_start) / 2;
    }
#    print "splitpoint=$split_point\n$left_start-$left_end\n"; 
    return $split_point;
}
