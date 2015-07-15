#!/usr/bin/perl
# this script works on top of the 'halTools' toolbox
# and exports a hal alignment to maf by splitting the alignment into several smaller
# alignment chunks of a certain size. An overlap between two consecutive alignment chunks can
# be specified. Furthermore, a list of 'genic' regions/intervals can be passed.
# In this case, the splitting is restricted to intergenic regions, e.g. position
# outside of the given intervals.
# Stefanie Koenig 2.06.2015

use strict;
use warnings;

use Getopt::Long; # for parameter specification on the command line

my $usage = <<'ENDUSAGE';
hal2maf_split.pl                        this script works on top of the 'halTools' toolbox
                                        and exports a hal alignment to maf by splitting the alignment into several smaller      
                                        alignment chunks of a certain size. An overlap between two consecutive alignment chunks can
                                        be specified. Furthermore, a list of 'genic' regions/intervals can be passed.    
                                        In this case, the splitting is restricted to intergenic regions, e.g. position
                                        outside of the given intervals.  	
SYNOPSIS

hal2maf_split.pl --halfile aln.hal --refGenome genome

    --halfile F                         F is the input hal file
    --refGenome S                       S is the name of the reference genome

OPTIONS

    --help             			output this help message
    --keepDupes                         keep duplicates, i.e. alignments of a sequence with itself (default: off)
    --keepAncestors                     export ancestral sequences (default: off)
    --refSequence S                     S is the name of the reference sequence within the reference genome
                                        (default: all sequences in the reference genome)
    --chunksize N                       size of the aligment chunk. N is the number of bases in the reference
                                        genome that are covered by the alignment chunks (default: 2500000)
    --overlap N                         overlap between to consecutive alignment chunks. N is the nunber of overlapping
                                        bases in the reference genome (default: 500000)
    --cpus N                            number of cpus (default: 1)
    --hal_exec_dir D                    D is the path to the hal executables. If not specified it must be in \$PATH environment variable.
    --no_split_list L                   list of 'genic' intervals. The splitting of the alignment is not allowed
                                        within these regions.  L is a file with the following format:
                                        seqname <tab> start <tab> end <newline>. Example:

                                        chr2 120567671 120601255
                                        chr2 120604238 120609520
                                        chr5 65261850  65335670
                                        chr5 56530780  865308994
                                        ...

DESCRIPTION
      
  Example:
    hal2maf_split.pl --halfile flies.hal --refGenome dmel --refSequence 3L --cpus 8 --hal_exec_dir /home/stefanie/tools/progressiveCactus/submodules/hal/bin

ENDUSAGE

sub cannot_overlap;

my ($halfile, $refGenome, $keepDupes, $keepAncestors, $refSequence, $no_split_list, $help); # options
my $h2m_param = "";
my $hal_exec_dir;
my $chunksize = 2500000;
my $overlap = 500000;
my $padding = 10000; 
my $min_intergenic = 2000; # when a list of genic intervals is given,
                           # all intervals that are within this distance are combined to a single interval
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

# reading in the name and length of each sequence in the reference genome ...
my @seqlist = qx($halStats --chromSizes $refGenome $halfile);
my $r_c = $?; # return code
if($r_c != 0){
    print "terminated after an error in halStats.\n";
    exit(1);
}
# ... and store it in a hash
my %seqHash = ();
foreach(@seqlist){
    chomp;
    my($seqname,$len)=split("\t",$_);
    next if(defined($refSequence) && $refSequence ne $seqname); # if a reference sequence is specified, skip all others
    $seqHash{$seqname}=$len;
}

if(defined($refSequence) && !exists($seqHash{$refSequence})){
    print "Reference sequence $refSequence not found in reference genome $refGenome\n";
    exit(1);
}

# reading in no_split_list
my @intervals = ();
if(defined($no_split_list)){
    open(IN, "<$no_split_list") or die ("Could not open $no_split_list");
    while(<IN>){
	chomp;
	my($chr,$start,$end)=split("\t",$_);
	next if(defined($refSequence) && $refSequence ne $chr); # if a reference sequence is specified, skip all others
	push @intervals,[$chr, $start, $end];
    }
    close(IN);
}

# join neighboring genic intervals within a distance of atmost $min_intergenic bases
# sort intervals by 1. chromosome, 2. start
@intervals = sort {$a->[0] cmp $b->[0] || $a->[1] <=> $b->[1]} @intervals;

my @joined=(); # array of sorted, joined genic intervals: Elements: <seq> <start> <end> <end of preceeding interval> <start of succeeding interval>
my ($chr, $start, $end);
my $pred_end = 0;

foreach (@intervals){
    if(defined($chr) && $_->[0] eq $chr && $_->[1] - $end <= $min_intergenic ) {
        if($end < $_->[2]){
            $end = $_->[2];
        }
    } else {
        if (defined($chr)){
	    my $succ_start = $seqHash{$chr} - 1; # start position of succeeding interval on this sequence
	    if($_->[0] eq $chr){ 
		$succ_start = $_->[1];
	    }
            push @joined,[$chr, $start, $end, $pred_end, $succ_start];
	    $pred_end = $end;  # end position of preceeding interval on this sequence
	    if($_->[0] ne $chr){
		$pred_end = 0;
	    }
        }
        ($chr, $start, $end) = ($_->[0], $_->[1], $_->[2]);
    }
}
if (defined($chr)){
    my $succ_start = $seqHash{$chr} - 1;
    push @joined,[$chr, $start, $end, $pred_end, $succ_start];
}


# calculate start and end positions of alignment chunks
#
# if a list of genic intervals is specified, the start/end of an alignment chunk is shifted to
# an intergenic region if it is contained in a genic interval:
#
# case 1)
# the start position s of an alignment chunk A is contained in a genic interval 
#
#                                                                 s
#                                                                 |----------------------|  alignment chunk A
#                               |---------|                  |---------|    genic intervals I1,...,Im (all disjunct) 
#                                         y                  x
#
# => s is shifted to the upstream intergenic region: s=max{(x+y)/2, x-10kb}
#
# case 2)
# the end position e of an alignment chunk A is contained in a genic interval 
#
#                                    e
#             |----------------------|  alignment chunk A
#                               |---------|                  |---------|    genic intervals I1,...,Im (all disjunct) 
#                                         x                 y
#
# => e is shifted to the downstream intergenic region: e=min{(x+y)/2, x+10kb}

my @aln_chunks = (); # each element is one alignment chunk, e.g. sequence segment (seqname:start-end) in the reference

foreach my $seq (sort {$a cmp $b} keys %seqHash ){
    my $seqlen = $seqHash{$seq};
    my $start = 0;
    while($start + $chunksize <= $seqlen){
	my $end = $start + $chunksize - 1;
        # go to first genic interval that may contain the start position of the alignment chunk 
	while (@joined && cannot_overlap($joined[0], [$seq, $start, $end])){
	    shift @joined;
	}	
	foreach(@joined) {
	    last if($_->[0] gt $seq || ($_->[0] eq $seq && $_->[1] > $start));
	    if($_->[0] eq $seq && $_->[1] <= $start && $_->[2] >= $start){ # start position of alignment chunk within a genic interval
		$start = int(($_->[1]+$_->[3])/2); #shift start position fo upstream intergenic region
		if($_->[1]-$padding > $start){
		    $start=$_->[1]-$padding;
		}
		last;
	    }
	}
	foreach(@joined) {
	    last if($_->[0] gt $seq || ($_->[0] eq $seq && $_->[1] > $end));
	    if($_->[0] eq $seq && $_->[1] <= $end && $_->[2] >= $end){ # end position of alignment chunk within a genic interval
		$end = int(($_->[2]+$_->[4])/2); # shift end position to downstream intergenic region
		if($_->[2]+$padding<$end){
		    $end=$_->[2]+$padding;
		}
		last;
	    }
	}
	push @aln_chunks,[$seq, $start, $end]; # add alignment chunk
	$start = $end + 1 - $overlap; # start position of next alignment chunk
    }
    # last alignment chunk
    my $end = $seqlen - 1;
    foreach(@joined) {
	last if($_->[0] gt $seq || ($_->[0] eq $seq && $_->[1] > $start));
	if($_->[0] eq $seq && $_->[1] <= $start && $_->[2] >= $start){
	    $start = int(($_->[1]+$_->[3])/2); 
	    if($_->[1]-$padding > $start){
		$start=$_->[1]-$padding;
	    }
	    last;
	}
    }
    push @aln_chunks,[$seq, $start, $end]; # add alignment chunk
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

# returns true if interval a
# - comes before interval b in the sorting order and
# - does not overlap with b
sub cannot_overlap { 
    my $a = shift;
    my $b = shift;
    if ($a->[0] lt $b->[0] or ( $a->[0] eq $b->[0] and $a->[2] < $b->[1])){
        return 1;
    }
    return 0;
}
