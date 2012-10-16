#!/usr/bin/perl
#
# filter a psl file (from BLAT,GMAP)
#
# Mario Stanke, September 2009
#
use strict;
use Getopt::Long;

my $usage = "$0 -- filter a psl file (e.g. BLAT or GMAP)\n";
$usage .= "\n";
$usage .= "Usage: $0 <in.psl >out.f.psl\n";
$usage .= "  PREREQUISITE: input psl file must be sorted by query name (standard with BLAT and GMAP)\n";
$usage .= "                Do a sort -k 10,10 but be aware: LC_ALL may have to be set to C because sort ignores characters like \":\"\n";
$usage .= "  if option 'paired' is used then it expects .f,.r or /1,/2 suffixes of mate pairs\n";
$usage .= "  \n";
$usage .= "  options:\n";
$usage .= "  --pairbed=s        file name of pairedness coverage:\n";
$usage .= "                     a .bed format file in which for each position the number of filtered\n";
$usage .= "                     read pairs is reported that contain the position in or between the reads\n";
$usage .= "  --minId=n          minimal percentage of identity (default 92)\n";
$usage .= "  --minCover=n       minimal percentage of coverage of the query read (default 80)\n";
$usage .= "  --uniq             take only best match and only, when second best is much worse (default false)\n";
$usage .= "  --uniqthresh       threshold % for uniq, second best must be at most this fraction of best (default .96)\n";
$usage .= "  --best             output all best matches that satisfy minId and minCover\n";
$usage .= "  --commongenefile=s file name in which to write cases where one read maps to several different genes\n";
$usage .= "  --nointrons        do not allow longer gaps (for RNA-RNA alignments)\n";
$usage .= "  --paired           require that paired reads are on opposite strands of same target(default false)\n";
$usage .= "  --maxintronlen=n   maximal separation of paired reads (default 500000\n";
$usage .= "  --verbose          output debugging info (default false)\n";


my $maxintronlen = 500000;
my $minintronlen = 35;
my $maxSortesTest = 100000; # check sortedness only for this many of the first lines to save memory
my $minId = 92;
my $minCover = 80;
my $uniqthresh = 0.96; # a match is considered unique if the second best match has less than this percentage of the best
my $uniq = 0;
my $nointrons = 0;
my $best = 0;
my $commongenefile;
my $pairbedfile;
my $paired = 0;
my $verbose = 0;
my $help = 0;
my %qnamestems = (); # for checking sortedness.
my $cmdline = join(" ", @ARGV);
my $maxCountInsert = 1000000;

GetOptions(
    'help!'=>\$help,
    'maxintronlen:i'=>\$maxintronlen,
    'minId:i'=>\$minId,
    'minCover:i'=>\$minCover,
    'uniqthresh:f'=>\$uniqthresh,
    'paired!'=>\$paired,
    'uniq!'=>\$uniq,
    'nointrons!'=>\$nointrons,
    'best!'=>\$best,
    'commongenefile:s'=>\$commongenefile,
    'pairbed:s'=>\$pairbedfile,
    'verbose!'=>\$verbose);

if ($help) {
    print "$usage";
    exit(0);
}


my ($match,$TgapCount,$strand,$qname,$qstart,$qend,$qsize,$targetname,$tstart,$tend,$blockSizes,$qStarts,$tStarts,$qBaseInsert,$tBaseInsert);
my ($qnamestem,$qsuffix);
my $skiplines=0;
my ($line,$lastCompactifyLine)=(0,0);
my $oldqnamestem = "";
my (@f,@b,@t,@q,@insertlen);
my ($outMinId,$outMinCover,$outPaired,$outUniq,$outBest,$outIntrons) = (0,0,0,0,0,0); # number of reasons for filtering (nested, this order)
my @qali = (); # array of array references: lines for each query (pair)
my %paircovsteps = (); # for pairedness coverage
                       # keys: target names (chromosomes)
                       # values: array references
                       #         elements: [pos,diff]
                       #                   position is 0-based, diff is +1 at start or -1 at end of mate pair
                       #                XXXXXXXXXXXX---------------------XXXXXXXXX
                       #               +1                                         -1


open (COMMON, ">$commongenefile") or die ("Could not open $commongenefile for writing.") if (defined($commongenefile));

while (<>) {
    $skiplines=5 if (/psLayout/);
    if ($skiplines>0) {
        $skiplines--;
        next;
    }

    $line++;
    s/^#.*//;
    next unless /\S/;
    if ($line%100000==1){
	$| = 1;
	print STDERR "\r"."processed line $line";
    }
    if (defined($pairbedfile) && $line % 10000000 == 0 && $line >= $lastCompactifyLine * 1.5){
	print STDERR "\ncompactifing coverage after $line lines ...";
	compactifyBed();
	$lastCompactifyLine = $line;
	print STDERR "done\n";
    }

    @f = split /\t/, $_, 21;
    if (@f < 20) { warn "Not PSL format"; next }
    
    $match       = $f[0];
    $qBaseInsert = $f[5];
    $TgapCount   = $f[6];
    $tBaseInsert = $f[7];
    $strand      = $f[8];
    $qname       = $f[9];
    $qsize       = $f[10];
    $qstart      = $f[11];
    $qend        = $f[12];    
    $targetname  = $f[13];
    $tstart      = $f[15];
    $tend        = $f[16];
    $blockSizes  = $f[18];
    $qStarts     = $f[19];
    $tStarts     = $f[20];
    
    $blockSizes =~ s/[, ]$//;
    $tStarts =~ s/[, ]$//;
    @b = split /,/, $blockSizes; #
    @t = split /,/, $tStarts;
    @q = split /,/, $qStarts;
    
    $qnamestem = $qname;
    if ($paired){
        $qnamestem =~ s/[\.\/]([fr12])$//;
        $qsuffix = $1; # f,1: forward mate, r,2: reverse mate "": no mates
    }
    if ($oldqnamestem ne $qnamestem && $oldqnamestem ne ""){
	if ($line <= $maxSortesTest && $qnamestems{$qnamestem}){
	    print STDERR "Input file not sorted by query name! $qnamestem occurred previously. Set LC_ALL=C and sort -k 10,10\n";
	    exit 1;
	}
	processQuery() if (@qali);
    }

    # filter for minimum percentage of identity
    my $tgap = 0; # inserted bases in target sequence, excluding introns
    my $qgap = 0; # inserted bases in query sequence
    my $gaps = 0; # total amount of gaps
    
    
    for (my $i=0; $i<@b-1; $i++){
	$tgap = $t[$i+1]-$t[$i]-$b[$i];
	$qgap = $q[$i+1]-$q[$i]-$b[$i];
	$tgap = 0 if ($qgap ==0 && $tgap >= $minintronlen && $tgap <= $maxintronlen); #target gap is intron
	$gaps += ($tgap>$qgap)?  $tgap : $qgap; # count the larger gap if both seqs happen to have a gap
    }
    my $percid = sprintf("%.1f", 100*$match/((($qend - $qstart + $gaps)))); 
    if ($percid < $minId){
	$outMinId++;
	next;
    }

    # filter for minimum coverage
    my $coverage =  sprintf("%.1f", 100*($qend - $qstart)/$qsize);
    if ($coverage < $minCover){
	$outMinCover++;
	next;
    }

    # filter for introns
    if ($nointrons && $qBaseInsert + $tBaseInsert > 10){
	$outIntrons++;
	next;
    }

    push @qali, [$_, $targetname, $qsuffix, $strand, $tstart, $tend, $percid, $coverage];
    # print "$targetname, $qsuffix, $strand, $tstart, $tend, $percid\n";

    $oldqnamestem = $qnamestem;
    $qnamestems{$qnamestem} = 1 if ($line <= $maxSortesTest);
}
processQuery() if ($qnamestem ne "");

close COMMON if (defined($commongenefile));

#
# write pairedness coverage info into the pairbedfile
#
if (defined($pairbedfile)){
    open (PAIRBED, ">$pairbedfile") or die ("Could not open $pairbedfile for writing.");
    print PAIRBED "track type=bedGraph name=\"pairedness coverage\" description=\"pairedness coverage\"";
    print PAIRBED " visibility=full color=200,100,0 altColor=200,100,0\n";
    compactifyBed();
    foreach my $chr (sort keys %paircovsteps){
	my $cov = 0;
	my $pos = 0;
	next if (!@{$paircovsteps{$chr}});
	foreach my $step (@{$paircovsteps{$chr}}){
	    print PAIRBED "$chr\t$pos\t$step->[0]\t$cov\n" if ($pos<$step->[0] && $cov>0);
	    $pos = $step->[0];
	    $cov += $step->[1];
	}
	warn ("inconsistent") if ($cov!=0);
    }
    close PAIRBED;
}

@insertlen = sort {$a <=> $b} @insertlen;

#foreach(@insertlen){
#    print STDOUT "Insert\t".$_."\n";
#}

print STDERR "\n        filtered:\n";
print STDERR "----------------:\n";
print STDERR "percent identity: $outMinId\n";
print STDERR "coverage        : $outMinCover\n";
print STDERR "nointrons       : $outIntrons\n" if ($nointrons);
if ($paired) {
    print STDERR "not paired      : $outPaired\n" if ($paired);
    print STDERR "quantiles of unspliced insert lengths: ";
    for (my $i=1;$i<10;$i++){
	print STDERR "q[" . (10*$i) . "%]=" . ($insertlen[int($i*@insertlen/10)]) . ", ";
    }
    print STDERR "\n";
}
print STDERR "uniq            : $outUniq\n" if ($uniq);
print STDERR "best            : $outBest\n" if ($best);
print STDERR "command line: $cmdline\n";

sub processQuery(){
    # print "processing " . scalar(@qali) . " alignments\n";
    
    # filter @qali based on mate pair consistency
    # keep only alignments for which there is a possible mate:
    # 1) same chromosome
    # 2) different strand
    # 3) distance in genome < minintronlen
    if ($paired){
	@qali = sort {$a->[1] cmp $b->[1] || $a->[4] cmp $b->[4]} @qali;
	my @matepairs = ();
	my %mated = ();
	for (my $i=0;$i < @qali-1; $i++){
	    for (my $j=$i+1; $j < @qali && $qali[$i]->[1] eq $qali[$j]->[1]; $j++){ # only loop until leave chromosome
		#print "comparing $i,$qali[$i]->[1], with $j,$qali[$j]->[1]\n";
		if ($qali[$i]->[2] ne $qali[$j]->[2]){    # different mate: (f,r) or (1,2)
		    if ($qali[$i]->[3] ne $qali[$j]->[3]){ # different strand
			my $dist = $qali[$j]->[4] - $qali[$i]->[5] - 1;
			$dist = $qali[$i]->[4] - $qali[$j]->[5] - 1 if ($qali[$i]->[4] > $qali[$j]->[4]);
			if ($dist < $maxintronlen && $dist>=0){ # not too far apart, not overlapping either
			    # print "found mate pair $i,$j\n";
			    push @matepairs, [$i,$j,scoreMate($i,$j,$dist)];
			    $mated{$i}=0 if (!defined($mated{$i}));
			    $mated{$j}=0 if (!defined($mated{$j}));
			    $mated{$i}++;
			    $mated{$j}++;
			    my $inslen = $qali[$j]->[5] - $qali[$i]->[4] - 1;
			    $inslen = $qali[$i]->[5] - $qali[$j]->[4] - 1 if ($inslen<0);
			    push @insertlen, $inslen if (@insertlen < $maxCountInsert); # if not limited, this may use on huge files an unwarranted amount of RAM for a simple statistics
			} else {
			    #print "distance not right\n";
			}
		    }
		}
	    }
	}
	# print "found " . scalar(@matepairs) . " mate pairs, involving " . scalar(keys %mated) . " mates.\n" if ($verbose);
	$outPaired += @qali - scalar(keys %mated);
	if ((!$uniq && !$best) || @matepairs<2){# let pass all read alignments that are involved in mate pairs
	    foreach my $i (sort {$a <=> $b} keys %mated){
		print $qali[$i]->[0];
	    }
	} else {# uniq or best
	    @matepairs = sort {$b->[2] <=> $a->[2]} @matepairs;
	    if ($uniq){# let pass only best mate pair, and only if second is significantly worse
		my $second = 1;
		while ($second < @matepairs && similar($qali[$matepairs[0]->[0]], $qali[$matepairs[$second]->[0]], $qali[$matepairs[0]->[1]], $qali[$matepairs[$second]->[1]])){
		    $second++; 
		}
		if ($second < @matepairs){
		    my $ratio = $matepairs[$second]->[2] / $matepairs[0]->[2];
		    if ($verbose) {
			print "\nbest two mates\n";
			print "" . join (", ", @{$qali[$matepairs[0]->[0]]}) . "\npaired with\n"
			    .join (" ", @{$qali[$matepairs[0]->[1]]})
			    . "\nscore=$matepairs[0]->[2]\n";
			print "" . join (", ", @{$qali[$matepairs[1]->[0]]}) . "\npaired with\n"
			    .join (" ", @{$qali[$matepairs[1]->[1]]})
			    . "\nscore=$matepairs[1]->[2]\n";
			print "ratio = $ratio\n";
		    }
		    if ($ratio < $uniqthresh){
			# print the two alignments for best mate pair only
			print $qali[$matepairs[0]->[0]]->[0];
			print $qali[$matepairs[0]->[1]]->[0];
			$outUniq += @qali-1;
		    } else {
			@matepairs = ();
			$outUniq += @qali;
		    }
		} else {
		    print "suboptimal mate pairs are similar\n" if ($verbose);
		    print $qali[$matepairs[0]->[0]]->[0];
		    print $qali[$matepairs[0]->[1]]->[0];
		}
		splice @matepairs, 1; # keep only the best pair (if any)
	    } else { # best: take all best alignment pairs
		my $optscore = $matepairs[0]->[2];
                my @bestTnames = ();
		my $numbest = 0;
		my %haveOutput = (); # remember alignments that were printed, as some alignments can be part of several equally good pairs
		while ($numbest < @matepairs && $matepairs[$numbest]->[2] == $optscore){
		    my ($aliline1, $aliline2) = ($qali[$matepairs[$numbest]->[0]]->[0], $qali[$matepairs[$numbest]->[1]]->[0]);
		    print $aliline1 if (!exists($haveOutput{$aliline1}));
		    $haveOutput{$aliline1} = 1;
                    print $qali[$matepairs[$numbest]->[1]]->[0] if (!exists($haveOutput{$aliline2}));
                    $haveOutput{$aliline2} = 1;
                    push @bestTnames, $qali[$matepairs[$numbest]->[0]]->[1];
		    $numbest++;
		}
		$outBest += @matepairs - $numbest;
		splice @matepairs, $numbest; # keep only the first $numbest pairs
                if (@bestTnames>1){
                    my %genenames = ();
                    foreach my $Tname (@bestTnames) { $Tname =~ s/\.t\d+//; $genenames{$Tname}=1; }
                    print COMMON $oldqnamestem . "\t" . join(" ", keys %genenames) . "\n" if (%genenames > 1 && defined($commongenefile));
                }
	    }
	}
	# output pairedbed info: go through list of all mate pairs and store start and end position
	if (defined($pairbedfile)){
	    while (@matepairs>0){
		my $chr = $qali[$matepairs[0]->[0]]->[1];
		$paircovsteps{$chr} = [] if (!defined($paircovsteps{$chr}));
		my $pend = $qali[$matepairs[0]->[1]]->[5];
		my $pstart = $qali[$matepairs[0]->[0]]->[4];
		push @{$paircovsteps{$chr}}, [$pstart-1, 1];
		push @{$paircovsteps{$chr}}, [$pend, -1];
		shift @matepairs;
	    }
	}
    } else { # not paired, single read
	if (($uniq || $best) && @qali>1){
	    my %rscores;
	    foreach my $ali (@qali){
		$rscores{$ali} = scoreAli($ali); # store scores in hash so later sorting is faster
	    }
	    @qali = sort {-($rscores{$a} <=> $rscores{$b})} @qali;
	    if ($uniq) {
		my $second = 1;
                while ($second < @qali && similar($qali[0], $qali[$second])){
                    $second++;
                }
                if ($second < @qali){
		    # let pass only best mate pair, and only if second is significantly worse
		    my $ratio = $rscores{$qali[$second]}/$rscores{$qali[0]};
		    if ($verbose){
			print "best two alignments\n";
			print "" . $qali[0]->[0] ."\nscore=$rscores{$qali[0]}\n";
			print "" . $qali[1]->[0] ."\nscore=$rscores{$qali[1]}\n";
			print "ratio = $ratio\n";
		    }
		    if ($ratio < $uniqthresh){
			print $qali[0]->[0];
			$outUniq += @qali-1;
		    } else {
			$outUniq += @qali; # drop all
		    }
		} else {
		    print "suboptimal alignments are similar\n" if ($verbose);
		    print $qali[0]->[0];
		    $outUniq += @qali-1;
		}
	    } else { # take all best alignments, that share the maximum score
		my $optscore = $rscores{$qali[0]};
		my @bestTnames = ();
		while ($rscores{$qali[0]} == $optscore){
		    print $qali[0]->[0];
		    push @bestTnames, $qali[0]->[1];
		    shift @qali;
		}
		$outBest += @qali;
		if (@bestTnames>1){
		    my %genenames = ();
		    foreach my $Tname (@bestTnames) { $Tname =~ s/\.t\d+//; $genenames{$Tname}=1; }
		    print COMMON $oldqnamestem . "\t" . join(" ", keys %genenames) . "\n" if (%genenames>1 && defined($commongenefile));
		}
	    }
	} else {
	    foreach my $ali (@qali){
		print $ali->[0];
	    }
	}
    }
    @qali = ();
}

# for comparing quality of two alignments
sub scoreAli{
    my $ali = shift;
    return $ali->[6]/100 # percent identity
	+ $ali->[7]/100; # percent coverage
}

# for comparing quality of two mate pair read alignments (i1,j1), (i2,j2)
sub scoreMate{
    my $i = shift;
    my $j = shift;
    my $dist = shift;
    my $score = ($qali[$i]->[6] + $qali[$j]->[6])/100 # percent identity
	+ ($qali[$i]->[7] + $qali[$j]->[7])/100; # percent coverage
    $score -= $dist/$maxintronlen/10 if (!$best); # penalty for distance between mates. Do not use if option 'best' is chosen, otherwise a one base difference may cause a difference
    return $score;
}

#
# checking whether two alignments (or two alignment pairs) are similar
# Purpose: Due to separate handling of spliced and unspliced alignments it can happen
# that very similar alignments are reported, e.g. an unspliced read going approximately up to an intron
# and a spliced read with a few base pairs on one exon.
# These should not be considered ambiguous when --uniq is specified.
sub similar {
    return similar(@_[0],@_[1]) && similar(@_[2],@_[3]) if (@_ == 4);
    my $r1 = shift;
    my $s1 = shift;
    if ($verbose){
	print "\nchecking whether these alignments are approximately the same\n";
	print "" . join (", ", @$r1) . "\nand\n" . join (", ", @$s1) . "\n\n";
    }
    return 0 if ($r1->[5] <= $s1->[4] || $r1->[4] >= $s1->[5]); # here: similar = overlapping target range
    print "they are similar\n" if ($verbose);
    return 1;
}


#
# compactifyBed
# if several steps coincide then summarize them equivalently by one step in order to 
# 1) save memory or
# 2) output a bed file
#
sub compactifyBed {
    my $before=0;
    my $after=0;
    foreach my $chr (sort keys %paircovsteps){
        next if (!@{$paircovsteps{$chr}});
        @{$paircovsteps{$chr}} = sort {$a->[0] <=> $b->[0]} @{$paircovsteps{$chr}};
	$before += scalar(@{$paircovsteps{$chr}});
#	print "before compactify\n";
#	foreach my $step (@{$paircovsteps{$chr}}){
#	    print $step->[0] . "\t" . $step->[1] . "\n";
#	}
	my $i=0;
	while ($i<@{$paircovsteps{$chr}}-1){
	    if ($paircovsteps{$chr}->[$i]->[0] eq $paircovsteps{$chr}->[$i+1]->[0]){
		$paircovsteps{$chr}->[$i]->[1] += $paircovsteps{$chr}->[$i+1]->[1]; # add value from i+1 to i
		splice @{$paircovsteps{$chr}}, $i+1, 1; # remove element i+1
	    } else {
		$i++;
	    }
	}
	$after += scalar(@{$paircovsteps{$chr}});
    }
    print STDERR "\nbefore compactifying: $before  after: $after\n";
}
