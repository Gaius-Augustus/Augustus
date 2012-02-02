#!/usr/bin/perl
#
# filter a psl file (from BLAT,GMAP)
#
# NOTE(S)
# This script only works with previously SORTED files!!!
#
# Created: 5-November-2011
# Last modified: 2-February-2012

use strict;
use Getopt::Long;
require 'auxFunctions/printAoA.pl';
require 'auxFunctions/simpleFilter.pl';
require 'auxFunctions/printMatedPairsInfo.pl';
require 'auxFunctions/printMatedHash.pl';
require 'auxFunctions/printHash.pl';
require 'auxFunctions/printRscores.pl';
require 'auxFunctions/printPairCovSteps.pl';
require 'auxFunctions/printChrOfPairCovSteps.pl';

my $usage = "$0 -- filter a psl file (e.g. BLAT or GMAP)\n"; # my $usage is an "strict way" of declaring the variable
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
my $processQueryCalls = 0;
my $computeMatePairsCalls = 0;
my $qaliSize = 0;
my $maxQali = 0;
# Counters
my $outMinCover = 0;
my $outMinId = 0;
my $outIntrons = 0;
my $maxQaliSize = 0;

# Defining filtering options
GetOptions(
    'help!'=>\$help,
    'maxintronlen:i'=>\$maxintronlen,     # optional integer argument
    'minId:i'=>\$minId,                   # optional integer argument
    'minCover:i'=>\$minCover,             # optional integer argument
    'uniqthresh:f'=>\$uniqthresh,         # optional floating point argument
    'paired!'=>\$paired,                  # negatable option
    'uniq!'=>\$uniq,                      # negatable option
    'nointrons!'=>\$nointrons,            # negatable option
    'best!'=>\$best,                      # negatable option
    'commongenefile:s'=>\$commongenefile, # optional string argument
    'pairbed:s'=>\$pairbedfile,           # optional string argument
    'verbose!'=>\$verbose);               # negatable option

if ($help) 
{
    print "$usage";
    exit(0);
}


my ($match,$TgapCount,$strand,$qname,$qstart,$qend,$qsize,$targetname,$tstart,$tend,$blockSizes,$qStarts,$tStarts,$qBaseInsert,$tBaseInsert);
my ($qnamestem,$qsuffix);
my $skiplines=0;
my $line=0;
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
my @params = ($minintronlen, $maxintronlen, $minId, $outMinId, $minCover, $outMinCover, $nointrons, $outIntrons);
my @results = ();
my ($percid, $coverage);

# Opening $commongenefile in order to write in it. COMMON is file handle
open (COMMON, ">$commongenefile") or die ("Could not open $commongenefile for writing.") if (defined($commongenefile));


while (<>) 
{
    $skiplines=5 if (/psLayout/); # If alignments were generated by psLayout
    if ($skiplines>0) 
	{
        $skiplines--;
        next;
    }

    $line++;
    s/^#.*//;             # Get rid of commentaries
    next unless /\S/;     # /\S/ means match whitespace character

	   # # Checks whether the 100000 line limit has been reached
       # if ($line%100000==1)
	   # { 
	   # $| = 1;  # Special variable for debugging and buffering
	   # print STDERR "\r"."processed line $line";
       # }

    # See help (above) for definition of $pairbedfile
    if (defined($pairbedfile) && $line % 10000000 == 0)
	{
		print STDERR "\ncompactifing coverage after $line lines ...";
		compactifyBed();      # Call to COMPACTIFYBED, see below
		print STDERR "done\n";
    }

	# If no line limit has been reached, then split the current line, with 
	# \t as the pattern. Recall that PSL format has 21 fields
    @f = split /\t/, $_, 21;
    if (@f < 20) 
	{ 
		warn "Not PSL format"; 
		next;
	}
    
	# Check the http://genome.ucsc.edu web page for more information about the 
	# format
    $match       = $f[0];  # No. of bases that match that aren't repeats
    $qBaseInsert = $f[5];  # No. of bases inserted in query
    $TgapCount   = $f[6];  # No. of inserts in target
    $tBaseInsert = $f[7];  # No. of bases inserted in target
    $strand      = $f[8];  # +, - for query strand. 2nd +,- for genomic strand
    $qname       = $f[9];  # Query name, e.g. SRR027108.21/1
    $qsize       = $f[10]; # Query sequence size
    $qstart      = $f[11]; # Alignment start position in query
    $qend        = $f[12]; # Alignment end position in query
    $targetname  = $f[13]; # Target sequence name
    $tstart      = $f[15]; # Alignment start position in target
    $tend        = $f[16]; # Alignment end position in target
    $blockSizes  = $f[18]; # List of sizes of each block
    $qStarts     = $f[19]; # List of starting positions of each block in query
    $tStarts     = $f[20]; # List of starting positions of each block in target
    
    $qnamestem = $qname;  # Update qnamestem with the current line's query name
    # if option 'paired' used then it expects .f,.r or /1,/2 
	# suffixes of mate pairs;
    if ($paired)
	{ 
        # Line below takes out "/1,/2,f" or "r" from e.g. stem SRR027108.21/1
		# print STDERR "qnamestem=$qnamestem before trimming\n";
        $qnamestem =~ s/[\.\/]([fr12])$//;  # REPLACE every pattern */[fr12] by nothing 
		# print STDERR "qnamestem=$qnamestem after trimming\n";
        $qsuffix = $1; # f,1: forward mate, r,2: reverse mate "": no mates
    }	
	# print "$qname, $qnamestem, $qsuffix, $oldqnamestem\n";


	 # # Printing contents of qnamestem
	 # printHash(\%qnamestems);
	 # print STDERR "Where the qnamestem to be processed is:$qnamestem\n";

     if ($oldqnamestem ne $qnamestem && $oldqnamestem ne "") # Check whether queries are different
	 {
        # print STDERR "Is $qnamestem contained already in \%qnamestems? ". $qnamestems{$qnamestem}."\n";

	 	# The if line below checks whether 10th field is sorted in ascending 
	 	# order
	 	if ($line <= $maxSortesTest && $qnamestems{$qnamestem})
	 	{
	 		print STDERR "Input file not sorted by query name! $qnamestem occurred previously. Set LC_ALL=C and sort -k 10,10\n";
	 		exit 1;
	 	}
	 	# By definition (above): @quali is an array of array references: 
	 	# lines for each query (pair)
	 	if (@qali)
		  {
			  $qaliSize = processQuery(\@qali, \@params) if (@qali); 
		  }
		if ($qaliSize > $maxQali)
		{$maxQali = $qaliSize;}
     }
	
	# Computing $coverage and $percid
      my @results = simpleFilter(\@f, \@params);
	  my $percid = $results[0];
	  my $coverage = $results[1];

 	# filter for percent identity
    if ($percid < $minId)
 	{
		printf STDERR "$qnamestem filtered out by percId=$percid\n";
 		$outMinId++;
 		next;
    }

    # filter for minimum coverage
    if ($coverage < $minCover)
 	{
		printf STDERR "$qnamestem filtered out by coverage=$coverage\n";
 		$outMinCover++;
 		next;
    }

    # filter for introns
	my $gapSize = $qBaseInsert + $tBaseInsert;
    if ($nointrons && $qBaseInsert + $tBaseInsert > 10)
 	{	
		printf STDERR "$qnamestem filtered out by gapSize=$gapSize (nointrons)\n";
 		$outIntrons++;
 		next;
    }


     push @qali, [$_, $targetname, $qsuffix, $strand, $tstart, $tend, $percid, $coverage];
	# push @qali, [$_, "targetname", $qsuffix, "strand", "tstart", "tend", "percid", "coverage"]; 

    $oldqnamestem = $qnamestem;
    $qnamestems{$qnamestem} = 1 if ($line <= $maxSortesTest);

} # END WHILE


if ($qnamestem ne "")
{
  print STDERR "processQuery() called because $qnamestem != ''\n";
  my ($refPairCovSteps, $maxQaliSize) = processQuery(\@qali, \@params); 
  %paircovsteps = %{$refPairCovSteps};
  printPairCovSteps(\%paircovsteps);
}

close COMMON if (defined($commongenefile));


 #
 # write pairedness coverage info into the pairbedfile
 #
 if (defined($pairbedfile))
 {
     open (PAIRBED, ">$pairbedfile") or die ("Could not open $pairbedfile for writing.");
     print PAIRBED "track type=bedGraph name=\"pairedness coverage\" description=\"pairedness coverage\"";
     print PAIRBED " visibility=full color=200,100,0 altColor=200,100,0\n";
	 printSizeOfCoverInfo(\%paircovsteps);
     compactifyBed();
	 printSizeOfCoverInfo(\%paircovsteps);

     foreach my $chr (sort keys %paircovsteps)
 	{
 		my $cov = 0;
 		my $pos = 0;
 		next if (!@{$paircovsteps{$chr}});
 		foreach my $step (@{$paircovsteps{$chr}})
 		{
 			print PAIRBED "$chr\t$pos\t$step->[0]\t$cov\n" if ($pos<$step->[0] && $cov>0);
 			$pos = $step->[0];
 			$cov += $step->[1];
 		}
 		warn ("inconsistent") if ($cov!=0);
     }
     close PAIRBED;
 } # end if


#Sorting insert lengths in Ascending order
print STDERR "----------------------------------------\n";
print STDERR "Insert lengths before sorting:\n";
foreach my $ins (@insertlen) { print STDERR "$ins\n";}
 @insertlen = sort {$a <=> $b} @insertlen;
print STDERR "Insert lengths after sorting:\n";
foreach my $ins (@insertlen) { print STDERR "$ins\n";}
print STDERR "----------------------------------------\n";


# Printing summary statistics of the filter
 print STDERR "\n        filtered:\n";
 print STDERR "----------------:\n";
 print STDERR "percent identity: $outMinId\n";
 print STDERR "coverage        : $outMinCover\n";
 print STDERR "nointrons       : $outIntrons\n" if ($nointrons);
 if ($paired) 
 {
     print STDERR "not paired      : $outPaired\n" if ($paired);
     print STDERR "quantiles of unspliced insert lengths: ";
     for (my $i=1;$i<10;$i++)
 	{
 		print STDERR "q[" . (10*$i) . "%]=" . ($insertlen[int($i*@insertlen/10)]) . ", ";
     }
     print STDERR "\n";
 } # end if
 
 print STDERR "unique          : $outUniq\n" if ($uniq);
 print STDERR "best            : $outBest\n" if ($best);
 print STDERR "command line: $cmdline\n";
 print STDERR "processQuery() was called: $processQueryCalls times\n";
 print STDERR "computeMatePairs() was called: $computeMatePairsCalls times\n";
 print STDERR "processed $line reads\n";
 print STDERR "maxQaliSize is: $maxQali\n";


###################
# Functions
###################


sub processQuery()
{
	$processQueryCalls++;
	my $maxQaliSize = 0;

	print STDERR "------------------------------------------------------\n";
	print STDERR "Calling processQuery: Size of \@qali=".scalar(@qali)." corresponding to $oldqnamestem. \n";
	print STDERR "------------------------------------------------------\n";

	my ($qaliRef, $paramsRef) = @_;
	# De-reference the array list
	my @f = @$qaliRef;
	my @params = @$paramsRef;
	my $minintronlen = $params[0];
	my $maxintronlen = $params[1];
	my $minId        = $params[2];
	my $outMinId     = $params[3];
	my $minCover     = $params[4];
	my $outMinCover  = $params[5];
	my $nointrons    = $params[6];
	my $outIntrons   = $params[7];
	my (@bestTnames, %genenames);


    # print "processing " . scalar(@qali) . " alignments\n";
    # filter @qali based on mate pair consistency
    # keep only alignments for which there is a possible mate:
    # 1) same chromosome
    # 2) different strand
    # 3) distance in genome < minintronlen
    if ($paired)
	{   
	 	$computeMatePairsCalls++;

	    # Computing and printing mated pairs
		my ($refmatepairs, $refinsertlen, $refmated, $maxQaliSize) = computeMatePairs();
		my @matepairs = @$refmatepairs;
		my @insertlen = @$refinsertlen;
		my %mated = %{$refmated}; # keeps track of mates (keys) and their number of matings (value)

		# Printing by-products of computeMatePairs()
		print STDERR "------------------------------------------------------\n";
		printMatedPairsInfo(\@matepairs, \@qali);
		printMatedHash(\%mated);
		print STDERR "------------------------------------------------------\n";
		print STDERR "(paired,uniq,best)=($paired,$uniq,$best)\n";

		$outPaired += @qali - scalar(keys %mated); # Taking out all alignments that were not paired
		# print STDERR "size(qali)=". scalar(@qali). " No.mated alignments= ". 
		#       scalar(keys %mated) . ", outPaired=$outPaired\n";
		if ((!$uniq && !$best) || @matepairs<2)
 	    { # let pass all read alignments that are involved in mate pairs
			print STDERR "------------------------------------------------\n";
			print STDERR "Letting pass all mated-paired alignments= " . scalar(keys %mated) 
			.", listed below:\n";
			print STDERR "Size of \@matepairs=".scalar(@matepairs)."\n";
			print STDERR "------------------------------------------------\n";
			my @mymate; 
			foreach my $i (sort {$a <=> $b} keys %mated)
 	 	    {
				@mymate = split /\t/, $qali[$i]->[0], 21;
				print STDERR "Letting pass paired-alignment ($i): $mymate[9]\n";
				print $qali[$i]->[0] . "percid=$qali[$i]->[6], coverage=$qali[$i]->[7]\n";
 	 	    }
 	    } else {# ($uniq or $best) selected

			print STDERR "------------------------------------------------\n";
			print STDERR "Sort (descending) mated pairs by scoreMate\n";
			print STDERR "BEFORE sorting\n";
			printMatedPairsInfo(\@matepairs, \@qali);
		    @matepairs = sort {$b->[2] <=> $a->[2]} @matepairs; # descending order
		    # @matepairs = sort {-($b->[2] <=> $a->[2])} @matepairs; # ascending order
			print STDERR "AFTER sorting\n";
			printMatedPairsInfo(\@matepairs, \@qali);
			print STDERR "------------------------------------------------\n";

			if ($uniq)
 		      {# let pass only best mate pair, and only if second is significantly worse
				print STDERR "Selecting a unique mate pair\n";
			    my $second = 1;
				print STDERR "------------------------------------------------\n";
				print STDERR "Comparing similarity between pairs:\n";
				print STDERR "$matepairs[0]->[0] and $matepairs[$second]->[0]\n"; 
				print STDERR "$matepairs[0]->[1] and $matepairs[$second]->[1]\n"; 
				print STDERR "------------------------------------------------\n";
				# Sweep through all mate-pairs and stop until a di-similar one is found
				while ($second < @matepairs && similar($qali[$matepairs[0]->[0]],
													   $qali[$matepairs[$second]->[0]],
													   $qali[$matepairs[0]->[1]],
													   $qali[$matepairs[$second]->[1]]))
 	 		     {
 	 			   $second++; 
			     }
				print STDERR "second=$second\n";
				if ($second < @matepairs)
 	 		      {
					  # let pass only best mate-pair, and only if $second is significantly worse
					  my $ratio = $matepairs[$second]->[2] / $matepairs[0]->[2];
					  # if ($verbose) 
					  #  {
						   print STDERR "\nbest two mates\n";
					       # Pair at position 0 (top position)
					       my @mate00 = split /\t/, $qali[$matepairs[0]->[0]]->[0], 21;
						   my @mate01 = split /\t/, $qali[$matepairs[0]->[1]]->[0], 21;
					       # Pair at position 1 (top-1 position)
					       my @mate20 = split /\t/, $qali[$matepairs[1]->[0]]->[0], 21;
					       my @mate21 = split /\t/, $qali[$matepairs[1]->[1]]->[0], 21;

						   print STDERR "" . $mate00[9] . " paired with " . @mate01[9] . 
						                ", score=$matepairs[0]->[2]\n";
						   print STDERR "" . $mate20[9] . " paired with " . @mate21[9] . 
						                ", score=$matepairs[1]->[2]\n";
					       print STDERR "Ratio between these two pairs: $ratio; uniqthresh=$uniqthresh\n";
					   # } 
					  if ($ratio < $uniqthresh)
					  {
						  # Let pass the two alignments corresponding to best-scored pair
						  print STDERR "ratio=$ratio<uniqthresh=$uniqthresh, passing unique " .
						               "mate pair: $mate00[9] and $mate01[9].\n";
						  print $qali[$matepairs[0]->[0]]->[0];
						  print $qali[$matepairs[0]->[1]]->[0];
					  } else {
						  # print STDERR "All mate-paired alignments dropped if the score between highest 
						  # and lowest mate-pairs doesn't satisfy $uniqthresh
						  print STDERR "Clearing the contents of \@matepairs()\n";
						  @matepairs = ();
					  } 
				  } else {# Implies: ($second == @matepairs), implies all pairs in @matepairs are similar
					  # if ($verbose)
					  #   {
						  print STDERR "suboptimal mate pairs are similar\n";
					    # }
					  # Letting pass only best mate pair
					  print $qali[$matepairs[0]->[0]]->[0];
					  print $qali[$matepairs[0]->[1]]->[0];
				  }
				print STDERR "------------------------------------------------\n";
				print STDERR "\@matepairs before SPLICING\n";
				printMatedPairsInfo(\@matepairs, \@qali);
				splice @matepairs, 1; # keep only the best pair (if any)
				print STDERR "\@matepairs after SPLICING\n";
				printMatedPairsInfo(\@matepairs, \@qali);
				print STDERR "------------------------------------------------\n";

			} else { # ($uniq, $best) = (0,1)
				my $optscore = $matepairs[0]->[2];
				my @bestTnames = ();
				my $numbest = 0;
				print STDERR "--------------------------------------------\n";
				print STDERR "(uniq,best)=($uniq,$best) selected\n";
				print STDERR "Letting pass all mate pairs that share maximum scoreMate\n";
				while ($numbest < @matepairs && $matepairs[$numbest]->[2] == $optscore)
				{
          			# Pair at position 0 (top position)
					my @mate00 = split /\t/, $qali[$matepairs[$numbest]->[0]]->[0], 21;
					my @mate01 = split /\t/, $qali[$matepairs[$numbest]->[1]]->[0], 21;

					print STDERR "Letting pass mate-pair (best): " . $mate00[9] . 
					             " paired with " . @mate01[9] . ", score=$matepairs[0]->[2]\n";
					print $qali[$matepairs[$numbest]->[0]]->[0];
					print $qali[$matepairs[$numbest]->[1]]->[0];
					print STDERR "Storing at \@bestTnames: $qali[$matepairs[$numbest]->[0]]->[1]\n";
					push @bestTnames, $qali[$matepairs[$numbest]->[0]]->[1];
					$numbest++;
				}
				$outBest += @matepairs - $numbest;
				print STDERR "------------------------------------------------\n";
				print STDERR "\@matepairs before SPLICING\n";
				printMatedPairsInfo(\@matepairs, \@qali);
				splice @matepairs, $numbest; # keep only the first $numbest pairs
				print STDERR "\@matepairs after SPLICING\n";
				printMatedPairsInfo(\@matepairs, \@qali);
				print STDERR "------------------------------------------------\n";

				if (@bestTnames>1)
				{
					my %genenames = ();
					foreach my $Tname (@bestTnames) 
 				    { 
						$Tname =~ s/\.t\d+//; $genenames{$Tname}=1; 
 				    }
					print COMMON $oldqnamestem . "\t" . join(" ", keys %genenames) . "\n" if (%genenames > 1 && defined($commongenefile));
 			    } 
			}#if ($uniq)
		}# if ((!$uniq && !$best) || @matepairs<2)

		# Whether "uniq, best" or "nothing" has been selected along with "paired", the thing is that 
		# the code saves to file the "pairbed" information
		# output pairedbed info: go through list of all mate pairs and store start and end position
		if (defined($pairbedfile))
		{
			print STDERR "---------------------Start: Pair bed section ------------------------\n";
			print STDERR "The size of \@matepairs=". scalar(@matepairs)."\n";
			print STDERR "---------------------End: Pair bed section ------------------------\n";
			while (@matepairs>0)
			{
				my $chr = $qali[$matepairs[0]->[0]]->[1]; #Extract target name from first mate pair
				$paircovsteps{$chr} = [] if (!defined($paircovsteps{$chr}));
				my $pend = $qali[$matepairs[0]->[1]]->[5]; # Extract tend position of Second mate
				my $pstart = $qali[$matepairs[0]->[0]]->[4]; # Extract tstart position fo first mate
				push @{$paircovsteps{$chr}}, [$pstart-1, 1];
				push @{$paircovsteps{$chr}}, [$pend, -1];
				shift @matepairs;
			}
		}

		
    } else { # IF NOT PAIRED, single read
		if (($uniq || $best) && @qali>1)
		{
			my %rscores;
			foreach my $ali (@qali)
			{
				$rscores{$ali} = scoreAli($ali); # store scores in hash so later sorting is faster
			}
			#### Printing alignment info before sorting
			print STDERR "------------------------------------------------------\n";
			print STDERR "\@qali BEFORE sorting by score:\n";
			printQaliWithScores(\@qali, \%rscores);
			@qali = sort {-($rscores{$a} <=> $rscores{$b})} @qali; # Sorting alignemtns in descending order
			#### Printing alignment info after sorting
			print STDERR "------------------------------------------------------\n";
			print STDERR "\@qali AFTER sorting by score:\n";
			printQaliWithScores(\@qali, \%rscores);
			print STDERR "------------------------------------------------------\n";

			if ($uniq) 
			{
				my $second = 1;
				# Sweep through all alignments and stop until a di-similar one is found
                while ($second < @qali && similar($qali[0], $qali[$second]))
				{
                    $second++;
                }
				print STDERR ">>> Second = $second\n";
				# The $second alignment is the one with lowest $score, but that is still similar
                if ($second < @qali) 
				{
					# let pass only best alignment, and only if $second is significantly worse
					my $ratio = $rscores{$qali[$second]}/$rscores{$qali[0]};
					# if ($verbose)
					# {
						print STDERR "Comparing scores between:\n";
						print STDERR $qali[0]->[0]."percid=$qali[0]->[6],coverage=$qali[0]->[7],score=$rscores{$qali[0]} AND \n";
						print STDERR $qali[1]->[0]."percid=$qali[1]->[6],coverage=$qali[1]->[7],score=$rscores{$qali[1]}\n";
						print STDERR "ratio between these best two alignments = $ratio\n";
					# }
					if ($ratio < $uniqthresh) # ... significantly worse in terms of "scoreAli"
					{
						# Letting pass only the alignment with the top score
						print $qali[0]->[0];
						$outUniq += @qali-1;
					} else {
						# drop out all alignments if the score between highest and lowest 
						# alignment doesn't satisfy $uniqthresh
						$outUniq += @qali;  
						print STDERR "All queries (". scalar(@qali) . " alignments): $oldqnamestem " . 
						             "filtered out by uniqueness ratio= $ratio \> uniqthresh=" .
									 "$uniqthresh\n";
					}
				} else {# Implies: ($second == @qali), which implies all alignments in @qali are similar
					print "suboptimal alignments are similar\n" if ($verbose);
					# Letting pass only best alignment
					print $qali[0]->[0];  
					$outUniq += @qali-1;
				}
			} else { # ($unique, $best) = (0,1) 
				print "\%rscores has contents\n";
				printRscores(\%rscores);
				my $optscore = $rscores{$qali[0]}; 
				print "-------------------------------------------------------------------------\n";
				print "optscore = $optscore:@qali[0]->[0]:percid=@qali[0]->[6]:coverage=@qali[0]->[7]\n";
				print "-------------------------------------------------------------------------\n";
				my @bestTnames = ();
				while ($rscores{$qali[0]} == $optscore)
				{
					print $qali[0]->[0]; # Let pass all alignments that share maximum score
					print STDERR "Let pass: @qali[0]->[0]:percid=@qali[0]->[6]:coverage=@qali[0]->[7]\n";
					push @bestTnames, $qali[0]->[1]; # Store target names in an array, for future use
					shift @qali; # "Delete" passed alignment and move to the next one
				}
				$outBest += @qali; # Number of filtered out alignments is the size of 
				                   # remaining alignments in @qali
				print "Size of \@qali is " . scalar(@qali) . "\n";
				my @myQuery = ();
				for (my $it=0; $it<@qali; $it++)
				{
					@myQuery = split /\t/, $qali[$it]->[0], 21;
					print STDERR "$myQuery[9] filtered out by best criterion: (percId & coverage)\n";
				}
				my @myQuery = ();
				for (my $it=0; $it<@qali; $it++)
				{
					@myQuery = split /\t/, $qali[$it]->[0], 21;
					print STDERR "$myQuery[9] filtered out by best, alignment score=" . 
					"$rscores{$qali[$it]} \< optimalScore=$optscore\n";
				}

				if (@bestTnames>1)
				{
					my %genenames = ();
					foreach my $Tname (@bestTnames) 
					{ 
						#$Tname = $Tname."t123";
						#print STDERR "Values Tname= $Tname, before regexp\n";
						$Tname =~ s/\.t\d+//; #### REPLACE anything of the sort "chr12.t123" by "chr12" 
						#print STDERR "Values Tname= $Tname, after regexp\n";
						$genenames{$Tname}=1; 
					}
					# Prints to commongenefile the names of the common genes
					if (%genenames>1 && defined($commongenefile))
					{
						print COMMON $oldqnamestem . "\t" . join(" ", keys %genenames) . "\n";
					}
				} # end if (@bestTnames>1)

			} # end of if ($unique)

		} else { # This else corresponds to: "if !(($uniq || $best) && @qali>1)"
			  print STDERR "------------------------------------------------\n";
			  print STDERR "(paired,uniq,best)=($paired,$uniq,$best)\n";
			  print STDERR "Letting pass all alignments= " . scalar(@qali) .", listed below:\n";
			  foreach my $ali (@qali)
			  {
			  	print $ali->[0] . "percid=$ali->[6], coverage=$ali->[7]\n"; 
			  }
			  print STDERR "------------------------------------------------\n";
		}
    } # End if ($paired)
    @qali = ();

	return (\%paircovsteps, $maxQaliSize);


}  # End processQuery()



sub computeMatePairs()
{
	#### Printing alignment info before sorting
	print STDERR "------------------------------------------------------\n";
	print STDERR "\@qali BEFORE sorting by \$tname and \$tStart:\n";
	printQali(\@qali);

	# Sorts @qali by $targetname (first criterion) and then by $tstart (second criterion)
	@qali = sort {$a->[1] cmp $b->[1] || $a->[4] cmp $b->[4]} @qali; 
	if (@qali > $maxQaliSize){$maxQaliSize = @qali;}

	#### Printing alignment info after sorting
	print STDERR "------------------------------------------------------\n";
	print STDERR "\@qali AFTER sorting by \$tname and \$tStart \n";
	printQali(\@qali);
	print STDERR "------------------------------------------------------\n";


	my @matepairs = ();
	my %mated = ();

	for (my $i=0;$i < @qali-1; $i++)
	{
		for (my $j=$i+1; $j < @qali && $qali[$i]->[1] eq $qali[$j]->[1]; $j++) 
		{			# only loop until leave chromosome
			my @mate1 = split /\t/, $qali[$i]->[0], 21;
			my @mate2 = split /\t/, $qali[$j]->[0], 21;
			# print STDERR "comparing [$i,$mate1[9],$qali[$i]->[1]], with [$j,$qali[$j]->[1], $mate2[9]]\n";
			print STDERR "comparing pair: $i,$j:[$mate1[9],$mate2[9]]=" . 
			"[$qali[$i]->[1] ,mate $qali[$i]->[2], strand $qali[$i]->[3]; ".
			"$qali[$j]->[1] ,mate $qali[$j]->[2], strand $qali[$j]->[3]]\n";
			if ($qali[$i]->[2] ne $qali[$j]->[2])
	 		       {# different mate: (f,r) or (1,2); given by $qsuffix
	 		         	   if ($qali[$i]->[3] ne $qali[$j]->[3])
	 		         	   {# different strand
			     		    # $qali[$j]->[4] is $tstart
			     		    # $qali[$i]->[5] is $tend
						      my $dist = $qali[$j]->[4] - $qali[$i]->[5] - 1;
						      $dist = $qali[$i]->[4]- $qali[$j]->[5]- 1 if ($qali[$i]->[4] > $qali[$j]->[4]);

							  print STDERR "Alignments $i,$j (above) with different mates, " . 
							        " different strands, dist=$dist, " . 
							        "maxintronlen=$maxintronlen\n";

						      if ($dist < $maxintronlen && $dist>=0)
							  { # not too far apart, not overlapping either
								  # print "found mate pair $i,$j\n";

								  # Storing mate-pair score (a function of "coverage and percid"
								  push @matepairs, [$i,$j,scoreMate($i,$j,$dist)];

								  $mated{$i}=0 if (!defined($mated{$i}));
								  $mated{$j}=0 if (!defined($mated{$j}));
								  $mated{$i}++;
								  $mated{$j}++;
								  my $inslen = $qali[$j]->[5] - $qali[$i]->[4] - 1;
								  $inslen = $qali[$i]->[5] - $qali[$j]->[4] - 1 if ($inslen<0);

								  print STDERR ">>>found mate pair $i,$j (above alignment)!!!\n"; 
								  # print STDERR "[$mate1[9],$mate2[9]]=" . 
								  # "[$qali[$i]->[1] ,mate $qali[$i]->[2], strand $qali[$i]->[3]; ".
								  # "$qali[$j]->[1] ,mate $qali[$j]->[2], strand $qali[$j]->[3]]\n";
								  print STDERR "mated[$i]=$mated{$i},mated[$j]=$mated{$j},inslen=$inslen,".
								  "dist=$dist,scoreMates=".scoreMate($i,$j,$dist)."\n";

								  # Storing insert length between mated pairs
								  push @insertlen, $inslen;

							  } else {
								  print STDERR "dist=$dist>maxintronlen=$maxintronlen || dist<0, " . 
								        "not right between [$i,$mate1[9],$qali[$i]->[1]],".
										" with [$j,$qali[$j]->[1], $mate2[9]]\n";
							  }

				      	   }		# end mid if
				      }			# end outer if

		} # end inner for
	} # end middle for


	print STDERR "Summary evaluation of mate-pairs\n";
	print STDERR "found " . scalar(@matepairs) . " mate pairs, involving " . 
	              scalar(keys %mated) . " mates\n";
	print STDERR "------------------------------------------------------\n";


	return(\@matepairs, \@insertlen, \%mated, $maxQaliSize);


} # end computeMatePairs()


# for comparing quality of two alignments
sub scoreAli
{
    my $ali = shift;
    return $ali->[6]/100 # percent identity
		+ $ali->[7]/100; # percent coverage
}


# for comparing quality of two mate pair read alignments (i1,j1), (i2,j2)
sub scoreMate()
{
	my $i = shift();
	my $j = shift();
	my $dist = shift();
    my $score = ($qali[$i]->[6] + $qali[$j]->[6])/100; # percent identity
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
sub similar 
{
    return similar(@_[0],@_[1]) && similar(@_[2],@_[3]) if (@_ == 4);
	# Alignment r
    my $r1 = shift;
	my @rali = @$r1;
	my @rArray = split /\t/, $rali[0], 21;
	my $rname = $rArray[9];
	my $rstart = $r1->[4];
	my $rend = $r1->[5];

	# Alignment s
    my $s1 = shift;
	my @sali = @$s1;
	my @sArray = split /\t/, $sali[0], 21;
	my $sname = $sArray[9];
	my $sstart = $s1->[4];
	my $ssend = $s1->[5];

    # if ($verbose)
	# {
	    print STDERR "\n[SIMILAR]: checking whether $rname and $sname are approx. the same:";
	    print STDERR "\n[SIMILAR]: Overlap: ($rname):$rend <= ($sname):$sstart";
	    print STDERR "\n[SIMILAR]: Overlap: ($rname):$rstart >= ($sname):$ssend\n";
    # }
##    if ($r1->[5] <= $s1->[4] || $r1->[4] >= $s1->[5])
    if ($rend <= $sstart || $rstart >= $ssend)
	{
		print STDERR "[SIMILAR]: Alignments are NOT SIMILAR\n";
		return 0; # here: similar = overlapping target range
	}
	print STDERR "[SIMILAR]: Alignments are SIMILAR\n";
    return 1;

}


#
# compactifyBed
# if several steps coincide then summarize them equivalently by one step in order to 
# 1) save memory or
# 2) output a bed file
#
sub compactifyBed 
{
    my $before=0;
    my $after=0;
    foreach my $chr (sort keys %paircovsteps)
	{
        next if (!@{$paircovsteps{$chr}});

		print STDERR "------------------------------------------------------\n";
		print STDERR "Pairedness coverage BEFORE sorting:\n";
		printChrOfPairCovSteps(\%paircovsteps, \$chr);

		# Foreach "$chr" in %paircovsteps, sort the contents of first column (target coords) in ascending order
        # {-1} corresponds to $pend, and {1} corresponds to $pstart
		@{$paircovsteps{$chr}} = sort {$a->[0] <=> $b->[0]} @{$paircovsteps{$chr}}; 

		print STDERR "------------------------------------------------------\n";
		print STDERR "Pairedness coverage AFTER sorting:\n";
		printChrOfPairCovSteps(\%paircovsteps, \$chr);
		print STDERR "------------------------------------------------------\n";

		$before += scalar(@{$paircovsteps{$chr}});
		my $i=0;
		print STDERR "-----------------------------------------\n";
		print STDERR "Before compactifying...\n";
		printChrOfPairCovSteps(\%paircovsteps, \$chr);
		while ($i<@{$paircovsteps{$chr}}-1)
		{
			if ($paircovsteps{$chr}->[$i]->[0] eq $paircovsteps{$chr}->[$i+1]->[0]) #[0] is the "position"
			{
				$paircovsteps{$chr}->[$i]->[1] += $paircovsteps{$chr}->[$i+1]->[1]; # add value from i+1 to i
				splice @{$paircovsteps{$chr}}, $i+1, 1; # remove element i+1

			} else {
				$i++;
			}
		}
		print STDERR "After compactifying...\n";
		printChrOfPairCovSteps(\%paircovsteps, \$chr);
		print STDERR "-----------------------------------------\n";
		$after += scalar(@{$paircovsteps{$chr}});
    }
    print STDERR "\nbefore compactifying: $before  after: $after\n";
}


# Auxiliary function that prints the contents of @qali arrray long with the associated 
# coverage and percId of each alignment
# Move to auxfunctions directory???
sub printQaliWithScores()
{
my $qaliRef = shift();
my $refRscores = shift();
my @qali = @$qaliRef;
my %rscores = %{$refRscores};

for (my $it=0; $it<@qali; $it++) 
  {
	my @align = split /\t/, $qali[$it]->[0], 21;
	print STDERR "$qali[$it]->[1]; $oldqnamestem; @align->[15]; @align->[16]; " .
	"percid=$qali[$it]->[6]; coverage=$qali[$it]->[7]; score=$rscores{$qali[$it]}\n";
  }

return 1;
}


sub printQali()
{
my $qaliRef = shift();
my @qali = @$qaliRef;

for (my $it=0; $it<@qali; $it++) 
  {
	my @align = split /\t/, $qali[$it]->[0], 21;
	print STDERR "$qali[$it]->[1]; $oldqnamestem; @align->[15]; @align->[16]; " . 
	"percid=$qali[$it]->[6]; coverage=$qali[$it]->[7];\n";
  }

}

sub printSizeOfCoverInfo()
{
	my $refPairCovSteps = shift();
	my %paircovsteps = %{$refPairCovSteps};

	foreach my $chr (sort keys %paircovsteps)
 	{
		print STDERR "Size of cover. info of chr=$chr is " . scalar(@{$paircovsteps{$chr}})."\n";
	}
}
