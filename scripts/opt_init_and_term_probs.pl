#!/usr/bin/perl

#############################################################
# opt_init_and_term_probs.pl
# optimally train the initial and terminal probabilities
#
# usage: opt_init_and_term_probs.pl
#
#
#
# Mario Stanke, 25.04.2007
#############################################################

use strict;
use IO::File;

# parameters
my $opt_trans_matrix = "/home2/mstanke/augustus/config/species/nt/nt_trans_shadow_partial_utr.pbl";
#my $fafile = "/home2/mstanke/trunc/traindata/test.fa";
my $fafile = "/home2/mstanke/trunc/traindata/train/mopt.train.and.random.fa";
my $gfffile = "/home2/mstanke/trunc/traindata/train/mopt.train.CDS.gff";

my $n=71; # number of states
my $rounds = 80;
my (@initProbs, @termProbs);

##############################################################
# read in initProbs and termProbs
##############################################################

open(TRANS, "<$opt_trans_matrix") or die ("Could not open transition matrix file $opt_trans_matrix");
my $initFlag=0;
my $termFlag=0;
my @transfilelines;
while (<TRANS>) {
    push @transfilelines, $_;
    if (/^\s*\[Initial\]\s*/){
	$initFlag=1;
	$termFlag=0;
    } elsif (/^\s*\[Terminal\]\s*/){
	$initFlag=0;
	$termFlag=1;
    }
    if (/^\s*\[Transition\]\s*/){
	$initFlag=$termFlag=0;
    }
    if ($initFlag && /^\s*(\d+)\s+(\S+)/){
	@initProbs[$1]=$2;
    } 
    if ($termFlag && /^\s*(\d+)\s+(\S+)/){
	@termProbs[$1]=$2;
    }
}
close TRANS;
print "Initial probabilities: " . (join " ", @initProbs) . "\n";
print "Terminal probabilities: " . (join " ", @termProbs) . "\n";



#######################################################################################
# initialization
#######################################################################################
my $optscore = evalscore();
print "Initial run: score=$optscore\n";


#######################################################################################
# optimization loop for initial and terminal probabilities
#######################################################################################

my $v; #reference to the initial or terminal probability vector
my (@tryi,@tryt);
my $initRound;
for (my $r=0; $r<$rounds; $r++) {
    my $found_improvement = 0;
    foreach $v (\@initProbs, \@termProbs) {
	$initRound = ($v == \@initProbs);
	print "Improving " . ($initRound? "initial":"terminal") . 
	    " probabilities. Currently: " . join (" ", @{$v}) . "\n";
	# make a list with all the varied probability vectors to try
	my @tryvectors = getVariedTransVectors($v, 1.0, 1);
	print "** Trying " . scalar(@tryvectors) . " variations of vector.\n";

	foreach my $varieddist (@tryvectors) {
	    print "Try varied distribution " . join (" ", @{$varieddist}) . " ";
	    if ($initRound) {
		@tryi = @{$varieddist};
		@tryt = @termProbs;
	    } else {
		@tryt = @{$varieddist};
		@tryi = @initProbs;	
	    }
	    save_vectors(\@tryi, \@tryt, $opt_trans_matrix);
	    my $score = evalscore();
	    print "current score=$score\n";
	    if ($score > $optscore) {
		print "************ Found improvement\n";
		$optscore = $score;
		# score current value
		if ($initRound) {
		    @initProbs = @tryi;
		    print "New initProbs: " . join (" ", @tryi) . "\n";
		} else {
		    @termProbs = @tryt;
		    print "New termProbs: " . join (" ", @tryt) . "\n";
		}
		save_vectors(\@initProbs, \@termProbs, $opt_trans_matrix . ".opt");
	    }
	}
#	if (!$found_improvement && $r < $rounds-1) {
#	    print "Could not further improve. Skipping last " . ($rounds-$r-1) ." rounds\n";
#	    last;
#	}
    }
}


#######################################################################################
# subroutines
#######################################################################################

################################################
# evalscore: determine the score
# belonging to the current parameters
################################################

my %storedsnsp = {}; # hash with the stored sn and sp array references

sub evalscore {
    my $score;
    
    my $output = qx(truncateAndTest.pl --seq=$fafile --gff=$gfffile --gffOut=temp.out.gff);
    $output =~ /score =\s*(\S+)/;
    $score = $1;
    return $score;
}


####################################################################
# save_vectors: replace existing start and end probabilities in the 
# transition probabilities file
####################################################################
sub save_vectors {
    my $iv = shift;
    my $tv = shift;
    my $filename = shift;
    my $initFlag=0;
    my $termFlag=0;
    open (TRANS, ">$filename") or die ("Could not open $filename for writing.");
    foreach my $line (@transfilelines) {
	if ($line =~ /^\s*\[Initial\]\s*/){
	    $initFlag=1;
	    $termFlag=0;
	} elsif ($line =~ /^\s*\[Terminal\]\s*/){
	    $initFlag=0;
	    $termFlag=1;
	} elsif ($line =~ /^\s*\[Transition\]\s*/){
	    $initFlag=$termFlag=0;
	}
	if ($initFlag && $line =~ /^(\s*)(\d+)(\s+)(\S+)(.*)/){
	    print TRANS $1 . $2 . $3 . $iv->[$2] . $5 . "\n";
	} elsif ($termFlag && $line =~ /^(\s*)(\d+)(\s+)(\S+)(.*)/){
	    print TRANS $1 . $2 . $3 . $tv->[$2] . $5 . "\n";
	} else {
	    print TRANS $line;
	}
    }
    close TRANS;
}


#
# norm a vector to sum up to $normsum
# if $normed is true
sub norm {
    my ($vecref, $normsum, $normed) = @_;
    if ($normed) {
	my $sum=0;
	foreach my $item (@{$vecref}){
	    $sum += $item;
	}
	if ($sum > 0 && $normsum > 0) {
	    foreach my $item (@{$vecref}){
		$item *= $normsum/$sum;
	    }
	}
    }
    # round to 6 places
    foreach my $item (@{$vecref}){
	$item = sprintf("%.6f", $item);
	$item =~ s/(\.\d*[1-9])0+$/$1/; # remove trailing zeros
    }
}

#
# getVariedTransVectors
# parameters: \@transvec, $normsum, $normed
# returns a list of references to varied transition probability vectors
# (this could be more fancy later, possibly with randomness)
#
sub getVariedTransVectors {
    my $transvec = shift;
    my $normsum = shift;
    my $normed = shift;
    my @tryvectors = (); #list of references to transition probability vectors
    my $factor;

    # randomly change a random state
    my $k= int(rand(@{$transvec}));
    if (rand(6) < 1) {
	$k=0;
    }
    for (my $v=0; $v<5; $v++){
	$factor = 1+3*rand();
	if (rand(2)<1){
	    $factor=1/$factor;
	    print "decreasing ";
	} else {
	    print "increasing ";
	}
	print " prob of state $k by a factor of $factor\n";
	# modify element k
	my @vartransvec = @{$transvec};
	$vartransvec[$k] = $transvec->[$k] * $factor;
	norm(\@vartransvec, $normsum, $normed);
	push @tryvectors, \@vartransvec;
    }

    return @tryvectors;
}

#
# roundVector
#
sub roundVector {
    my $sum=0;
    foreach my $v (@_) {
	$sum += $v;
	$v = sprintf("%.6f", $v);
    }
    my $roundedsum = sprintf("%.6f", $sum); # is usually 1
    my $newsum=0;
    my $maxel=0;
    for (my $i=0; $i < @_; $i++) {
	$newsum += $_[$i];
	if ($_[$i]> $_[$maxel]) {
	    $maxel = $i;
	}
    }
    # in case the rounding changed the sum a little, adjust the largest element only
    if ($newsum != $roundedsum){
	@_[$maxel] += $roundedsum-$newsum;
    }
}
