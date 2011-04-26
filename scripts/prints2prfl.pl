#!/usr/bin/perl -w
################################################################################
#
#   Name:    prints2prfl.pl
#   Project: Gene Prediction with Protein Family Patterns 
#   Author:  Oliver Keller
#   Date:    2011-01-03
#   Version: 0.1
#
#
#   prints2prfl.pl [<options>] <MSA> 
#
#   Converts a PRINTS database flat file into a set of block profiles
#   that can be used as input for AUGUSTUS-PPX.
#
#   Uses a BLOSUM matrix to calculate pseudocounts. Ensure that
#   the variable AUGUSTUS_CONFIG_PATH is set in order to make that work.
#

use strict;
use Getopt::Long;
use List::Util "sum";

my %data;

my $MIN_WIDTH = 6;
my $REG_MATR_FILE;
my $WEIGHTING_MODE="none";
my $DEFAULT_GLOB_WEIGHT=20;
my $GLOB_WEIGHT;
my $AMINO_ACIDS = "GDERKNQSTAVLIFYWHMCP";
my $QIJ_ORDER =   "ARNDCQEGHILKMFPSTWYV";	
my @AA_LIST = split "", $AMINO_ACIDS;
my $MAX_ENTROPY;
my $MIN_CONSERVE;
my $DISABLE_WEIGHTS;
my $AA_COUNT = length($AMINO_ACIDS);
my $LOG20 = log($AA_COUNT);
my $INITARGS = join(" ", @ARGV);


my $regmatrix;


### subroutines

#
# read in the qij (regularisation) matrix
#
sub read_qij {
    my $fh = shift;
    my @SM_AAS = ();
    my $matrix = {};
    my $sm_count=0;
    while (<$fh>) {
	chomp;
	s/\#.*$//;
	s/^\s*|\s*$//g;
	next if (/^$/);
	unless (@SM_AAS) {
	    my $FILE_AAS=uc;
	    $FILE_AAS =~ s/[^A-Z]//g;
	    my $patt = '^['.$FILE_AAS.']*$';
	    if (length($FILE_AAS)==20 && $AMINO_ACIDS  =~ /$patt/) {
		@SM_AAS = split "", $FILE_AAS;
		next;
	    } else {
		@SM_AAS = split "", $QIJ_ORDER;
	    }
	}
	my @line = split(/\s+/);
	if (@line <= $sm_count || @line > 20 || grep { !/^\d+\.\d*$|^\.\d+$/ } @line) {
	    print STDERR "Not a valid frequency file!\n";
	    return;
	}
	@{$matrix->{$SM_AAS[$sm_count]}}{@SM_AAS} = @line;
	$sm_count++;
	last if ($sm_count == 20);
    }

    # normalizing - the substitution matrix (SM) is scaled the following way:
    # - it is made symmetric by substituting SM_ij and SM_ji with (SM_ij+SM_ji/2)
    # - every entry is divided by the total sum of all entries
    my $sm_sum=0;
    foreach (0..$#SM_AAS) {
	my $aa1 = $SM_AAS[$_];
	foreach my $aa2 (@SM_AAS[0..$_-1]) {
	    my $ref1 = \$matrix->{$aa1}{$aa2};
	    unless (defined $matrix->{$aa2}{$aa1}) {
		$matrix->{$aa2}{$aa1} = $$ref1;
	    } else {
		my $ref2 = \$matrix->{$aa2}{$aa1};
		my $val = ($$ref1 + $$ref2)/2;
		$$ref1 = $$ref2 = $val;
	    }
	    $sm_sum += 2 * $$ref1;
	}
	$sm_sum += $matrix->{$aa1}{$aa1};
    }
    foreach (@SM_AAS) {
	foreach my $aa2 (@SM_AAS) {
	    $matrix->{$_}{$aa2} /= $sm_sum;
	}
	$matrix->{"*"}{$_} = sum(@{$matrix->{$_}}{@SM_AAS});
    }
    
    $regmatrix=$matrix;
}

sub calculate_pssm {
    #### change this according to msa2prfl
    my $blockref = shift;
    my ($size, $height) = 
	@$blockref{"width", "height"};
    my $result = { "lines" => [] };
    
    unless ($DISABLE_WEIGHTS) {
	$blockref->{weights} = [ (0) x $height ];
	foreach (@{$blockref->{columns}}) {
	    my $frq = $_->{frq};
	    my @res = grep { /^[$AMINO_ACIDS]$/ } keys %$frq;
	    my %colweights = map { ( $_ => 0 ) } keys %$frq;
	    @colweights{@res} = map { 1/$frq->{$_}/(scalar @res) } @res;
	    my $i=0;
	    foreach my $c (@{$_->{str}}) {
		$blockref->{weights}[$i++] += $colweights{$c};
	    }
	}
	$_ /= $blockref->{width} foreach (@{$blockref->{weights}});
    } 
    
    my @entropies;
    foreach (@{$blockref->{columns}}) {
	my $freq = $_->{frq};
	unless ($DISABLE_WEIGHTS) {
	    %$freq = ();
	    foreach my $i (0..$height-1) {
		$freq->{$_->{str}[$i]} += $blockref->{weights}[$i];
	    }
	} else {
	    $_ /= $height foreach values %$freq;
	    delete $freq->{$_} foreach grep {  /[^$AMINO_ACIDS]/ } keys %$freq;
	}	
 
	# substitute non-standard keys 
	# U -> [C]     O->[K]
	# B -> [DN]    J->[IL]
	# Z -> [EQ]    X->(any)
	foreach ( ["U","C"], ["B","D","N"], ["Z","E","Q"], 
		  ["O", "K"], ["J", "I", "L"], ["X", @AA_LIST] ) {
	    next unless exists $freq->{$_->[0]};
	    my $rem = $freq->{$_->[0]};
	    delete $freq->{shift @$_};
	    my @freqs = @{$regmatrix->{"*"}}{@$_};
	    my $factor = $rem / sum(@freqs);
	    foreach (@$_) {
		$freq->{$_} += (shift @freqs)*$factor;
	    }
	}
### DEBUG
	if (keys %$freq > $AA_COUNT) {
	    die "Have keys   '".join("", sort keys %$freq,)."'\nshould have '".join("", sort @AA_LIST)."'";
	}
	if (abs(sum (values %$freq) - 1) > 1e-4) {
	    die sprintf("col %d: sum is not 1 but %f ", $_->{num}, sum(values %$freq));
	} 
### 

	# calculate regularised counts (used as pseudo counts)
	my $res_count = scalar grep { $freq->{$_} != 0 } keys %$freq;
	my $regweight = $GLOB_WEIGHT / $height;
	$regweight *= $res_count if $WEIGHTING_MODE eq "blimps";

	my %reg = ();
	foreach my $aa (@AA_LIST) {
	    $reg{$aa} = 
		sum ( map { $freq->{$_} * $regmatrix->{$_}{$aa} / $regmatrix->{"*"}{$_} }
		      keys %$freq );
	}
	foreach my $aa (@AA_LIST) {
	    $freq->{$aa} += $regweight * $reg{$aa};
	    $freq->{$aa} /= (1 + $regweight);
	}

	if (defined $MAX_ENTROPY) {
	    my $entropy=0;
	    my $total=1;
	    foreach my $val (values %$freq) {
		$entropy -= $val * log($val);
		$total -= $val;
	    }
	    $entropy /= $LOG20;
	    if (abs($total) > 1e-4) {
		print STDERR "Warning: column $_->{num} does not add up to 1!";
		print STDERR "Deviation is $total.\n";
	    }
	    push @entropies, $entropy;
	}
    } # end for each column	
	    
    # mark entropy based cutoffs
    @{$blockref}{"from", "to"} = (0,$blockref->{width});
    if (defined $MAX_ENTROPY) {
	while ($blockref->{to} - $blockref->{from} >= $MIN_WIDTH) {
	    if ($entropies[$blockref->{from}]>$MAX_ENTROPY) {
		$blockref->{from}++;
	    } elsif ($entropies[$blockref->{to}-1]>$MAX_ENTROPY) {
		$blockref->{to}--;
	    } else {
		last;
	    }
	}

	# high-entropy blocks to be removed
	if (sum(@entropies[$blockref->{from}..$blockref->{to}-1]) > $MAX_ENTROPY*($blockref->{to}-$blockref->{from})) {
	    $blockref->{from} = $blockref->{to};
	}
    }
}	


############ main ###############

GetOptions("width=i" => \$MIN_WIDTH,
	   "qij=s" => \$REG_MATR_FILE,
	   "mode=s" => \$WEIGHTING_MODE,
	   "regweight=f" => \$GLOB_WEIGHT,
	   "max_entropy=f" => \$MAX_ENTROPY,
	   "min_conserve=f" => \$MIN_CONSERVE,
	   "noweights" => \$DISABLE_WEIGHTS
    );

unless (defined $GLOB_WEIGHT) {
    $GLOB_WEIGHT = $DEFAULT_GLOB_WEIGHT;
    $GLOB_WEIGHT /= 4 if $WEIGHTING_MODE eq "blimps";
}

if (!defined $MAX_ENTROPY && defined $MIN_CONSERVE) {
    $MAX_ENTROPY = 1 - $MIN_CONSERVE;
}


# directories for config files
# this is the method to be used for all AUGUSTUS-PPX scripts
my @CONFIG_DIRS = ("");
if (defined $ENV{"AUGUSTUS_CONFIG_PATH"}) {
    push (@CONFIG_DIRS, $ENV{"AUGUSTUS_CONFIG_PATH"}."/profile/");
}
if ($0 =~ /^(.*\/)/ && $1 ne "./") {
    push (@CONFIG_DIRS, $1);
}


# 
# look for regularisation matrix file
#
unless (defined $REG_MATR_FILE) {
    $REG_MATR_FILE="default.qij";
} elsif ($REG_MATR_FILE =~ s/^((~.*?|)\/)//) {
    # absolute directory given: do not check other dirs
    @CONFIG_DIRS=($1);
}
foreach (@CONFIG_DIRS) {
    my $fh;
    if (-e "$_$REG_MATR_FILE" && open $fh, "$_$REG_MATR_FILE") {
	&read_qij($fh);
	last;
    }
}
unless (defined $regmatrix) {
    print STDERR "Qij file '$REG_MATR_FILE' could not be read.\nCheck the parameter --qij=<filename>.\n\n";
    die;
}	

#
# read in PRINTS database flat file
#
my $name;
while (<>) {
    s/\s*$//;
    next unless /^(..)[;] (.*)/;
    my ($key, $value) = ($1, $2);
    if ($key eq "gc") {
	print STDERR "Found new entry '$value'\n";
	$name = $value; 
	$data{$name} = { blocks => [], ID => $name };
    } elsif ($key =~ /^g([xt])$/) {
	my $out = { "x" => "AC", "t" => "DE" }->{$1};
	$data{$name}{$out} = $value;
    } elsif ($key =~ /^f/) {
	my $blocklist = $data{$name}{blocks};
	my $currblock = $blocklist->[-1];
	if ($key eq "fc") {
	    push @$blocklist, {
		name => $value, 
		columns => [],
		height => 0
	    }
	} elsif ($key eq "fl")  {
	    $currblock->{width} = $value;
	    foreach (0..($value-1)) {
		push @{$currblock->{columns}}, { str=>[], num=>$_, frq=>{} };
	    }
	} elsif ($key eq "fd") {
	    my ($motif, $seqname, $pos, $dist) = split /\s+/, $value;
	    unless (defined $currblock->{width}) {
		$currblock->{width} = length($motif);
	    } elsif (length($motif) != $currblock->{width}) {
		die "ERROR: motif width mismatch. Aborting";
	    }
	    push @{$currblock->{ibr}}, $dist;
	    for (my $i=0; $i<$currblock->{width}; $i++) {
		my $c = substr($motif, $i, 1);
		my $col = $currblock->{columns}[$i];
		if (@{$col->{str}} != $currblock->{height}) {
		    die "ERROR: motif height mismatch. Aborting";
		}
		push @{$col->{str}}, $c;
		$col->{frq}{$c}++;
	    }
	    $currblock->{height}++;
	}
    }
}

foreach (keys %data) {
    my $dataset = $data{$_};
    my $filename = "$dataset->{AC}.prfl";
    open OUTFILE, ">$filename";
    my $blocklist = $dataset->{blocks};
    my $blcount = 0;
    my @ibr = ();
    my $ibradd = 0;
    print OUTFILE "[name]\n$dataset->{ID}\n";
    print OUTFILE "# $dataset->{DE}\n";
    foreach my $block (@$blocklist) {
	### HERE: define columns and weights (see msa2prfl)
	if (@ibr) {
	    foreach (0..$#ibr) {
		$ibr[$_] += $block->{ibr}[$_];
	    }
	} else {
	    @ibr = @{$block->{ibr}};
	}
	&calculate_pssm($block);
	my ($minprev, $maxprev);
	$minprev = $maxprev = $ibr[0];
	foreach (@ibr) {
	    if ($minprev > $_) {
		$minprev = $_;
	    } elsif ($maxprev < $_) {
		$maxprev = $_;
	    }
	}
	if ($minprev < -$block->{from}) {
	    $block->{from}=-$minprev;
	}
	if ($block->{to} - $block->{from} < $MIN_WIDTH) {
	    $ibradd += $block->{width};
	    next;
	}
	$_ += $block->{from} + $ibradd foreach ($minprev, $maxprev);
	$ibradd = $block->{width}-$block->{to};
	@ibr=();
	print OUTFILE "\n[dist]\n# distance from previous block\n";
	print OUTFILE "# <min> <max>\n";
	print OUTFILE "$minprev\t$maxprev\n\n";
	print OUTFILE "[block]\n# block no. $blcount follows, $block->{height} sequences, length $block->{width}\nname=$block->{name}\n";
	print OUTFILE "# note: PSSM column 0 corresponds to original block column $block->{from}\n" if $block->{from}>0;
	$blcount++;
	print OUTFILE "#\n# <colnr> <probs for $AMINO_ACIDS>\n";
	print OUTFILE "#\t".join("\t", @AA_LIST)."\n";
	for (my $colno = $block->{from}; $colno < $block->{to}; $colno++) {
	    my $freq = $block->{columns}[$colno]{frq};
	    print OUTFILE "".($colno-$block->{from})."\t".join("\t", map { sprintf(($_ < 1e-4 ? "%.2g" : "%7.5f"), $_) } @{$freq}{@AA_LIST})."\n";
	}
    }
    print OUTFILE "\n# created by:\n";
    print OUTFILE "# $0 $INITARGS\n";
    close OUTFILE;
    print "Wrote '$filename'.\n";
}
