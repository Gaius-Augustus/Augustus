#!/usr/bin/perl -w
################################################################################
#
#   Name:    block2prfl.pl
#   Project: Gene Prediction with Protein Family Patterns 
#   Author:  Oliver Keller
#   Date:    2011-01-07
#   Version: 0.2 (rev.374)
#
#
#   block2prfl.pl [<options>] <MSA> 
#
#   Converts a BLOCKS database flat file into a protein block profile
#   that can be used as input for AUGUSTUS-PPX.
#
#   Uses a BLOSUM matrix to calculate pseudocounts. Ensure that
#   the variable AUGUSTUS_CONFIG_PATH is set in order to make that work.
#

use strict;
use File::Basename;
use List::Util ("min","max","sum");
use POSIX;
use Getopt::Long;

my $DISABLE_WEIGHTS = 0;

my $AMINO_ACIDS = "GDERKNQSTAVLIFYWHMCP";
my $EXT_AAS = $AMINO_ACIDS."BZU";
my @AA_LIST = split "", $AMINO_ACIDS;

my $SHOW_ENTROPY = 0;
my $MIN_BLOCK_WIDTH = 6;

my ($MIN_CONSERVE, $MAX_ENTROPY);
my $LOG20 = log(scalar @AA_LIST);
my ($SKIPSTR, $USESTR, $name, $desc);
my $RELAX_AROUND_SKIPPED = 0; # whether to allow arbitrary distances around explicitly skipped blocks
my $RELAX_TERM; # allow arbitrary distance to sequence start and end

my $INITARGS = join(" ",@ARGV);
my $OWNSCOREFILE;

my $GLOB_WEIGHT = 20;
my $WEIGHTING_MODE = "none";

GetOptions("show_entropy" => \$SHOW_ENTROPY,
	   "max_entropy=f" => \$MAX_ENTROPY,
	   "min_conserve=f" => \$MIN_CONSERVE,
	   "skip=s" => \$SKIPSTR,
	   "useonly=s" => \$USESTR,
	   "relax_around_skipped" => \$RELAX_AROUND_SKIPPED,
	   "relax_terminal" => \$RELAX_TERM,
	   "setname=s" => \$name,
	   "setdesc=s" => \$desc,
	   "min_blockwidth=i" => \$MIN_BLOCK_WIDTH,
	   "ownscorefile=s" => \$OWNSCOREFILE,
	   "mode=s" => \$WEIGHTING_MODE,
	   "regweight=f" => \$GLOB_WEIGHT);

if (!defined $MAX_ENTROPY && defined $MIN_CONSERVE) {
    $MAX_ENTROPY = 1 - $MIN_CONSERVE;
}
if ($WEIGHTING_MODE eq "blimps") {
    $GLOB_WEIGHT /= 4;
}

my $SMfilename = $ENV{"AUGUSTUS_CONFIG_PATH"}."/profile/default.qij" if defined $ENV{"AUGUSTUS_CONFIG_PATH"};
$SMfilename = $ENV{"QIJ"} if (defined $ENV{"QIJ"});
# try "./default.qij" if QIJ is not specified 
$SMfilename = "default.qij" unless defined $SMfilename;
if (! -e $SMfilename && $0 =~ /^(.*\/)/ && -e "$1/$SMfilename") {
    $SMfilename="$1/$SMfilename";
}

my $fh;

open QIJFILE, $SMfilename or die "Could not open '$SMfilename'\n";

my %matrix;
my %marg;
my @SM_AAS=();
my $sm_count=0;

sub cut_sides { # rightcut MUST be <= 0 !!!
    my ($array, $leftcut, $rightcut) = @_; 
    splice(@$array, 0, $leftcut);
    splice(@$array, $rightcut) if $rightcut < 0;
}

sub block_to_pssm {
    my $blockref = shift;
    my ($size, $seq_count) = 
	@$blockref{"width", "num_sequences"};
    my $result = { "lines" => [] };

    my @probs;
    my $total_weight = 0;
    foreach (0..$size-1) {
	$probs[$_] = { map { ($_ => 0) } split("",  $EXT_AAS) };
    }
  
    foreach my $arrayref (values %{$blockref->{sequences}})
    {
	my ($seq, $weight) = @$arrayref[1,2];
#	print STDERR "sequence: $seq ($weight)\n";
	$weight = 1 / $seq_count if  (!defined $weight || $DISABLE_WEIGHTS);
	$total_weight += $weight;
	for (0..$size-1)
	{
	    my $key = substr($seq,0,1,"");
	    $probs[$_]{$key} += $weight;
	}
    }	

    my @entropies = ();
    foreach my $pos (0..$size-1) {
	my $total=1;
	foreach ( ["U","C"], ["B","D"], ["Z","E"] ) {
	    $probs[$pos]{$_->[1]} += $probs[$pos]{$_->[0]};
	}
	my $rem_weight = 0;
	foreach (grep(/[^$EXT_AAS]/,keys %{$probs[$pos]})) {
	    print  STDERR "Unknown char '$_' found in column $pos in block $blockref->{number}.\n";
	    $rem_weight += $probs[$pos]{$_};
	}
	if ($rem_weight > 0) {
	    printf STDERR "Equally distributing a fraction of %5.3f\n", $rem_weight / $total_weight;
	    $_ += $rem_weight/20   foreach (@{$probs[$pos]}{@AA_LIST}) 
	}
	
	# calculate regularised counts (used as pseudo counts)
	my %reg = ();
	foreach my $aa (@AA_LIST) {
	    $reg{$aa} = 
		sum ( map { $probs[$pos]{$_} * $matrix{$_}{$aa} / $marg{$_} }
		      @AA_LIST );
	}
	my $res_count = grep { $probs[$pos]{$_} != 0 } @AA_LIST;

	my $reg_weight = $GLOB_WEIGHT / $seq_count;
	$reg_weight *= $res_count if $WEIGHTING_MODE eq "blimps";


	my $entropy = 0;
	my $line = "";
	foreach (@AA_LIST)
	{
	    my $val = ($probs[$pos]{$_} + $reg_weight*$reg{$_}) / $total_weight / (1+$reg_weight);
	    $entropy -= $val * log($val);
	    $total-=$val;
	    $line .= sprintf ($val < 1e-4 ? "\t%.2g" : "\t%.5f", $val);
	}
	$entropy /= $LOG20;
	if (abs($total) > 1e-4) {
	    $line .= " # Rounding error: $total" ;
	    print STDERR "Warning: line $pos in block $blockref->{number} does not add up to 1!\n";
	    print STDERR "Deviation is $total\n";
	}
	push @entropies, $entropy;
	push @{$result->{lines}}, $line;
    } 
    my ($leftcut, $rightcut) = (0, 0);
    if (defined $MAX_ENTROPY) {
	while ($leftcut < @entropies && $entropies[$leftcut] > $MAX_ENTROPY) {
	    $leftcut++;
	}
	while ($rightcut > $leftcut-@entropies && $entropies[$rightcut-1] > $MAX_ENTROPY) {
	    $rightcut--;
	}
	if (($leftcut - $rightcut) > 0) {
	    &cut_sides(\@entropies, $leftcut, $rightcut);
	    if (@entropies < $MIN_BLOCK_WIDTH || sum(@entropies)/@entropies > $MAX_ENTROPY) {
		$leftcut = $size;
		$size = $rightcut = 0;
		@entropies = ();
		$result->{lines} = [];
		print STDERR "Removed block $blockref->{number} from profile (entropy too high)\n";
	    } else {
		$size -= ($leftcut - $rightcut);
		&cut_sides($result->{lines}, $leftcut, $rightcut);
		printf STDERR ("Removed %d columns (%d+%d) from block $blockref->{number} (entropy too high)\n", $leftcut - $rightcut, $leftcut, -$rightcut);
	    }
	}
    }
    $result->{entropies} = [@entropies];
    $result->{size} = $size;
    @$result{"left_cutoff", "right_cutoff"} = ($leftcut, -$rightcut);
    return $result;
}

sub read_a_block {
    my $result = {};
    my $status = 0;
    while (<>) {
	last if /^\/\/$/;
	if ($status == 0) {
	    if (/^ID   (.+); BLOCK$/) {
		$result->{id}=$1;
		$status = 1;
	    }
	} elsif ($status == 1) {
	    if (s/; distance from previous block=\((-?\d+),(\d+)\).*$//) {
		@{$result}{qw/min_prev max_prev/} = ($1, $2);
	    }
	    if (/^AC   (.+)/) {
		$result->{number} = $1;
	    } elsif (/^DE   (.+)$/) {
		$result->{de}=$1;
	    } elsif (/^BL   (.+)width=(\d+);? seqs=(\d+)/)  {
		@{$result}{qw/motif width num_sequences/} = ($1, $2, $3);
		$status = 2;
	    }
	} elsif ($status == 2) {
	    next unless /^ *(\S+) +\( *(-?\d+)\) (\S*) *([\d|\.]+)?.*$/;
	    my $seq_name = "$1:$2";
	    next if (exists $result->{sequences}{$seq_name});
	    push (@{$result->{sequences}{$seq_name}}, $2, $3) ;
# save sequence weight if it was found
	    if (defined $4){ 
		push (@{$result->{sequences}{$seq_name}}, $4) ;
	    }
	}
    }
    return $result;
}
     
    
    
# Reading in substitution matrix Q_ij
while (<QIJFILE>) {
    chomp;
    s/\#.*$//;
    s/^\s*|\s*$//g;
    next if (/^$/);
    unless (@SM_AAS) {
	@SM_AAS=split (/\s+/, uc, 20);
	if (@SM_AAS == 20 ) {
	    my $patt = '^['.join("",@SM_AAS).']*$';
	    next if ($AMINO_ACIDS  =~ /$patt/);
	}
	@SM_AAS = @AA_LIST;
    }
    my @line = split(/\s+/);
    if (@line <= $sm_count || @line > 20 || grep { !/^\d+\.\d*$|^\.\d+$/ } @line) {
	die "Not a valid frequency file";
    }
    @{$matrix{$SM_AAS[$sm_count]}}{@SM_AAS} = @line;
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
	my $ref1 = \$matrix{$aa1}{$aa2};
	unless (defined $matrix{$aa2}{$aa1}) {
	    $matrix{$aa2}{$aa1} = $$ref1;
	} else {
	    my $ref2 = \$matrix{$aa2}{$aa1};
	    my $val = ($$ref1 + $$ref2)/2;
	    $$ref1 = $$ref2 = $val;
	}
	$sm_sum += 2 * $$ref1;
    }
    $sm_sum += $matrix{$aa1}{$aa1};
}
foreach (@SM_AAS) {
    foreach my $aa2 (@SM_AAS) {
	$matrix{$_}{$aa2} /= $sm_sum;
    }
    $marg{$_} = sum(@{$matrix{$_}}{@SM_AAS});
}
#@{$marg}{@SM_AAS} = map { sum(@{$matrix{$_}}{@SM_AAS}) } @SM_AAS;
 

my $blkcount=0;
my ($totalminp, $totalmaxp) = (0,0);
my $endwithempty;

my $blset_inverted = !defined $USESTR;
my %blset;
if (!$blset_inverted) {
    $blset{$_} = 1 foreach (split(/,/, $USESTR));
}
if (defined $SKIPSTR) {
    my @skiplist = split (/,/, $SKIPSTR);
    if ($blset_inverted) {
	$blset{$_} = 1 foreach (@skiplist);
    } else {
	delete $blset{$_} foreach (@skiplist);
    }
}

my $relax = $RELAX_TERM;
until (eof()) {
    my $blockref = read_a_block();
    next unless defined $blockref;
    my ($minp, $maxp)
	= @$blockref{"min_prev", "max_prev"};
    $minp = max (0, $minp);
    if (defined $maxp) {
	$maxp = max ($minp, $maxp);
    } else {
	$maxp="*";
    }

    # if the last block is empty, we assume it to be the sequence end
    # in this case we finish after printing the distance to the previous block
    $endwithempty = $blockref->{width} == 0;

    $totalminp += $minp;
    if ($totalmaxp ne "*") {
	if ($maxp eq "*" ) {
	    $totalmaxp = "*";
	} else {
	    $totalmaxp += $maxp;
	}
    }
    
    if ($blset_inverted == exists $blset{$blockref->{number}}) {
	# this block is skipped explicitly
	# in case we remove the last block, we treat it as empty
	$totalminp += $blockref->{width};
	$totalmaxp += $blockref->{width} unless $totalmaxp eq "*";
	$relax = $RELAX_AROUND_SKIPPED unless $relax;
	next;
    }

    my $pssm = &block_to_pssm($blockref);
    $totalminp += $pssm->{left_cutoff};
    $totalmaxp += $pssm->{left_cutoff} unless $totalmaxp eq "*";

    unless (defined $name) {
	$name = $blockref->{id};
	$name = "" unless defined $name;
    } elsif ($name && defined $blockref->{id} && $blockref->{id} ne $name) {
	print STDERR "Block name is '".$blockref->{id}."', while pattern will be named '$name'\n";
    }

    # do not print anything if the block was deleted due to high entropy
    next if ($pssm->{size} == 0);

    if ($name) {
	$desc = $blockref->{de} unless defined $desc;
	if (defined $desc) {
	    $desc = "# $desc\n";
	} else {
	    $desc = "";
	}
	print "[name]\n$name\n$desc";
	$name = "";
    }
    
    print "\n[dist]\n# distance from previous block\n# <min> <max>\n";
    print ($relax ? "0\t*\n" : "$totalminp\t$totalmaxp\n");
    print
	"\n[block]\n# block no. ".($blkcount++)." follows, $blockref->{num_sequences} sequences, length $pssm->{size}\n";
    my $blockname = $blockref->{number};

    print "name=$blockname\n" if defined $blockname;
    if ($pssm->{left_cutoff}) {
	print "# note: PSSM column 0 corresponds to original block column $pssm->{left_cutoff}\n";
    }
    
    print "#\n".
	"# <colnr> <probs for $AMINO_ACIDS>\n".
	"#	G	D	E	R	K	N	Q	S	T	".
	"A	V	L	I	F	Y	W	H	M	C	P\n";
    foreach my $pos (0..$pssm->{size}-1) {
	print $pos.$pssm->{lines}[$pos]."\n";
	if ($SHOW_ENTROPY) {
	    printf ("# conservedness in previous line (1-normed entropy): %.3f\n", 1-$pssm->{entropies}[$pos]);
	}
    }
    
    ($totalmaxp, $totalminp) = ($pssm->{right_cutoff}) x 2;
    $relax=0;
}

if (($totalminp>0 || $endwithempty) && !$RELAX_TERM) {
    $totalmaxp = "*" unless $endwithempty;
    print "\n[dist]\n# distance from previous block\n# <min> <max>\n$totalminp\t$totalmaxp\n";
}

print "\n# created by:\n";
print "# $0 $INITARGS\n";
