#!/usr/bin/perl -w
################################################################################
#
#   Name:    msa2prfl.pl
#   Project: Gene Prediction with Protein Family Patterns 
#   Author:  Oliver Keller
#   Date:    2011-01-07
#   Version: 0.2 (rev.374)
#
#
#   msa2prfl.pl [<options>] <MSA> 
#
#   Converts a Multiple Sequence Alignment in FASTA or CLUSTAL format into
#   a protein block profile that can be used as input for AUGUSTUS-PPX.
#
#   Uses a BLOSUM q_ij matrix to calculate pseudocounts which is searched for
#   in a subdirectory of $AUGUSTUS_CONFIG_PATH. Ensure that the environment 
#   variable AUGUSTUS_CONFIG_PATH is set.
#
#   Options:
#   --width=i              minimal block width (default: 6)
#   --qij=s                alternative file containing BLOSUM q_ij matrix
#   --max_entropy=f        maximal entropy of a block column (default: disabled)
#   --keep_empty           do not remove empty columns from alignment
#   --prefix_from_seqnames if MSA contains partial sequences, determine
#                          sequence offset from sequence name as in XXXX/<from>-<to>
#   --relax                consider MSA as partial (allow arbitrary distance at 
#                          beginning and end)
#   --blockscorefile=s     create a log file with blocks from the MSA and their scores
#   --setname=s            set a name for the profile
#   --setdesc=s            set a description for the profile
#   --setacc=s             set an accession id for the profile
#   --info=s               provide a text file containing the meta information

use strict;
use Getopt::Long;
use List::Util qw/sum min max/;

my @seqnames;
my %sequences;

my $MIN_WIDTH=6;
my $REG_MATR_FILE;
my $WEIGHTING_MODE="none";
my $DEFAULT_GLOB_WEIGHT=20;
my $GLOB_WEIGHT;
my $AMINO_ACIDS = "GDERKNQSTAVLIFYWHMCP";
my $QIJ_ORDER =   "ARNDCQEGHILKMFPSTWYV";	
my @AA_LIST = split "", $AMINO_ACIDS;
my $SETNAME;
my $SETDESC;
my $SETACC;
my $INFOFILE;
my $KEEP_EMPTY;
my $USE_MSA_DIST;
my $MAX_ENTROPY;
my $MIN_CONSERVE;
my $DISABLE_WEIGHTS;
my $PREFIX_FROM_SEQNAMES;
my $RELAX;
my $AA_COUNT = length($AMINO_ACIDS);
my $LOG20 = log($AA_COUNT);
my $INITARGS = join(" ",@ARGV);
my $blockscorefile;

my %backcol;
@backcol{@AA_LIST} = ( 0.07088, 0.05268, 0.06270, 0.05256, 0.05807, 
		       0.04439, 0.04037, 0.07068, 0.05837, 0.07689, 
		       0.06538, 0.09200, 0.05527, 0.03995, 0.03224,
		       0.01312, 0.02253, 0.02353, 0.01793, 0.05046 );
my $minFreq = 0.0001;
if (defined $ENV{"minFreq"}) {
    $minFreq = $ENV{"minFreq"};
}

# set this to 1 if the same behavious as block2prfl.pl is desired
my $COMPATIBILITY_MODE = 0;
if ($COMPATIBILITY_MODE) {
    $WEIGHTING_MODE="blimps";
}

# note: if this is changed, 
#       block dist calculation has to be modified
my $MAX_GAP_COUNT = 0;

### subroutines


sub usage {
    print STDERR "
Usage:

   msa2prfl.pl [<options>] <MSA> 

   Converts a Multiple Sequence Alignment in FASTA or CLUSTAL format into
   a protein block profile that can be used as input for AUGUSTUS-PPX.

   Uses a BLOSUM q_ij matrix to calculate pseudocounts which is searched for
   in a subdirectory of \$AUGUSTUS_CONFIG_PATH. Ensure that the environment 
   variable AUGUSTUS_CONFIG_PATH is set.

Options:
   --width=i              minimal block width (default: 6)
   --qij=s                alternative file containing BLOSUM q_ij matrix
   --max_entropy=f        maximal entropy of a block column (default: disabled)
   --keep_empty           do not remove empty columns from alignment
   --prefix_from_seqnames if MSA contains partial sequences, determine
                          sequence offset from sequence name as in XXXX/<from>-<to>
   --relax                consider MSA as partial (allow arbitrary distance at 
                          beginning and end)
   --blockscorefile=s     create a log file with blocks from the MSA and their scores
   --setname=s            set a name for the profile
   --setdesc=s            set a description for the profile
   --setacc=s             set an accession id for the profile
   --info=s               provide a text file containing the meta information
   --help                 show this message

";
}

#
# merge a block into a interblock region
#
sub merge_ibr {
    my ($ibr, $bl) = @_;
    
    foreach (0..$#$ibr) {
	$ibr->[$_] += 
	    $bl->{ibr}[$_] + scalar @{$bl->{columns}};
    }
}
    


#
# read in a clustal file
#
sub read_clustal {
    my ($seqref, $seqnameref) = @_;
    my $partlen=0;
    my $newpartlen=0;
    my $seqcount=0;
    my $status=0;
    while (<>) {
	last if /^\/\/$/;
	if (/^(\s*\S+)\s+([A-Za-z. -]+)/) {
	    $status=1;
	    my ($seqname, $newseq) = ($1, $2);
	    # Be tolerant with sequence names, but require
	    # at least one character not used in the conservedness line
	    # ignore line if it only contains conservedness characters 
	    if ($seqname=~s/^\s+//) {
		next unless $seqname=~/[^:.* ]/;
	    }
	    $newseq =~ s/ //g;
	    if ($newpartlen == 0) {
		$newpartlen = length($newseq);
	    } elsif ($newpartlen != length($newseq)) {
		die "CLUSTAL format error: Sequence length mismatch";
	    }
	    unless (exists $seqref->{$seqname}) { 
		push @$seqnameref, $seqname;
		$seqref->{$seqname} = "-" x $partlen;
	    } 
	    my $to_add = $partlen - length($seqref->{$seqname});
	    if ($to_add < 0) {
		die "CLUSTAL format error: Sequence names not unique";
	    } else {
		$seqref->{$seqname} .= $newseq;
	    }
	} elsif (/^\s*$/ && $status==1) {
	    $partlen += $newpartlen;
	    $newpartlen=0;
	    foreach (values %$seqref) {
### DEBUG part
		if (length($_) > $partlen) {
		    die "Internal bug in CLUSTAL parser. Please report";
		}
### 
		$_ .= "-" x ($partlen - length($_));
	    }
	    $status=0;
	}
    }
}

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
	    return undef;
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
    
    return $matrix;
}

sub output_format {
    return sprintf(($_[0] < 1e-4 ?  "%.2g" : $COMPATIBILITY_MODE ? "%.5f" : "%7.5f" ), @_);
}

sub format_collist {
    my @pos = map { $_->{num} } @{$_[0]};
    my $firstpos = shift @pos;
    my @output = ( [$firstpos, $firstpos] );
    foreach (@pos) {
	if ($_ == $output[-1][1] + 1) {
	    $output[-1][1]++;
	} else {
	    push @output, [ $_, $_ ];
	}
    }
    foreach (@output) {
	if ($_->[0] == $_->[1]) {
	    $_ = $_->[0];
	} elsif ($_->[0]+1 == $_->[1]) {
	    $_ = join(",",@$_);
	} else {
	    $_ = join("-", @$_);
	}
    }
    return join(",",@output);
}


############ main ###############
my $HELP;
GetOptions("width=i" => \$MIN_WIDTH,
	   "qij=s" => \$REG_MATR_FILE,
	   "mode=s" => \$WEIGHTING_MODE,
	   "regweight=f" => \$GLOB_WEIGHT,
	   "max_entropy=f" => \$MAX_ENTROPY,
	   "min_conserve=f" => \$MIN_CONSERVE,
	   "keep_empty" => \$KEEP_EMPTY,
	   "use_msa_dist" => \$USE_MSA_DIST,
	   "prefix_from_seqnames" => \$PREFIX_FROM_SEQNAMES,
	   "relax" => \$RELAX,
	   "setname=s" => \$SETNAME,
	   "setdesc=s" => \$SETDESC,
	   "setacc=s" => \$SETACC,
	   "noweights" => \$DISABLE_WEIGHTS,
	   "blockscorefile=s" => \$blockscorefile,
	   "info=s" => \$INFOFILE,
	   "help" => \$HELP,
    );

if ($HELP) {
    &usage();
    exit(0);
}

unless (defined $GLOB_WEIGHT) {
    $GLOB_WEIGHT = $DEFAULT_GLOB_WEIGHT;
    $GLOB_WEIGHT /= 4 if $WEIGHTING_MODE eq "blimps";
}

if (!defined $MAX_ENTROPY && defined $MIN_CONSERVE) {
    $MAX_ENTROPY = 1 - $MIN_CONSERVE;
}

if (defined $INFOFILE && -f $INFOFILE && open (INFOFILE, $INFOFILE)) {
    while (<INFOFILE>) {
	$SETACC = $1 if /^Accession number:\s*(\S+)/ && !defined $SETACC;
	$SETNAME = $1 if /^Name:\s*(\S+)/ && !defined $SETNAME;
	$SETDESC = $1 if /^Description:\s*(\S.*)/ && !defined $SETDESC;
    }
    close INFOFILE;
}

foreach ($SETNAME, $SETACC) {
    $_ = "unknown" unless defined;
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
my $regmatrix;
unless (defined $REG_MATR_FILE) {
    $REG_MATR_FILE="default.qij";
} elsif ($REG_MATR_FILE =~ s/^((~.*?|)\/)//) {
    # absolute directory given: do not check other dirs
    @CONFIG_DIRS=($1);
}
foreach (@CONFIG_DIRS) {
    my $fh;
    if (-e "$_$REG_MATR_FILE" && open $fh, "$_$REG_MATR_FILE") {
	$regmatrix = read_qij($fh);
	last;
    }
}
unless (defined $regmatrix) {
    print STDERR "Qij file '$REG_MATR_FILE' could not be read.\nCheck the parameter --qij=<filename>.\n\n";
    die;
}	


#
# read in fasta/clustal file
#
while (<>) {
    if ($. == 1 && /^CLUSTAL/) {
	&read_clustal(\%sequences, \@seqnames);
	last;
    }
    chomp;
    if (s/^>//) {
	chomp;
	while (exists $sequences{$_}) {
	    if (s/_\((\d+)\)$//) {
		$_ .= sprintf("_(%d)", $1+1);
	    } else {
		$_ .= "_(1)";
	    }
	}
	push @seqnames, $_;
    } elsif (@seqnames) {
	$sequences{$seqnames[-1]} .= $_;
    }
}

my $width=0;
my $height = @seqnames;

foreach (values %sequences) {
    s/[^A-Za-z.-]//g;
    $width=length if $width<length;
}

foreach (values %sequences) {
    my $rem = $width - length;
    if ($rem) {
	print STDERR "Warning: sequence lengths differ. Filled up with gaps.\n";
	$_ .=  "-" x $rem;
    }
}

my @blocks; 
my $last_ibr = $PREFIX_FROM_SEQNAMES ?
    [ map { /\/(\d+)-(\d+)$/ && $1 > 0 ? $1-1 : 0 } @seqnames ] :
    [ (0) x ($height+1) ];

push @$last_ibr, max @$last_ibr if ($PREFIX_FROM_SEQNAMES);

print STDERR "Determining block columns...\n";
foreach my $i (0..($width-1)) {
    my $val={ str=>[], num=>$i, frq=>{ "-" => 0, "." => 0} };
    my $percentage = int(20*$i/$width)*5;
    foreach (@seqnames) {
	my $c = substr($sequences{$_}, $i, 1);
	push @{$val->{str}}, $c;
	$val->{frq}{$c}++;
    }
    my $gapcount = sum(@{$val->{frq}}{"-","."});
    # if ($i == 180) {
    # 	print STDERR "\n".join("", @{$val->{str}})."\n";
    # }
    if ($gapcount <= $MAX_GAP_COUNT && !grep { /[a-z]/ } keys %{$val->{frq}}) {
#	print STDERR "*";
	# a new block column
	if (defined $last_ibr || @blocks == 0) {
	    # begin new block
	    push @blocks, { columns => [ $val ], ibr => $last_ibr };
	    undef $last_ibr;
	} else {
	    # add column to existing block
	    push @{$blocks[-1]{columns}}, $val;
	}
    } elsif ($gapcount < $height || $KEEP_EMPTY) {
#	print STDERR ".";
	unless (defined $last_ibr) {
	    $last_ibr = [ (0) x ($height+1) ];
	    if (@blocks && scalar @{$blocks[-1]{columns}} < $MIN_WIDTH) {
		merge_ibr ($last_ibr, pop @blocks);
	    }
	}
	foreach (0..($height-1)) {
	    $last_ibr->[$_]++ if $val->{str}[$_] !~ /[.-]/;
	}
	$last_ibr->[$height]++;
    } else {
	if ($gapcount > $height) {
	    die "gaps=$gapcount height=$height";
	}
#	print STDERR " ";
    }
#    print STDERR "\n" unless ($i % 80);
    printf STDERR "\r%2d%%", $percentage;
}
print STDERR "\rdone\n";

unless (defined $last_ibr) { 
    # profile ends with a block
    my $lastbl = $blocks[-1];
    my $blw = scalar @{$lastbl->{columns}};
    $last_ibr= [(0) x ($height+1)];
    if ($blw < $MIN_WIDTH) {
	merge_ibr($last_ibr, $lastbl);
	pop @blocks;
    } 
}

foreach (@blocks) {
    $_->{width} = scalar @{$_->{columns}};
}

unless ($DISABLE_WEIGHTS) {
    foreach my $bl (@blocks) {
	next if $bl->{width}==0;
	# calculate the weights
	my $weights = [ (0) x $height ];
	foreach (@{$bl->{columns}}){
	    my $frq = $_->{frq};
	    my @res = grep { /^[$AMINO_ACIDS]$/ } keys %$frq;
	    my %colweights = map { ( $_ => 0 ) } keys %$frq;
	    @colweights{@res} = map { 1/$frq->{$_}/(scalar @res) } @res;
	    my $i=0;
	    foreach my $c (@{$_->{str}}) {
		$weights->[$i++] += $colweights{$c};
	    }
	}
### DEBUG part
	# 
	if ($COMPATIBILITY_MODE) {
	    $bl->{compweights} = [ map { sprintf("%.3f", $_ * $height / $bl->{width}) } @$weights ];
	}
###
	$_ /= $bl->{width} foreach (@$weights);
	$bl->{weights} = $weights;
    }
}
    
my @suffices = ("A".."Z", "AA".."ZZ");

foreach my $bl (@blocks) {
    my @entropies;
    foreach (@{$bl->{columns}}) {
	my $freq = $_->{frq};
	# apply weights
	unless ($DISABLE_WEIGHTS) {
	    %$freq = ();
	    my $comptotal = sum (@{$bl->{compweights}}) if $COMPATIBILITY_MODE;
	    foreach my $i (0..$height-1) {
		$freq->{$_->{str}[$i]} += ($COMPATIBILITY_MODE ? 
					   $bl->{compweights}[$i]/$comptotal : 
					   $bl->{weights}[$i]);
	    }
	} else {
	    $_ /= $height foreach values %$freq;
	    delete $freq->{$_} foreach grep {  /[^$AMINO_ACIDS]/ } keys %$freq;
	}
### DEBUG
	if ($COMPATIBILITY_MODE) {
	    foreach ( ["U","C"], ["B","D"], ["Z","E"] ) {
		if (exists $freq->{$_->[0]}) {
		    $freq->{$_->[1]} += $freq->{$_->[0]};
		    delete $freq->{$_->[0]};
		}
	    }
	    my $rem_weight = 0;
	    foreach (grep { /[^$AMINO_ACIDS]/ } keys %$freq ) {
		$rem_weight += $freq->{$_};
		delete $freq->{$_};
	    }
	    if ($rem_weight>0) {
		$freq->{$_} += $rem_weight/$AA_COUNT foreach @AA_LIST;
	    }
	} else {
###
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
	}
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
    }  # end for each column

    # mark entropy based cutoffs
    last unless @{$bl->{columns}};
    @{$bl}{"from", "to"} = (0,$bl->{width});
    if (defined $MAX_ENTROPY) {
	while ($bl->{to} - $bl->{from} >= $MIN_WIDTH) {
	    if ($entropies[$bl->{from}]>$MAX_ENTROPY) {
		$bl->{from}++;
	    } elsif ($entropies[$bl->{to}-1]>$MAX_ENTROPY) {
		$bl->{to}--;
	    } else {
		last;
	    }
	}

	# high-entropy blocks to be removed
	if (sum(@entropies[$bl->{from}..$bl->{to}-1]) > $MAX_ENTROPY*($bl->{to}-$bl->{from})) {
	    $bl->{from} = $bl->{to};
	}
    }
}


push @blocks, { width=>0, ibr => $last_ibr };
for (my $i = 0; $i < $#blocks; $i++) {
    my ($bl, $nextbl) = @blocks[$i, $i+1];

    # set name to AC, plus consecutive suffix
    $bl->{name} = $SETACC."_".$suffices[$i];

    # remove cut off and high-entropy blocks
    my $cutoff = $bl->{width} - $bl->{to};
    $bl->{width} = $bl->{to} - $bl->{from};
    if ($bl->{width} < $MIN_WIDTH) {
	merge_ibr($nextbl->{ibr}, $bl);
	splice(@blocks, $i, 1);
	last if $i >= $#blocks;
	redo;
    }

    # apply right cutoff
    splice(@{$bl->{columns}}, $bl->{to});
    $_ += $cutoff foreach (@{$nextbl->{ibr}});

    # apply left cutoff
    splice(@{$bl->{columns}}, 0, $bl->{from});
    $_ += $bl->{from} foreach (@{$bl->{ibr}});
}

if (sum (map{$_->{width}} @blocks)==0) {
    print STDERR "No blocks found in MSA. Use \"prepareAlign\" to eliminate sequences.\n";
    exit;
}


# calculate interblock ranges
foreach (map { $_->{ibr} } @blocks) {
    my $maxdist=pop @$_;
    $maxdist=max @$_ unless $USE_MSA_DIST;
    @$_ = ( min(@$_), $maxdist );
}
if ($RELAX) {
    $blocks[0]{ibr}[1] = $blocks[-1]{ibr}[1] = "*";
} elsif ($PREFIX_FROM_SEQNAMES && grep { /\/\d+-\d+$/ } @seqnames) {
   $blocks[-1]{ibr}[1] = "*";
}

### output the profile
print "[name]\n$SETNAME\n";
print "# $SETDESC\n" if defined $SETDESC;

my $blcount=0;
my $colcount=0;
my $minlen=0;
foreach (@blocks) {
    print "\n[dist]\n# distance from previous block\n";
    print "# <min> <max>\n";
    print join("\t", @{$_->{ibr}})."\n";
    $minlen += $_->{ibr}[0];
    if (defined $_->{columns} && @{$_->{columns}}) {
	print "\n[block]\n# block no. $blcount follows, $height sequences, length $_->{width}\n";
	print "# corresponding to MSA columns:\n";
	print "# ".&format_collist($_->{columns})."\n";
	print "name=$_->{name}\n";
	$blcount++;
	$colcount += $_->{width};
	$minlen += $_->{width};
	print "#\n# <colnr> <probs for $AMINO_ACIDS>\n";
	print "#\t".join("\t", @AA_LIST)."\n";
	$_->{mu} = $_->{var} = 0;
	for (my $colno = 0; $colno < $_->{width}; $colno++) {
	    my $freq = $_->{columns}[$colno]{frq};
	    my $mu = 0;
	    foreach my $aa (@AA_LIST) {
		my $oddlog = log($freq->{$aa} / $backcol{$aa});
		$mu += $oddlog * $backcol{$aa};
		$_->{var} += $oddlog * $oddlog * $backcol{$aa};
	    }
	    $_->{mu} += $mu;
	    $_->{var} -= $mu * $mu;
	    print "$colno\t".join("\t", map { output_format($freq->{$_}) } @AA_LIST)."\n";
	}
    }
}

print "\n# created by:\n";
print "# $0 $INITARGS\n";

if (defined $blockscorefile && open OUTFILE, ">$blockscorefile") {
    my %output;
    my @outseqs;
    foreach my $seqname (@seqnames) {
	my $seq = $sequences{$seqname};
	my $pos=0;
	my $netcount=0;
	my $offset = ($PREFIX_FROM_SEQNAMES && $seqname=~s/\/(\d+)-\d+$//) ? $1-1 : 0;
	unless (exists $output{$seqname}) {
	    push @outseqs, $seqname ;
	    $output{$seqname} = "$seqname:\n";
	}
	foreach (@blocks) {
	    next unless $_->{width};
	    my $startpos=$_->{columns}[0]{num};
	    my $endpos=$_->{columns}[-1]{num};
	    my $netdiff = $offset;
	    $offset=0;
	    while ($pos < $startpos) {
		$netdiff++ if substr($seq, $pos, 1) =~ /^[A-Za-z]$/;
		$pos++;
	    }
	    $netcount += $netdiff;
	    my $score=1;
	    my $outseq="";
	    foreach my $col  (@{$_->{columns}}) {
		my $c=substr($seq, $col->{num}, 1);
		$outseq .= $c;
		$c =~ tr/a-z/A-Z/;
		if ($c =~ /^[$AMINO_ACIDS]$/) {
		    $score *= $col->{frq}{$c} / $backcol{$c};
		} else {
		    $score = 0;
		    last;
		}
	    }
	    my $spec;
	    if ($score==0) { $spec=0 } else { 
		$score=log($score) ;
		$spec = ($score - $_->{mu})/sqrt($_->{var});
	    }
	    $output{$seqname}.="$netdiff\t$_->{name}\t";
	    $output{$seqname}.=sprintf("%8.5f\t%8.5f\t%5d %s", exp($score/$_->{width}), $spec, $netcount, $outseq);
	    $pos = $endpos+1;
	    $netcount+=$_->{width};
	    $output{$seqname}.=sprintf(" %d\n", $netcount);
	}
    }
    print OUTFILE $output{$_}."--\n" foreach (@outseqs);
    close OUTFILE;
}
		

print STDERR "Profile has $blcount blocks with $colcount columns.\n";
print STDERR "Minimum admissible sequence length: $minlen\n\n";
	

