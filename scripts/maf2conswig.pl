#!/usr/bin/perl
#
# Convert maf alignment format to a conservation graph in wig format
# Use simply average percentage of identity in a sliding window as VISTA does.
#
# input file format:
# 
# a score=5037
# s chlamy4.chromosome10    5424 36280 + 6579462 ATCAGCCACAG--ACC...
# s    volvox.scaffold_9 2188403 51128 + 9999999 ACCAGCCACGGGCACC...
# 
#
# 13.06.2009, Mario Stanke, mstanke@gwdg.de

use strict;
use Getopt::Long;

my $targetspecies=1; # use first species in alignment as default target species
my @names;           # array of species.seqnames
my @begins;          # array of all begin positions
my @strands;         # array of all strands
my @lengths;         # array of all lengths
my @alirows;         # array of all alignment rows (including gaps)

my $help = 0;
my $winsize = 31;
my @cumids;             #array of number of identical bases
my $len; 

GetOptions('help!'=>\$help,
	   'winsize:i' =>\$winsize);

exec("perldoc $0") if ($help);

my $radius = int($winsize/2);


while (<>) {
    chomp;
    if (!/^s / && @alirows){
	compute_wig();
	@names = @begins = @strands = @alirows = ();
    }
    if (/^s /){
	my @f = split (/\s+/, $_);
	push @names, $f[1];
	push @begins, $f[2];
	push @strands, $f[4];
	push @lengths, $f[5];
	push @alirows, $f[6];
    }
}

if (@alirows){
    compute_wig();
}

sub compute_wig{
    die ("This script cannot deal with multiple aligmnments yet.\n") if (@alirows > 2); 
    $len = length($alirows[0]);
    my $chrom = $names[$targetspecies-1];
    $chrom =~ s/.*?\.//;
  
    compute_identical(); # number of rows identical to reference for each position in the alignment
    my @scores = ();
    my $seqpos = 0;
    for (my $i=0; $i < $len; $i++){
	next if (substr($alirows[$targetspecies-1], $i,1) eq '-'); # skip when gaps are in the reference sequence
	my $a = $i - $radius;
	my $b = $i + $radius;
	$b = $len-1 if($b > $len-1);
	my $numid;
	if ($a>0){
	    $numid = $cumids[$b] - $cumids[$a-1];
	} else {
	    $numid = $cumids[$b];
	}
	my $percentid = sprintf ("%.3f", $numid/(@alirows-1)/(2*$radius+1));
	$seqpos++;
	push @scores, $percentid;
	# print "i=$i is pos $seqpos\n";
    }

    my $start;
    if ($strands[$targetspecies-1] eq '-') {
	@scores = reverse @scores;
	$start = $lengths[$targetspecies-1] - $begins[$targetspecies-1] - scalar(@scores);
    } else {
	$start = $begins[$targetspecies-1]+1;
    }
    print "fixedStep chrom=$chrom start=$start step=1\n";
    print "" . join("\n", @scores) . "\n";
}

# number of rows identical to target species in columns 1 through 'pos'
# assume pairwise alignment for now
# In case of multiple alignment this script should probably 
# use the pairwise induced alignments (delete '-' aligned against '-'). 
sub compute_identical {
    @cumids = ();
    my $id=0;
    for (my $pos = 0; $pos< $len; $pos++){
	my $tchar = substr($alirows[$targetspecies-1], $pos,1);
	if ($tchar ne '-'){
	    for (my $k=0; $k < @alirows; $k++){
		next if ($k == $targetspecies-1); # skip comparison with target species itself
		$id++ if ($tchar eq substr($alirows[$k], $pos, 1));
	    }
	}
	# print "pos=$pos id=$id\n";
	push @cumids, $id;
    }
}


__END__

=head1 NAME

maf2conswig.pl convert maf alignment format to a conservation graph in wig format 
               Use simply average percentage of identity in a sliding window as VISTA does.

=head1 SYNOPSIS

maf2conswig.pl < alignment.maf > consscore.wig
options
     winsize=n      size of averaging window, default 40
=cut
