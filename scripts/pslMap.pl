#!/usr/bin/perl
# combine two psl files transitively
# 
#     map target (genome)
#  \           \      /           /
#   \           \    /           /
#    \           \  /           /
#     \           \/           /
#      ------------------------   map query = in target (mRNA)
#        -- ----  -----   ---    
#       / / |  |  \   \   | |
#      / /  |  |   \   \  | |
#     / /   |  |    \   \ | |
#    / /    |  |     \   \| |
#    --     ----      -------    
#    in query (transcript reads)
#
#
#
# Mario Stanke, 7.1.2010

use strict;
use Getopt::Long;

my $filterunspliced = 0;
my $help = 0;
my $verbose = 0;
my $multerr = 0;
my $minmatch = 0;
my $maperr = 0;
my $numunspliced = 0;
my ($infile, $mapfile, $outfile);

GetOptions('in=s'=>\$infile,
           'out=s'=>\$outfile,
           'map=s'=>\$mapfile,
	   'minmatch:n'=>\$minmatch,
           'filterunspliced!'=>\$filterunspliced,
           'help!'=>\$help); 

exec("perldoc $0") if ($help || !defined($infile) || !defined($outfile) || !(defined($mapfile)));

die ("$infile does not exist.") if (! -e $infile);
open(MAP, "<$mapfile") || die "Couldn't open $mapfile.\n";

#
# First read in and store the map
#
my %map = (); # keys: query names: values psl array refs

while (<MAP>) {
    my @f = split /\t/, $_;
    if (@f < 21) { warn @f,"Not PSL format"; next }
    if (exists($map{$f[9]})){
	$multerr++;
	if ($multerr <= 10){
	    print STDERR "$f[9] occurs more than once in map file. Will ignore all but first.\n";
	} elsif ($multerr == 10) {
	    print STDERR "Suppressing this error message from now on.\n";
	}
    }
    die ("Map file is required to map the whole query in this version. Violating line is $_;") 
	if ($f[10] != $f[12] || $f[11] != 0);
    $map{$f[9]} = \@f;
}
close MAP;


#
# Now read in.psl sequentially and convert to out.psl
#
open(IN, "<$infile") || die "Couldn't open $infile.\n";
open(OUT, ">$outfile") || die "Couldn't open $outfile for writing\.n";

while(<IN>){
    my @f = split /\t/, $_;
    if (@f < 21) { warn @f,"Not PSL format"; next }
    if (exists($map{$f[13]})){
	my $m = $map{$f[13]};
	print "mapping\n" . join("\t", @f) . "\nalong\n"  . join("\t", @$m) . "\n" if ($verbose);
	next if ($f[0] < $minmatch); # filter out bad alignments with few matches
	if ($filterunspliced && !spliced(\@f, $m)){
	    $numunspliced++;
	    next;
	}
	my $o = mapali(\@f, $m);
	print OUT "" . join("\t", @$o) . "\n";
    } else {
	$maperr++;
	if ($maperr <= 10){
	    print STDERR "The target $f[13] is not contained in $mapfile. Will not map alignment.\n";
	} elsif ($maperr == 10) {
	    print STDERR "Suppressing this error message from now on.\n";
	}
    }
}

close IN;
close OUT;

print STDERR "$maperr alignments could not be mapped because the target was not a query in $mapfile.\n" if ($maperr);
print STDERR "$multerr lines in $mapfile had a duplicated query.\n" if ($multerr);
print STDERR "$numunspliced unspliced read alignments were filtered out.\n" if ($filterunspliced);

sub mapali {
    my $f = shift;
    my $m = shift;
    
    # copy these values and assume that map alignment is perfect here
    my @o = ($f->[0],$f->[1],$f->[2],$f->[3],$f->[4],$f->[5]);
    
    push @o, $f->[6]+$m->[6]; # number of inserts in target
    push @o, $f->[7]+$m->[7]; # number of bases in inserts in target
    die ("Negative strand alignments in map are not implemented yet.") if ($m->[8] eq "-");
    my $strand = ($f->[8] eq $m->[8])? "+" : "-";
    push @o, $strand;
    push @o, $f->[9]; # query name
    push @o, $f->[10]; # query size
    push @o, $f->[11]; # query start
    push @o, $f->[12]; # query end
    push @o, $m->[13]; # target name
    push @o, $m->[14]; # target size

    # b=block sizes, q=query starts, t=target starts
    my @ob; my @oq; my @ot; # out
    # 'in' alignment
    my @fb = split(",", $f->[18]);
    my @fq = split(",", $f->[19]);
    my @ft = split(",", $f->[20]);   
    # 'map' alignment
    my @mb = split(",", $m->[18]);
    my @mq = split(",", $m->[19]);
    my @mt = split(",", $m->[20]);

    # TODO strand treatment

    # go over the map blocks and 'translate' the 'in' blocks falling in each 'map' block
    while(@mb){
	while(@fb && $ft[0]+$fb[0] <= $mq[0]+$mb[0]){ # 'in' blocks falling completely in 'map; block
	    push @ob, $fb[0]; # block size unchanged
	    push @oq, $fq[0]; # query start stays the same
	    # target start is sum of 'in' and 'map' target starts minux 'map' query start
	    push @ot, $mt[0] + $ft[0] - $mq[0];
	    shift @fb; shift @fq; shift @ft; # remove another 'in' block
	}
	if (@fb && $ft[0] < $mq[0]+$mb[0]){ # have an 'in' block that is intron split
	    my $fracsize =  $mq[0]+$mb[0] - $ft[0]; # block fraction size that fits into first 'map' block
	    push @ob, $fracsize;
	    push @oq, $fq[0]; # query start is 'in' target start
	    push @ot, $mt[0] + $ft[0]; # target start is sum of 'in' and 'map' target starts
	    # 'remove' fraction from first 'in' block
	    $fb[0] -= $fracsize;
	    $fq[0] += $fracsize;
	    $ft[0] += $fracsize;
	}
	shift @mb; shift @mq; shift @mt; # remove another 'map' block
    }
    

    push @o, $ot[0]; # target start
    push @o, $ot[-1] + $ob[-1]; # target end
    push @o, scalar(@ob); # block count
    push @o, "" . join (",", @ob) . ","; # block sizes
    push @o, "" . join (",", @oq) . ","; # query starts
    push @o, "" . join (",", @ot) . ","; # target starts
    print "result:\n" . join("\t", @o) . "\n\n" if ($verbose);
    return \@o;
}

# spliced
# check whether target range of read to mRNA alignment contains an intron in the map file
#
sub spliced {
    my $f = shift;
    my $m = shift;
    my @mq = split(",", $m->[19]);

    shift @mq;
    foreach my $mapqstart (@mq){
	return 1 if ($mapqstart > $f->[15] && $mapqstart < $f->[16]);
    }
    print "not spliced\n" if ($verbose);
    return 0;
}


__END__

=pod

=head1 NAME

pslMap.pl       combine two psl files transitively    

=head1 SYNOPSIS

pslMap.pl [ --filterunspliced ] --in=in.psl --map=map.psl --out=out.psl

    If in and map are considered to be mappings from the query to the target, then
    out is the consecutive application of in and map. The map is assumed to not contain mismatches.
    In contrast to the UCSC pslMap C program, this program is reporting the mismatches of in.psl in
    the output alignments.
    
=head1 OPTIONS

  in               psl input file, e.g. from an alignment of transcript reads against mRNA sequences
  map              psl input file, e.g. from an alignment of mRNA sequences against a genome
                   The queries of map must be the targets of in.
  out              psl output file, e.g. the inferred alignment of transcript reads against the genome
  filterunspliced  alignments in in.psl that do not cover an intron in map.psl are not reported in out.psl
  minmatch         filter out alignments with a number of matches (first column) that is below this (default: 0)

=head1 DESCRIPTION
    
    The planned typical usage of this script is for RNA-Seq alignment.
    Then in.psl can contain alignments of reads against small exon-exon windows (or candidates thereof)
    and can be done with a fast almost exact mapper that does no spliced alignments.
    out.pl will then be the spliced alignments of the reads against the genome.
    The option --filterunspliced can then be used to ensure that read alignments are not reported twice.
    The unspliced alignments can be done directly to the genome and the spliced only through this indirect approach.
    Example input visualized:
    
        map target (genome)
     \           \      /           /
      \           \    /           /
       \           \  /           /
        \           \/           /
         ------------------------   map query = in target (mRNA)
           -- ----  -----   ---    
          / / |  |  \   \   | |
         / /  |  |   \   \  | |
        / /   |  |    \   \ | |
       / /    |  |     \   \| |
       --     ----      -------    
       in query (transcript reads)


=cut
