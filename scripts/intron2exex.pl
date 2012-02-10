#!/usr/bin/perl
# make a fasta file with exon-exon sequences
# 
# Mario Stanke, 7.1.2010

use strict;
use Getopt::Long;

my $flank = 75; # amount of flanking DNA from both exons
my $help = 0;
my $N = 1; # counter for the id of the sequences
my ($gfffile, $seqfile, $exexfile, $mapfile, $seq);

GetOptions('introns=s'=>\$gfffile,
           'seq=s'=>\$seqfile,
           'exex=s'=>\$exexfile,
           'map=s'=>\$mapfile,
	   'flank:n'=>\$flank,
           'help!'=>\$help); 

exec("perldoc $0") if ($help || !defined($gfffile) || !defined($seqfile) || !(defined($exexfile)));

die ("$seqfile does not exist.") if (! -e $seqfile);
open(GFFFILE, "<$gfffile") || die "Couldn't open $gfffile\n";

#
# First read in and store all intron annotations
#

#
# data structure:
# annos: hash of annotations
# annotation: hash of gbfkeys, keys: CDS:genename or mRNA:genename
# gbfkey: list:
#         {name, strand, exon1start, exon1start, ..., exonXstart, exonXstart}
#         name is "mRNA" or "CDS"
#
my %allintrons; # keys: chromosomes, reference to array of hash refs
                #                                 keys: "start" "end" "strand"
my $intron;     # reference to a hash
my ($chr, $type, $begin, $end, $strand);
my @f;

while (<GFFFILE>) {
    s/#.*//;
    next unless /\S/;
    if (/^(\S+):(\d+)-(\d+)$/){
	$chr = $1;
	$begin = $2;
	$end = $3;
	$allintrons{$chr} = [] if (!defined($allintrons{$chr}));
	push @{$allintrons{$chr}}, {"start"=> $begin, "end" => $end, "strand" => "."}; # create hash ref
    } else {
	split /\t/, $_, 9;
	if (@f < 8) { warn @f,"Neither simple nor GFF format"; next }
	$chr = $f[0];
	$type = $f[2];
	$begin = $f[3];
	$end = $f[4];
	$strand = $f[6];
	next if ($type ne "intron");
	$allintrons{$chr} = [] if (!defined($allintrons{$chr}));
	push @{$allintrons{$chr}}, {"start"=> $begin, "end" => $end, "strand" => $strand}; # create hash ref
    }
}

#
# Now read in the genome sequentially, chromosome by chromosome
#
open(GENOME, "<$seqfile") || die "Couldn't open $seqfile.\n";
open(EXEX, ">$exexfile") || die "Couldn't open $exexfile for writing\.n";
open(MAP, ">$mapfile") || die "Could not open $mapfile for writing.\n" if (defined($mapfile));

$/="\n>";
while(<GENOME>) {
    /[>]*(\S+)/;
    $chr = $1;
    $seq = $'; #'
    $seq =~ s/>//;
    $seq =~ s/\n//g;

    my $length = length $seq;
    my $introns = $allintrons{$chr};

    foreach my $intron (@{$introns}){
	next if ($intron->{"start"} < 1 || $intron->{"end"} > $length);
	my $a = $intron->{"start"} - 1 - $flank;
	my $d = $intron->{"end"} + $flank;
	$a = 0 if ($a < 0);
	$d = $length-1 if ($d >= $length);
	my $L1 = $intron->{"start"} - 1 - $a;
	my $L2 = $d - $intron->{"end"};
	my $seqwin = substr($seq, $a, $L1) . substr($seq, $intron->{"end"}, $L2);
	my $id = "exex$N";
	$N++;
	print EXEX ">$id" 
	    . " $L1-" . substr($seq, $intron->{"start"}-1, 2) . substr($seq, $intron->{"end"}-2, 2) . "-$L2"
	    . " $chr" . $intron->{"strand"} . ":". $intron->{"start"} . "-" . $intron->{"end"}
	    . "\n$seqwin\n";
	if (defined($mapfile)){
	    print MAP ($L1+$L2) . "\t0\t0\t0\t0\t0\t1\t" . ($intron->{"end"}- $intron->{"start"} + 1)
		. "\t" . (($strand ne "-")? "+":"-") . "\t$id\t" 
		. ($L1+$L2) . "\t0\t" . ($L1+$L2) . "\t$chr\t$length\t$a\t$d\t2\t$L1,$L2,\t"
		. "0,$L1,\t" . "$a," . $intron->{"end"} . ",\n";
	}
   }
}
close EXEX;
close GENOME;
close MAP if (defined($mapfile));

sub reversecomplement {
    my $s = shift;
    $s = reverse $s;
    $s =~ tr/acgtACGT/tgcaTGCA/;
    return $s;
}

__END__

=pod

=head1 NAME

intron2exex.pl   make a fasta file with exon-exon sequences

=head1 SYNOPSIS

intron2exex.pl --introns=introns.gff --seq=genome.fa --exex=exex.fa --map=map.psl

    The input coordinates are used to make small artificial transcript sequences by splicing out the intron
    and writing a fixed flanking region on both sides of the intron. This script can be used as a help for
    spliced short read alignment. See also pslMap.pl.
    
=head1 OPTIONS

   introns input file with introns in GFF format or in this format:
           chr10:100008749-100010821
           chr10:100010934-100011322
           ...
           coordinates are the 1-based start and end of intron candiates
   seq     input genome sequence
   exex    ouput fasta file with intron-flanking sequences
   map     psl file with mapping information between the exon-exon sequences and the genome
   flank   amount of flanking 'exon' sequence on both sides of the intron (default: 75)

=head1 DESCRIPTION
 

=cut
