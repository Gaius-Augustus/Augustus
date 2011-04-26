#!/usr/bin/perl
# gffGetmRNA.pl
# Creates mRNA fasta sequence files from a GFF with gene annotation
# 
# Mario Stanke, January, 2010
#

use strict;
use Getopt::Long;

my $usage = "cat genes.gff | gffGetmRNA.pl --genome=s --mrna=s\n";
$usage .= "   Makes a fasta file with mRMA sequences\n";
$usage .= "   --genome=s   Input a fasta file with the genomic sequences.\n";
$usage .= "   --mrna=s     Output fasta file with mRNA sequences.\n";

my ($seqfile, $mrnafile);
my %sequence = ();
my %ignorechr = ();

GetOptions('genome=s'=>\$seqfile, 
	   'mrnae=s'=>\$mrnafile);

if (!defined($seqfile) || !defined($mrnafile)) {
    print $usage;
    exit;
} 

my @gfflines; # complete input file

# hash of transcripts
# keys: chromosome names
# values: hash reference
#         keys: transcript ids
#         values: gff array reference
my %txs = ();


readGenome();
parseAndStoreGFF();
sortAndCheck();
open MRNA, ">$mrnafile" or die ("Could not open mrnafile for writing.");
printmrnas();
close MRNA;
if (keys %ignorechr > 0){
    print "Annotation from these chromosomes was ignored because they were not contained in $seqfile.:\n"
	. join (", ", sort keys %ignorechr) . "\n";
}

# Read in the sequence file in one chunk.
# And sort it in the sequence hash.
# Yes, this requires a lot of memory for large genomes.
sub readGenome {
    my $seqname;
    if ($seqfile){
	open (SEQ, "<$seqfile") or die ("Could not open sequence file $seqfile\n");
	$/=">";
	while (<SEQ>){
	    s/>$//;
	    next unless /\S+/;
	    /(.*)\n/;
	    
	    $seqname = $1;
	    my $sequencepart = $'; #'
	    $seqname =~ s/\s.*//; # seqname only up to first white space
	    $sequencepart =~ s/\s//g;
	    $sequence{$seqname} = $sequencepart;
	}
	print "Read in " . (scalar keys %sequence) . " sequence(s) from $seqfile.\n";
	$/="\n";
    }
}

#
# Read in the GFF file and store the exons by their transcript id
#
sub parseAndStoreGFF{
    my ($txid, $type, $chr);
    while(<STDIN>){
        my @f = split /\t/;
        next if (@f<8);
	if ($f[8] =~ /(transcript_id|Transcript)."([^"]+)"/){ #"
	    $txid = $2;
	} elsif ($f[8] =~ /gene_id."([^"]+)"/){ #"
	    $txid = $1;
	} elsif ($f[8] =~ /Parent=Transcript:([^;]+)/){
	    $txid = $1;
	} elsif ($f[8] =~ m/Parent=([^;]+)/){
	    $txid = $1;
	} else {
	    $txid = $f[8];
	}
	$txid =~ s/\n//;
	my ($chr,$type) = ($f[0],$f[2]);
	if ($type eq "exon"){
	    $txs{$chr}{$txid} = [] if (!exists($txs{$chr}{$txid}));
	    push @{$txs{$chr}{$txid}}, \@f;
	}
    }
}

sub sortAndCheck {
    foreach my $chr (sort keys %txs){
	foreach my $txid (sort keys %{$txs{$chr}}){
	    @{$txs{$chr}{$txid}} = sort {$a->[3] <=> $b->[3]} @{$txs{$chr}{$txid}}; # sort by start position
	    my $strand;
	    my $end;
	    foreach my $exon (@{$txs{$chr}{$txid}}){
		if (!defined($strand)){
		    $strand = $exon->[6];
		} else {
		    die ("Strand not consistent in transcript $txid.") if ($strand ne $exon->[6]);
		}
		if (defined($end) && $end >= $exon->[3]){
		    die ("Exons overlapping in $txid.");
		}
	    }
	}
    }
}

sub printmrnas() {
    foreach my $chr (sort keys %txs){
	if ($sequence{$chr}){
	    foreach my $txid (sort keys %{$txs{$chr}}){
		my $seq = "";
		foreach my $exon (@{$txs{$chr}{$txid}}){
		    die ("Coordinate of $txid out of sequence range. $exon->[4] >= " . length($sequence{$chr}) . "\n") 
			if ($exon->[4] >= length($sequence{$chr}));
		    my $seqpart = lc(substr($sequence{$chr}, $exon->[3]-1, $exon->[4] - $exon->[3] + 1));
		    $seq .= $seqpart;
		}
		$seq = rc($seq) if ($txs{$chr}{$txid}->[0]->[6] eq "-");
		print MRNA ">$txid\n" . getFa($seq, 100);
	    }
	} else {
	    $ignorechr{$chr}++;
	}
    }
}


# reverse complement
sub rc{
    my $s = shift;
    $s = reverse $s;
    $s =~ tr/acgtACGT/tgcaTGCA/;
    return $s;
}

sub getFa{
    my $seq = shift;
    my $cols = 100;
    $cols = shift if (@_);
    
    my $start = 0;
    my $ret = "";
    while (length($seq)-$start >= $cols) {
	my $shortline = substr($seq, $start, $cols);
	$ret .= "$shortline\n";
	$start += $cols;
    }
    $ret .= substr($seq, $start, $cols) . "\n";
}
