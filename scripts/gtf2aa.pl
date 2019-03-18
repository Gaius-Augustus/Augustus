#!/usr/bin/env perl

# Convert a gtf file and corresponding genome fasta file to protein file with standard translational code
# Katharina J. Hoff
# March 15th 2018
# uses functions from braker.pl

# usage: gtf2aa.pl genome.fa genes.gtf prot.fa

#use strict;
use warnings;

my $genome = $ARGV[0];
my $gtf = $ARGV[1];
my $prot = $ARGV[2];

$usage = "gtf2aa.pl genome.fa genes.gtf prot.fa\n";

if(scalar(@ARGV)!=3){
    print $usage;
    exit(1);
}

gtf2fasta($genome, $gtf, $prot);

###################################################################################################
# extract DNA sequence of CDS in gtf from genome fasta file, write to CDS fasta file               #
####################################################################################################
sub gtf2fasta {
    my $genome_file = shift;
    my $gtf_file = shift;
    my $fasta_file = shift;
    my %gtf;
    my %genome;
    open (GTF, "<", $gtf_file ) or die ("ERROR: in file " . __FILE__ ." at line ". __LINE__
                    . "\nCould not close file $gtf_file!\n");
    while ( <GTF> ) {
    if ( $_ =~ m/\tCDS\t/ ) {
        $_ =~ m/transcript_id \"(\S+)\"/;
        my $txid = $1;
        my @line = split(/\t/);
        push @{$gtf{$line[0]}{$txid}}, $_;
    }
    }
    close (GTF) or die ("ERROR: in file " . __FILE__ ." at line ". __LINE__
            . "\nCould not close file $gtf_file!\n");
    open (GENOME, "<", $genome_file ) or die ("ERROR: in file " . __FILE__
                          . " at line ". __LINE__ ."\nCould not close file $genome_file!\n");
    my $seq = "";
    my $locus;
    my %cds_seq;
    while (<GENOME>) {
    chomp;
    if( not ($_ =~ m/^>/ ) ) {
        $seq .= $_;
    } elsif ( $_ =~ m/^>(\S+)/ ) {
        if (defined ($locus) ) {
        while ( my ( $txid, $txgtf ) = each %{$gtf{$locus}} ) {
            foreach ( @{$txgtf} ) {
            my @line = split(/\t/);
            if ( not( defined( $cds_seq{$txid} ) ) ) {
                $cds_seq{$txid} = substr ( $seq, ( $line[3] -1 ),
                               ( $line[4] - $line[3] + 1 ) );
                if ( $line[6] eq '-' ) {
                $cds_seq{$txid} = reverse_complement ( $cds_seq{$txid} );
                }
            }else {
                if ( $line[6] eq '+') {
                $cds_seq{$txid} .= substr ( $seq, ( $line[3] -1 ),
                                ( $line[4] - $line[3] + 1 ) );
                } else {
                $cds_seq{$txid} = reverse_complement(
                    substr ( $seq, ( $line[3] -1 ), ( $line[4] - $line[3] + 1 ) ) )
                    . $cds_seq{$txid};
                }
            }
            }
        }
        }
        $locus = $1;
        $seq = "";
    }
    }
    # excise seqs for last contig:
    if (defined ($locus) ) {
    while ( my ( $txid, $txgtf ) = each %{$gtf{$locus}} ) {
        foreach ( @{$txgtf} ) {
        my @line = split(/\t/);
        if ( not( defined( $cds_seq{$txid} ) ) ) {
            $cds_seq{$txid} = substr ( $seq, ( $line[3] -1 ), ( $line[4] - $line[3] + 1 ) );
            if ( $line[6] eq '-' ) {
            $cds_seq{$txid} = reverse_complement ( $cds_seq{$txid} );
            }
        }else {
            if ( $line[6] eq '+') {
            $cds_seq{$txid} .= substr ( $seq, ( $line[3] -1 ),
                            ( $line[4] - $line[3] + 1 ) );
            } else {
            $cds_seq{$txid} = reverse_complement(
                substr ( $seq, ( $line[3] -1 ), ( $line[4] - $line[3] + 1 ) ) )
                . $cds_seq{$txid}
            }
        }
        }
    }
    }
    close(GENOME) or die ("ERROR: in file " . __FILE__ ." at line ". __LINE__
              . "\nCould not close file $genome_file!\n");
    open (FASTA, ">", $fasta_file) or die ("ERROR: in file " . __FILE__ ." at line " . __LINE__
                       . "\nCould not close file $fasta_file!\n");
    while ( my ( $txid, $dna ) = each %cds_seq ) {
    print FASTA ">$txid\n".dna2aa($dna)."\n";
    }
    close (FASTA) or die ("ERROR: in file " . __FILE__ ." at line ". __LINE__
              . "\nCould not close file $fasta_file!\n");
}

####################################################################################################
# reverse complement of DNA sequence                                                               #
####################################################################################################

sub reverse_complement {
    my $in = shift;
    $in =~ tr/ACGTacgt/TGCAtgca/;
    $in = reverse ($in);
    return $in;
}

####################################################################################################
# Translate DNA to protein sequence                                                                #
####################################################################################################

sub dna2aa {
    my %genetic_code = (
    'TCA' => 'S', # Serine
    'TCC' => 'S', # Serine
    'TCG' => 'S', # Serine
    'TCT' => 'S', # Serine
    'TTC' => 'F', # Phenylalanine
    'TTT' => 'F', # Phenylalanine
    'TTA' => 'L', # Leucine
    'TTG' => 'L', # Leucine
    'TAC' => 'Y', # Tyrosine
    'TAT' => 'Y', # Tyrosine
    'TAA' => '*', # Stop
    'TAG' => '*', # Stop
    'TGC' => 'C', # Cysteine
    'TGT' => 'C', # Cysteine
    'TGA' => '*', # Stop
    'TGG' => 'W', # Tryptophan
    'CTA' => 'L', # Leucine
    'CTC' => 'L', # Leucine
    'CTG' => 'L', # Leucine
    'CTT' => 'L', # Leucine
    'CCA' => 'P', # Proline
    'CAT' => 'H', # Histidine
    'CAA' => 'Q', # Glutamine
    'CAG' => 'Q', # Glutamine
    'CGA' => 'R', # Arginine
    'CGC' => 'R', # Arginine
    'CGG' => 'R', # Arginine
    'CGT' => 'R', # Arginine
    'ATA' => 'I', # Isoleucine
    'ATC' => 'I', # Isoleucine
    'ATT' => 'I', # Isoleucine
    'ATG' => 'M', # Methionine
    'ACA' => 'T', # Threonine
    'ACC' => 'T', # Threonine
    'ACG' => 'T', # Threonine
    'ACT' => 'T', # Threonine
    'AAC' => 'N', # Asparagine
    'AAT' => 'N', # Asparagine
    'AAA' => 'K', # Lysine
    'AAG' => 'K', # Lysine
    'AGC' => 'S', # Serine
    'AGT' => 'S', # Serine
    'AGA' => 'R', # Arginine
    'AGG' => 'R', # Arginine
    'CCC' => 'P', # Proline
    'CCG' => 'P', # Proline
    'CCT' => 'P', # Proline
    'CAC' => 'H', # Histidine
    'GTA' => 'V', # Valine
    'GTC' => 'V', # Valine
    'GTG' => 'V', # Valine
    'GTT' => 'V', # Valine
    'GCA' => 'A', # Alanine
    'GCC' => 'A', # Alanine
    'GCG' => 'A', # Alanine
    'GCT' => 'A', # Alanine
    'GAC' => 'D', # Aspartic Acid
    'GAT' => 'D', # Aspartic Acid
    'GAA' => 'E', # Glutamic Acid
    'GAG' => 'E', # Glutamic Acid
    'GGA' => 'G', # Glycine
    'GGC' => 'G', # Glycine
    'GGG' => 'G', # Glycine
    'GGT' => 'G'  # Glycine
    );
    my $seq = shift;
    $seq = uc($seq);
    my @codons = $seq =~ /(.{1,3})/g;
    my $aa = "";
    foreach ( @codons ) {
        if($_ =~ m/N/i){
            $aa .= "X";
        }else{
            $aa .= $genetic_code{$_};
        }
    }
    return $aa;
}
