#!/usr/bin/perl
# convert old gbrowse format of gene predictions to gff3 format
# Mario Stanke, 3/8/2011, mario.stanke@uni-greifswald.de

use strict;
use Getopt::Long;

my $help = 0;
my $hints = 0;
my $id = 0;
my ($genename, $txname,$source);

GetOptions('help!'=>\$help,
    'hints!'=>\$hints,
    'source=s'=>\$source);

exec("perldoc $0") if ($help);

while(<>){
    chomp;
    my @f = split /\t/;
    if (@f == 9){
	if (!$hints){
	    if ($f[2] eq "gene"){
		@f[8] =~ /Gene (\S+)/;
		$genename = $1;
		@f[8] =~ s/Gene (\S+)/ID=$genename/;
	    } elsif ($f[2] eq "mRNA"){
		@f[8] =~ s/mRNA (\S+)\s*[;]?\s*(Gene \S+)?/ID=$1;Parent=$genename/;
		$txname = $1;
	    } else { 
		@f[8] =~ s/mRNA (\S+);?/Parent=$1/;
	    }
	    if (@f[2] =~ /5\S*UTR/i){
		@f[2] = "five_prime_utr";
	    }
	    if (@f[2] =~ /3\S*UTR/i){
		@f[2] = "three_prime_utr";
	    }
	    next if (@f[2] eq "internal" || @f[2] eq "single" ||  @f[2] eq "initial" || @f[2] eq "terminal" || @f[2] eq "tts" || @f[2] eq "tss");
	} else {
	    $id++;
	    @f[8] =~ s/$/ID=$id/;
	    @f[8] =~ s/mult=/Note=/;
	    @f[8] =~ s/Note /Note=/;
	}
	@f[1] = $source if (defined($source));
    }

    print "" . join("\t", @f) . "\n";
}


__END__

=pod

=head1 NAME

gbrowseold2gff3.pl    format convert a GFF file

=head1 SYNOPSIS

gbrowseold2gff3.pl < in.gbrowse > out.gff3

    This is a simple conversion, almost line by line, from a GBrowse 1.6 compatible to a GBrowse 2.0 compatible format for genes.
    
=head1 OPTIONS

    hints       convert gff format for intron hints to gff3
    source      fill this into the second column

=head1 DESCRIPTION
    
    example input:

   chromosome_1    AUGUSTUS      gene    12656   14013   0.05    +       .       Gene g1
   chromosome_1    AUGUSTUS      mRNA    12656   14013   0.05    +       .       mRNA g1.t1 ; Gene g1
   chromosome_1    AUGUSTUS      tss     12656   12656   .       +       .       mRNA g1.t1
   chromosome_1    AUGUSTUS      5'-UTR  12656   12867   0.21    +       .       mRNA g1.t1
   chromosome_1    AUGUSTUS      start_codon     12868   12870   .       +       0       mRNA g1.t1
   chromosome_1    AUGUSTUS      initial 12868   12993   0.85    +       0       mRNA g1.t1
   chromosome_1    AUGUSTUS      terminal        13249   13479   0.99    +       0       mRNA g1.t1
   chromosome_1    AUGUSTUS      intron  12994   13248   0.99    +       .       mRNA g1.t1
   chromosome_1    AUGUSTUS      CDS     12868   12993   0.85    +       0       mRNA g1.t1
   chromosome_1    AUGUSTUS      CDS     13249   13479   0.99    +       0       mRNA g1.t1
   chromosome_1    AUGUSTUS      stop_codon      13477   13479   .       +       0       mRNA g1.t1
   chromosome_1    AUGUSTUS      3'-UTR  13480   14013   0.22    +       .       mRNA g1.t1
   chromosome_1    AUGUSTUS      tts     14013   14013   .       +       .       mRNA g1.t1

    example output:

   chromosome_1    AUGUSTUS      gene    12656   14013   0.05    +       .       ID=g1
   chromosome_1    AUGUSTUS      mRNA    12656   14013   0.05    +       .       ID=g1.t1;Parent=g1
   chromosome_1    AUGUSTUS      five_prime_utr  12656   12867   0.21    +       .       Parent=g1.t1
   chromosome_1    AUGUSTUS      start_codon     12868   12870   .       +       0       Parent=g1.t1
   chromosome_1    AUGUSTUS      intron  12994   13248   0.99    +       .       Parent=g1.t1
   chromosome_1    AUGUSTUS      CDS     12868   12993   0.85    +       0       Parent=g1.t1
   chromosome_1    AUGUSTUS      CDS     13249   13479   0.99    +       0       Parent=g1.t1
   chromosome_1    AUGUSTUS      stop_codon      13477   13479   .       +       0       Parent=g1.t1
   chromosome_1    AUGUSTUS      three_prime_utr 13480   14013   0.22    +       .       Parent=g1.t1

=cut
