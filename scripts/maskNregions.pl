#!/usr/bin/perl

# create a gff-like file that marks regions in a fasta file
# that contain "N" or "n".

# Katharina Hoff, 18.7.2011

my $usage = "maskNregions.pl genome.fa > out.file\n";

if (@ARGV != 1) {
    print $usage;
    exit;
}

my $genome = $ARGV[0];
my $fastaContent = "";
my @seqLetters = ();

sub procSeq;

open (FASTA, $genome) or die ("\n\ncould not open file $genome!\n");
LINE: while ($line = <FASTA>){
    next LINE if $line =~ m/^#/; #discard comments
    if ($line =~ /^>/){
	procSeq($fastaContent);
	$fastaContent = "";
    }else{
    $line =~ s/[\x0A\x0D]+//g; 
    $line =~ s/(\s+)(\n)(\r)//g;
    $fastaContent = $fastaContent.$line;
    }
}
procSeq($fastaContent) if ($fastaContent ne "");

close(FASTA) or die "Could not close file $genome!\n";

sub procSeq{
    my $seq = shift; # nucleotide sequence
    if(length($seq)>0){
	@seqLetters = split(//, $seq);
	$n = scalar(@seqLetters);
	$pos = 1;
	foreach(@seqLetters){
	    if($_ eq "n" or $_ eq "N"){
		if($pos == 1 || ($seqLetters[$pos-2] ne "n" and $seqLetters[$pos-2] ne "N")){
		    #print "$pos on $fastaHeader is the beginning of N-region\n";
		    print $fastaHeader."\t"."maskN"."\t"."N-region"."\t".$pos."\t";
		}
		if($pos == $n || ($seqLetters[$pos] ne "n" and $seqLetters[$pos] ne "N")){
		    #print "$pos on $fastaHeader is the end of N-region\n";
		    print $pos."\t"."."."\t"."."."\t"."."."\t"."N-region\n";
		}
	    }
	    $pos = $pos+1;
	}
    }
    $line =~ s/[\x0A\x0D]+//g; #removing those ugly whitelines
    $line =~ s/(\n)(\r)//g; #remove them alllll!
    $line =~ m/(^>\w+)/i; #matches a word starting with > (Fasta)
    $fastaHeader = substr($1, 1);
}


