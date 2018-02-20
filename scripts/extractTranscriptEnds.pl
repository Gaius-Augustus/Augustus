#!/usr/bin/perl

# extract transcript ends of a given length from a gene fasta file (no UTR)
# Katharina Hoff, 9.8.2011, katharina.hoff@gmail.com

my $n = 30;

my $seq;
my $seqLen;
my $i = 1;
while(<STDIN>){
    if($_=~m/^>/){
	$lenLen = length($seq);
	print ">seq".$i."\n";
	print substr($seq, $seqLen-30, 27)."\n";
	$seq = "";
	$i = $i+1;
    }else{
	chomp;
	$seq = $seq.$_;
    }
}

if(length($seq)>1){
    $lenLen = length($seq);
    print ">seq".$i."\n";
    print substr($seq, $seqLen-30, 30)."\n";
}
