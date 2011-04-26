#!/usr/bin/perl
use strict;
use warnings;

sub calculateACGT {
  my $seq=shift(@_);
  my $contents=shift(@_);
  my $length= length($seq);
  my $char="";
  for (my $i=0;$i<$length;$i++) {
    $char=uc(substr($seq,$i,1));
    ${$contents}{$char}++;
  }
  ${$contents}{'length'}=$length;
}

sub output {
  my $data=shift(@_);
  if ((shift(@_))==1) {
    print("summary: BASE COUNT     ${$data}{'A'} a   ${$data}{'C'} c  ${$data}{'G'} g   ${$data}{'T'} t");  
    if (${$data}{'N'} > 0) {
      print("   ${$data}{'N'} n");
    }
    if (${$data}{'rest'} > 0) {
      print("   ${$data}{'rest'} ?");
    }
    print("\n");
    print("total ${$data}{'length'}bp\n");
    print("gc: ${$data}{'gc'}%\n"); 
   }
  else {
    print("${$data}{'length'} bases.\t${$data}{'name'} ");
    print("BASE COUNT     ${$data}{'A'} a   ${$data}{'C'} c  ${$data}{'G'} g   ${$data}{'T'} t");
    if (${$data}{'N'} > 0) {
      print("   ${$data}{'N'} n");
    }
    if (${$data}{'rest'} > 0) {
      print("   ${$data}{'rest'} ?");
    }
    print("\n");
  }
}

sub main {
  my $usage="";
  my $seqfilename="";
  my $seq="";
  +my %totaldata=(A=>0,C=>0,G=>0,T=>0,N=>0,rest=>0,gc=>0);
  $usage .= "$0 -- count the a,c,g,t \n";
  $usage .= "\n";
  $usage .= "Usage: $0 seq-file\n";
  $usage .= "\n";

  if (scalar(@ARGV) < 1) {
      die "\n$usage";
  }

  $seqfilename = $ARGV[0];
  open(FASTA, "<$seqfilename") || die "Couldn't open $seqfilename\n";

  $/="\n>";

  while(<FASTA>) {
    my %data=(A=>0,C=>0,G=>0,T=>0,N=>0); # initialize counts to 0
    /[>]*(.*)\n/;
    
    $data{'name'}= $1;
    $seq = $';
    $data{'name'}=~ s/\s+.*//;
    $seq =~ s/>//;
    $seq =~ s/\n//g;
    &calculateACGT($seq,\%data);
    $data{'rest'}=$data{'length'}-$data{'A'}-$data{'T'}-$data{'G'}-$data{'C'}-$data{'N'};
    foreach my $key ("A","C","T","G","N","rest") {
      $totaldata{$key} += $data{$key};
    }
    &output(\%data,0); 
  }
  my $sum=$totaldata{'T'}+$totaldata{'G'}+$totaldata{'A'}+$totaldata{'C'};
  $totaldata{'length'}=$totaldata{'N'}+$totaldata{'rest'}+$sum;
  $totaldata{'gc'}=($totaldata{'C'}+$totaldata{'G'})/$sum  if ($sum>0);
  &output(\%totaldata,1);
}

&main;





