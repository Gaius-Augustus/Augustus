#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

# a gff file looks like this:
#<seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
sub generateGff {
  my $data=shift;
  my $returnstring = "";
  $returnstring .= "$$data{'seqname'}\t";
  $returnstring .= "$$data{'source'}\t";
  $returnstring .= "$$data{'feature'}\t";
  $returnstring .= "$$data{'start'}\t";
  $returnstring .= "$$data{'end'}\t";
  $returnstring .= "$$data{'score'}\t";
  $returnstring .= "$$data{'strand'}\t";
  $returnstring .= "$$data{'frame'}\t";
  $returnstring .= "$$data{'additional'}"; 
  return $returnstring;
}


#this sub takes 3 parameters, 
# $fastastring: the fasta string we want to parse, 
# $data: a reference to a hash where we save all the results in
# $radius: we are calculating the start and the end values with this variable. 
sub parseFasta {
  my $fastaString=shift;
  my $data=shift;
  my $radius=shift;
  if ($fastaString =~ m/^(.*)-([\d]*)_([\+-])\s([\d]*)\stranscripts:\s([^,]*).*$/) {
    $$data{'seqname'} = $1;
    $$data{'source'} = "PASA";
    $$data{'feature'} = "tts";
    $$data{'start'} = $2-$radius;
    $$data{'end'} = $2+$radius;
    $$data{'score'} ="$4";
    $$data{'strand'} ="$3";
    $$data{'frame'} =".";
    my $temp=$5;
    if ($temp =~ m/(gi\|\d+\|).*/) {
      $temp = $1;
    }
    $$data{'additional'} ="pri=4;EST=$temp;src=E";
    return $data;
  }
  else {
    return undef;
  }
}

sub main {
  #carries the information from one format to the other
  my %shuttle=();   
  #setting the USAGE variable
  my $USAGE="USAGE: pasapolyA2hints.pl <PASA-Output file in fasta-format> [--radius=integer]\n";
  my $rad=5;
  my $pasapolyAfile="";
  
  GetOptions("radius:i"=> \$rad);


  if (scalar(@ARGV) < 1) {
    die $USAGE;
  } 
  else {
    $pasapolyAfile=$ARGV[0];
    open(FASTA, "<$pasapolyAfile") || die "Couldn't open $ARGV[0]\n";
    $/="\n>";
    while(<FASTA>) {
      /[>]*(.*)\n/;
      my $temp = $1;
      chomp($temp);
      &parseFasta($temp,\%shuttle,$rad);
      print(&generateGff(\%shuttle),"\n");
    }
  }

  

}

&main;

