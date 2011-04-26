#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;


	sub main {
	  my $usage="USAGE:  ratewig /path/to/refseqfile/ /path/to/wigfile\n\tratewig --refseqfile=<filename> --wigfile=<filename>\n"; 
	  my $refseqfile;
	  my $wigfile;
	  my $help;
	  
	  if (scalar(@ARGV)!=2) { print "$usage\n"; exit(0); }
	  $refseqfile = $ARGV[0];
	  $wigfile = $ARGV[1];

	  if (@ARGV==0) {print "$usage\n"; exit(0);}

	  GetOptions('refseqfile=s' => \$refseqfile,
		     'wigfile=s' => \$wigfile,
		     'help!' => \$help );
	  if ($help) { print "$usage\n"; exit(0); }
	   
	   

	}

&main;
