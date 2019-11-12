#!/usr/bin/env perl

#############################################################################
# filterGenesIn_mRNAname.pl
# filter genes from a genbank flat file database
# those genes whose mRNA names are given in a gtf file
# are printed to STDOUT.
#
# This script is used by the braker.pl pipeline.
# Please be extremely careful when changing this script because the braker.pl
# pipeline may fail upon custom modification of this script.
# In case of doubt, contact katharina.hoff@uni-greifswald.de 
#
# usage: fileterGenesIn_mRNAname.pl gtffile dbfile
#
#
# Mario Stanke, Simone Lange, Katharina Hoff; 20.02.2018
############################################################################

use strict;
use warnings;

if ( $#ARGV != 1 ) {
    print "usage: filterGenesIn_mRNAname.pl gtffile dbfile\n\n";
    print "gtffile         genes to be kept in genbank file in gtf format\n";
    print "                (transcript_id "...") is essential.\n";
    print "dbfile          genbank file\n\n";
    print "Only the the first of identically named RNA loci is kept\n";
    exit;
}
my $origfilename = $ARGV[1];
my $goodfilename = $ARGV[0];

my %goodids;
open( GOODFILE, "<", "$goodfilename" )
    || die "Couldn't open goodfile $goodfilename\n";
while (<GOODFILE>) {
    if ( $_ =~ m/transcript_id \"([^"]*)\"/ ) {
        if (not (defined $goodids{$1})){
	    $goodids{$1} = 1;
	}
    }
}
close(GOODFILE) || die ( "Couldn't close goodfile $goodfilename!\n" );

open( my $ORIGFILE, "$origfilename" ) || die ( "Couldn't open dbfile $origfilename!\n" );
my @data = <$ORIGFILE>;
close($ORIGFILE) || die ( "Couldn't close dbfile $origfilename!\n" );

$/ = "\n//\n";

my $head;
my $mRNAflag = 0;
my $cdsFlag  = 0;
my $genename;
my $printFlag      = 0;
my $firstPrintFlag = 0;

foreach (@data) {
    if ( $_ =~ m/^LOCUS/ ) {
        $head      = "";
        $printFlag = 0;
        $genename  = "";
        $head      = $head . $_;
    }
    if ( $_ =~ m/^FEATURES/ ) {
        $head = $head . $_;
    }
    if ( $_ =~ m/\s+source/ ) {
        $head = $head . $_;
    }
    if ( $mRNAflag == 1 and not( $_ =~ m/CDS/ ) ) {
        $head = $head . $_;
    }
    if ( $_ =~ m/\s+mRNA/ ) {
        $mRNAflag = 1;
        $head     = $head . $_;
    }
    if ( $cdsFlag == 1 ) {
        if ( $_ =~ m/\s+\/gene="/ ) {
            my @tmp = split(/\"/);
            $genename       = $tmp[1];
            $cdsFlag        = 0;
            $firstPrintFlag = 1;
        }
        else {
            $head = $head . $_;
        }
    }
    if ( $_ =~ m/\s+CDS/ ) {
        $mRNAflag = 0;
        $head     = $head . $_;
        $cdsFlag  = 1;
    }
    if ( $firstPrintFlag == 1 and length($head) >= 2 ) {
        if ( $goodids{$genename} ) {
            print $head;
            $head      = "";
            $printFlag = 1;
        }
        $firstPrintFlag = 0;
    }
    if ( $printFlag == 1 ) {
        print $_;
    }
}

