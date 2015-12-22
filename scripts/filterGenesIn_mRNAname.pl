#!/usr/bin/perl

#############################################################
# filterGenes
# filter genes from a genbank flat file database
# usage: fileterGenesIn_mRNAname.pl namefile dbfile
#
#
# Mario Stanke, Simone Lange, Katharina Hoff; 21.12.2015
#############################################################

use strict;
use warnings;

if ($#ARGV != 1) {
    print "usage:filterGenes namefile dbfile\n";
    print "names of the loci to be kept come from\n";
    print "the first parameter. Only the the first of identical loci is kept\n";
    exit;
} 
my $origfilename = $ARGV[1];
my $goodfilename = $ARGV[0];

my %goodids;
open(GOODFILE, "<", "$goodfilename") || die "Couldn't open goodfile $goodfilename\n";
while(<GOODFILE>) {
            if($_ =~ m/transcript_id \"(.*)\"/) {                                                         
                $goodids{$1} = 1;                                                                         
            }                                                                                             
}
close(GOODFILE) || die "Couldn't close goodfile $goodfilename!\n";

open(my $ORIGFILE, "$origfilename") || die "Couldn't open dbfile\n";
my @data = <$ORIGFILE>;
close($ORIGFILE);


$/="\n//\n";

my $head;
my $mRNAflag = 0;
my $cdsFlag = 0;
my $genename;
my $printFlag = 0;
my $firstPrintFlag = 0;

foreach(@data) {
    if($_=~m/^LOCUS/){
	$head = "";
	$printFlag = 0;
	$genename = "";
	$head = $head.$_;
    }
    if($_=~m/FEATURES/){
	$head = $head.$_;
    }
    if($_=~m/source/){
	$head = $head.$_;
    }
    if($mRNAflag==1 and not($_=~m/CDS/)){
	$head = $head.$_;
    }
    if($_=~m/mRNA/){
	$mRNAflag = 1;
	$head = $head.$_;
    }
    if($cdsFlag==1){
        if($_=~m/gene="/){
		my @tmp = split(/\"/);
		$genename = $tmp[1];
		$cdsFlag = 0;
		$firstPrintFlag = 1;
        }else{
		$head = $head.$_;
        }
    }
    if($_=~m/CDS/){
	$mRNAflag = 0;
        $head = $head.$_;
	$cdsFlag = 1;
    }
    if($firstPrintFlag==1 and length($head)>=2){
	if($goodids{$genename}){
		print $head;
		$head = "";
		$printFlag = 1;
        }
	$firstPrintFlag = 0;
    }
    if($printFlag==1){
	print $_;
    }
}




	        
		

