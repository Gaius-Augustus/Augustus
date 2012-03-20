#!/usr/bin/perl

#############################################################
# filterGenes
# filter genes from a genbank flat file database
# usage: filterGenes namefile dbfile
#
#
# Mario Stanke, 13.08.2002
#############################################################

if ($#ARGV != 1) {
    print "usage:filterGenes.pl namefile dbfile\n";
    print "names of the genes to be filtered out come from\n";
    print "the first parameter.";
    exit;
} 
$origfilename = $ARGV[1];
$badfilename = $ARGV[0];
open(badfile, "<$badfilename") || die "Couldn't open name file";

while(<badfile>){
   /.*/; 
   $badlist{$&}=1; 
}


open(origfile, "<$origfilename") || die "Couldn't open dbfile\n";


$/="\n//\n";
while(<origfile>) {
    $gendaten=$_;
    m/^LOCUS +(\S+) .*/;
    $genname=$1;
    if (!exists $badlist{$genname}) {
	print "$gendaten";
    }
}
