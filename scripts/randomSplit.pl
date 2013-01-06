#!/usr/bin/perl

#############################################################
# randomSplit
# randomly split a genbank file
#
# usage: randomSplit dbfile size 
#
# dbfile: genbank file containing the genes
# size: size of the test set
# output: two files are created with names ending in .train 
#         and .test, being a training set and test set, respectively.
#
#
# Mario Stanke, 24.06.2002
#############################################################

srand 4;

if ($#ARGV != 1) {
    print "$0: randomly split a genbank file in two subsets of given sizes\n";
    print "usage: randomSplit dbfile size\n";
    exit;
}
$dbfilename = $ARGV[0];
$size = $ARGV[1];

open (STDIN, $dbfilename);

@list = <STDIN>;

@namelines = grep /^LOCUS   +/, @list;

if ($size > @namelines) {
    print "size $size is greater than the number of genes in file\n",
	"$dbfilename. Aborting.\n\n";
    exit;
}

my %unique = ();

foreach (@namelines) {
    /LOCUS +([^ ]+) */;
    #print "$1\n";
    push @names, "$1";
    if(not(defined($unique{$1}))){$unique{$1} = 1;}else{
	die( "ERROR in randomSplit.pl line 47: LOCUS names in genbank file are not unique!\n");
    }
}

%testnames=();
while ($size--) {
  $rand = rand (@names);

  $testnames{$names[$rand]}=1; 
  splice @names, $rand, 1;          # delete array element
}

open (TRAINFILE, ">${dbfilename}.train");
open (TESTFILE, ">${dbfilename}.test");

open (STDIN, $dbfilename);

#print "random test set:\n";

$/="\n//\n";
while(<STDIN>) {
    $gendaten=$_;
    m/^LOCUS +(\S+) .*/;
    $genname=$1;
    
    if (exists($testnames{$genname})) {
	#print "$genname\n";
	print TESTFILE "$gendaten";
    } else {
	print TRAINFILE "$gendaten";
    }
}
