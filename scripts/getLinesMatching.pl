#!/usr/bin/perl

# Author: Mario Stanke
# Contact: mario.stanke@uni-greifswald.de
# Last modification: May 1st 2018

use Getopt::Long;
use strict;

my $complement = 0;
my $whitespace = 0;
my $patfrom    = "";
my $patto      = "";

GetOptions(
    'v!'          => \$complement,
    'patfrom:s'   => \$patfrom,
    'patto:s'     => \$patto,
    'whitespace!' => \$whitespace
);

if ( @ARGV != 2 ) {
    print "Efficient script for printing only those lines where the n-th "
        . " column is a word in the match.lst file\n\n";
    print "Usage:\n\n";
    print "cat input | $0 match.lst n > output\n\n";
    print "Columns are based on tab-separation and are 1-based.\n";
    print "Options:\n";
    print "  --v               Use complement. Print all lines NOT matching.\n";
    print "  --patfrom --patto Apply query replace regular expression to the \n"
        . "                    entry in the n-th column first before \n"
        . "                    checking for membership in the list.\n"
        . "                    Will use the pattern s/\$patfrom/\$patto/\n"
        . "                    Useful for removing modifications, e.g.\n"
        . "                    tripping a trailing -1 before performing the\n"
        . "                    check.\n"
        . "                    Both default to empty patterns.\n";
    print "  --whitespace      Split columns at whitespace rather than tab.\n";
    exit(1);
}

my $matchfile = $ARGV[0];
my $n         = $ARGV[1];
open( MATCH, "<$matchfile" ) or die("Could not open $matchfile\n");

my %include = ();

while (<MATCH>) {
    chomp;
    $include{$_}++;
}
close MATCH;

my @f;
my $field;
while (<STDIN>) {
    chomp;
    if ($whitespace) {
        @f = split( /\s+/, $_, $n + 1 );
    } else {
        @f = split( /\t/, $_, $n + 1 );
    }
    $field = $f[ $n - 1 ];
    $field =~ s/$patfrom/$patto/oee if ( $patto ne "" || $patfrom ne "" );
    print $_ . "\n"
        if ( ( !$complement && $include{$field} )
        || ( $complement && !$include{$field} ) );
}
