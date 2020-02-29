#!/usr/bin/env perl

use strict;
use warnings;

my $offset=0;
my @output;

if (@ARGV) {
    $offset = shift @ARGV;
}

my $target;
my ($from, $to) = ($offset, $offset+10000);
my $name="AUGUSTUS";

while (<>) {
    chomp;
    my $comment = "";
    if (/^# .* --predictionStart=(\d+) --predictionEnd=(\d+)/) {
	($from, $to) = ($1 + $offset -1000, $2 + $offset +1000);
	$from = 0 if $from < 0;
    }
    if (/^# .* --proteinprofile/) {
	$name .= "-PPX";
    }
    $comment = $1 if s/\s*\#.*$//;
    my @line = split /\t/;
    next unless @line == 9;
    next if $line[2] eq "interblock_region";
    if ($line[8] =~ /ID=([^;]*)/) {
	$line[8]=$1;
    } else {
	next;
    }
    $target = $line[0] unless defined $target;
    push @output,join("\t", @line)."\n";
}

print "browser position $target:$from-$to hide all\n";
print "track offset=$offset name=\"$name\" description=\"$name prediction\" visibility=3\n";
print foreach @output;










