#!/usr/bin/perl

# "Clean" a MS/DOS fasta file from weird characters that screw with e.g. Perl

my $usage = "cleanDOSfasta.pl inFile > outFile\n";

if(@ARGV != 1){
	print $usage;
	exit;
}

my $inFile = $ARGV[0];

open(IN, "<", $inFile) or die "Could not open file $inFile\n!";
LINE: while($line = <IN>){
	next LINE if $line =~ m/^#/; #discard comments
	$line =~ s/[\x0A\x0D]+//g; # remove all those ugly whitespaces
	$line =~ s/(\n)(\r)//g; #remove them alllll!
	print $line."\n";
}
close(IN) or die "Could not close file $inFile!\n";
