#!/usr/bin/perl
#
# summarize multiple identical hints to one with mult=n
#
# This script is used by the braker.pl pipeline.
# Please be extremely careful when changing this script because the braker.pl
# pipeline may fail upon custom modification of this script.
# In case of doubt, contact katharina.hoff@uni-greifswald.de
#
# Mario Stanke, 4.1.2010, last modification by Katharina J. Hoff on March 1st 2018

use strict;
use Getopt::Long;

my $usage = "$0 -- summarize multiple identical hints to one with mult=n\n"
	 	. "\n"
		. "Usage: $0 <in.hints >joined.hints\n"
		. "  PREREQUISITE: input GFF file must be sorted so that hints that should be summarized are below each other\n"
		. "  e.g. do a cat hints.gff | sort -n -k 4,4 | sort -s -n -k 5,5 | sort -s -k 3,3 | sort -s -k 1,1 | join_mult_hints.pl\n"
		. "  hints are only joined if they are from the same src (e.g. src=E)\n"
		. "  score in column 6 is adapted to mult= value for usage with GeneMark\n";

my $help = 0;
GetOptions( 'help!' => \$help );
if ($help) {
    print "$usage";
    exit(0);
}

my @lf;
my ( $lm, $m );
my %identical; # holds all hints with identical coordinates, locus, strand & frame

while ( <STDIN> ) {
    my @f = split( /\t/, $_ );
    my $srcKey;
    if($f[8] =~ m/so?u?rce?=(\w)/){
    	$srcKey = $1;
    }else{
    	$srcKey = 'no_src';
    }
    print $srcKey ." is $srcKey\n";

    if ( !(@lf) ) {
        @lf = @f;
        print "Source key when pushing is $srcKey\n";
        print "Pushing:\n";
        print join ("\t", @f);
        my @to_be_pushed = @f;
    	push(@{$identical{$srcKey}}, \@to_be_pushed);
    }
    elsif (
        !(     ( $f[0] eq $lf[0] )
            && ( $f[2] eq $lf[2] )
            && ( $f[3] == $lf[3] )
            && ( $f[4] == $lf[4] )
            && ( $f[6] eq $lf[6] )
            && ( $f[7] eq $lf[7] )
        )
        )
    {
    	#print "The script believes that hint\n";
    	#print join ("\t", @f);
    	#print "and\n";
    	#print join ("\t", @lf);
    	#print "are not identical and proceeds to summarize\n";
    	#print "At this point, the keys are:\n";
    	my @s = keys %identical;
    	foreach(@s) {
    		print "key $_\n";
    	}
        summarizeHint(\%identical);
        undef %identical;
        @lf = @f;
    }
    else {
    	#print "script believes that hint\n";
    	#print join ("\t", @f);
    	#print "and\n";
    	#print join ("\t", @lf);
    	#print "are identical and pushes according to key $srcKey\n";
    	print "Source key when pushing is $srcKey\n";
    	print "Pushing:\n";
        print join ("\t", @f);
    	my @to_be_pushed = @f;
    	push(@{$identical{$srcKey}}, \@to_be_pushed);
    }
}

summarizeHint(\%identical);

sub summarizeHint {
	my $hints = shift;
	print "Script is in summarizeHints\n";
	print "Number of keys is ". scalar(keys %{$hints})."\n";
	foreach my $src (keys %{$hints}) {
		print "Source key when processing is $src\n";
		my @h = @{${$hints->{$src}}[0]};
		# if there is only one hint, print exactly as it was
		if(scalar(@{$hints->{$src}}) == 1 ) {
			print join ("\t", @h);
		}else{
			my $mult = 1;
			for (my $i = 1; $i < scalar (@{$hints->{$src}}); $i++ ) {
				my @l = @{${$hints->{$src}}[$i]};
				print "Field 8 is $l[8]\n";
				if($l[8] =~ m/mult=(\d+)/) {
					$mult += $1;
				}else{
					$mult++;
				}
			}
			$h[8] =~ s/gro?u?p=[^;]*;//;
			$h[8] =~ s/mult=\d+;//;
			for (my $i = 0; $i < 9; $i++) {
				if( not($i==5) && not ($i==8) ) {
					print $h[$i]."\t";
				}elsif($i==5){
					print $mult."\t";
				}elsif($i==8){
					print "mult=$mult;$h[$i]";
				}
			}
		}
	}
}