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
    if ( !(@lf) ) {
        @lf = @f;
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
        summarizeHint(\%identical);
        undef %identical;
        @lf = @f;
        my @to_be_pushed = @f;
        push(@{$identical{$srcKey}}, \@to_be_pushed);
    }
    else {
    	my @to_be_pushed = @f;
    	push(@{$identical{$srcKey}}, \@to_be_pushed);
    }
}

summarizeHint(\%identical);

sub summarizeHint {
	my $hints = shift;
	foreach my $src (keys %{$hints}) {
		my @h = @{${$hints->{$src}}[0]};
		if(scalar(@{$hints->{$src}}) == 1 ) {
			print join ("\t", @h);
		}else{
			my $mult = 0;
			for (my $i = 0; $i < scalar (@{$hints->{$src}}); $i++ ) {
				my @l = @{${$hints->{$src}}[$i]};
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
