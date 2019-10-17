#!/usr/bin/env perl

####################################################################################################
#                                                                                                  #
# fix_joingenes_gtf.pl                                                                             #
# Script that adds missing gene feature line to the output of joingenes (gtf format)               #
# This script should be removed when joingenes output format has been fixed!                       #
#                                                                                                  #
# Author: Katharina Hoff                                                                           #
#                                                                                                  #
# Contact: katharina.hoff@uni-greifswald.de                                                        #
#                                                                                                  #
# Acknowledgements: I thank Martijn Holterman for suggesting to print introns in the output file.  #
#                                                                                                  #
# This script is under the Artistic Licence                                                        #
# (http://www.opensource.org/licenses/artistic-license.php)                                        #
####################################################################################################

use Getopt::Long;
use strict;
use warnings;

my $usage = <<'ENDUSAGE';

fix_joingenes_gtf.pl     add missing gene features lines to joingenes output file

SYNOPSIS

fix_joingenes_gtf.pl [OPTIONS] < in.gtf > out.gtf

OPTIONS

--help	Display this help message

ENDUSAGE

my $help;

GetOptions(
    'help!' => \$help
);

if ($help) {
	print $usage;
    exit(0);
}

my %genes;
my $comments = "";
my $tx_id;
my @previous;
while(<STDIN>){
	if(m/^\#/){
		$comments .= $_;
		@previous = ();
	}else{
		chomp;
		my @t = split(/\t/);
		my @t2 = split(/\./, $t[8]);
		if(m/\ttranscript\t/){
			$tx_id = $t2[0];
			if(not(defined($genes{$tx_id}))){
				$genes{$tx_id}{'seq'} = $t[0];
				$genes{$tx_id}{'src'} = $t[1];
				$genes{$tx_id}{'strand'} = $t[6];
				$genes{$tx_id}{'comments'} = $comments;
				$comments = "";
			}
			if(defined($genes{$tx_id}{'start'})){
				if($genes{$tx_id}{'start'} > $t[3]){
					$genes{$tx_id}{'start'} = $t[3];
				}
			}else{
				$genes{$tx_id}{'start'} = $t[3];
			}
			if(defined($genes{$tx_id}{'stop'})){
				if($genes{$tx_id}{'stop'} < $t[4]){
					$genes{$tx_id}{'stop'} = $t[4];
				}
			}else{
				$genes{$tx_id}{'stop'} = $t[4];
			}
			push(@{$genes{$tx_id}{'lines'}}, $t[0]."\t".$t[1]."\t".$t[2]."\t".$t[3]."\t".$t[4]."\t".$t[5]."\t".$t[6]."\t".$t[7]."\ttranscript_id \"".$t[8]."\"; gene_id \"".$tx_id."\"");
		}elsif(m/\tCDS\t/){
			if ($previous[2] eq "CDS"){
				my $start = $previous[4] + 1;
				my $end = $t[3] - 1;
				push(@{$genes{$tx_id}{'lines'}}, $t[0]."\t".$t[1]."\tintron\t".$start."\t".$end."\t.\t".$t[6]."\t.\t".$t[8]);
			}
			push(@{$genes{$tx_id}{'lines'}}, $_);
			push(@{$genes{$tx_id}{'lines'}}, $t[0]."\t".$t[1]."\texon\t".$t[3]."\t".$t[4]."\t.\t".$t[6]."\t.\t".$t[8]);
		}else{
			push(@{$genes{$tx_id}{'lines'}}, $_);
		}
		@previous = @t;
	}
}

foreach my $key (keys(%genes)) {
	print $genes{$key}{'comments'};
    print $genes{$key}{'seq'}."\t".$genes{$key}{'src'}."\tgene\t".$genes{$key}{'start'}."\t".$genes{$key}{'stop'}."\t.\t".$genes{$key}{'strand'}."\t.\t".$key."\n";
    foreach(@{$genes{$key}{'lines'}}){
    	print $_."\n";
    }
}