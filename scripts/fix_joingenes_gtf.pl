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
while(<STDIN>){
	if(m/^\#/){
		$comments .= $_;
	}else{
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
		}
		push(@{$genes{$tx_id}{'lines'}}, $_);
	}
}

foreach my $key (keys(%genes)) {
	print $genes{$key}{'comments'};
    print $genes{$key}{'seq'}."\t".$genes{$key}{'src'}."\tgene\t".$genes{$key}{'start'}."\t".$genes{$key}{'stop'}."\t.\t".$genes{$key}{'strand'}."\t.\t".$key."\n";
    foreach(@{$genes{$key}{'lines'}}){
    	print $_;
    }
}