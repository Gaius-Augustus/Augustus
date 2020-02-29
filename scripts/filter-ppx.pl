#!/usr/bin/env perl

use strict;
use warnings;
my $mode = 0;
my $pending = "";

while (<>) {
    if (/^# start gene /) {
	$mode=1;
    } elsif ($mode==3) {
	$mode=0;
	next if /^###/;
    }
    if ($mode > 0) {
	$pending .= $_;
	if ($mode == 1 && /\tprotein_match\t/) {
	    $mode=2;
	}
	if (/^# end gene /) {
	    if ($mode == 2) {
		print $pending;
		$mode=0;
	    } else {
		$mode=3;
	    }
	    $pending="";
	}
    } else {
	print;
    }
}
    
	
