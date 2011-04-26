#!/usr/bin/perl -w
use strict;

die "Need two arguments" unless @ARGV >= 2;

my %todel;
foreach (split(",", $ARGV[1])) {
    $todel{$_}=1;
}
splice (@ARGV, 1);

my  @input;
while (<>) {
    if (/^\[/ || @input == 0) {
	push @input, $_;
    } else {
	$input[-1].=$_;
    }
}

my $index = 0;
my $blockno = 0;
my @output = ();
my $deleted_lines = 0;
my $ADDMODE = 1;


sub infsum {
    my $result = 0;
    foreach (@_) {
	return "*" if /^\*$/;
	print STDERR "Unexpected value: '$_'\n" unless /^\d+$/;
	$result += $_ if /^\d+$/;
    }
    return $result;
}
 
foreach (@input) {
    if (/^\[block\]/m) {
	my $name;
	$name=$1 if (/^name=(.*?)$/m);
	if ($ADDMODE && $deleted_lines > 0) {
	    # if deleted_lines > 0 => block was deleted but no [dist] followed
	    # if there is [dist] in output: add deleted_lines to <min> and set <max> to inf
	    # else push new output deleted_lines <max>

	    if (@output && $output[-1] =~ /^\[dist\]/m) {
		$output[-1] =~ s/^(\d+)(\s+)(\d+|\*)/($1 + $deleted_lines)."$2*"/me;
	    } else {
		push @output, "[dist]\n# distance to previous block\n# <min> <max>\n$deleted_lines\t*\n\n";
	    }
	}
	if ((defined $name && exists $todel{$name}) || exists $todel{$blockno}) {
	    my $blsize = 0;
	    while (/^(\d+)\s/mg && $blsize == $1) {
		$blsize++;
	    }
	    ### DEBUG part
	    if ($blsize == 0) {
		die "This file is corrupted (trying to delete empty block). Gave up";
	    }
	    $deleted_lines = $blsize;
	    ### end of DEBUG part
	    unless ($ADDMODE) {
		# ADDMODE == false
		while (@output && $output[-1] =~ /^\[dist\]/m) {
		    pop @output;
		}
	    }
	} else {
	    push @output, $_;
	    $deleted_lines = 0;
	}
	$blockno++;
    } elsif (/^\[dist\]/)  {
	if ($deleted_lines>0)  {
	    if ($ADDMODE) {
		# if deleted_lines > 0 && there is [dist] in output (must be
		# from before deleted block): add output and deleted_lines to input, replace output by input
		# if deleted_lines > 0 && no [dist] in output:
		# push current on output add deleted_lines to min and set max to inf
		my ($outmin, $outmax) = 
		    (@output && $output[-1] =~ /^\[dist\]/m) ?
		    ((pop @output) =~ /^(\d+)\s+(\d+|\*)/m) : (0, "*");
		### DEBUG part
		die "Undefined" unless defined $_;
		die "Bad syntax: $_" unless defined $outmax;
		s/^(.*\n)(\d+)(\s+)(\d+|\*)//s;
		my ($prefix, $inmin, $tab, $inmax) = ($1, $2, $3, $4);
		### DEBUG part
		die "Bad syntax: $_" unless (defined $inmax);
		push @output, $prefix.($outmin + $deleted_lines + $inmin).$tab.&infsum($outmax, $inmax, $deleted_lines).$_;
		$deleted_lines = 0;
	    }
	    # if ADDMODE==false: skip [dist] in output
	} else {
	    push @output, $_;
	}
    } else {
	push @output, $_;
    }
}

# this is needed if last block was deleted
if ($ADDMODE && $deleted_lines > 0) {
    # if deleted_lines > 0 => block was deleted but no [dist] followed
    # if there is [dist] in output: add deleted_lines to <min> and set <max> to inf
    # else push new output deleted_lines <max>
    
    if (@output && $output[-1] =~ /^\[dist\]/m) {
	$output[-1] =~ s/^(\d+)(\s+)(\d+|\*)/($1 + $deleted_lines)."$2*"/me;
    } else {
	push @output, "[dist]\n# distance to previous block\n# <min> <max>\n$deleted_lines\t*\n\n";
    }
}

print foreach @output;
	
    
