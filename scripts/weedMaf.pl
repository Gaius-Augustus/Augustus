#!/usr/bin/perl
#
# weedMaf remove only-gap columns from a multiple alignment in MAF format
#
# input file format:
# 
# a score=5037
# s chlamy4.chromosome10    5424 36280 + 6579462 ATCA-CCACAG--ACC...
# s    volvox.scaffold_9 2188403 51128 + 9999999 ACCA-CCACGGGCACC...
# 
#
# 10.03.2013, Mario Stanke, mario.stanke@uni-greifswald.de

use strict;
use Getopt::Long;

my @alirows;         # array of all alignment rows (including gaps)
my @alidata;           # array of everying up to alignment rows

my $help = 0;
my $len;

GetOptions('help!'=>\$help);

exec("perldoc $0") if ($help);

while (<>) {
    if (!/^s\s/){
	my $line = $_;
	if (@alirows){
	    weed();
	    @alidata = @alirows = ();
	}
	print $line;
    } elsif (/(^s\s.*\s)(\S+)$/){
	push @alidata, $1;
	push @alirows, $2;
    }
}

if (@alirows){
    weed();
}

sub weed{
    $len = length($alirows[0]);
    my $k = @alirows;
    my @newalirows = ();
    foreach my $alirow (@alirows){
	if (length($alirow) != $len){
	    die ("Inconsistent alignment lengths");
	}
	push @newalirows, "";
    }
    for (my $i=0; $i<$len; $i++){
	my $onlygaps = 1;
	for (my $j=0; $j<@alirows && $onlygaps; $j++){
	    if (substr($alirows[$j], $i, 1) ne '-'){
		$onlygaps = 0;
	    }
	}
	if (!$onlygaps){
	    for (my $j=0; $j<@alirows; $j++){
		$newalirows[$j] .= substr($alirows[$j], $i, 1);
	    }
	}
    }
    for (my $j=0; $j<@alirows; $j++){
	print $alidata[$j] . $newalirows[$j] . "\n";
    }
    
}


__END__

=head1 NAME

weedMaf.pl remove only-gap columns from a multiple alignment in MAF format

=head1 SYNOPSIS

weedMaf.pl < in.maf > out.maf

=cut
