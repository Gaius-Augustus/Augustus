#!/usr/bin/perl
#
# make for each matching est contig a list of matching genomic contigs in the same line
#

$oldestcontigname = "---begin---";
while (<>) {
    chomp;
    @f = split /\t/, $_;
    $estcontigname = $f[0];
    $genomiccontigname = $f[1];
    if ($estcontigname ne $oldestcontigname) {
	if ($oldestcontigname ne "---begin---") {
	    print "\n";
	}
	print $estcontigname;
    }
    print "\t$genomiccontigname";
    $oldestcontigname = $estcontigname;
}
