#!/usr/bin/perl

# Katharina J. Hoff, May 26th 2015
#
# Convert a cegma output.gff file to a gff file that is compatible with AUGUSTUS training
#
# Format example CEGMA (v.2.4.010312):
# CHROMOSOME_I    cegma   First   636237  636337  38.40   -       0       KOG0328.7
# CHROMOSOME_I    cegma   Exon    636237  636337  38.40   -       .       KOG0328.7
# CHROMOSOME_I    cegma   Internal        633894  634800  641.47  -       1       KOG0328.7
# CHROMOSOME_I    cegma   Exon    633894  634800  641.47  -       .       KOG0328.7
# CHROMOSOME_I    cegma   Internal        632564  632587  -9.15   -       0       KOG0328.7
# CHROMOSOME_I    cegma   Exon    632564  632587  -9.15   -       .       KOG0328.7
# CHROMOSOME_I    cegma   Terminal        631097  631288  138.28  -       0       KOG0328.7
# CHROMOSOME_I    cegma   Exon    631097  631288  138.28  -       .       KOG0328.7
#
# Format example CEGMA (v.2.5):
# CHROMOSOME_I    cegma   First    636237  636337  38.40   -       .       KOG0328.7
# CHROMOSOME_I    cegma   Internal    633894  634800  641.47  -       .       KOG0328.7
# CHROMOSOME_I    cegma   Internal    632564  632587  -9.15   -       .       KOG0328.7
# CHROMOSOME_I    cegma   Terminal    631097  631288  138.28  -       .       KOG0328.7  

# Format example for AUGUSTUS:
# CHROMOSOME_I    cegma   CDS    636237  636337  38.40   -       .       transcript_id "g1.KOG0328.7"
# CHROMOSOME_I    cegma   CDS    633894  634800  641.47  -       .       transcript_id "g1.KOG0328.7"
# CHROMOSOME_I    cegma   CDS    632564  632587  -9.15   -       .       transcript_id "g1.KOG0328.7"
# CHROMOSOME_I    cegma   CDS    631097  631288  138.28  -       .       transcript_id "g1.KOG0328.7"

my $usage = "cegma2gff.pl cegma.gff > augustus.gff\n";

if(@ARGV!=1){print STDERR $usage; exit -1;}

my $geneCounter = 0;
my @t;

open(CEGMA, "<", $ARGV[0]) or die("Could not open CEGMA file $ARGV[0]!\n");
while(<CEGMA>){
    if(($_=~m/First/) or ($_=~m/Single/)){$geneCounter = $geneCounter + 1;}
    if(not($_=~m/Exon/)){
	chomp;
	@t = split(/\t/);
	print "$t[0]\t$t[1]\tCDS\t$t[3]\t$t[4]\t$t[5]\t$t[6]\t$t[7]\ttranscript_id \"g$geneCounter.$t[8]\"\n";
    }
}
close(CEGMA) or die("Could not close CEGMA file $ARGV[0]!\n");
