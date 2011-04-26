#!/usr/bin/perl
#
# run sim4 once for each line in the matchlist file 
# 
#
my @contigs = ();

system ("rm -f sim4.results");
while (<>) {
    chomp;
    @f = split /\t/, $_;
    $estfilename = shift @f;
    system ("echo '$estfilename' | cdbyank ~/galdieria/galdi3/Newbler_all_454AllContigs.fna.long.uc.cidx > mrna.fa");
    system ("echo '@f' | cdbyank ~/galdieria/galdi3/build3.fa.cidx > genomic.fa");
    system ("perl -i.orig -p -e 's/^Incorrectly.*fixed.\n//' mrna.fa genomic.fa");
    system ("keepOnlyLonger.pl mrna.fa genomic.fa");
    if (open(LONG, "<genomic.fa.long")) {
	close LONG;
	system ("echo '--------\n>$estfilename' >> sim4.results");
	system ("sim4 mrna.fa genomic.fa.long >> sim4.results");
    } 
}
