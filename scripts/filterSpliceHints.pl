#!/usr/bin/perl

# filter augustus splice hints file (e.g. from RNA-seq data)
# for splice hints that contain a particular splicing motif,
# usually GT-AG
#
# Usage: filterSpliceHints.pl genome.fa hints.gff > filtered.hints.gff

# Katharina Hoff, 9.6.2011

my $usage = "filterSpliceHints.pl genome.fa hints.gff splice-pattern> out.gff\n\nThe most typical splice pattern should be GTAG! Sometimes, GCAG is also used.\n";

if (@ARGV != 3) {
    print $usage;
    exit;
}

my $genome = $ARGV[0];
my $hints = $ARGV[1];
my $splice = $ARGV[2];

open (FASTA, $genome) or die ("\n\ncould not open file $genome!\n");
LINE: while ($line = <FASTA>){
    next LINE if $line =~ m/^#/; #discard comments
    if ($line =~ /^>/){
	   chomp($line);
           #$line =~ s/[\x0A\x0D]+//g; #removing those ugly whitelines
           #$line =~ s/(\n)(\r)//g; #remove them alllll!
           #$line =~ m/(^>\w+)/i; #matches a word starting with > (Fasta)
           
	   $hash_key = substr($line, 1, length($line)-1)
    }else{
    $line =~ s/[\x0A\x0D]+//g; 
    $line =~ s/(\s+)(\n)(\r)//g;
    $line = uc($line);
#    print "Hash key: $hash_key\n";
#    print "Content: $line\n";
    $fasta_hash{$hash_key}.=$line; 
}
}
close(FASTA) or die "Could not close file $genome!\n";

open (HINTS, $hints) or die ("Could not open file $hints!\n");
LINE: while($line = <HINTS>){
    @gff = split(/\t/, $line);
#    print "Scaffold: ".$gff[0]."\n";
#    print "Length: ".length($fasta_hash{$gff[0]})."\n";
    $siteA = substr($fasta_hash{$gff[0]}, ($gff[3]-1), 2);   
    $siteB = substr($fasta_hash{$gff[0]}, ($gff[4]-2), 2);
    $given = $siteA.$siteB;
    #print "Splice site: $given\n";
    if($given =~ m/$splice/){
	print $gff[0]."\t".$gff[1]."\t".$gff[2]."\t".$gff[3]."\t".$gff[4]."\t".$gff[5]."\t+\t".$gff[7]."\t".$gff[8];
    }else{
	$given = reverse $given;
	$given =~ tr/ACGTacgt/TGCAtgca/;
	if($given =~ m/$splice/){
	    print $gff[0]."\t".$gff[1]."\t".$gff[2]."\t".$gff[3]."\t".$gff[4]."\t".$gff[5]."\t-\t".$gff[7]."\t".$gff[8];
	}
    }
}
close(HINTS) or die "Could not close file $hints!\n";


