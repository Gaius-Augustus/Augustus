#!/usr/bin/perl

# filter augustus splice hints file (e.g. from RNA-seq data)
# for splice hints that contain a particular splicing motif,
# usually GT-AG
#
# Usage: filterSpliceHints.pl genome.fa hints.gff > filtered.hints.gff

# Katharina Hoff, 9.6.2011

my $usage = "filerIntrons.pl genome.fa hints.gff > out.gff\n";

if (@ARGV != 2) {
    print $usage;
    exit;
}

my $genome = $ARGV[0];
my $hints = $ARGV[1];
my $splice = "GTAG";

open (FASTA, $genome) or die ("\n\ncould not open file $genome!\n");
LINE: while ($line = <FASTA>){
    next LINE if $line =~ m/^#/; #discard comments
    if ($line =~ /^>/){
           $line =~ s/[\x0A\x0D]+//g; #removing those ugly whitelines
           $line =~ s/(\n)(\r)//g; #remove them alllll!
           $line =~ m/(^>\w+)/i; #matches a word starting with > (Fasta)
           $hash_key = substr($1, 1)
    }else{
    $line =~ s/[\x0A\x0D]+//g; 
    $line =~ s/(\s+)(\n)(\r)//g;
    $fasta_hash{$hash_key}.=$line; 
}
}
close(FASTA) or die "Could not close file $genome!\n";

open (HINTS, $hints) or die ("Could not open file $hints!\n");
LINE: while($line = <HINTS>){
    @gff = split(/\t/, $line);
    $siteA = substr($fasta_hash{$gff[0]}, ($gff[3]-1), 2);   
    $siteB = substr($fasta_hash{$gff[0]}, ($gff[4]-2), 2);
    $given = $siteA.$siteB;
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


