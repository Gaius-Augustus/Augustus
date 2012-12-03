#!/usr/bin/perl

# Katharina J. Hoff, Dec 3rd 2012
#
# Simply fasta headers for usage with web service or autoAug.pl

my $usage = "simplifyFastaHeaders.pl in.fa nameStem out.fa header.map\n\nin.fa is the input fasta sequence\nnameStem is the new beginning of each fasta entry\nout.fa is the output fasta file with simplified headers\nheader.map contains the new and old headers in tabulator separated format\n\n";

if(@ARGV!=4){print STDERR $usage; exit -1;}

my $inputFile = $ARGV[0];
my $nameStem = $ARGV[1];
my $outFile = $ARGV[2];
my $mapFile = $ARGV[3];

my $c = 0;
my $emptyC = 0;
my $wrongNL = 0;
my $prot = 0;
my $dna = 0;

open(INPUT, "<", $inputFile) or die("Could not open input fasta file $inputFile!\n");
open(OUTPUT, ">", $outFile) or die("Could not open output fasta file $outFile, check whether directory exists and check writing permissions!\n");
open(MAP, ">", $mapFile) or die("Could not open map file $mapFile, check whether directory exists and check writing permission!\n");
while(<INPUT>){
    if(not($_=~m/\n$/)){
        if($wrongNL < 1){
            print STDERR "Warning: something seems to be wrong with the newline character! This is likely to cause problems with the autoAug.pl pipeline and the AUGUSTUS web service! Please adapt your file to UTF8! This warning will be supressed from now on!\n";
            $wrongNL++;
        }
    }
    if(m/^>/){
	$c++;
	print OUTPUT ">$nameStem$c\n";
	print MAP ">$nameStem$c\t$_";
    }else{
	if(length($_)>1){	    
	    chomp;
	    print OUTPUT $_."\n";
	    if($_=!m/[ATGCNatgcn]/){
		if($dna==0){
		    print STDOUT "Assuming that this is not a DNA fasta file because other characters than A, T, G, C, N, a, t, g, c, n were contained. If this is supposed to be a DNA fasta file, check the content of your file! If this is supposed to be a protein fasta file, please ignore this message.\n";
		    $dna++;
		}
	    }
	    if($_=!m/[AaRrNnDdCcEeQqGgHhIiLlKkMmFfPpSsTtWwYyVvBbZzJjXx]/){
		if($prot==0){
		    print STDOUT "Assuming that this is not a protein fasta file because other characters than AaRrNnDdCcEeQqGgHhIiLlKkMmFfPpSsTtWwYyVvBbZzJjXx were contained. If this is supposed to be DNA fasta file, please ignore this message.\n";
		    $prot++;
		}
	    }
	}else{
	    if($emptyC < 1){
		print STDERR "Warning: empty line was removed! This warning will supressed from now on!\n";
	    } 
	    $emptyC++;
	}
    }
}
close(INPUT) or die("Could not close input fasta file $inputFile!\n");
close(OUTPUT) or die("Could not close output fasta file $outFile!\n");
close(MAP) or die("Could not close map file $mapFile!\n");
