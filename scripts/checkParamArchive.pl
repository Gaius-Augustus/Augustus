#!/usr/bin/perl

# This script checks whether a parameter *.tar.gz archive contains valid AUGUSTUS parameter files. Errors due to missing files are written to STDERR. If UTR parameters are not present, this is written to STDOUT.
# param-archive.tar.gz is the archive in question
# projectDir is a directory to which the archive may be extracted (this script does not clean up after archive extraction!)
# This script assumes that it will produce the first directory in projectDir! Make sure there ain't any other files around!

my $usage = "checkParamArchive.pl param-archive.tar.gz projectDir > utrPossible 2> errors\n";

if(@ARGV != 2){
        print $usage;
        exit;
}

my $parArch = $ARGV[0];
my $projectDir = $ARGV[1];

# extract archive
`cd $projectDir; tar xzvf $parArch`;
my $speciesName = `ls $projectDir`;
chomp($speciesName);

# required files
my $prefix = "$projectDir/$speciesName/$speciesName";
my $metaparsCfg = $prefix."_metapars.cfg";
my $parsCfg = $prefix."_parameters.cfg";
my $weights = $prefix."_weightmatrix.txt";
my $intronProbs = $prefix."_intron_probs.pbl";
my $intronNoCRF = $prefix."_intron_probs.pbl.withoutCRF";
my $exonProbs = $prefix."_exon_probs.pbl";
my $exonProbsNoCRF = $prefix."_exon_probs.pbl.withoutCRF";
my $igenicProbs = $prefix."_igenic_probs.pbl";
my $igenicProbsNoCRF = $prefix."_igenic_probs.pbl.withoutCRF";
my $utrProbs = $prefix."_utr_probs.pbl";
my $utrCfg = $prefix."_metapars.utr.cfg";

my %shouldBeThere;
$shouldBeThere{$speciesName."_metapars.cfg"} = 1;
$shouldBeThere{$speciesName."_weightmatrix.txt"} = 1;
$shouldBeThere{$speciesName."_parameters.cfg"} = 1;
$shouldBeThere{$speciesName."_intron_probs.pbl"} = 1;
$shouldBeThere{$speciesName."_intron_probs.pbl.withoutCRF"} = 1;
$shouldBeThere{$speciesName."_exon_probs.pbl"} = 1;
$shouldBeThere{$speciesName."_exon_probs.pbl.withoutCRF"} = 1;
$shouldBeThere{$speciesName."_igenic_probs.pbl"} = 1;
$shouldBeThere{$speciesName."_igenic_probs.pbl.withoutCRF"} = 1;
$shouldBeThere{$speciesName."_utr_probs.pbl"} = 1;
$shouldBeThere{$speciesName."_metapars.utr.cfg"} = 1; 
$shouldBeThere{"."} = 1;
$shouldBeThere{".."} = 1;

# check whether all files exist
if(not(-e $metaparsCfg)){
	print STDERR "Metaparameter config file $metaparsCfg is missing!\n";
}
if(not(-e $parsCfg)){
	print STDERR "Parameter config file $parsCfg is missing!\n";
}
if(not(-e $weights)){
	print STDERR "Weight file $weights is missing!\n";
}
if(not(-e $intronProbs)){
	print STDERR "Intron probability file $intronProbs is missing!\n";
}
if(not(-e $intronNoCRF)){
	print STDERR "Intron probability file (without CRF) $intronNoCRF is missing!\n";
}
if(not(-e $exonProbs)){
	print STDERR "Exon probability file $exonProbs is missing!\n";
}
if(not(-e $exonProbsNoCRF)){
	print STDERR "Exon probability file (without CRF) $exonProbsNoCRF is missing!\n";
}
if(not(-e $igenicProbs)){
	print STDERR "Intergenic probability file $igenicProbs is missing!\n";
}
if(not(-e $igenicProbsNoCRF)){
	print STDERR "Intergenic probabiliy file (without CRF) $igenicProbsNoCRF is missing!\n";
}
# The following file is not required, because if UTR training was not performed, it does not exist!
#if(not(-e $utrProbs)){
#	print STDOUT "UTR probability file $utrProbs is missing!\n";
#}

# Check whether there are any other files that SHOULD NOT be there:

opendir(D, "$projectDir/$speciesName") || die "Can't open directory $projectDir/$speciesName!\n";
my @list = readdir(D);
closedir(D);

foreach my $f (@list) {
    if(not(defined($shouldBeThere{$f}))){
	print STDERR "File $f should not be in the archive!\n";
    }
}

# Check whether UTR predictions are impossible due to missing parameter files:
if(not(-e $utrProbs)){
    print STDOUT "$utrProbs is missing, disable UTR prediction!\n";
}elsif(not(-e $utrCfg)){
    print STDOUT "$utrCfg is missing, disable UTR prediction!\n";
}
