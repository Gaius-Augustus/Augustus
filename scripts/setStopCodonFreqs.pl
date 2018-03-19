#!/usr/bin/perl
# Author: Katharina J. Hoff
# Date: February 19th 2018
# set stop codon frequencies in parameters.cfg of a particular species

use strict;
use Getopt::Long;

my $usage = "setStopCodonFreqs.pl --AUGUSTUS_CONFIG_PATH=path --species=someName --etrainOut=file\n";

my $AUGUSTUS_CONFIG_PATH;
my $species;
my $etrainOut;
my $help;

GetOptions( 'AUGUSTUS_CONFIG_PATH=s'  => \$AUGUSTUS_CONFIG_PATH,
	    'species=s'              => \$species,
	    'etrainOut=s'              => \$etrainOut,
	    'help!'                  => \$help);
	    
if($help){
    print $usage;
    exit 0;
}

if(not(defined($AUGUSTUS_CONFIG_PATH))){
    if(-e $ENV{'AUGUSTUS_CONFIG_PATH'}){
	$AUGUSTUS_CONFIG_PATH=$ENV{'AUGUSTUS_CONFIG_PATH'};
    }else{
	print $usage;
	die("$AUGUSTUS_CONFIG_PATH missing!\n");	
    }		  
}

if(not(defined($species))){
    print $usage;
    die("No species defined!\n");
}

if(not(defined($etrainOut))){
    print $usage;
    die("No etraining output file defined!\n");
}

my $freqOfTag;
my $freqOfTaa;
my $freqOfTga;
open(ETRAIN, "<", $etrainOut) or die ("Could not open file $etrainOut!\n");
while(<ETRAIN>){
    if(/tag:\s*.*\((.*)\)/){     
	$freqOfTag = $1;         
    }elsif(/taa:\s*.*\((.*)\)/){ 
	$freqOfTaa = $1;         
    }elsif(/tga:\s*.*\((.*)\)/){ 
	$freqOfTga = $1;         
    }
}
close(ETRAIN) or die ("Could not close file $etrainOut!\n");

my $cfgFile = "$AUGUSTUS_CONFIG_PATH/species/$species/$species\_parameters.cfg";
if(-e "$cfgFile"){
    setParInConfig($cfgFile, "/Constant/amberprob", $freqOfTag);
    setParInConfig($cfgFile, "/Constant/ochreprob", $freqOfTaa);
    setParInConfig($cfgFile, "/Constant/opalprob", $freqOfTga);
}else{
    die("Species config file $cfgFile does not seem to exist!\n");
}

sub setParInConfig{
    my $configFileName = shift;
    my $parName = shift;
    my $value = shift;
    open(CFGFILE, "+<$configFileName") or die ("Could not read config file $configFileName\n");
    my @lines = <CFGFILE>;
    foreach my $line (@lines){
	$line =~ s/(\s*$parName +)(\S+?)(\s|\#|$)/$1$value$3/;
    }
    seek(CFGFILE, 0,0);
    print CFGFILE @lines or die ("Could not write $configFileName");
    truncate(CFGFILE, tell(CFGFILE));
    close(CFGFILE);
}
