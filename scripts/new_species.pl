#!/usr/bin/perl
#
# newspecies.pl
# Create the parameter files necessary for training AUGUSTUS for a new species.
#
# Mario Stanke, 20.12.06

use strict;
use Getopt::Long;

my $usage = "$0 -- create the parameter files necessary for training AUGUSTUS for a new species.\nUsage:\n";
$usage .= "$0 --species=myspecies\n";
$usage .= "myspecies                  species name, prefix of the parameter files\n";
$usage .= "options:\n";
$usage .= "--AUGUSTUS_CONFIG_PATH=dir full path to augustus/config directory\n";
$usage .= "--silent                   suppress help messages\n";
$usage .= "--ignore                   don't do anything if species already exists (default: off)\n";
$usage .= "--prokaryotic              use prokaryotic template instead of eukaryotic template\n";

my $species;
my $silent;
my $AUGUSTUS_CONFIG_PATH;
my $ignore;
my $prokaryotic;

##############################################################
# Check the command line
##############################################################

if ($#ARGV<0) {
    print "$usage";
    exit;
}

GetOptions( 'species=s' => \$species, 'AUGUSTUS_CONFIG_PATH=s' => \$AUGUSTUS_CONFIG_PATH, 'silent!' => \$silent,
    'ignore!' => \$ignore, 'prokaryotic!' => \$prokaryotic);


if (!defined($species) || $species eq ""){
    print "no species specified\n$usage";
    exit;
}

my $configdir;
if (length($AUGUSTUS_CONFIG_PATH) > 0) {
    $configdir = $AUGUSTUS_CONFIG_PATH;
} else {
    exists($ENV{AUGUSTUS_CONFIG_PATH}) or die("Environment variable AUGUSTUS_CONFIG_PATH not set.");
    $configdir = $ENV{AUGUSTUS_CONFIG_PATH};
}
if ($configdir !~ /\/$/){
    $configdir .= "/";
}

my $speciesdir = $configdir . "species/" . $species . "/";
my $cfgfilename = $speciesdir . $species . "_parameters.cfg";
my $weightfilename = $speciesdir . $species . "_weightmatrix.txt";
my $metafilename = $speciesdir . $species . "_metapars.cfg";
my $metautrfilename = $speciesdir . $species . "_metapars.utr.cfg";
my $metacgpfilename = $speciesdir . $species . "_metapars.cgp.cfg";
my $ovlpfilename = $speciesdir . $species . "_ovlp_len.pbl";
my $transshadowfilename = $speciesdir . $species ."_trans_shadow_bacterium.pbl";


##############################################################
# check whether the config files already exist
##############################################################
my $exists =0;
if (stat $cfgfilename > 0){ 
    print "$cfgfilename already exists.\n" if (!$ignore);
    $exists = 1;
}
if (stat $weightfilename > 0){ 
    print "$weightfilename already exists.\n" if (!$ignore);
    $exists = 1;
}
if (stat $metafilename > 0){ 
    print "$weightfilename already exists.\n" if (!$ignore);
    $exists = 1;
}
if ($exists) {
    if ($ignore){
	print "Already have files for species $species. Ignoring.\n" if (!$silent);
	exit(0);
    } else {
	print STDERR "If you are sure you want to delete the existing parameters for $species, delete the files ";
	print STDERR "for $species and run again.\nOtherwise use another name for that species.\n";
	exit(1);
    }
}

##############################################################
# create the new files
##############################################################

# info prokaryotic/eukaryotic
if(!$prokaryotic){
    print "Will create parameters for a EUKARYOTIC species!\n";
}else{
    print "Will create parameters for a PROKARYOTIC species!\n";
}

# directory

if (stat $speciesdir == 0){
    print "creating directory $speciesdir ...\n";
    mkdir $speciesdir;
}

# config file

print "creating $cfgfilename ...\n";
if(!$prokaryotic){
    open (GENERIC, "<$configdir/species/generic/generic_parameters.cfg") or die ("Could not open $configdir/species/generic/generic_parameters.cfg.");
    open (CFG, ">$cfgfilename") or die ("Could not write $cfgfilename.");
    while (<GENERIC>) {
	s/generic/$species/;
	s/use as template for your own species//;
	print CFG;
    }
}else{
    open (PROK, "<$configdir/species/template_prokaryotic/template_prokaryotic_parameters.cfg") or die ("Could not open $configdir/species/template_prokaryotic/template_prokaryotic_parameters.cfg.");
    open (CFG, ">$cfgfilename") or die ("Could not write $cfgfilename.");
    while (<PROK>) {
        s/template_prokaryotic/$species/;
        s/use as template for your own species//;
        print CFG;
    }

}

# weightmatrix
print "creating $weightfilename ...\n";
if(!$prokaryotic){
    if (stat "$configdir/species/generic/generic_weightmatrix.txt" == 0){
	die ("$configdir/species/generic/generic_weightmatrix.txt doesn't exist.");
    }
    system ("cp $configdir/species/generic/generic_weightmatrix.txt $weightfilename");
}else{
    if (stat "$configdir/species/template_prokaryotic/template_prokaryotic_weightmatrix.txt" == 0){
        die ("$configdir/species/template_prokaryotic/template_prokaryotic_weightmatrix.txt doesn't exist.");
    }
    system ("cp $configdir/species/template_prokaryotic/template_prokaryotic_weightmatrix.txt $weightfilename");
}

# meta parameters for training
print "creating $metafilename ...\n";
if(!$prokaryotic){
    if (stat "$configdir/species/generic/generic_metapars.cfg" == 0){
	die ("$configdir/species/generic/generic_metapars.cfg doesn't exist.");
    }
    system ("cp $configdir/species/generic/generic_metapars.cfg $metafilename");
    system ("cp $configdir/species/generic/generic_metapars.utr.cfg $metautrfilename");
    system ("cp $configdir/species/generic/generic_metapars.cgp.cfg $metacgpfilename");
}else{
    if (stat "$configdir/species/template_prokaryotic/template_prokaryotic_metapars.cfg" == 0){
        die ("$configdir/species/template_prokaryotic/template_prokaryotic_metapars.cfg doesn't exist.");
    }
    system ("cp $configdir/species/template_prokaryotic/template_prokaryotic_metapars.cfg $metafilename");
    system ("cp $configdir/species/template_prokaryotic/template_prokaryotic_metapars.utr.cfg $metautrfilename");
}

# files that are only created for prokaryotes
if($prokaryotic){
    print "creating $ovlpfilename ...\n";
    if (stat "$configdir/species/template_prokaryotic/template_prokaryotic_ovlp_len.pbl"== 0){
	die ("$configdir/species/template_prokaryotic/template_prokaryotic_ovlp_len.pbl does not exist.");
    }
    system ("cp $configdir/species/template_prokaryotic/template_prokaryotic_ovlp_len.pbl $ovlpfilename");
    print "creating $transshadowfilename ..\n";
    if (stat "$configdir/species/template_prokaryotic/template_prokaryotic_trans_shadow_bacterium.pbl"== 0){
	die ("$configdir/species/template_prokaryotic/template_prokaryotic_trans_shadow_bacterium.pbl does not exist.");
    }
    system ("cp $configdir/species/template_prokaryotic/template_prokaryotic_trans_shadow_bacterium.pbl $transshadowfilename");
}


print "The necessary files for training $species have been created.\n";
if (!$silent){
    print "Now, either run etraining or optimize_parameters.pl with --species=$species.\n";
    print "etraining quickly estimates the parameters from a file with training genes.\n";
    print "optimize_augustus.pl alternates running etraining and augustus to find optimal metaparameters.\n\n";
}
