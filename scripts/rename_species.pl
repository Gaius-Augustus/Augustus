#!/usr/bin/env perl

# Katharina J. Hoff
# last modified on September 25th 2022
# purpose of this script is to rename AUGUSTUS species parameters (e.g. from fly to Drosophila_melanogaster)

use File::Path qw(make_path);
use File::Copy;

my $usage = "renameParameters.pl oldname newname augustus_config_path\n";

# oldname could be fly
# newname could be Drosophila_melanogaster
# augustus_config_path could be e.g. ~/Augustus/config

if(@ARGV != 3){
    print $usage;
    exit;
}

# build original directory name
my $oldDir = $ARGV[2]."/species/".$ARGV[0];

# build new directory name
my $newDir = $ARGV[2]."/species/".$ARGV[1];


# find string that needs to be replaced, i.e. species name
opendir(my($dir), $oldDir) || die "Could not open directory $oldDir!\n";
my @dirs = readdir $dir;
closedir $dir;

# create the target directory
eval { make_path($newDir) };
if($@){
    print STDERR "Could not create directory $newDir: @$\n";
}


# retrieve files that need to be copied:
my %filehash;
opendir my($specdir), $oldDir or die "Could not open directory $oldDir!\n";
while( my $file = readdir($specdir)){
    if(not($file =~m/^\./)){
        $fstr = $file;
        $fstr =~ s/$ARGV[0]/$ARGV[1]/;
        if($file =~ m/_parameters.cfg/){
            open(INFILE, "<", "$oldDir/$file") or die ("Could not open file $oldDir/$file!\n");
            open(OUTFILE, ">", "$newDir/$fstr") or die ("Could not open output file $newDir/$fstr!\n");
            while(<INFILE>){
                $_=~s/(.+[Ff]ile\s+)$ARGV[0](_.+)/$1$ARGV[1]$2/;
                if($_=~m/^#/){
                    $_=~s/$ARGV[0]/$ARGV[1]/;
                }
                print OUTFILE $_;
            }
            close(INFILE);
            close(OUTFILE);
        }else{
            copy("$oldDir/$file","$newDir/$fstr") or die "Copying file $oldDir/$file to $newDir/$fstr failed: $!";
        }       
    }
}
closedir $specdir;

