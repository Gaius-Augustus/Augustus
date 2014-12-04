#!/usr/bin/perl

# Katharina J. Hoff
# November 14th 2014
# This script is part of WebAUGUSTUS.
# The purpose is to copy user-uploaded parameters to the AUGUSTUS_CONFIG_PATH and rename the parameters according to the job id

use File::Path qw(make_path);

my $usage = "moveParameters.pl parameterdir newid copytodir\n";

# parameterdir is the directory in which the parameter file originally reside, without / at the end
# newid is the new species identifier of the parameter set, e.g. the WebAUGUSTUS job id
# copytodir is the folder in which AUGUSTUS looks for parameter files, without / at the end

if(@ARGV != 3){
    print $usage;
    exit;
}

# find string that needs to be replaced, i.e. species name
opendir my($dir), $ARGV[0] or die "Could not open directory $ARGV[0]!\n";
my @dirs = readdir $dir;
closedir $dir;
# params in WebAUGUSTUS contains three folders: ., .., and the parameter archive. Determine which is the archive:
my $archC = 0;
foreach(@dirs){
    if(not($_ eq ".") and not($_ eq "..")){
	$species = $_;
    }
    $archC++;
}
if($archC > 3){
    print STDERR "The parameter archive appears to have contained more than one folder, not sure which one to use!\n";
}

print "Species is $species\n";

my $targetdir = $ARGV[2]."/".$ARGV[1];
# create the target directory
eval { make_path($targetdir) };
if($@){
    print STDERR "Could note create directory $targetdir: @$\n";
}

# retrieve files that need to be copied:
my $speciesdir = "$ARGV[0]/$species";
my %filehash;
opendir my($specdir), $speciesdir or die "Could not open directory $ARGV[0]!\n";
while( my $file = readdir($specdir)){
    if(not($file =~m/^\./)){
	open(INFILE, "<", "$speciesdir/$file") or die ("Could not open file $speciesdir/$file!\n");
	$fstr = $file;
	$fstr =~ s/$species/$ARGV[1]/;
	open(OUTFILE, ">", "$targetdir/$fstr") or die ("Could not open output file $targetdir/$fstr!\n");
	while(<INFILE>){
	    $_=~s/$species/$ARGV[1]/;
	    print OUTFILE $_;
	}
	close(INFILE);
	close(OUTFILE);
    }

}
closedir $specdir;


