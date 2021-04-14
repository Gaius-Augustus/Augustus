#!/usr/bin/env perl
#
# parse psl format, fiter the alignments and
# output them to stdout
#
# Katharina J. Hoff, April 14th 2021, katharina.hoff@uni-greifswald.de

use File::Basename;
use lib dirname (__FILE__);
use SplicedAlignment;
use strict;
use Getopt::Long;
use File::Which qw(which);

my $usage = "$0 -- parse psl format to orf\n";
$usage .= "Usage: $0 in.psl seqfile.fa\n";
$usage .= "output to stdout. first messages. then complete genes. then genes complete only at the 5' end\n\t--CDBTOOLS_PATH=/path/tocdbfasta/cdbyank can be set if not in \$PATH\n";
  
if ($#ARGV != 1) {
    die "Unknown option\n\n$usage";
}
my $pslFile = $ARGV[0];
my $seqFile = $ARGV[1];
my @alignmentList=();
my $CDBTOOLS_PATH;
my $cdbtools_path;

GetOptions(
        'CDBTOOLS_PATH=s' => \$cdbtools_path
);

set_CDBTOOLS_PATH();

###########################################################################################
#
# read in the psl alignments
#
###########################################################################################

open(PSL, "<$pslFile") or die ("Could not open $pslFile\n");

my ($estName, $estLength, $genomicName, $genomicLength, $estExonBegin, $estExonEnd, $exonBegin, $exonEnd, $intronFlag, $quality, $complement);
my @genomicMatches;

while (<PSL>) {
    my @pslLine = split(/\t/, $_);
    $estName = $pslLine[9];
    $estLength = $pslLine[12];
    $genomicName = $pslLine[13];
    $genomicLength = $pslLine[14]; # it is unclear whether the genome sequence size is meant or the range of the alignment on the genome sequence 
    $complement = ($pslLine[8] =~ m/-/);
    my $sa = new SplicedAlignment($estName, $estLength, $genomicName, $genomicLength, $complement);
    my $nBlocks = $pslLine[17];
    my @blockSizes = split(/,/, $pslLine[18]);
    my @qStarts = split(/,/, $pslLine[19]); # zero based, we want one based in the end, I think
    my @tStarts = split(/,/, $pslLine[20]); # zero based
    for(my $i = 0; $i < $nBlocks; $i++){
        $estExonBegin = $qStarts[$i] + 1;
        $estExonEnd = $estExonBegin + $blockSizes[$i] - 1;
        $exonBegin = $tStarts[$i] + 1;
        $exonEnd = $exonBegin + $blockSizes[$i] - 1;;
        $quality = 1; # don't know how to compute quality for single exons from psl
        if($i < ($nBlocks-1)){
            $intronFlag = 1;
        }else{
            $intronFlag = 0;
        }
        $sa->addExon($estExonBegin-1, $estExonEnd-1, $exonBegin-1, $exonEnd-1 , $quality, $intronFlag);
    }
	push @alignmentList, $sa;
    print("Pushed one line\n");
}

close(PSL) or die ("Could not close $pslFile!\n");

print "Read ".scalar(@alignmentList)." alignments\n";

my @geneList=();

###########################################################################################
#
# filter and enrich the data
#
###########################################################################################

foreach my $sa (@alignmentList){
    $sa->checkAlignment();
    print STDERR $sa->get_status(), "\n";
    if ($sa->get_status() eq "OK"){
	    my $seq = getSequence($sa->get_contigname);
	    $sa->makeGene($seq);
	$sa->findSplicedUTR();
	print STDERR "\t\t", $sa->get_status(), "\n"; 
    if ($sa->get_status() eq "OK"){
	        push @geneList, $sa;
	    }
    }
}

###########################################################################################
#
# output the data
#
###########################################################################################

print "### complete genes\n";
foreach my $sa (@geneList){
    if ($sa->get_complete5prime() && $sa->get_complete3prime()){
	print $sa->output;
    }
}
print "### genes incomplete at the 3' end and complete at the 5' end\n";
foreach my $sa (@geneList){
    if ($sa->get_complete5prime() && !$sa->get_complete3prime()){
	print $sa->output;
    }
}



###########################################################################################
#
# subroutines
#
###########################################################################################

sub getSequence {
    my $seqname = shift;
    my $seq;
    if (not(-e $seqFile.".idx")){
        # create fasta idx file
        print("$CDBTOOLS_PATH/cdbfasta $seqFile\n");
        system("$CDBTOOLS_PATH/cdbfasta $seqFile");
    }
    # let cdbyank use the index to search the sequence
    system("echo '$seqname' | $CDBTOOLS_PATH/cdbyank ${seqFile}.cidx -d $seqFile > genomic.fa");
    system ("perl -i.orig -p -e 's/^Incorrectly.*fixed.\n//' genomic.fa");
    open (SEQ, "<genomic.fa") or die ("Could not open genomic.fa");
    $/="\n";
    my @lines=<SEQ>;
    shift @lines;
    $seq = join ("", @lines);
    $seq =~ s/\n//g;
    return $seq;
}

####################### set_CDBTOOLS_PATH #######################################
# * set path to cdbfasta/cdbyank
################################################################################

sub set_CDBTOOLS_PATH {
    # try to get path from ENV
    if ( defined( $ENV{'CDBTOOLS_PATH'} ) && not (defined($cdbtools_path)) ) {
        if ( -e $ENV{'CDBTOOLS_PATH'} ) {
            print  "\# " . (localtime) . ": Found environment variable \$CDBTOOLS_PATH. Setting "
                  . "\$CDBTOOLS_PATH to ".$ENV{'CDBTOOLS_PATH'}."\n";
            $CDBTOOLS_PATH = $ENV{'CDBTOOLS_PATH'};
        }
    }elsif(not(defined($cdbtools_path))) {
        print "\# " . (localtime) . ": Did not find environment variable \$CDBTOOLS_PATH\n";
    }

    # try to get path from command line
    if ( defined($cdbtools_path) ) {
        my $last_char = substr( $cdbtools_path, -1 );
        if ( $last_char eq "\/" ) {
            chop($cdbtools_path);
        }
        if ( -d $cdbtools_path ) {
            print "\# " . (localtime)
                . ": Setting \$CDBTOOLS_PATH to command line argument "
                . "--CDBTOOLS_PATH value $cdbtools_path.\n";
            $CDBTOOLS_PATH = $cdbtools_path;
        } else {
            print "#*********\n"
                . "# WARNING: Command line argument --CDBTOOLS_PATH was "
                . "supplied but value $cdbtools_path is not a directory. Will not "
                . "set \$CDBTOOLS_PATH to $cdbtools_path!\n"
                . "#*********\n";
        }
    }

    # try to guess
    if ( not( defined($CDBTOOLS_PATH) )
        || length($CDBTOOLS_PATH) == 0 )
    {
        print "\# " . (localtime)
            . ": Trying to guess \$CDBTOOLS_PATH from location of cdbfasta"
            . " executable that is available in your \$PATH\n";
        my $epath = which 'cdbfasta';
        if(defined($epath)){
            if ( -d dirname($epath) ) {
                print "\# " . (localtime) . ": Setting \$CDBTOOLS_PATH to "
                    . dirname($epath) . "\n";
                $CDBTOOLS_PATH = dirname($epath);
            }
        } else {
            print "#*********\n"
                . "# WARNING: Guessing the location of \$CDBTOOLS_PATH "
                . "failed with \"which cdbfasta\"!\n"
                . "#*********\n";
        }
    }

    if ( not( defined($CDBTOOLS_PATH) ) ) {
        print "cdbfasta and cdbyank are required for running psl2orf.pl\n"
           .  "You have 3 different "
           .  "options to provide a path to cdbfasta/cdbyank:\n"
                    .  "   a) provide command-line argument\n"
                    .  "      --CDBTOOLS_PATH=/your/path\n"
                    .  "   b) use an existing environment variable\n"
                    . "       \$CDBTOOLS_PATH\n"
                    .  "      for setting the environment variable, run\n"
                    .  "           export CDBTOOLS_PATH=/your/path\n"
                    .  "      in your shell. You may append this to your "
                    .  ".bashrc or .profile file in\n"
                    .  "      order to make the variable available to all your\n"
                    .  "      bash sessions.\n"
                    .  "   c) braker.pl can try guessing the "
                    .  "      \$CDBTOOLS_PATH from the location of a cdbfasta\n"
                    .  "      executable that is available in your \$PATH\n"
                    .  "      variable. If you try to rely on this option, you\n"
                    . "       can check by typing\n"
                    .  "           which cdbfasta\n"
                    .  "      in your shell, whether there is a cdbfasta\n"
                    .  "      executable in your \$PATH\n";
        print STDERR "\# " . (localtime) . " ERROR: in file " . __FILE__
                   . " at line ". __LINE__ . "\n" . "\$CDBTOOLS_PATH not set!\n";
        exit(1);
    }
    if ( not ( -x "$CDBTOOLS_PATH/cdbfasta" ) ) {
        print STDERR "\# " . (localtime) . " ERROR: in file " . __FILE__
             ." at line ". __LINE__ ."\n"
             . "$CDBTOOLS_PATH/cdbfasta is not an executable file!\n";
        exit(1);
    }elsif ( not ( -x "$CDBTOOLS_PATH/cdbyank" ) ) {
        print STDERR "\# " . (localtime) . " ERROR: in file " . __FILE__
                ." at line ". __LINE__ ."\n"
                . "$CDBTOOLS_PATH/cdbyank is not an executable file!\n";
        exit(1);
    }
}