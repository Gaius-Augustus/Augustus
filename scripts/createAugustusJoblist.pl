#!/usr/bin/perl
#
# createAugustusJoblist.pl
# Create a joblist with overlapping sequence chunks from multiple fasta files.
# Mario Stanke, 22.01.2007
#
#

use strict;
use Getopt::Long;

my $usage = "$0 -- Create a joblist with sequence chunks from multiple fasta files.\n\n";
$usage .= "Usage: $0\n\n";
$usage .= "parameters:\n";
$usage .= "--sequences seqs.lst input sequences, format: each line contains one sequence including the full path and its size, e.g.\n";
$usage .= "                     /cluster/data/panTro2/1/chr1.fa\t1\t229974691\n";
$usage .= "                     /cluster/data/panTro2/1/chr1_random\t1\t9420409\n";
$usage .= "                     /cluster/data/panTro2/2/chr2a\t1\t114460064\n";
$usage .= "                     or\n";
$usage .= "                     /cluster/data/panTro2/1/chr1_random\t/hints/chr1_random\t1\t9420409\n";
$usage .= "                     /cluster/data/panTro2/2/chr2a\t/hints/chr2a\t1\t114460064\n";
$usage .= "--outputdir s        directory, in which later the AUGUSTUS output will be written.\n";
$usage .= "--command s          AUGUSTUS command, e.g. \"augustus --species=human --maxDNAPieceSize=600000\".\n";
$usage .= "--joblist job.lst    filename with list of jobs as given to parasol.\n";
$usage .= "--chunksize n        chunk size. Each sequence is (imaginarily) cut into chunks of this size\n\n";
$usage .= "options:\n";
$usage .= "--overlap n          overlap. Neighboring chunks overlap by this number of bases.\n";
$usage .= "--padding n          padding on both sides (default 0).\n";
$usage .= "--errordir errdir    directory, in which later the AUGUSTUS error messages will be written.\n";
$usage .= "--check              insert parasol input/output checks.\n";
$usage .= "--wrap=s             have each job in a separate file, preceded by command s.\n";
$usage .= "--jobprefix=s        prefix of job name (default: \"job.\")\n";

if (@ARGV < 2) {
    print $usage;
    exit;
}

my ($seqfilename, $chunksize, $overlap, $outputdir, $errordir, $command, $joblist, $seqnr, $name, $predStart, $predEnd, $pE, $chunkid, $padding, $check, $wrap, $wholecommand);
my $jobnr = 0;
$padding = 0;
$check = 0;
$wrap = "";
my $jobprefix = "job.";

GetOptions('sequences=s'=>\$seqfilename,
	   'chunksize=i'=>\$chunksize,
	   'overlap:i'=>\$overlap,
	   'outputdir=s'=>\$outputdir,
	   'errordir=s'=>\$errordir,
	   'command=s'=>\$command,
	   'padding:i'=>\$padding,
	   'joblist=s'=>\$joblist,
	   'check!'=>\$check,
	   'wrap=s'=>\$wrap,
	   'jobprefix:s'=>\$jobprefix);

die ("Need to specify chunksize.\n") unless (defined $chunksize);
die ("Need to specify command.\n") unless (defined $command);
die ("Need to specify joblist.\n") unless (defined $joblist);
die ("Need to specify chunksize.\n") unless (defined $chunksize);
die ("Need to specify outputdir.\n") unless (defined $outputdir);

open(SEQ, "<$seqfilename") or die ("Could not open $seqfilename");
open (BATCH, ">$joblist") or die "Couldn't open $joblist\n";

$outputdir =~ s/\/$//;
$seqnr = 1;
while (<SEQ>) {
    my ($path, $start, $end, $hints);
    if (/(.+)\s+(.+)\s+(.+)\s+(.+)/) {
	$path = $1;
	$hints = $2;
	$start = $3;
	$end = $4;
    } elsif (/(.+)\s+(.+)\s+(.+)/){
	$path = $1;
	$start = $2;
	$end = $3;
	undef $hints;
    }
    $start -= $padding;
    $end += $padding;

    my $ovlp;
    if (!$overlap) {
	$ovlp = (($chunksize > 3000000)? 500000 : $$chunksize/6);
    } else {
	$ovlp = $overlap;
    }
    my $chunknr = 0;
    $name = $path;
    $name =~ s/.*\///;
    $predEnd = -1;
    for ($predStart = $start; $predStart <= $end && $predEnd < $end; $predStart = $predEnd+1-$ovlp) {
	$chunknr++;
	$predEnd = $predStart + $chunksize-1;
	$pE = $predEnd;
	if ($pE > $end){
	    $pE = $end;
	}
	$chunkid = sprintf ("%03d", $chunknr);
	my $gfffilename = "$outputdir/$seqnr.$chunkid.${name}.$predStart..$pE.gff";
	my $errorfilename;
	if (defined $errordir) {
	    $errorfilename = "$errordir/$seqnr.$chunkid.${name}.$predStart..$predEnd.err";
	} else {
	    $errorfilename = "$outputdir/$seqnr.$chunkid.${name}.$predStart..$pE.err";
	}
	$wholecommand = $command;
	if (defined $hints) {
	    if ($check){
		$wholecommand .= " --hintsfile={check in exists $hints}";
	    } else {
		$wholecommand .= " --hintsfile=$hints";
	    }
	}
	if ($check) {
	    $wholecommand .=  " --predictionStart=$predStart --predictionEnd=$pE {check in line+ $path} --outfile={check out line+ $gfffilename} --errfile=$errorfilename\n";
	} else {
	    $wholecommand .=  " --predictionStart=$predStart --predictionEnd=$pE $path --outfile=$gfffilename --errfile=$errorfilename\n";
	}
	if ($wrap eq "") {
	    print BATCH $wholecommand;
	} else {
	    $jobnr++;
	    my $jobname = "$jobprefix$jobnr";
	    open (JOB, ">$jobname");
	    print JOB "$wrap\n" . $wholecommand;
	    print BATCH "$jobname\n";
	    close JOB;
	    system("chmod +x $jobname");
	}
    }
    $seqnr++;
}
