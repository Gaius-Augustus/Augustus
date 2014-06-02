#!/usr/bin/perl
#
# Convert a psl format output file (from BLAT or GMAP) from cDNA alignments
# to a hints file
#
# Mario Stanke, 26.6.2008

use strict;
use Getopt::Long;

my $usage = "$0 -- convert blat file with mRNA or EST alignments to hints file for AUGUSTUS\n";
$usage .= "\n";
$usage .= "Usage: $0 --in=blatfile --out=hintsfile\n";
$usage .= "  PREREQUISITE: input psl file must be sorted by target (=genome) sequence names\n";
$usage .= "  and should be sorted within the sequences by begin coordinates for efficiency\n";
$usage .= "  e.g. do\n";
$usage .= "  cat blat.psl | sort -n -k 16,16 | sort -s -k 14,14 > blat.sorted.psl\n";
$usage .= "  for output from the UCSC table browser do:\n";
$usage .= "  cat blat.psl | sort -n -k 17,17 | sort -s -k 15,15 > blat.sorted.psl\n";
$usage .= "  when the 17th field is the begin coordinate and the 15th field is the sequence name\n\n";
$usage .= "  options:\n";
$usage .= "  --priority=n       priority of hint group (default 4)\n";
$usage .= "  --maxgaplen=n      gaps at most this length are simply closed (default 14)\n";
$usage .= "  --minintronlen=n   alignments with gaps shorter than this and longer than maxgaplen are discarded (default 41)\n";
$usage .= "  --maxintronlen=n   alignments with longer gaps are discarded (default 350000)\n";
$usage .= "  --maxQgaplen=n     maximum length of gap in query (cDNA) sequence (default 5)\n";
$usage .= "  --ep_cutoff=n      this many bp are cut off of each exonpart hint at end of alignment (default 10)\n";
$usage .= "  --source=s         source identifier (default 'E')\n";
$usage .= "  --intronsonly      only retrieve intron hints (e.g. because the exon(part) hints are retreived by converting to a wig track, default: off)\n";
$usage .= "  --nomult           do not summarize multiple identical intron hints to a single one\n";
$usage .= "  --remove_redundant only keep the strongest hint for a region (default false)\n";
$usage .= "  --maxcoverage=n    maximal number of hints at a given position. A high value causes long running time of\n";
$usage .= "                     AUGUSTUS in regions with thousands of cDNA alignments. (default 3000)\n";
$usage .= "  --ssOn             include splice site (dss, ass) hints in output (default false)\n";
$usage .= "  --trunkSS          include splice sites hints from the ends of a truncated alignment (contig too short)\n";
$usage .= "  --coloffset=n      column offset, 0 for direct blat output, 1 for psl format from UCSC database (default 0)\n";
$usage .= "  --score=f          fill this number in in the score column (default 0)\n";
$usage .= "  --clonefile=s      provide a file with clone names so close alignments from the same clone can be grouped.\n";
$usage .= "                     AUGUSTUS will try to put those hints into a single transcripts even if different ends of\n";
$usage .= "                     the clones were sequenced. File format (tab delimited):\n";
$usage .= "                     cloneA\tread1\tread2\n";
$usage .= "                     cloneA\tread3\n";
$usage .= "                     cloneB\tread4\tread5\n";
$usage .= "  --terminusfile=s   provide a file with EST terminus information to infer tss/tts hints.\n";
$usage .= "                     AUGUSTUS will use tss/tts hints to predict transcription start/termination sites\n";
$usage .= "                     File format (tab delimited):\n";
$usage .= "                     # ESTname    EstDir    Type FrontTerminus  EndTerminus\n";
$usage .= "                     CACW5781.b1     5       A2      5TSS       Unknown\n";
$usage .= "                     CACW6759.g1     3       F23     5TNS       3TNS\n";
$usage .= "                     CACW14459.g2    3       D2      Unknown    3TNS\n";
$usage .= "                     CACW21662.g1    3       C2      5TNS       Unknown\n";
$usage .= "                     CACW25491.g1    3       F21     5TNS       3TNS-NP\n";
$usage .= "                     \n";
$usage .= "                     cloneB\tread4\tread5\n";
$usage .= "  --maxgenelen=n     alignments of the same clone are considered to be of the same gene if not separeted by more than this (400000)\n";
$usage .= "                     Alignments that span more than this are ignored, but better filter long introns through alignment program.\n";

my $blatfilename;
my $hintsfilename;
my $minintronlen = 41;
my $maxintronlen = 350000;
my $maxgaplen = 14;
my $maxQgaplen = 5;
my $ep_cutoff = 10;
my $min_endblock_len = 8;
my $termhintradius = 15;
my $source="E";
my $priority = 4;
my $prgsrc = "b2h";
my $line=0;
my $coloffset=0;
my $mult=1;
my $remove_redundant=0;
my @coverage=();
my $maxcoverage = 3000;
my $ssOn=0;
my $trunkSS=0;
my $exonpartid = "ep"; # abbreviate to decrease file size
my $minESToverhang = 40; # for trunkSS
my $maxDNAoverhang = 10000; # for trunkSS
my $score = 0;
my ($clonefile, $terminusfile);
my $maxgenelen = 400000;
my $intronsonly=0;
my %seqseen = (); # for testing whether it is sorted by target sequences
my $warnslow = 0;

if ($#ARGV < 1) {
    print "$usage";
    exit;
}
GetOptions(
	   'in=s'=>\$blatfilename,
	   'out=s'=>\$hintsfilename,
	   'minintronlen:i'=>\$minintronlen,
	   'maxintronlen:i'=>\$maxintronlen,
	   'maxgaplen:i'=>\$maxgaplen,
	   'maxQgaplen:i'=>\$maxQgaplen,
	   'ep_cutoff:i'=>\$ep_cutoff,
	   'source:s'=>\$source,
	   'priority:i'=>\$priority,
	   'coloffset:i'=>\$coloffset,
	   'intronsonly!'=>\$intronsonly,
	   'mult!'=>\$mult,
	   'remove_redundant!'=>\$remove_redundant,
	   'maxcoverage:i'=>\$maxcoverage,
	   'ssOn!'=>\$ssOn,
	   'trunkSS!'=>\$trunkSS,
           'score:f'=>\$score,
	   'clonefile:s'=>\$clonefile,
	   'terminusfile:s'=>\$terminusfile,
	   'maxgenelen:i'=>\$maxgenelen);

if (!defined($hintsfilename)) {
    print "Missing output file name.\n$usage";
    exit;
}

my %grp = (); # keys: "tname-qname-tstart-tend", values: grpname
my ($match,$TgapCount,$strand,$qname,$qsize,$blockSizes,$tStarts, $qStarts, $tstart, $tend, $qstart, $qend, $tsize, $targetname);
my %clone = (); # keys: seq names, values: clone names
my %grpnr = ();
my $skiplines=0;

if ($clonefile){
    # have a file with clone names that tell us that some cDNAs belong to the same transcript
    # use this to group alignments of cDNAs of the same clone
    open(CLONE, "<$clonefile") || die "Could not open $clonefile\n";

    while (<CLONE>){
	# store a hash that gives the clone name for each seqname
	chomp;
	my @f = split /\t/;
	print STDERR "Wrong format in clone file. Line was:\n$_" if (@f<2);
	for(my $i=1; $i<@f; $i++){
	    $clone{$f[$i]}=$f[0];
	}
    }
    close CLONE;

    my @hits = (); # a hit consists of: contig name, queryname, tbegin, tend
    open(BLAT, "<$blatfilename") || die "Couldn't open $blatfilename\n";

    while(<BLAT>){
	$skiplines=5 if (/psLayout/);
	if ($skiplines>0) {
	    $skiplines--;
	    next;
	}
	s/#.*//;
	next unless /\S/;
	my @f = split /\t/, $_, $coloffset+21;
	if (@f < $coloffset+20) { warn "Not BLAT format"; next } # blat format from the GenomeBrowser has an additional first column    
	$qname       = $f[$coloffset+9];
	$targetname  = $f[$coloffset+13];
	$tstart      = $f[$coloffset+15];
	$tend        = $f[$coloffset+16]-1;
	
	push @hits, [$targetname,$qname,$tstart,$tend];
    }
    close BLAT;
    @hits = sort {$a->[0] cmp $b->[0] or $clone{$a->[1]} cmp $clone{$b->[1]} or $a->[2] cmp $b->[2] or $a->[3] cmp $b->[3] or $a <=> $b} @hits; # sort by contig name, clone name and then by position


    my @samegrp = ();
    foreach my $hit (@hits){
	if (@samegrp<1 || !exists($clone{$hit->[1]}) || $clone{$samegrp[0]->[1]} ne $clone{$hit->[1]} || $samegrp[0]->[0] ne $hit->[0] || $hit->[3]-$samegrp[0]->[2]+1 > $maxgenelen){
	    # belongs to a different group, save previous group
	    addgrp(\@samegrp);
	    @samegrp=();
	} else {
	    
	}
	push @samegrp, $hit;
    }
    # add very last group
    addgrp(\@samegrp);
}# if $clonefile

# addgrp
# A group of two or more neighboring alignments of the same clone is stored for later referene.
sub addgrp {
    my $samegrpref = shift;
    if (@{$samegrpref} > 1 && exists($clone{$samegrpref->[0]->[1]})){ # if two or more alignments are neighbors and from the same clone, then give them a group name
	my $cname = $clone{$samegrpref->[0]->[1]};
	# check whether the same cDNA has been aligned twice in this group
	# if yes, don't consider this a group, because this happens when a genomic region is repeated and then there should be rather two than one gene
	my $duplicate=0;
	my %cdnanames=();
	foreach my $h (@{$samegrpref}){
	    if (exists($cdnanames{$h->[1]})) {
		$duplicate=1;
	    }
	    $cdnanames{$h->[1]} = 1;
	}
	if (!$duplicate) {
	    my $nr = 1;
	    if (exists($grpnr{$cname})){
		$nr = ++$grpnr{$cname};
	    } else {
		$nr = $grpnr{$cname} = 1;
	    }
	    foreach my $h (@{$samegrpref}){
		$grp{"$h->[0]-$h->[1]-$h->[2]-$h->[3]"} = $clone{$h->[1]} . "-" . $nr; #($nr>1)? $clone{$h->[1]} : $clone{$h->[1]} . "-" . $nr;
	    }
	}
    }
}

my %termini = (); # keys: qname, values: [FrontTermius, EndTerminus, 5Termpos, 3Termpos]
                  # FrontTerminus in {5TSS, 5TNS}, FrontTerminus in {3TSS, 3TNS}, 5/3Termpos is the position of the first/last transcript position in the cDNA, right after/before the adapter/vector sequence

if ($terminusfile){
    # have a file with terminus information that tells us about tss/tss sites and directions of ESTs
    open(TERM, "<$terminusfile") || die "Could not open $terminusfile\n";

    while (<TERM>){
	chomp;
	next if (/^\s*\#/);
	my @f = split /\t/;
	if (@f==5) { # old file format with 5 columns
	    if (exists($termini{$f[0]})){
		print STDERR "Duplicate entry in terminus file for $f[0]\n";
	    } else {
		$termini{$f[0]} = [$f[3], $f[4]]; # assume in this case that termpos is right at end of sequence
	    }
	} elsif (@f==6) { # new file format with 6 columns
	    if (!exists($termini{$f[0]})){
		$termini{$f[0]} = ["", "", -1, -1];
	    }
	    if ($f[1] eq "5TSS" || $f[1] eq "5TNS"){
		if ($termini{$f[0]}->[0] eq ""){
		    $termini{$f[0]}->[0] = $f[1];
		    $termini{$f[0]}->[2] = $f[4]+1;
		} else {
		    #   print STDERR "Duplicate 5' terminus for same cDNA $f[0].\n";
		}
	    } elsif ($f[1] eq "3TSS" || $f[1] eq "3TNS"){
		if ($termini{$f[0]}->[1] eq ""){
		    $termini{$f[0]}->[1] = $f[1];
		    $termini{$f[0]}->[3] = $f[3]-1;
		} else {
		    #print STDERR "Duplicate 3' terminus for same cDNA $f[0].\n";
		}
	    } else {
		print STDERR "Wrong format in second field in terminus file. Line was:\n$_";	
	    }
	} else {
	    print STDERR "Wrong format in terminus file. Line was:\n$_";
	}
    }
    close TERM;
}



open(BLAT, "<$blatfilename") || die "Couldn't open $blatfilename\n";
open(HINTS, ">$hintsfilename") || die "Could not open $hintsfilename";

my ($i, $j, $mstart, $mend, $badalignment, $gaplen);
my (@dsshints, @asshints, @exonhints, @exonparthints, @intronhints, @tsshints, @ttshints);
my (@f,@b,@t,@q);
my (@blockbegins, @blockends, @folintronok);
my $numBlocks;

# hint lists are sorted by by increasing begin position
my @hint; # (begin, end, strand, tname, qname)
my $hintref;
my $oldtargetname;
$oldtargetname = "no name yet";
$skiplines=0;
my %qnamefreqs = ();

while (<BLAT>) {
    if (/psLayout/){
	$skiplines=5;
    }
    if ($skiplines>0) {
        $skiplines--;
        next;
    }
 
    $line++;
    s/^#.*//;
    next unless /\S/;
    if ($line%1000==1){
	$| = 1;
	print "\r"."processed line $line";
    }
    
    @f = split /\t/, $_, $coloffset+21;
    if (@f < $coloffset+20) { warn "Not BLAT format"; next } # blat format from the GenomeBrowser has an additional first column
    
    $match       = $f[$coloffset+0];
    $TgapCount   = $f[$coloffset+6];
    $strand      = $f[$coloffset+8];
    $qname       = $f[$coloffset+9];
    $qsize       = $f[$coloffset+10];
    $qstart      = $f[$coloffset+11];
    $qend        = $f[$coloffset+12];
    $targetname  = $f[$coloffset+13];
    $tsize       = $f[$coloffset+14];
    $tstart      = $f[$coloffset+15];
    $tend        = $f[$coloffset+16]-1; # BLAT reports one base after the end
    $blockSizes  = $f[$coloffset+18];
    $qStarts     = $f[$coloffset+19];
    $tStarts     = $f[$coloffset+20];
    
    #print "match=", $match;
    #print " TgapCount=", $TgapCount;
    #print " strand=", $strand;
    #print " qname=", $qname;
    #print " blockSizes=", $blockSizes;
    #print " tStarts=", $tStarts, "\n";
 
    next if ($tend - $tstart + 1 > $maxgenelen);
 
    my $qidx = ++$qnamefreqs{$qname};
    my $grpname = $qname;
    my $cdnaname = "";
    if ($clonefile){
	my $key = "$targetname-$qname-$tstart-$tend";
	if (exists($grp{$key})){
	    $grpname = $grp{$key}; # this cdna and other cdnas of the same clone align nearby, use one group name for all hints
	    $cdnaname= $qname;
	} elsif ($qidx > 1){
	    $grpname = "$qname:$qidx"; # same cdna aligned to several places, distinguish different groups, one per alignment
	    $cdnaname= $qname;
	}
    }
 
    if ($targetname ne $oldtargetname) {
	if ($seqseen{$oldtargetname}){
	    if ($intronsonly && $mult){
		die ("\nPSL file MUST be sorted by target sequence names when 'intronsonly' and 'mult' options are acive (do: cat my.psl | sort -n -k 16,16 | sort -s -k 14,14).\n");
	    } elsif (!$warnslow) {
		print STDERR "WARNING: PSL file not sorted by target sequence names. Will be slow. Do cat my.psl | sort -n -k 16,16 | sort -s -k 14,14 for better performance.\n";
		$warnslow = 1;
	    }
	}
	printHints();
	$seqseen{$oldtargetname} = 1;
	$#coverage = -1;
    }
    $oldtargetname = $targetname;

    # filter hits
    # 
    #
    #if ($TgapCount < 0) {next}
    my $filterout=0;
    for (my $i = int($tstart/10); $i <= int($tend/10) && !$filterout;$i++) {
	if ($coverage[$i] >= $maxcoverage) {
	    #print "omitting " , $qname, " too much coverage at ", 10*$i , "\n";
	    $filterout=1;
	}
    }
    if ($filterout) {
	next;
    }
    for (my $i=int($tstart/10); $i <= int($tend/10); $i++) {
	if (defined $coverage[$i]) {
	    $coverage[$i]++;
	} else {
	    $coverage[$i]=1;
	}
    }

    $blockSizes =~ s/[, ]$//;
    $tStarts =~ s/[, ]$//;
    @b = split (/,/, $blockSizes);
    @t = split (/,/, $tStarts);
    @q = split (/,/, $qStarts);
    
    $numBlocks = scalar @t;
    # Go throught the line and decide for each gap in the target sequence whether
    # it is an intron or a deletion. If neither seems likely, ignore that line.
    # Close the deletion gaps by joining the adjacent blocks.
    #
    $badalignment = 0;
    @blockbegins=();
    @blockends=();
    @folintronok=();
    for ($i=0; $i<$numBlocks; $i++) {
	$mstart = $t[$i]+1; # blat is 0-based
	$mend = $mstart + $b[$i] - 1;
	if ($#blockends>=0) {
	    $gaplen = $mstart - $blockends[$#blockends] - 1;
	} else {
	    $gaplen = $minintronlen; # so the block is added below
	}
	#print " gaplen ", $gaplen; 
	if ($gaplen >= $minintronlen && $gaplen <= $maxintronlen) {
	    push @blockbegins, $mstart;
	    push @blockends, $mend;
	    if ($i+1 < $numBlocks && ($q[$i] + $b[$i] >= $q[$i+1] - $maxQgaplen)) {
		push @folintronok , 1;
	    } else {
		push @folintronok , 0;
	    }
#	    print "intron qgap from " . ($q[$i] + $b[$i]) . " to " . ($q[$i+1]) . " len= "  . ($q[$i+1] - $q[$i] - $b[$i]) . " ok=" . $folintronok[$#folintronok] . "\n" if ($i+1<$numBlocks);
	} elsif ($gaplen <= $maxgaplen){
	    # close the gap and igore it
	    $blockends[$#blockends] = $mend;
	    if ($i+1 < $numBlocks && ($q[$i] + $b[$i] >= $q[$i+1] - $maxQgaplen)) {
		$folintronok[$#folintronok] = 1;
	    } else {
		$folintronok[$#folintronok] = 0;
	    }
#	    print "new qgap from " . ($q[$i] + $b[$i]) . " to " . ($q[$i+1]) . " len= "  . ($q[$i+1] - $q[$i] - $b[$i]) . " ok=" . $folintronok[$#folintronok] . "\n" if ($i+1<$numBlocks);

	} else {
	    $badalignment = 1;
	} 
    }
    next unless ($badalignment == 0);
#    print "\n Gesamt:\n";
#    for ($i=0; $i<=$#blockbegins; $i++) {
# 	print $blockbegins[$i], "-", $blockends[$i], ",ok=", $folintronok[$i], "\n";
#    }
    
    # now add the hints
    $numBlocks = scalar @blockbegins;
    for ($i=0; $i<$numBlocks; $i++) {
	if ($i==0 && $i==$numBlocks-1 && !$intronsonly) {
	    # just one exonpart, should not happen when spliced EST
	    if ($blockbegins[$i] + 2*$ep_cutoff <= $blockends[$i]){
		@hint = ($blockbegins[$i]+$ep_cutoff, $blockends[$i]-$ep_cutoff, '.', $grpname, $cdnaname, 1);
		addExonpartHint([@hint]);
	    }
	} elsif ($i==0) {
	    # first block
	    if ($blockbegins[$i] + $min_endblock_len-1 <= $blockends[$i]){
		if ($blockbegins[$i] + $ep_cutoff <= $blockends[$i] && !$intronsonly){
		    @hint = ($blockbegins[$i]+$ep_cutoff, $blockends[$i], '.', $grpname, $cdnaname, 1);
		    addExonpartHint([@hint], $grpname);
		}
		if ($ssOn && !$intronsonly) {
		    # add both a dss hint on the plus strand and a ass hint on the reverse strand as 
		    # blat does not label the strand reliably
		    @hint = ($blockends[$i]+1, $blockends[$i]+1, '.', $grpname, $cdnaname, 1);
		    addSignalHint(\@dsshints, [@hint]);
		    @hint = ($blockends[$i]+1, $blockends[$i]+1, '.', $grpname, $cdnaname, 1);
		    addSignalHint(\@asshints, [@hint]);
		    @hint = ($blockbegins[$i+1]-1, $blockbegins[$i+1]-1, '.', $grpname, $cdnaname, 1);
		    addSignalHint(\@dsshints, [@hint]);
		    @hint = ($blockbegins[$i+1]-1, $blockbegins[$i+1]-1, '.', $grpname, $cdnaname, 1);
		    addSignalHint(\@asshints, [@hint]);
		}
		if ($folintronok[$i] && ($i<$numBlocks-2 || $blockends[$i+1]-$blockbegins[$i+1]+1 > $min_endblock_len)) {
		    # add following intron hint
		    @hint = ($blockends[$i]+1, $blockbegins[$i+1]-1, '.', $grpname, $cdnaname, 1);
		    addIntervalHint(\@intronhints, [@hint]);
		}
	    }
	} elsif ($i==$numBlocks-1 && !$intronsonly) {
	    # last block
	    if ($blockends[$i] - $min_endblock_len + 1 >= $blockbegins[$i]){
		if ($blockbegins[$i] <= $blockends[$i]-$ep_cutoff){
		    @hint = ($blockbegins[$i], $blockends[$i]-$ep_cutoff, '.', $grpname, $cdnaname, 1);
		    addExonpartHint([@hint]);
		}
	    }
	} else { 
	    # internal block, add following intron hint 
	    if (!$intronsonly){
		@hint = ($blockbegins[$i], $blockends[$i], '.', $grpname, $cdnaname, 1);
		insertIntervalHint(\@exonhints, [@hint]);
	    }
	    if ($folintronok[$i] && ($i<$numBlocks-2 || $blockends[$i+1]-$blockbegins[$i+1]+1 > $min_endblock_len)) {
		@hint = ($blockends[$i]+1, $blockbegins[$i+1]-1, '.', $grpname, $cdnaname, 1);
		addIntervalHint(\@intronhints, [@hint]);
		if ($ssOn && !$intronsonly){
		    @hint = ($blockends[$i]+1, $blockends[$i]+1, '.', $grpname, $cdnaname, 1);
		    addSignalHint(\@dsshints, [@hint]);
		    @hint = ($blockends[$i]+1, $blockends[$i]+1, '.', $grpname, $cdnaname, 1);
		    addSignalHint(\@asshints, [@hint]);
		    @hint = ($blockbegins[$i+1]-1, $blockbegins[$i+1]-1, '.', $grpname, $cdnaname, 1);
		    addSignalHint(\@asshints, [@hint]);
		    @hint = ($blockbegins[$i+1]-1, $blockbegins[$i+1]-1, '.', $grpname, $cdnaname, 1);
		    addSignalHint(\@dsshints, [@hint]);
		}
	    }
	}
    }
    # transcription start and termination hints (tss,tts)
    # 
    #
    if ($terminusfile){
	my $tref = $termini{$qname};
	if ($tref) {
	    my ($FrontTerminus, $EndTerminus) = ($tref->[0], $tref->[1]);
	    my ($vector5len, $vector3len) = (0,0); # length of vector sequence before real transcript starts
	    if (@{$tref}==4){
		$vector5len = $tref->[2]-1 if ($tref->[2] != -1);
		$vector3len = $qsize - $tref->[3] if ($tref->[3] != -1);
	    }
	    # print "$FrontTerminus, $EndTerminus, $vector5len, $vector3len, " . ($qstart - $vector5len) . " " . ($qsize-$qend - $vector3len)  . "\n";
	    my ($termcenter, $termstart, $termend, $termstrand);
	    if (($FrontTerminus eq "5TSS" || $FrontTerminus eq "5TNS") && $qstart - $vector5len < 3) { # terminus at left cDNA end
		$termcenter = ($strand eq "+")? $tstart - $qstart + $vector5len : $tend + $qstart - $vector5len ;
		$termcenter++; # as gff is 1-based
		$termstart = $termcenter - $termhintradius;
		$termend = $termcenter + $termhintradius;
		$termstart = 1 if ($termstart < 1);
		$termend = $tsize if ($termend > $tsize);
		$termstrand = (($strand eq '+' && $FrontTerminus eq "5TSS") || ($strand eq '-' && $FrontTerminus eq "5TNS"))? '+' : '-';
		@hint = ($termstart, $termend, $termstrand, $grpname, $cdnaname, 1);
		if ($FrontTerminus eq "5TSS"){
		    addSignalHint(\@tsshints, [@hint]);
		} else {
		    addSignalHint(\@ttshints, [@hint]);
		}
	    }
	    if (($EndTerminus eq "3TSS" || $EndTerminus eq "3TNS") && $qsize-$qend - $vector3len < 3) { # terminus at right cDNA end
		$termcenter = ($strand eq "+")? $tend + $qsize - $qend - $vector3len : $tstart - $qsize + $qend + $vector3len;
		$termcenter++; # as gff is 1-based
		$termstart = $termcenter - $termhintradius;
		$termend = $termcenter + $termhintradius;
		$termstart = 1 if ($termstart < 1);
		$termend = $tsize if ($termend > $tsize);
		$termstrand = (($strand eq '+' && $EndTerminus eq "3TSS") || ($strand  eq '-' && $EndTerminus eq "3TNS"))? '+' : '-';
		@hint = ($termstart, $termend, $termstrand, $grpname, $cdnaname, 1);
		if ($EndTerminus eq "3TSS"){
		    addSignalHint(\@ttshints, [@hint]);
		} else {
		    addSignalHint(\@tsshints, [@hint]);
		}
	    }
#	    print HINTS "$targetname\t$prgsrc\tterminus\t$tstart\t$tend\t.\t$strand\t.\tNote \"$qname:" . $strand . (join ",", @{$tref}) . "\";\n";
	} else {
#	    print "no terminus info for $qname\n";
	}
    }

    # splice site hints, when the "cDNA is longer than the gDNA":
    # - there is unaligned cDNA left at the end
    # - the genomic DNA ends soon after the aligned part
    if ($trunkSS && !$intronsonly){
	# left genomic end
	if ($tstart < $maxDNAoverhang && $tstart > 5 &&
	    (($strand eq '+' && $qstart > $minESToverhang) || 
	     ($strand eq '-' && $qsize-$qend > $minESToverhang))){
		# small window [-4,+4] for up to 4 random matches, or up to 4 mismatches
	    @hint = ($tstart-4, $tstart+4, '.', $grpname, $cdnaname, 1);
		addSignalHint(\@dsshints, [@hint]);
		addSignalHint(\@asshints, [@hint]);
	    }
	# right genomic end
	if ($tsize-$tend < $maxDNAoverhang && $tsize-$tend > 5 &&
	    (($strand eq '+' && $qsize-$qend > $minESToverhang) || 
	     ($strand eq '-' && $qstart > $minESToverhang))){
		# small window [-4,+4] for up to 4 random matches, or up to 4 mismatches
		@hint = ($tend+2-4, $tend+2+4, '.', $grpname, $cdnaname, 1);
		addSignalHint(\@dsshints, [@hint]);
		addSignalHint(\@asshints, [@hint]);
	    }
	 
    }
}
close BLAT;

print "\n";
printHints();

close HINTS;

###########################################################################################
#
# subroutines
#
###########################################################################################

#
# printHints
# print and delete the hints
#
sub printHints {
    # finished computing the hints for the old sequence. output them.
    if($remove_redundant){
	# delete all exonpart hints that are contained in an exon hint
	my $startidx=0;
	my $curidx;
	foreach my $exon (@exonhints){
	    my $start = $exon->[0];
	    my $end = $exon->[1];
	    my $strand = $exon->[2];
	    # increase $startidx until the exonpart does not start to the left of start
	    while ($startidx <= $#exonparthints && @exonparthints[$startidx]->[0] < $start){
		$startidx++;
	    }
	    $curidx = $startidx;
	    while ($curidx <= $#exonparthints && @exonparthints[$curidx]->[0] <= $end){
		if (@exonparthints[$curidx]->[0] >= $start && @exonparthints[$curidx]->[1] <= $end && @exonparthints[$curidx]->[2] eq $strand ) {
		    #redundant, delete it
		    #print "deleting " , (join " ", @{$exonparthints[$curidx]}), " as it is contained in " , (join " ", @{$exon}), "\n";
		    splice @exonparthints, $curidx, 1;
		} else {
		    $curidx++;
		}
	    }
	}
    }
    # sort introns by start and then by end
    @intronhints = sort {$a->[0] <=> $b->[0] || $a->[1] <=> $b->[1]} @intronhints;
    if ($mult) {
	# summarize multiple identical intron hints to a single one
	my @sumintronhints = ();
	foreach my $href (@intronhints){
	    if (@sumintronhints>0 && $href->[0] == $sumintronhints[-1]->[0] && $href->[1] == $sumintronhints[-1]->[1]){
		$sumintronhints[-1]->[5]++;
	    } else {
		push @sumintronhints, $href;
	    }
	}
	@intronhints = @sumintronhints;
    }

    printTypeHints("tss", \@tsshints);
    printTypeHints("tts", \@ttshints);
    printTypeHints("ass", \@asshints);
    printTypeHints("dss", \@dsshints);
    printTypeHints($exonpartid, \@exonparthints);
    printTypeHints("exon", \@exonhints);
    printTypeHints("intron", \@intronhints);

    # delete all hints as the new sequence starts
    @tsshints = ();
    @ttshints = ();
    @dsshints = ();
    @asshints = ();
    @exonhints = ();
    @exonparthints = ();
    @intronhints = ();
}

sub printTypeHints {
    my $type = shift;
    my $hintsarray = shift;
    foreach my $href (@{$hintsarray}) {
	print HINTS "$oldtargetname\t$prgsrc\t$type\t$href->[0]\t$href->[1]\t$score\t$href->[2]\t.\t";
	print HINTS "grp=$href->[3];" if ($href->[5]==1);
	print HINTS "cdna=$href->[4];" if ($href->[4] ne "");
	print HINTS "mult=$href->[5];" if ($href->[5]>1); # multiplicity of hint
	print HINTS "pri=$priority;src=$source\n";
    }
}

#
# addExonpartHint(hintref)
# search in the list of existing exonpart hints for an including or included one
# if no such hint exists, sort the parameter hint into this list
# if there is a stronger one, do nothing, if there are weaker ones, replace them with this one
#

sub addExonpartHint {
    my $href = shift;
    my $begin = $href->[0];
    my $end = $href->[1];
    my $strand = $href->[2];
    my $k;
    #print (join " ", @{$href});
    #print "\n";
    
    my $redundant = 0;
    if ($remove_redundant) {
	# check whether the exonpart hint is contained in one of the exon hints.
	$k = $#exonhints;
	# check the list of previous exonpart hints
	if ($#exonparthints>=0) {
	    $k = $#exonparthints;
	    while ($k>=0 && $exonparthints[$k]->[0]> $begin - 10000 && !$redundant) { #assume exonpart hints are less than 10000bp
		if ($exonparthints[$k]->[0]<= $begin && $exonparthints[$k]->[1] >= $end && $exonparthints[$k]->[2] eq $strand){
		    #print "found including hint: ", (join " ", @{$exonparthints[$k]}), "\n";
		    $redundant=1;
		} elsif ($exonparthints[$k]->[0] >= $begin && $exonparthints[$k]->[1] <= $end && $exonparthints[$k]->[2] eq $strand){
		    #print "found included hint: ", (join " ", @{$exonparthints[$k]}), " delete it now.\n";
		    splice @exonparthints, $k, 1; #delete k-th element
		}
		$k--;
	    }
	}
    }
    if (!$redundant) {
	#insert hint at the right position
	#print "found no redundant hint\n";	
	$k = $#exonparthints;
	if ($remove_redundant) {
	    while ($k>=0 && $exonparthints[$k]->[0] > $begin) {
		$k--;
	    }
	}
	my @temparray = ($href);
	if ($k == $#exonparthints) {
	    #print "insert at end\n";
	    push @exonparthints, @temparray;
	} else {
	    #print "*** splicing list ***\n";
	    splice (@exonparthints, $k+1, 0, @temparray);
	}
    }
#    print "\n----------------------------\nexonpart hints\n";
#    foreach $hintref (@exonparthints) {
#	print (join ", ", @{$hintref});
#	print "\n"; 
#    } 
#    print "\n----------------------------\n";
}

#
# addSignalHint(hintref)
# add the hint if it is not already there

sub addSignalHint {
    my $hintlistref = shift;
    my $href = shift;
    my $begin = $href->[0];
    my $end = $href->[1];
    my $strand = $href->[2];
    #print (join " ", @{$href});

    # add it if the same hint does not exist already
    if (@{$hintlistref}<1 || !$remove_redundant) {
	push @{$hintlistref}, $href;
    } else {
	my $k = @{$hintlistref}-1;
	# sort by 1. begin 2. end 3. strand
	while ($k>=0 &&
	       ($hintlistref->[$k]->[0] > $begin ||
		($hintlistref->[$k]->[0] == $begin && $hintlistref->[$k]->[1] > $end) ||
		($hintlistref->[$k]->[0] == $begin && $hintlistref->[$k]->[1] == $end && $hintlistref->[$k]->[2] ge $strand))) {
	    $k--;
	}
	my @temparray = ($href);
	if (!(($k+1 <= @{$hintlistref}-1 && $hintlistref->[$k+1]->[0] == $begin && $hintlistref->[$k+1]->[1] == $end && $hintlistref->[$k+1]->[2] eq $strand))) {
	    # hint does not previously exist, insert it
	    if ($k== @{$hintlistref}-1) {
		push @{$hintlistref}, @temparray;
	    } else {
		splice (@{$hintlistref}, $k+1, 0, @temparray);
	    }
	}
    }
#    print "\n----------------------------\n hints\n";
#    foreach $hintref (@{$hintlistref}) {
#	print (join ", ", @{$hintref});
#	print "\n"; 
#    } 
#    print "\n----------------------------\n";
}

#
# insertIntervalHint(hintref)
# for exon and intron hints (not exonpart)
# add the hint if it is not already there

sub insertIntervalHint {
    my $hintlistref = shift;
    my $href = shift;
    my $begin = $href->[0];
    my $end = $href->[1];
    my $strand = $href->[2];
    #print (join " ", @{$href});
    #print "\n";

    # shortcut to add all hints, regardless whether they are mutiples
    #push @{$hintlistref}, $href;
    #return;

    # add it if the same hint does not exist already
    if (@{$hintlistref}<1) {
	push @{$hintlistref}, $href;
    } else {
	my $k = @{$hintlistref}-1;
	if ($remove_redundant) {
	    while ($k>=0 && ($hintlistref->[$k]->[0]> $begin || ($hintlistref->[$k]->[0]== $begin && $hintlistref->[$k]->[1]>= $end))) {
		$k--;
	    }
	}
	my @temparray = ($href);
	if ( !$remove_redundant || 
	    !(($k+1 <= @{$hintlistref}-1 && $hintlistref->[$k+1]->[0]== $begin && $hintlistref->[$k+1]->[1]== $end && $hintlistref->[$k+1]->[2] eq $strand) ||
	    ($k+2 <= @{$hintlistref}-1 && $hintlistref->[$k+2]->[0]== $begin && $hintlistref->[$k+2]->[1]== $end && $hintlistref->[$k+2]->[2] eq $strand))
	    ) {
	    # hint does not previously exist, insert it
	    if ($k== @{$hintlistref}-1) {
		push @{$hintlistref}, @temparray;
	    } else {
		splice (@{$hintlistref}, $k+1, 0, @temparray);
	    }
	}
    }
#    print "\n----------------------------\n interval hints\n";
#    foreach $hintref (@{$hintlistref}) {
#	print (join ", ", @{$hintref});
#	print "\n"; 
#    } 
#    print "\n----------------------------\n";
}


#
# addIntervalHint(hintref)
# for exon and intron hints (not exonpart)
# add the hint if it is not already there

sub addIntervalHint {
    my $hintlistref = shift;
    my $href = shift;
    push @{$hintlistref}, $href;
}

