#!/usr/bin/env perl

# Author: Katharina J. Hoff, December 18th 2012
#
# This script runs Marias UTR finding tool with a set of given parameters (aim: find optimal parameters)

if(scalar(@ARGV) != 0){
    die("\nThis script runs Maria's UTR finding tool with a set of given parameters and files. Beware: parameters and files are hard-coded in the script!\n")
}

# specify executables
my $findUtr = "/home/katharina/SVN/utrrnaseq/trunks/Debug/main";
my $gff2lst = "/home/katharina/SVN/utrrnaseq/trunks/scripts/gff2lst.pl";
my $overlapStat = "/usr/local/bin/overlapStat.pl";

# specify input files
my $genomeFile = "/home/katharina/SVN/utrrnaseq/trunks/input/chr19.fa";
my $startStopFile = "/home/katharina/SVN/utrrnaseq/trunks/input/final.start_stop.gff";
my $intronFile = "/home/katharina/SVN/utrrnaseq/trunks/input/lib1.introns.ff.gff";
my $repeatFile = "/home/katharina/SVN/utrrnaseq/trunks/input/repeats.gff";
my $wiggleFile = "/home/katharina/SVN/utrrnaseq/trunks/input/trimmed05.blat.sf.wig";

# specify reference file
my $utrRef = "/home/katharina/SVN/utrrnaseq/trunks/input/final.ref.utrs.lst";

# specify output file
my $output = "/home/katharina/utrrnaseq.out";

# specify parameters that are to be tested:
my @w = (1, 20, 40, 60, 80);
my @v = (10, 50, 100, 150, 200);
my @c = (1, 10, 50, 100, 200);
my @p = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6);
my @i = (0.3, 0.5, 0.7);
my @m = (0.1, 0.2, 0.3);
my @z = ("false", "true");

# test whether exectuables exist
if(not(-e $findUtr)){
    die("findUtr executable $findUtr does not exist! Please check script source code and set a correct executable!\n");
}
if(not(-e $gff2lst)){
    die("gff2lst executable $gff2lst does not exist! Please check script source code and set a correct exectuable!\n");
}
if(not(-e $overlapStat)){
    die("overlapStat exectuable $overlapStat does not exist! Please check script source code and set a correct exectuable!\n");
}

# test whether input files exist
if(not(-e $genomeFile)){
    die("Genome file $genomeFile does not exist!\n");
}
if(not(-e $startStopFile)){
    die("StartStopFile $startStopFile does not exist!\n");
}
if(not(-e $intronFile)){
    die("intronFile $intronFile does not exist!\n");
}
if(not(-e $repeatFile)){
    die("repeatFile $repeatFile does not exist!\n");
}
if(not(-e $wiggleFile)){
    die("wiggleFile $wiggleFile does not exist!\n");
}
if(not(-e $utrRef)){
    die("reference file $utrRef does not exist!\n");
}

# count number of annotated UTRs
my %seen;
my $annoUTRs;
my @t;
open(REF, "<", $utrRef) or die("Could not open reference UTR file $utrRef!\n");
while(<REF>){
    @t = split(/\t/);
    if(not(defined($seen{$t[0]}))){
	$seen{$t[0]} = 1;
    }
}
$annoUTRs = scalar(keys %seen);
close(REF) or die("Could not close reference UTR file $utrRef!\n");

for (keys %seen){
    delete $seen{$_};
}

# test whether parameters are set (may not be empty)
my $runs;
my $len;
$len = @w;
$runs = $len;
if($len == 0){die("Parameter set w must be defined! Please set at least one possible parameter!\n");}
$len = @v;
$runs = $runs * $len;
if($len == 0){die("Parameter set v must be defined! Please set at least one possible parameter!\n");}
$len = @c;
$runs = $runs * $len;
if($len== 0){die("Parameter set c must be defined! Please set at least one possible parameter!\n");}
$len = @p;
$runs = $runs * $len;
if($len== 0){die("Parameter set p must be defined! Please set at least one possible parameter!\n");}
$len = @i;
$runs = $runs * $len;
if($len== 0){die("Parameter set i must be defined! Please set at least one possible parameter!\n");}
$len = @m;
$runs = $runs * $len;
if($len== 0){die("Parameter set m must be defined! Please set at least one possible parameter!\n");}
$len = @z;
$runs = $runs * $len;
if($len== 0){die("Parameter set z must be defined! Please set at least one possible parameter!\n");}

# setting variables
my $this_w;
my $this_v;
my $this_c;
my $this_p;
my $this_m;
my $this_z;

# create a temporary directory
my $stemDir = "/tmp/";
if(not(-w $stemDir)){die("No writing permissions for temporary stem directory $stemDir!\n");}
my @chars = ( "A" .. "Z", "a" .. "z", 0 .. 9);
my $randSubDir = join("", @chars[ map { rand @chars } ( 1 .. 8 ) ]);
my $tempDir = $stemDir.$randSubDir;
unless(-d $tempDir){
    mkdir $tempDir or die("Failed in creating the temporary directory $tempDir. You need to set the variable stemDir $stemDir to a place where you have writing permissions!\n");
}

# temporary prediction gff file
my $tempOutGff = $tempDir."/utr.gff";
# temporary prediction lst file
my $tempOutLst = $tempDir."/utr.lst";
# temporary overlapOutput
my $overlapOut = $tempDir."/ov.out";

# output variables
my $nuc_tp_tot;
my $nuc_fn_tot;
my $nuc_fp_tot;
my $nuc_tp;
my $nuc_fn;
my $nuc_fp;
my $nuc_sens;
my $nuc_spec;
my $nuc_hm;
my $nuc_sens_tot;
my $nuc_spec_tot;
my $nuc_hm_tot;
my $utr_sens;

open(OUTPUT, ">", $output) or die ("Could not open output file $output!\n");
print OUTPUT "-w\t-v\t-c\t-p\t-i\t-m\t-z\tnAnnoUTRs\tnFoundUTRs\tUTRsensitivity\tnuc_tp_total\tnuc_fp_total\tnuc_fn_total\tnuc_sensitivity_total\tnuc_specificity_total\tnuc_harmonic_mean_total\tnuc_tp_detected\tnuc_fp_detected\tnuc_fn_detected\tnuc_sensitivity_detected\tnuc_specificity_detected\tnuc_harmonic_mean_detected\n";

# test whether the program call itself does run at all, no loop, yet
my $call;
my $alreadyRun = 0;
foreach(@w){
    $this_w = $_;
    foreach(@v){
	$this_v = $_;
	foreach(@c){
	    $this_c = $_;
	    foreach(@p){
		$this_p = $_;
		foreach(@i){
		    $this_i = $_;
		    foreach(@m){
			$this_m = $_;
			foreach(@z){
			    $alreadyRun = $alreadyRun + 1;
			    print "### Processing $alreadyRun out of $runs ... that means ".(100*$alreadyRun/$runs)."% are soon to be finished!\n";
			    $this_z = $_;
			    print OUTPUT "$this_w\t$this_v\t$this_c\t$this_p\t$this_i\t$this_m\t$this_z\t";
			    # execute the tool
			    $call  = "$findUtr -G $genomeFile -W $wiggleFile -I $intronFile -O $startStopFile -R $repeatFile -o $tempOutGff -w $this_w -v $this_v -c $this_c -p $this_p -i $this_i -m $this_m -z $this_z 1> $tempDir/findUtr.out 2> $tempDir/findUtr.err";
			    print "executing: $call\n";
			    system $call;
			    print "Done!\n";
			    # convert output
			    $call = "$gff2lst $tempOutGff 1> $tempOutLst 2> $tempDir/gff2lst.err";
			    print "executing: $call\n";
			    system $call;
			    # analyse pure nucleotide accuracy
			    $call = "$overlapStat $utrRef $tempOutLst 1> $tempDir/overlapStat.out 2> $tempDir/overlapStat.err";
			    print "executing: $call\n";
			    system $call;
			    open(NUC, "<", "$tempDir/overlapStat.out") or die("Could not open overlapStat output file $tempDir/overlapStat.out!\n");
			    while(<NUC>){
				if(m/combination  1 0/){
				    $_=~m/\| (\d+)/;
				    $nuc_fn_tot = $1;
				}elsif(m/combination  0 1/){
				    $_=~m/\| (\d+)/;
				    $nuc_fp_tot = $1;
				}elsif(m/combination  1 1/){
				    $_=~m/\| (\d+)/;
				    $nuc_tp_tot = $1;
				}
			    }
			    close(NUC) or die("Could not close overlapStat output file $tempDir/overlapStat.out!\n");
			    $nuc_sens_tot = $nuc_tp_tot/($nuc_tp_tot+$nuc_fn_tot);
			    $nuc_spec_tot = $nuc_tp_tot/($nuc_tp_tot+$nuc_fp_tot);
			    $nuc_hm_tot = 2*$nuc_sens_tot*$nuc_spec_tot/($nuc_sens_tot+$nuc_spec_tot);
			    print "Total nucleotide sensitivity: $nuc_sens_tot\n";
			    print "Total nucleotide specificity: $nuc_spec_tot\n";
			    print "Total nucleotide harmonic mean: $nuc_hm_tot\n";
			    # analyse nucleotide accuracy in detected UTRs, only
			    open(LST, "<", $tempOutLst) or die("Could not open temporary list output file $tempOutLst!\n");
			    open(FOUND, ">", "$tempDir/foundUTRs") or die("Could not open temporary found UTRs file $tempDir/foundUTRs!\n");
			    for (keys %seen)
			    {
				delete $seen{$_};
			    }
			    while(<LST>){
				@t = split(/\t/);
				if(not(defined($seen{$t[0]}))){
				    $seen{$t[0]} = 1;
				    print FOUND $t[0]."\n";
				}
				$len = scalar(keys %seen);
			    }
			    close(FOUND) or die("Could not close temporary found UTRs file $tempDir/foundUTRs!\n");
			    close(LST) or die("Could not close temporary list output file $tempOutLst!\n");
			    $call = "grep -f $tempDir/foundUTRs -F $utrRef 1> $tempDir/temp.ref 2> $tempDir/grep1.err";
			    print "executing: $call\n";
			    system $call;

			    $call = "$overlapStat $tempDir/temp.ref $tempOutLst 1> $tempDir/overlapStat2.out 2> $tempDir/overlapStat2.err";
			    print "executing: $call\n";
			    system $call;
			    open(NUC, "<", "$tempDir/overlapStat2.out") or die("Could not open overlapStat2 output file $tempDir/overlapStat2.out!\n");
			    while(<NUC>){
				if(m/combination  1 0/){
				    $_=~m/\| (\d+)/;
				    $nuc_fn = $1;
				}elsif(m/combination  0 1/){
				    $_=~m/\| (\d+)/;
				    $nuc_fp = $1;
				}elsif(m/combination  1 1/){
				    $_=~m/\| (\d+)/;
				    $nuc_tp = $1;
				}
			    }
			    close(NUC) or die("Could not close overlapStat2 output file $tempDir/overlapStat2.out!\n");
			    print "Detected $len UTRs out of $annoUTRs !\n";
			    $nuc_sens = $nuc_tp/($nuc_tp+$nuc_fn);
			    $nuc_spec = $nuc_tp/($nuc_tp+$nuc_fp);
			    $nuc_hm = 2 * $nuc_sens * $nuc_spec /( $nuc_sens + $nuc_spec);
			    print "Nucleotide sensitivity in detected UTRs: $nuc_sens\n";
			    print "Nucleotide specificity in detected UTRs: $nuc_spec\n";
			    print "Nucleotide harmonic mean in detected UTRs: $nuc_hm\n";
                            $utr_sens = $len/$annoUTRs;
			    print OUTPUT "$annoUTRs\t$len\t$utr_sens\t";
                            print OUTPUT "$nuc_tp_tot\t$nuc_fp_tot\t$nuc_fn_tot\t$nuc_sens_tot\t$nuc_spec_tot\t$nuc_hm_tot\t";
			    print OUTPUT "$nuc_tp\t$nuc_fp\t$nuc_fn\t$nuc_sens\t$nuc_spec\t$nuc_hm\n";
			    print "UTR detection sensitivity: $utr_sens\n";
			}

		    }
		}
	    }
	}
    }
}

# clean up
`cd $stemDir; rm -r $randSubDir;`;
close OUTPUT or die ("Could not close output file $output!\n");
