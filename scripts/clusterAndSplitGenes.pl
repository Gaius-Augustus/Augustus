#!/usr/bin/perl
#
# read in a file with inter-gene dependencies and
# create clusters without inter-cluster dependencies
# optionally split transcript-based files to smaller chunks
#
# Mario Stanke, 18.9.2009, mstanke@gwdg.de
use strict;
use Getopt::Long;
use FileHandle;

my $help = 0;
my $verbose = 0;
my $genefile;
my $genesprefix;
my $mapfile;
my $mapprefix;
my $seqfile;
my $clusterfile;
my $quantfile;
my $thresh = 10;
my $maxclustersize = 50;
my $chunksize = 3000;
my $numC = 0;
my $txformat = 0;
my $numClusterChunks = 0;


my $usage = "$0 -- read in a file with inter-gene dependencies and create clusters without\n";
$usage .= "inter-cluster dependencies. Optionally split transcript-based files to smaller chunks\n";
$usage .= "\n";
$usage .= "Usage: $0 < commonfile\n";
$usage .= "Options\n";
$usage .= " --verbose  \n";
$usage .= "   -t                     threshold: minimal number of dependencies between two genes\n";
$usage .= "                          in order to force them into same cluster (default: $thresh)\n";
$usage .= "   --txformat             commonfile holds transcript ids instead of gene ids\n";
$usage .= "                          the mapping given by the --genes file is used to map transcript ids to gene ids\n";
$usage .= "   --chunksize=n          number of genes per chunk (default: $chunksize)\n";
$usage .= "   --genes=genes.gtf      gene set in GTF format that is split accoring to clustering\n";
$usage .= "   --genesprefix=s        prefix for GTF output, can include directory (default: value of --genes)\n";
$usage .= "   --filtermap=ali.psl  \n";
$usage .= "   --mapprefix=s          prefix for short alignment output (default: value of --filtermap)\n";
#$usage .= "   --mrna=mrna.fa  \n";
$usage .= "   --clusters=s           output file with clusters \n";
$usage .= "   --maxclustersize=n     maximal size of cluster (default: $maxclustersize)\n";



GetOptions(
    'help!'=>\$help,
    'verbose!'=>\$verbose,
    'txformat!'=>\$txformat,
    'genes=s'=>\$genefile,
    'genesprefix:s'=>\$genesprefix,
    'clusters:s'=>\$clusterfile,
    't:n'=>\$thresh,
    'maxclustersize:n'=>\$maxclustersize,
    'chunksize:n'=>\$chunksize,
    'filtermap=s'=>\$mapfile,
    'mapprefix:s'=>\$mapprefix,
    'mrna=s'=>\$seqfile);

if ($help) {
    print "$usage";
    exit(0);
}

$genesprefix = $genefile if (!defined($genesprefix));
my %deppairs = (); # keys seqnameAseqnameB (lexicographic)
                   # values: multiplicities

my %cbyg = ();  # keys: geneid, values: clusters
my %clusters = ();
# cluster: list of geneids
my %geneByTxid = (); # keys: transcript ids, values gene ids

#
# read assignment from txids to gene ids
#
if (defined($genefile) || $txformat){
    die ("No --genes specified but option --txformat. No tx2gene mapping given.\n") if (!defined($genefile));
    open (GENE, "<$genefile") or die ("Could not open $genefile.");
    my ($geneid, $txid, $line);
    while (<GENE>){
	$line = $_;
	if ($line =~ /gene_id."([^"]+)"/){
	    $geneid = $1;
	    # assign gene to transcript id
	    if ($line =~ /transcript_id."([^"]+)"/){
		$txid = $1;
		if (defined($geneByTxid{$txid}) && ($geneid ne $geneByTxid{$txid})){
		    die ("Transcript $txid was previously assigned to $geneid, now to $geneByTxid{$txid}.\n" .
			"First offending line:\n$line");
		}
		$geneByTxid{$txid} = $geneid;
	    }
	}
    }
    close GENE;
}


# parse dependency pairs
while(<>){
    chomp;
    my @f = split /\s+/;
    shift @f; # ignore query name
    if ($txformat) {
	@f = map{$geneByTxid{$_}} @f; # map txids to gene ids
	# make result list unique
	my %seen = ();
	@f = grep { ! $seen{$_} ++ } @f;
    }
    next if (@f<2); # single gene, no clustering necessary
    for (my $i=0; $i<@f;$i++){
	$cbyg{$f[$i]} = [$f[$i]] if (!defined($cbyg{$f[$i]})); # initialize cluster with single gene
	for (my $j=$i+1; $j<@f;$j++){
	    if ($f[$i] < $f[$j]){
		$deppairs{"$f[$i] <=> $f[$j]"}++;
	    } else {
		$deppairs{"$f[$j] <=> $f[$i]"}++;
	    }
	}
    }
}

# create clusters
foreach my $pair (sort {$deppairs{$b} <=> $deppairs{$a}} keys %deppairs){ # sort by decreasign number of observations of pair
    next if $deppairs{$pair} < $thresh;
    # join the two clusters if not already joined
    my ($g1, $g2) = split / <=> /, $pair;
    my ($c1, $c2) = ($cbyg{$g1}, $cbyg{$g2});
    if ($c1 != $c2 && @{$c1} + @{$c2} < $maxclustersize){
	# join clusters
	push @{$c1}, @{$c2};
	foreach my $g (@$c2){
	    $cbyg{$g} = $c1;
	}
    }
}

foreach my $g (keys %cbyg){
    if (@{$cbyg{$g}} < 2){
	delete $cbyg{$g} ;
	next;
    }
    $clusters{$cbyg{$g}} = $cbyg{$g} if (!defined($clusters{$cbyg{$g}}));
}

# print clusters
if (defined($clusterfile)){
    open (CLUSTER, ">$clusterfile") or die ("Could not open $clusterfile for writing.");
    foreach my $cref (keys %clusters){
	print CLUSTER "" . join("\t", @{$clusters{$cref}}) . "\n";
    }
    close (CLUSTER);
}

# make chunks from clusters
my %chunkbyc =  (); # keys: clusterrefs
my %chunkbyg =  (); # keys: geneids
my $ng = 0;
foreach my $cref (keys %clusters){
    if ($numC == 0 || $ng > $chunksize/3){ # make chunks from clusers only one third full, they are harder to process later
	$ng = 0;
	$numC++;
    }
    $chunkbyc{$cref} = $numC;
    $ng += scalar(@{$clusters{$cref}});
}
$numClusterChunks = $numC;
print "Creating $numClusterChunks chunks with clusters.\n" if ($verbose);

# now split input files into chunks
if (defined($genefile)){
    my $chunkid = 1;
    my $warnMissing = 0;
    $ng = 0;
    my ($line,$geneid, $txid);
    my @fhs = (); # file handles, one for each chunk
    open (GENE, "<$genefile") or die ("Could not open $genefile.");
    while (<GENE>){
	next if (/^#/);
	$line = $_;
	if ($line =~ /gene_id."([^"]+)"/){
	    $geneid = $1;
	} else {
	    next; #die ("Could not find gene id in the following line:\n$line\n");
	}
	if (defined($cbyg{$geneid})){
	    $chunkid = $chunkbyc{$cbyg{$geneid}};
	} elsif (defined($chunkbyg{$geneid})) {
	    $chunkid = $chunkbyg{$geneid};
	} else  {
	    if ($numC == $numClusterChunks || $ng >= $chunksize){
		$ng = 0;
		$numC++;
	    }
	    $chunkid = $chunkbyg{$geneid} = $numC;
	    $ng += 1;
	}
	if (!defined($fhs[$chunkid])){
	    my $fh = FileHandle->new(">$genesprefix$chunkid");
	    $fhs[$chunkid] = $fh;
	}
	my $fh = $fhs[$chunkid];
	print $fh $line;
    }
    close (GENE);
    my $i=0;
    foreach my $fh (@fhs){
	if (defined($fh)){
	    $fh->close;
	} else {
	    if (!$warnMissing && $fh != $fhs[0]){ # 0th element does not exist
		print STDERR "Warning: Not all clusters have genes. If you used transcript names in commonfile\n";
		print STDERR "then invoke with --txformat flag.\n";
		$warnMissing = 1;
	    }
	}
    }

    # split map file into chunks
    @fhs = ();
    if (defined($mapfile)){
	open (MAP, "<$mapfile") or die ("Could not open $mapfile.");
	while (<MAP>){
	    $line = $_;
	    my @f = split /\t/, $line;
            $txid = $f[1];
	    die ("Transcript $txid has not been assigned a gene id\n") if (!defined($geneByTxid{$txid}));
	    $geneid = $geneByTxid{$txid};
	    if (defined($cbyg{$geneid})){
		$chunkid = $chunkbyc{$cbyg{$geneid}};
	    } elsif (defined($chunkbyg{$geneid})) {
		$chunkid = $chunkbyg{$geneid};
	    } else  {
		next; # only consider reads mapped to transcripts in GTF
	    }
	    if (!defined($fhs[$chunkid])){
		my $fh = FileHandle->new(">$mapprefix$chunkid");
		$fhs[$chunkid] = $fh;
	    }
	    my $fh = $fhs[$chunkid];
	    print $fh $line;
	}
    }
    close (MAP);
    foreach my $fh (@fhs){
        $fh->close if (defined($fh));
    }
}
