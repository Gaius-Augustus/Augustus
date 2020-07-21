# training and prediction of single-genome augustus
# Tests are longrunning because of CRF training and evaluations on larger chunks of genome:
# 9 regions of 2.5 MB each (human chr1, hg38)
# execute this script (single-threaded) with
# ~/Augustus/longrunning_examples/single_genome$ bash -v test_single.sh -e path/to/eval

# set given eval directory
while getopts e: option
do
case "${option}"
in
e) EVAL_DIR=${OPTARG};;
esac
done

export PERL5LIB=$EVAL_DIR # so Eval.pm is found
export AUGUSTUS_CONFIG_PATH=../../../config

########## create output folder
OUTDIR=output
EVAL_OUT_DIR=$OUTDIR/eval
mkdir -p $EVAL_OUT_DIR

########## download training data and sample genome data
DATADIR=data
TDIR=$DATADIR/training
mkdir -p $TDIR

# download training set
TRSET=$TDIR/train.1784.gb
if [ ! -e $TRSET ]; then
   wget http://bioinf.uni-greifswald.de/bioinf/downloads/data/aug-test/train.1784.gb.gz -O - | gunzip -c > $TRSET
fi

# download human chromosome 1 (hg38)
TESTSEQ=$DATADIR/chr1.fa.gz
if [ ! -e $TESTSEQ ]; then
   wget http://bioinf.uni-greifswald.de/bioinf/downloads/data/aug-test/chr1.fa.gz -O - > $TESTSEQ
fi

# download evaluation data
REFANNOFNAME=ensembl.ensembl_and_ensembl_havana.chr1.CDS.gtf.dupClean.FILTERED.gtf
REFANNO=$DATADIR/$REFANNOFNAME
if [ ! -e $REFANNO ]; then
   wget http://bioinf.uni-greifswald.de/bioinf/downloads/data/aug-test/$REFANNOFNAME -O - > $REFANNO
fi
   

########## HMM training
# create parameter set from template
HMMspecies=human_longrunningtest_hmm
if [ -d $AUGUSTUS_CONFIG_PATH/species/$HMMspecies ]; then
   rm -r $AUGUSTUS_CONFIG_PATH/species/$HMMspecies
fi
../../../scripts/new_species.pl --species=$HMMspecies

# train parameters
../../../bin/etraining --species=$HMMspecies $TRSET --UTR=on > $OUTDIR/etrain_hmm.out 2> $OUTDIR/etrain_hmm.err



########## CRF training
CRFspecies=human_longrunningtest_crf
if [ -d $AUGUSTUS_CONFIG_PATH/species/$CRFspecies ]; then
   rm -r $AUGUSTUS_CONFIG_PATH/species/$CRFspecies
fi   
../../../scripts/new_species.pl --species=$CRFspecies

# HMM-training to obtain UTR parameters
../../../bin/etraining --species=$CRFspecies $TRSET --UTR=on > /dev/null 2> /dev/null

# CRF training of CDS parameters
../../../bin/etraining --species=$CRFspecies $TRSET --CRF=on --UTR=off > $OUTDIR/etrain_crf.out 2> $OUTDIR/etrain_crf.err
# 46m



######## evaluation on long genomic regions
# using 9 regions on chr1, the same regions as used in multi-genome gene prediction (CGP)
#
# try and compare a few parameter sets (HMM, CRF, existing)
# on existing human parameters, try a few prediction options

numEvalRuns=6 # must match the size of the next three arrays
OPTIONSLIST=("--species=human --softmasking=0" \
	     "--species=human --softmasking=1" \
             "--species=human --softmasking=1 --UTR=1 --alternatives-from-sampling=1 --sample=100" \
             "--species=human --softmasking=1 --UTR=1" \
             "--species=$HMMspecies --softmasking=1 --UTR=1" \
	     "--species=$CRFspecies --softmasking=1 --UTR=1"\
	     )

OPTIONSDESCR=("standard human parameters, softmasking off" \
              "standard human parameters, softmasking" \
	      "standard human parameters, softmasking, UTR, alternatives-from-sampling" \
	      "standard human parameters, softmasking, UTR" \	      
      	      "HMM-trained parameters, softmasking, UTR" \
       	      "CRF-trained parameters, softmasking, UTR" \
	     )

OPTIONSNAMES=(human-nosm \
             human-sm \
	     human-sm-UTR-alt \
	     human-sm-UTR \
	     HMM-sm-UTR \
	     CRF-sm-UTR \
	     )


RFILE=$OUTDIR/pred-report.txt
rm -r $RFILE

for ((r=0; r<$numEvalRuns; r++)); do
   OPTIONS=${OPTIONSLIST[$r]}
   RUNNAME=${OPTIONSNAMES[$r]}
   DESCR=${OPTIONSDESCR[$r]}
   GFFFNAME=$OUTDIR/augustus-long-$RUNNAME.gff
   echo -e "$GFFFNAME\t$DESCR\trun options: $OPTIONS" | tee -a $RFILE

   DIR=`mktemp -d -p .`
   CHUNKLIST="27 30 47 54 57 80 86 101 118" # same as Giovanna uses in multi-genome longrunning test
   NUMCHUNKS=${#CHUNKLIST[@]}
   parallel -j 6 --will-cite ./aug_run_chunk.sh {} $TESTSEQ $DIR/augustus_{}.gff \"$OPTIONS\" ::: $CHUNKLIST
   # with sampling option this takes up to 4GB RAM per job
   
   for C in $CHUNKLIST; do
      cat $DIR/augustus_$C.gff
   done | ../../../scripts/join_aug_pred.pl > $GFFFNAME

   # count how many jobs finished
   NUM_FINISHED=`tail -n 2 $DIR/augustus_*.gff | grep "command line" | wc -l`
   if ((NUM_FINISHED < NUMCHUNKS)); then
       echo "Error when trying to compute $GFFFNAME\t$DESCR\trun options: $OPTIONS"
       echo "Please inspect $DIR"
       echo "ERROR: $RUNNAME.eval could not be computed."
   else
       rm -r $DIR
       $EVAL_DIR/evaluate_gtf.pl $REFANNO $GFFFNAME > $EVAL_OUT_DIR/$RUNNAME.eval
   fi
done

### cleanup intermediate results
# rm *.{gff,out,err}

# delete data with
# rm $REFANNO $TRSET $TESTSEQ 

# useful: head -n 14 *.eval
