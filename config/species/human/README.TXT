retrained on 1.1.2010
 see documentation in mario@giovanni:/data/human/training
# retraining human for rGASP 1.1 on 1.1.2010
# trainingset 1784 random RefSeq genes that are neither in an ENCODE region nor on chr3 or chr4
#
# see also /c1/scratch/mario/test/human
# see also mario@monty:~/human/hg18/refseq/gb
# using augustus revision r312

# HMM-train UTR, < 1min
etraining --species=human train.1784.gb --UTR=on

# CRF-train coding model, about 7h
etraining --species=human train.1784.gb --CRF=on --CRF_N=2 --UTR=off --CRFtrainSS=1 --CRFtrainIntron=1 --tieIgenicIntron=1

# performance on 10MB chr21 against RefSeq, using masked regions as hints:
# much better than previous version

/c1/project/rGASP/augustus/hs/hs.crftest/eval.24.out (with UTR)
/c1/project/rGASP/augustus/hs/hs.crftest/eval.21.out (only CDS)
/c1/project/rGASP/augustus/hs/hs.crftest/eval.0.out  (old human version, with UTR)
