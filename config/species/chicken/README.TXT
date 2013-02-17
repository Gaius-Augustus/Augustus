# chicken training from 17.02.2013
# genome downloaded from http://hgdownload.soe.ucsc.edu/goldenPath/galGal4/bigZips/galGal4.fa.gz
# training set: 3424 RefSeq genes (also from UCSC browser)

new_species.pl --species=chicken

# all parameter values in chicken_parameters.cfg are set to the human default parameters in human_parameters.cfg
# 10 different gc-content classes are used: /Constant/decomp_num_steps=10 (min range=0.37, max range=0.61 -> recommended values
# from etraining)

etraining --species=chicken genes.train.3423.gb | tee chicken_etrain
augustus --species=chicken genes.test.200.gb | tee chicken_pred

# estimated accuracy:
#####################################
#                SN           SP    #
#####################################
# gene level #  0.28        0.209   #
# exon level #  0.834       0.803   #
# nuc. level #  0.92        0.846   #
#####################################

# see also
stefanie@greifserv2:/fs1.data/species/birds/gallus/training_on_unmasked_genome/