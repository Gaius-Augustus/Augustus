# zebrafish trained on May 9th using the RefSeq annotation on assembly version danRer7
# for further info see /fs1.cm/species/zebrafish
# CDS model only

# performance on holdout set of 400 zebrafish genes
zebrafish      base level      0.895   0.711
zebrafish      exon level      0.802   0.71
zebrafish      gene level      0.275   0.193

# HMM training performance is significantly worse:
zebrafish-HMM  base level      0.937   0.554
zebrafish-HMM  exon level      0.807   0.564
zebrafish-HMM  gene level      0.268   0.109

# cross-species performance on same holdout set
chicken base level      0.881   0.669
chicken exon level      0.736   0.661
chicken gene level      0.182   0.125

human   base level      0.872   0.624
human   exon level      0.718   0.593
human   gene level      0.168   0.104

# conclusion: zebrafish model is better than chicken is better than human on zebrafish genome




