Contributed by Shishir K Gupta, Dept. Bioinformatik, Uni Wuerzburg in October 2013.

"I initially used scipio for training set. The protein set I choosed for scipio, was the proteins collected from the different GO subclasses. Then I used scipio + 107 highly conserved proteins in 7 sequenced ants. This parameters were best for me. Then also used Pasa and Cegma sets but it could not improve much more accuracy. I also generated hints from different sources (assembled RNAseq, raw RNASeq, EST and repeatmasker) and did the final predictions.  My training set had 330 genes and test set had 100 genes."

*******      Evaluation of gene prediction     *******

---------------------------------------------\
                 | sensitivity | specificity |
---------------------------------------------|
nucleotide level |       0.953 |       0.906 |
---------------------------------------------/

----------------------------------------------------------------------------------------------------------\
           |  #pred |  #anno |      |    FP = false pos. |    FN = false neg. |             |             |
           | total/ | total/ |   TP |--------------------|--------------------| sensitivity | specificity |
           | unique | unique |      | part | ovlp | wrng | part | ovlp | wrng |             |             |
----------------------------------------------------------------------------------------------------------|
           |        |        |      |                143 |                123 |             |             |
exon level |    702 |    682 |  559 | ------------------ | ------------------ |        0.82 |       0.796 |
           |    702 |    682 |      |   58 |    2 |   83 |   59 |    2 |   62 |             |             |
----------------------------------------------------------------------------------------------------------/

----------------------------------------------------------------------------\
transcript | #pred | #anno |   TP |   FP |   FN | sensitivity | specificity |
----------------------------------------------------------------------------|
gene level |   128 |   100 |   42 |   86 |   58 |        0.42 |       0.328 |
----------------------------------------------------------------------------/