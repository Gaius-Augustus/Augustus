#!/bin/sh
for species in hg38 rheMac3 mm10 rn6 oryCun2 bosTau8 canFam3 loxAfr3 echTel2 dasNov3 monDom5 galGal4; do
    "$1" getfasta -fi "$2"$species.fa -bed "$3"out$4prepare/$species.bed > "$5$species".MINIMAL.fasta
done

