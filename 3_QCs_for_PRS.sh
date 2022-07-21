#!/bin/bash

# Sex chromosomes #

# Sometimes sample mislabelling can occur, which may lead to invalid results. islabelled sample if difference in sex and gender identity.
# A sex check can be performed in PLINK, in which individuals are called as females if their X chromosome homozygosity estimate (F statistic) is < 0.2 and as males if the estimate is > 0.8.

plink \
    --bfile EUR \
    --extract EUR.QC.prune.in \
    --keep EUR.valid.sample \
    --check-sex \
    --out EUR.QC

# This will generate a file called EUR.QC.sexcheck containing the F-statistics for each individual.

#!/usr/bin/Rscript

# Read in file
valid <- read.table("EUR.valid.sample", header=T)
dat <- read.table("EUR.QC.sexcheck", header=T)
valid <- subset(dat, STATUS=="OK" & FID %in% valid$FID)
write.table(valid[,c("FID", "IID")], "EUR.QC.valid", row.names=F, col.names=F, sep="\t", quote=F) 
q() # exit R


# Sample overlap#
# not here because 2 different datasets but check otherwise
