#!/bin/bash

# Relatedness #

# Closely related individuals in the target data may lead to overfitted results, limiting the generalisability of the results.
# Individuals that have a first or second degree relative in the sample (π>0.125)

#!/bin/bash


plink \
    --bfile EUR \
    --extract EUR.QC.prune.in \
    --keep EUR.QC.valid \
    --rel-cutoff 0.125 \
    --out EUR.

# 503 people (240 males, 263 females) loaded from .fam.
# --extract: 268457 variants remaining.
# --keep: 483 people remaining.
# 268457 variants and 483 people pass filters and QC (before --rel-cutoff).
# Excluding 5229 variants on non-autosomes from relationship matrix calc.
# Remaining sample IDs written to EUR.QC.rel.id .

# A greedy algorithm is used to remove closely related individuals in a way that optimizes the size of the sample retained. 
# However, the algorithm is dependent on the random seed used, which can generate different results. 
# PLINK's algorithm for removing related individuals does not account for the phenotype under study. 
# To minimize the removal of cases of a disease, the following algorithm can be used instead: GreedyRelated.


# Generate final QC'ed target data file #

plink \
    --bfile EUR \
    --make-bed \
    --keep EUR.QC.rel.id \
    --out EUR.QC \
    --extract EUR.QC.snplist \
    --exclude EUR.mismatch \
    --a1-allele EUR.a1

# bfile	EUR	Informs plink that the input genotype files should have a prefix of EUR
# keep	EUR.QC.rel.id	Informs plink that we only want to keep samples in EUR.QC.rel.id
# extract	EUR.QC.snplist	Informs plink that we only want to use SNPs in EUR.QC.snplist in the analysis
# exclude	EUR.mismatch	Informs plink that we wish to remove any SNPs in EUR.mismatch
# a1-allele	EUR.a1	Fix all A1 alleles to those specified in EUR.a1
# out	EUR.QC	Informs plink that all output should have a prefix of EUR.QC




