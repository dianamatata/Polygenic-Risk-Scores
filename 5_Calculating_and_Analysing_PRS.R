#!/bin/bash

# after PRS QC of GWAS summary and target data

# Calculating and Analysing PRS

# In the previous sections, we have generated the following files:
# Height.QC.gz	The post-QCed summary statistic
# EUR.QC.bed	The genotype file after performing some basic filtering
# EUR.QC.bim	This file contains the SNPs that passed the basic filtering
# EUR.QC.fam	This file contains the samples that passed the basic filtering
# EUR.height	This file contains the phenotype of the samples
# EUR.cov	This file contains the covariates of the samples


# Update Effect Size #

# When the effect size relates to disease risk and is thus given as an odds ratio (OR), rather than BETA (for continuous traits), then the PRS is computed as a product of ORs. 
# To simplify this calculation, we take the natural logarithm of the OR so that the PRS can be computed using summation instead (which can be back-transformed afterwards). 
# We can obtain the transformed summary statistics with R

#!/usr/bin/Rscript

dat <- read.table(gzfile("Height.QC.gz"), header=T)
dat$BETA <- log(dat$OR)
write.table(dat, "Height.QC.Transformed", quote=F, row.names=F)
q() # exit R

#!/bin/bash


# Clumping
# Linkage disequilibrium, which corresponds to the correlation between the genotypes of genetic variants across the genome, 
# makes identifying the contribution from causal independent genetic variants extremely challenging. 
# One way of approximately capturing the right level of causal signal is to perform clumping, which removes SNPs in ways that only weakly correlated SNPs are retained 
# but preferentially retaining the SNPs most associated with the phenotype under study. 
# Clumping can be performed using the following command in plink:

plink \
    --bfile EUR.QC \
    --clump-p1 1 \
    --clump-r2 0.1 \
    --clump-kb 250 \
    --clump Height.QC.Transformed \
    --clump-snp-field SNP \
    --clump-field P \
    --out EUR

# 489805 variants loaded from .bim file.
# 483 people (232 males, 251 females) loaded from .fam.
# 489805 variants and 483 people pass filters and QC.
# Note: No phenotypes present.

# clump-p1	1	P-value threshold for a SNP to be included as an index SNP. 1 is selected such that all SNPs are include for clumping
# clump-r2	0.1	SNPs having r2 higher than 0.1 with the index SNPs will be removed
# clump-kb	250	SNPs within 250k of the index SNP are considered for clumping
# clump	Height.QC.Transformed	Base data (summary statistic) file containing the P-value information
# clump-snp-field	SNP	Specifies that the column SNP contains the SNP IDs
# clump-field	P	Specifies that the column P contains the P-value information

# This will generate EUR.clumped, containing the index SNPs after clumping is performed. We can extract the index SNP ID by performing the following command:

awk 'NR!=1{print $3}' EUR.clumped >  EUR.valid.snp
# $3 because the third column contains the SNP ID

# If your target data are small (e.g. N < 500) then you can use the 1000 Genomes Project samples for the LD calculation. 
# Make sure to use the population that most closely reflects represents the base sample.


# Generate PRS #

# plink provides a convenient function --score and --q-score-range for calculating polygenic scores.

# Here calculate PRS corresponding to a few thresholds for illustration purposes:
# For example, for the 0.05 threshold, we include all SNPs with P-value from 0 to 0.05, including any SNPs with P-value equal to 0.05.

echo "0.001 0 0.001" > range_list 
echo "0.05 0 0.05" >> range_list
echo "0.1 0 0.1" >> range_list
echo "0.2 0 0.2" >> range_list
echo "0.3 0 0.3" >> range_list
echo "0.4 0 0.4" >> range_list
echo "0.5 0 0.5" >> range_list

plink \
    --bfile EUR.QC \
    --score Height.QC.Transformed 3 4 12 header \
    --q-score-range range_list SNP.pvalue \
    --extract EUR.valid.snp \
    --out EUR


# score	Height.QC.Transformed 3 4 12 header	
#       We read from the Height.QC.Transformed file, assuming that the 3st is SNP ID; 4th  is the effective allele information; the 12th  is the effect size estimate; and that the file contains a header
# q-score-range	range_list SNP.pvalue	We want to calculate PRS based on the thresholds defined in range_list, where the threshold values (P-values) were stored in SNP.pvalue

# the above command and range_list will generate 7 files: of type EUR.0.5.profile
# https://choishingwan.github.io/PRS-Tutorial/plink/ for the exact plink formula used



# Accounting for Population Stratification #

# Population structure is the principal source of confounding in GWAS and is usually accounted for by incorporating principal components (PCs) as covariates. 
# We can incorporate PCs into our PRS analysis to account for population stratification.


# First, we need to perform prunning
plink \
    --bfile EUR.QC \
    --indep-pairwise 200 50 0.25 \
    --out EUR
    
# Then we calculate the first 6 PCs
plink \
    --bfile EUR.QC \
    --extract EUR.prune.in \
    --pca 6 \
    --out EUR






