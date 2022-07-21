#!/bin/bash

# PRS QC of GWAS summary and target data






# QC of Base Data ##############################################################


# Read GWAS summary statistics # 

gunzip -c Height.gwas.txt.gz | head


# Requirements # 

# - heritability check, greater then h^2 > 0.05
# 	The chip-heritability of a GWAS can be estimated using e.g. LD Score Regression (LDSC)
# - effect allele: which allele is the effect allele and which is the non-effect allele for PRS association results to be in the correct direction.
# 	Some GWAS results files do not make clear which allele is which. 
# - same genome built with target data


# File transfer # 

# if the file is intact, then md5sum generates a string of characters, which in this case should be: a2b15fb6a2bbbe7ef49f67959b43b160. If a different string is generated, then the file is corrupted.
md5sum Height.gwas.txt.gz


# Standard GWAS QC #

# If the base data have been obtained as summary statistics from a public source, then the typical QC steps that you will be able to perform on them are to filter the SNPs according to INFO score and MAF. 
# SNPs with low minor allele frequency (MAF) or imputation information score (INFO) are more likely to generate false positive results due to their lower statistical power (and higher probability of genotyping errors in the case of low MAF). 
# Therefore, SNPs with low MAF and INFO are typically removed before performing downstream analyses. 
# We recommend removing SNPs with MAF < 1% and INFO < 0.8 (with very large base sample sizes these thresholds could be reduced if sensitivity checks indicate reliable results). 

gunzip -c Height.gwas.txt.gz |\
awk 'NR==1 || ($11 > 0.01) && ($10 > 0.8) {print}' |\
gzip  > Height.gz

# Prints the header line (NR==1)
# Prints any line with MAF above 0.01 ($11 because the eleventh column of the file contains the MAF information)
# Prints any line with INFO above 0.8 ($10 because the tenth column of the file contains the INFO information)


# Duplicate SNPs#

# Most PRS software do not allow duplicated SNPs in the base data input and thus they should be removed

gunzip -c Height.gz |\
awk '{seen[$3]++; if(seen[$3]==1){ print}}' |\
gzip - > Height.nodup.gz

# Count number of time SNP ID was observed, assuming the third column contian the SNP ID (seen[$3]++). If this it the first time seeing this SNP ID, print it.


# Ambiguous SNPs#

# If the base and target data were generated using different genotyping chips and the chromosome strand (+/-) that was used for either is unknown, 
# then it is not possible to pair-up the alleles of ambiguous SNPs (i.e. those with complementary alleles, either C/G or A/T SNPs) across the data sets, 
# because it will be unknown whether the base and target data are referring to the same allele or not. 
# While allele frequencies could be used to infer which alleles are on the same strand, 
# the accuracy of this could be low for SNPs with MAF close to 50% or when the base and target data are from different populations. 
# Therefore, we recommend removing all ambiguous SNPs to avoid introducing this potential source of systematic error.

gunzip -c Height.nodup.gz |\
awk '!( ($4=="A" && $5=="T") || \
        ($4=="T" && $5=="A") || \
        ($4=="G" && $5=="C") || \
        ($4=="C" && $5=="G")) {print}' |\
    gzip > Height.QC.gz





# QC of Target Data ##############################################################

unzip EUR.zip

# Sample size: We recommend that users only perform PRS analyses on target data of at least 100 individuals.# 


# Check File transfer with md5sum # 


# Check same Genome Built # 


# Standard GWAS QC # 

# super important: Filtering should be performed on the control samples to avoid filtering SNPs that are causal 

plink \
    --bfile EUR \
    --maf 0.01 \
    --hwe 1e-6 \
    --geno 0.01 \
    --mind 0.01 \
    --write-snplist \
    --make-just-fam \
    --out EUR.QC

# bfile: Informs plink that the input genotype files should have a prefix of EUR
# maf: Removes all SNPs with minor allele frequency less than 0.01. Genotyping errors typically have a larger influence on SNPs with low MAF. Studies with large sample sizes could apply a lower MAF threshold
# hwe: Removes SNPs with low P-value from the Hardy-Weinberg Equilibrium Fisher's exact or chi-squared test. 
#      SNPs with significant P-values from the HWE test are more likely affected by genotyping error or the effects of natural selection. 
#      Filtering should be performed on the control samples to avoid filtering SNPs that are causal (under selection in cases). When phenotype information is included, plink will automatically perform the filtering in the controls.
# geno: Excludes SNPs that are missing in a high fraction of subjects. A two-stage filtering process is usually performed
#       here includes only SNPs with a 99% genotyping rate (1% missing)
# mind: Excludes individuals who have a high rate of genotype missingness, since this may indicate problems in the DNA sample or processing. 
#       here means exclude with more than 1% missing genotypes
# make-just-fam: informs plink to only generate the QC'ed sample name to avoid generating the .bed file: EUR.QC.fam
# write-snplist: Informs plink to only generate the QC'ed SNP list to avoid generating the .bed file.: EUR.QC.snplist


# Remove highly correlated SNPs #

# Very high or low heterozygosity rates in individuals could be due to DNA contamination or to high levels of inbreeding. 
# Therefore, samples with extreme heterozygosity are typically removed prior to downstream analyses.

plink \
    --bfile EUR \
    --keep EUR.QC.fam \
    --extract EUR.QC.snplist \
    --indep-pairwise 200 50 0.25 \
    --out EUR.QC

# bfile	EUR	Informs plink that the input genotype files should have a prefix of EUR
# keep	EUR.QC.fam	Informs plink that we only want to use samples in EUR.QC.fam in the analysis
# extract	EUR.QC.snplist	Informs plink that we only want to use SNPs in EUR.QC.snplist in the analysis
# indep-pairwise	200 50 0.25	Informs plink that we wish to perform pruning with a window size of 200 variants, sliding across the genome with step size of 50 variants at a time, and filter out any SNPs with LD higher than 0.25
# out	EUR.QC	Informs plink that all output should have a prefix of EUR.QC

# This will generate two files 1) EUR.QC.prune.in and 2) EUR.QC.prune.out. All SNPs within EUR.QC.prune.in have a pairwise r^2<0.25


# Compute Heterozygosity rates #

# This will generate the EUR.QC.het file, which contains F coefficient estimates for assessing heterozygosity. 
# We will remove individuals with F coefficients that are more than 3 standard deviation (SD) units from the mean
# F Inbreeding Coefficient?

plink \
    --bfile EUR \
    --extract EUR.QC.prune.in \
    --keep EUR.QC.fam \
    --het \
    --out EUR.QC

#!/usr/bin/env Rscript
#or
#!/usr/bin/Rscript
dat <- read.table("EUR.QC.het", header=T) # Read in the EUR.het file, specify it has header
m <- mean(dat$F) # Calculate the mean  
s <- sd(dat$F) # Calculate the SD
valid <- subset(dat, F <= m+3*s & F >= m-3*s) # Get any samples with F coefficient within 3 SD of the population mean
write.table(valid[,c(1,2)], "EUR.valid.sample", quote=F, row.names=F) # print FID and IID for valid samples
q() # exit R


# Remove Ambiguous SNPs #
# already done?

