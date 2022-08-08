############### AUTHOR: CINDY SARMENTO ######################

#!/usr/bin/env bash

### First argument: vcf file we want to analyse
### Second argument: prefix for output files
### USage: bash script.sh file.vcf.gz out_filename

VCF=$1 
OUT=$2

set -eux -o pipefail

#Calculate allele frequency
vcftools --gzvcf $VCF --freq2 --out $OUT --max-alleles 2

#Calculate mean depth per individual
vcftools --gzvcf $VCF --depth --out $OUT

#Calculate mean depth per site
vcftools --gzvcf $VCF --site-mean-depth --out $OUT

#Calculate site quality
vcftools --gzvcf $VCF --site-quality --out $OUT

#Calculate proportion of missing data per individual
vcftools --gzvcf $VCF --missing-indv --out $OUT

#Calculate proportion of missing data per site
vcftools --gzvcf $VCF --missing-site --out $OUT

#Calculate heterozygosity and inbreeding coefficient per individual
vcftools --gzvcf $VCF --het --out $OUT

#Reports  a p-value for each site from a Hardy-Weinberg Equilibrium test
vcftools --gzvcf $VCF --hardy --out $OUT

##Distribution of site quality, mapping quality and depth of coverage
bcftools query $VCF -f "%QUAL\t%MQ\t%INFO/DP\n" -o "$OUT"_QUAL_MQ_DP.txt

##Distribution of genotypes quality
bcftools query $VCF -f "%FORMAT/GQ\n" -o "$OUT"_GQ.txt
