#!/bin/bash
#Script to run entire analysis at once

#This script won't work if you have less than 5% of space available in your server.

display_usage() {
echo '
1st argument is the path to the vcf file 
2nd argument is the complete path to the directory where results must be saved
'
}


#check if required arguments are there and display usage message
if [ -z "$1" ] || [ -z "$2" ] ; then
        printf "Please provide the arguments required for the script.\n\n"
        display_usage
        exit 1
fi


#variantes
IN="$1"
OUT="$2"

################################################
##Identify prune sites (LD sites and Plink QC)##
################################################

#change values according to what's more adequate given the VCF statistics
geno=0.05	#Missingness per SNP
mind=0.1	#Missingness per individual
maf=0.01	#Minor allele frequency


cd $OUT

touch temp-cases.txt #create empty file

# --allow-extra-chr: avoid Error: Invalid chromosome code 'chrUn_...' 
# --make-pheno cases_ID.txt '*' & --allow-no-sex: avoid Error: "all individuals set ignore.  Likely input problem (col 6)" when converting to EIGENSTRAT format
# --keep-allele-order: force the original A1/A2 allele encoding to be preserved instead of the PLINK encoding
# --set-missing-var-ids @:#: gives an ID to each variant in the format chr:position
# --pca var-wts: asks for variants weight.
# --nonfounders: to include animals with unknown parents. the --maf option is not considered in founder animals, i.e. those with unknown parents

#Filter missingness

plink --vcf $IN --dog --allow-extra-chr --autosome --double-id --make-pheno temp-cases.txt '*' --allow-no-sex --nonfounders --keep-allele-order --set-missing-var-ids @:# --geno $geno --make-bed --vcf-half-call missing --out dog.cleaned.geno"$geno" 

plink --bfile dog.cleaned.geno"$geno" --dog --allow-extra-chr --autosome --double-id --make-pheno temp-cases.txt '*' --allow-no-sex --nonfounders --keep-allele-order --mind $mind --make-bed --out dog.cleaned.geno"$geno".mind"$mind" 

plink --bfile dog.cleaned.geno"$geno".mind"$mind" --dog --allow-extra-chr --autosome --double-id --make-pheno temp-cases.txt '*' --allow-no-sex --nonfounders --keep-allele-order --maf $maf --make-bed --out dog.cleaned.geno"$geno".mind"$mind".maf"$maf" 


# Perform Relatedness. 

plink2 --dog --threads 10 --memory 300000 --bfile dog.cleaned.geno"$geno".mind"$mind".maf"$maf" --make-king-table --out dog.kinship.table 

plink2 --dog --threads 10 --memory 300000 --bfile dog.cleaned.geno"$geno".mind"$mind".maf"$maf" --king-cutoff .177 --make-bed --out dog.kingship_clean 

plink --bfile dog.kingship_clean --keep-allele-order --recode --make-bed --dog --out dog.kingship_clean_recode


rm temp-cases.txt #remove temporary file

# Create missing genotype statistic individual and loci after filtering
plink --bfile dog.cleaned.geno"$geno".mind"$mind" --dog --keep-allele-order --allow-extra-chr --autosome --missing --out missing_genotype








