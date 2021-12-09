#!/bin/bash
#Script to run entire aDNA analysis at once

display_usage() {
echo '
1th argument is the path to the first bam file to merge
2th argument is the path to the second bam file to merge
3th argument is the path to the reference genome index
4th argument is the number of threads available to use. example "15"
5th argument is the complete path to the directory were results must be saved
6th argument is the sample name
7th argument is the minimum snp coverage. example 7.0
8th argument is the minimum snp quality. example 30.0
9th argument is the minimum genotype quality. example 20.0
10th argument is the name of the snpEff database (ex: canis_mt or CanFam3.1.95)'
}


#check if required arguments are there and display usage message
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ -z "$5" ] || [ -z "$6" ] || [ -z "$7" ] || [ -z "$8" ] || [ -z "$9" ] || [ -z "$10" ] ; then
        printf "Please provide the arguments required for the script.\n\n"
        display_usage
        exit 1
fi


#variantes 
bam1="$1"
bam2="$2"
reference="$3"
threads="$4"
base_output="$5"
base_name="$6"
snp_coverage="$7"
snp_quality="$8"
genotype_quality="$9"
snpEffDB="${10}"



#########################
######### MERGE #########
#########################

cd $base_output

samtools merge "$base_name".bam "$bam1" "$bam2"

samtools sort "$base_name".bam > "$base_name"_sorted.bam
samtools index "$base_name"_sorted.bam

qualimap bamqc -nt "$threads" -c -bam "$base_name"_sorted.bam --java-mem-size=32G



#########################
#### VARIANT CALLING ####
#########################

mkdir $base_output'/variant_calling'

gatk HaplotypeCaller -R "$reference" -I "$base_name"_sorted.bam --bam-output "$base_name"_GATK_out.bam -O "$base_output"/variant_calling/"$base_name".vcf

#vcf filter using GATK

cd $base_output'/variant_calling'
gatk VariantFiltration -R "$reference" -V "$base_output"/variant_calling/"$base_name".vcf --filter-name 'FAILED_qual' --filter-expression "QD < $snp_quality" --genotype-filter-name 'FAILED_DP' --genotype-filter-expression "DP < $snp_coverage" --genotype-filter-name "FAILED_GQ" --genotype-filter-expression "GQ < $genotype_quality" --set-filtered-genotype-to-no-call true -O "$base_output"/variant_calling/"$base_name"_filtered.vcf

gatk SelectVariants -V "$base_output"/variant_calling/"$base_name"_filtered.vcf --exclude-filtered true --select-type-to-exclude INDEL --exclude-non-variants true --remove-unused-alternates true --restrict-alleles-to BIALLELIC -O "$base_output"/variant_calling/"$base_name"_PASS_ONLY.vcf

########################
###### SNP effects #####
########################

java -jar /opt/anaconda3/share/snpeff-4.3.1t-1/snpEff.jar "$snpEffDB" "$base_output"/variant_calling/"$base_name"_PASS_ONLY.vcf > "$base_output"/variant_calling/"$base_name"_effects.vcf