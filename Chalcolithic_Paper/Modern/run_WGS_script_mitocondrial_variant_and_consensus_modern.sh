#!/bin/bash
#Script to create consensus file

display_usage() { 
echo '
1st argument is the sample name 
2nd argument is the complete path to the directory were results must be saved 
3rd argument is the path to the reference genome index
4th argument is the minimum snp coverage 
5th argument is the minimum snp quality
6th argument is the name of the snpEff database (ex: canis_mt or canfam3.1)'
}

#check if required arguments are there and display usage message
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ][ -z "$4" ] || [ -z "$5" ] || [ -z "$6" ] ; then
        printf "Please provide the arguments required for the script.\n\n"
        display_usage
        exit 1
fi

#variantes 

base_name="$1" 
base_output="$2" 
reference="$3" 
snp_coverage="$4"
snp_quality="$5"
snpEffDB="$6"


#########################
#### VARIANT CALLING ####
#########################

mkdir $base_output'/variant_calling'


gatk HaplotypeCaller -R "$reference" --sample-ploidy 1 -I "$base_name"_no_dups.bam --bam-output "$base_name"_GATK_out.bam -O "$base_output"/variant_calling/"$base_name".vcf 	#--sample-ploidy 1 tell the software that it is a haploid genome and when there are heterozygous positions, keep with the most frequent allele. 

#vcf filter using GATK

cd $base_output'/variant_calling'

gatk VariantFiltration -R "$reference" -V "$base_output"/variant_calling/"$base_name".vcf --filter-name 'FAILED_qual' --filter-expression "QD < $snp_quality" --genotype-filter-name 'FAILED_DP' --genotype-filter-expression "DP < $snp_coverage" --set-filtered-genotype-to-no-call true -O "$base_output"/variant_calling/"$base_name"_filtered.vcf

gatk SelectVariants -V "$base_output"/variant_calling/"$base_name"_filtered.vcf --select-type-to-exclude INDEL -O "$base_output"/variant_calling/"$base_name"_filtered_NO_INDELS.vcf

gatk SelectVariants -V "$base_output"/variant_calling/"$base_name"_filtered.vcf --exclude-filtered true --exclude-non-variants true --remove-unused-alternates true --restrict-alleles-to BIALLELIC -O "$base_output"/variant_calling/"$base_name"_PASS_ONLY.vcf

gatk SelectVariants -V "$base_output"/variant_calling/"$base_name"_filtered.vcf --exclude-filtered true --exclude-non-variants true --remove-unused-alternates true --select-type-to-include INDEL --restrict-alleles-to BIALLELIC -O "$base_output"/variant_calling/"$base_name"_INDELS.vcf 


########################
###### SNP effects #####
########################

java -jar /opt/anaconda3/share/snpeff-4.3.1t-1/snpEff.jar "$snpEffDB" "$base_output"/variant_calling/"$base_name"_PASS_ONLY.vcf > "$base_output"/variant_calling/"$base_name"_effects.vcf

##########################
####CONSENSUS SEQUENCE####
##########################

mkdir $base_output'/consensus_sequence'


bedtools genomecov -ibam "$base_output"/mapping/"$base_name"_no_dups.bam -bga > "$base_output"/consensus_sequence/"$base_name"_cov_regions.bed 

python /mnt/storage1/DATA/SCRIPTS/resolver_depth_nas_delecoes.py -v "$base_output"/variant_calling/"$base_name"_INDELS.vcf -b "$base_output"/consensus_sequence/"$base_name"_cov_regions.bed -o "$base_output"/consensus_sequence/"$base_name"_cov_regions_DELok.bed  #to solve the problem generated between the coverage difference generated by bedtools in the region of the deletions and the coverage value noted in the vcf for this deletion.

cat "$base_output"/consensus_sequence/"$base_name"_cov_regions_DELok.bed | awk '$4 < 10' > "$base_output"/consensus_sequence/"$base_name"_lowcov_regions.bed 

bedtools maskfasta -fi "$reference" -bed "$base_output"/consensus_sequence/"$base_name"_lowcov_regions.bed -fo "$base_output"/consensus_sequence/"$base_name"_MT_REFERENCE_with_Nscov.fa 

bgzip -c "$base_output"/variant_calling/"$base_name"_filtered_NO_INDELS.vcf > "$base_output"/variant_calling/"$base_name"_filtered_NO_INDELS.vcf.gz

tabix -p vcf "$base_output"/variant_calling/"$base_name"_filtered_NO_INDELS.vcf.gz

#################

zcat "$base_output"/variant_calling/"$base_name"_filtered_NO_INDELS.vcf.gz | grep "^#" > "$base_output"/variant_calling/"$base_name"_FAILED_ONLY.vcf #takes all lines starting with "#", i.e. header of vcf file.

zcat "$base_output"/variant_calling/"$base_name"_filtered_NO_INDELS.vcf.gz | grep -v "^#" | grep "FAILED" >> "$base_output"/variant_calling/"$base_name"_FAILED_ONLY.vcf   #takes all those that don't start with "#", then takes only those with the word "FAILED" (i.e. they didn't pass the filters) and adds it to the header that it had saved in the previous command

bedtools maskfasta -fi "$base_output"/consensus_sequence/"$base_name"_MT_REFERENCE_with_Nscov.fa -bed "$base_output"/variant_calling/"$base_name"_FAILED_ONLY.vcf -fo "$base_output"/consensus_sequence/"$base_name"_MT_REFERENCE_with_Ns.fa #take the fasta with previous Ns (low coverage) and replace it with N.

#################

bgzip -c "$base_output"/variant_calling/"$base_name"_PASS_ONLY.vcf > "$base_output"/variant_calling/"$base_name"_MT_PASS_ONLY.vcf.gz 

tabix -p vcf "$base_output"/variant_calling/"$base_name"_MT_PASS_ONLY.vcf.gz 

bcftools consensus -f "$base_output"/consensus_sequence/"$base_name"_MT_REFERENCE_with_Ns.fa -o "$base_output"/consensus_sequence/"$base_name"_MT_REFERENCE_with_Ns_and_SNPs.fa "$base_output"/variant_calling/"$base_name"_MT_PASS_ONLY.vcf.gz 
