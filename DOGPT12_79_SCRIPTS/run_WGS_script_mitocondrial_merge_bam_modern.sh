#!/bin/bash
#Script to create consensus file

display_usage() { 
echo '
1st argument is the path to the first bam file to merge
2nd argument is the path to the second bam file to merge
3rd argument is the number of threads available to use. example "15"
4th argument is the sample name 
5th argument is the complete path to the directory were results must be saved 
6th argument is the path to the reference genome index
7th argument is the minimum snp coverage 
8th argument is the minimum snp quality
9th argument is the name of the snpEff database (ex: canis_mt or canfam3.1)'
}

#check if required arguments are there and display usage message
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ][ -z "$4" ] || [ -z "$5" ] || [ -z "$6" ] || [ -z "$7" ] || [ -z "$8" ] || [ -z "$9" ]; then
        printf "Please provide the arguments required for the script.\n\n"
        display_usage
        exit 1
fi

#variantes 

bam1="$1"
bam2="$2"
threads="$3"
base_name="$4" 
base_output="$5" 
reference="$6" 
snp_coverage="$7"
snp_quality="$8"
snpEffDB="$9"

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


gatk HaplotypeCaller -R "$reference" --sample-ploidy 1 -I "$base_name"_sorted.bam --bam-output "$base_name"_GATK_out.bam -O "$base_output"/variant_calling/"$base_name".vcf 	#--sample-ploidy 1 para dizer ao software que é um genoma haploide e quando houver posições heterozigoticas, ficar com o alelo mais frequente.

#vcf filter using GATK

cd $base_output'/variant_calling'

gatk VariantFiltration -R "$reference" -V "$base_output"/variant_calling/"$base_name".vcf --filter-name 'FAILED_qual' --filter-expression "QD < $snp_quality" --genotype-filter-name 'FAILED_DP' --genotype-filter-expression "DP < $snp_coverage" --set-filtered-genotype-to-no-call true -O "$base_output"/variant_calling/"$base_name"_filtered.vcf

gatk SelectVariants -V "$base_output"/variant_calling/"$base_name"_filtered.vcf --select-type-to-exclude INDEL -O "$base_output"/variant_calling/"$base_name"_filtered_NO_INDELS.vcf

gatk SelectVariants -V "$base_output"/variant_calling/"$base_name"_filtered.vcf --exclude-filtered true --exclude-non-variants true --remove-unused-alternates true --restrict-alleles-to BIALLELIC -O "$base_output"/variant_calling/"$base_name"_PASS_ONLY.vcf

gatk SelectVariants -V "$base_output"/variant_calling/"$base_name"_filtered.vcf --exclude-filtered true --exclude-non-variants true --remove-unused-alternates true --select-type-to-include INDEL --restrict-alleles-to BIALLELIC -O "$base_output"/variant_calling/"$base_name"_INDELS.vcf #para criar um ficheiro só de INDELS


########################
###### SNP effects #####
########################

java -jar /opt/anaconda3/share/snpeff-4.3.1t-1/snpEff.jar "$snpEffDB" "$base_output"/variant_calling/"$base_name"_PASS_ONLY.vcf > "$base_output"/variant_calling/"$base_name"_effects.vcf

##########################
####CONSENSUS SEQUENCE####
##########################

mkdir $base_output'/consensus_sequence'


bedtools genomecov -ibam "$base_output"/mapping/"$base_name"_sorted.bam -bga > "$base_output"/consensus_sequence/"$base_name"_cov_regions.bed 

python /mnt/storage1/DATA/SCRIPTS/resolver_depth_nas_delecoes.py -v "$base_output"/variant_calling/"$base_name"_INDELS.vcf -b "$base_output"/consensus_sequence/"$base_name"_cov_regions.bed -o "$base_output"/consensus_sequence/"$base_name"_cov_regions_DELok.bed  #para resolver o problema gerado entre a diferença de coverage gerada pelo bedtools na região das deleções e o valor de coverage anotado no vcf para essa deleção

cat "$base_output"/consensus_sequence/"$base_name"_cov_regions_DELok.bed | awk '$4 < 10' > "$base_output"/consensus_sequence/"$base_name"_lowcov_regions.bed 

bedtools maskfasta -fi "$reference" -bed "$base_output"/consensus_sequence/"$base_name"_lowcov_regions.bed -fo "$base_output"/consensus_sequence/"$base_name"_MT_REFERENCE_with_Nscov.fa 

bgzip -c "$base_output"/variant_calling/"$base_name"_filtered_NO_INDELS.vcf > "$base_output"/variant_calling/"$base_name"_filtered_NO_INDELS.vcf.gz

tabix -p vcf "$base_output"/variant_calling/"$base_name"_filtered_NO_INDELS.vcf.gz

#################

zcat "$base_output"/variant_calling/"$base_name"_filtered_NO_INDELS.vcf.gz | grep "^#" > "$base_output"/variant_calling/"$base_name"_FAILED_ONLY.vcf # pega em todas as linhas que começam com "#", ou seja, o cabeçalho do ficheiro vcf

zcat "$base_output"/variant_calling/"$base_name"_filtered_NO_INDELS.vcf.gz | grep -v "^#" | grep "FAILED" >> "$base_output"/variant_calling/"$base_name"_FAILED_ONLY.vcf  #pega em todas as que não começam com "#", depois pega só nas que têm a palavra "FAILED" (ou seja. não passaram os filtros) e adiciona ao cabeçalho que tinha guardado no comando anterior

bedtools maskfasta -fi "$base_output"/consensus_sequence/"$base_name"_MT_REFERENCE_with_Nscov.fa -bed "$base_output"/variant_calling/"$base_name"_FAILED_ONLY.vcf -fo "$base_output"/consensus_sequence/"$base_name"_MT_REFERENCE_with_Ns.fa #aqui pego no fasta com Ns anterior (baixa cobertura) e substituo por N.

#################

bgzip -c "$base_output"/variant_calling/"$base_name"_PASS_ONLY.vcf > "$base_output"/variant_calling/"$base_name"_MT_PASS_ONLY.vcf.gz 

tabix -p vcf "$base_output"/variant_calling/"$base_name"_MT_PASS_ONLY.vcf.gz 

bcftools consensus -f "$base_output"/consensus_sequence/"$base_name"_MT_REFERENCE_with_Ns.fa -o "$base_output"/consensus_sequence/"$base_name"_MT_REFERENCE_with_Ns_and_SNPs.fa "$base_output"/variant_calling/"$base_name"_MT_PASS_ONLY.vcf.gz 
