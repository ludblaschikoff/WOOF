#!/bin/bash
#Script to run entire analysis at once

display_usage() {
echo '
1st argument is the path to the reference genome index
2nd argument is the complete path for the directory where bam files are located. usually ends with "2.mapping"
3rd argument is the complete path to the directory where results must be saved
4th argument is the minimum per read SNP quality to filter vcf variants
5th argument is the minimum depth of coverage to filter vcf variants
6th argument is the minimum genotype quality. example 20.0
7th argument is the snpEff database name (ex: canis_mt or CanFam3.1.95)
8th argument is the max RAM available for this job, in Gb (ex: 200)
'
}


#check if required arguments are there and display usage message
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ -z "$5" ] || [ -z "$6" ] || [ -z "$7" ] || [ -z "$8" ]; then
        printf "Please provide the arguments required for the script.\n\n"
        display_usage
        exit 1
fi


#variantes
reference="$1"
bam_files_dir="$2"
base_output="$3"
snp_quality="$4"
snp_coverage="$5"
genotype_quality="$6"
snpEffDB="$7"
maxRAM="$8"				

mkdir $base_output'/variant_calling'
cd $base_output'/variant_calling'

#CALCULAMOS A RAM DISPONIVEL PARA CADA RUN

threads=30

ind_ram=$(( $maxRAM/$threads ))		#run é o numero de amostras a serem processadas. é o mesmo que threads.

ind_ram=$ind_ram'g'					##atualizo o valor da variavel $ind_ram para $ind_ram + g para o nome da variável levar a letra g.

#### 1.VARIANT CALLING ####
#LOOP PARA FAZER GVCFs

ls $bam_files_dir/*.bam | xargs -n 1 -P $threads -I {} sh -c 'base_name="$(basename {} | cut -d_ -f2)"; echo "$base_name"; gatk --java-options "-Xmx'$ind_ram'" HaplotypeCaller -R '$reference' -I {} -O '"$base_output"'/variant_calling/"$base_name".g.vcf.gz -ERC GVCF' 		 # -P $runs é o numero de processadores que distribuo para cada amostra; #xargs (short for "eXtended ARGuments") serve para correr comandos em paralelo. is a command on Unix and most Unix-like operating systems used to build and execute commands from standard input. It converts input from standard input (ficheiros) into arguments to a command.
																																																																									#{} o input vem de trás.

#### 2.DATA AGGREGATION STEP ####

gatk --java-options "-Xmx"$maxRAM"g -Xms"$maxRAM"g" GenomicsDBImport --genomicsdb-workspace-path "$base_output"/variant_calling/genomicsDBworkspace/ --intervals "$base_output"/variant_calling/CanFam3.1_REF_Chr_names.list --sample-name-map "$base_output"/variant_calling/sample_names_map.tab --tmp-dir=/dev/shm --max-num-intervals-to-import-in-parallel 5			#--sample-name-map Provide sample GVCFs in a map file. It must contain sample name + tab separated + gvcf location of the sample		#Do not create the folder genomicsDBworkspace before, let the comand creates it. Otherwise it will produce an ERROR.
																																																																																																																			#--batch size controls the number of samples for which readers are open at once and therefore provides a way to minimize memory consumption.

#### 3.JOINT GENOTYPING ####

gatk --java-options "-Xmx"$maxRAM"g" GenotypeGVCFs -R $reference -V gendb://"$base_output"/variant_calling/genomicsDBworkspace/ -O "$base_output"/variant_calling/raw_variants.vcf.gz			#Perform joint genotyping on one or more samples pre-called with HaplotypeCaller ##gendb:// prefix to the database input directory path (created before); a GenomicsDB workspace created by GenomicsDBImport ; a final VCF in which all samples have been jointly genotyped

#### 4.VARIANT filtration ####

gatk VariantFiltration -R "$reference" -V "$base_output"/variant_calling/raw_variants.vcf.gz --filter-name 'FAILED_qual' --filter-expression "QD < $snp_quality" --genotype-filter-name 'FAILED_DP' --genotype-filter-expression "DP < $snp_coverage" --genotype-filter-name "FAILED_GQ" --genotype-filter-expression "GQ < $genotype_quality" --set-filtered-genotype-to-no-call true -O "$base_output"/variant_calling/All_variants_filter_mark.vcf   

#### 5. VARIANT SELECTION ####

gatk SelectVariants -V "$base_output"/variant_calling/All_variants_filter_mark.vcf --exclude-filtered true --restrict-alleles-to BIALLELIC --exclude-non-variants true --remove-unused-alternates true -O "$base_output"/variant_calling/All_variants_filtered_DP"$snp_coverage"QD"$snp_quality"GQ"$genotype_quality"_Biallelic.vcf						 

gatk SelectVariants -V "$base_output"/variant_calling/All_variants_filter_mark.vcf --exclude-filtered true --restrict-alleles-to BIALLELIC --exclude-non-variants true --select-type-to-include SNP --remove-unused-alternates true -O "$base_output"/variant_calling/SNP_variants_filtered_DP"$snp_coverage"QD"$snp_quality"GQ"$genotype_quality"_Biallelic.vcf		## Devemos fazer um vcf só para INDEL e outro só para SNPs

gatk SelectVariants -V "$base_output"/variant_calling/All_variants_filter_mark.vcf --exclude-filtered true --restrict-alleles-to BIALLELIC --exclude-non-variants true --select-type-to-include INDEL --remove-unused-alternates true -O "$base_output"/variant_calling/INDEL_variants_filtered_DP"$snp_coverage"QD"$snp_quality"GQ"$genotype_quality"_Biallelic.vcf			

###########################
######snpEFF###############
###########################
java -Xmx"$maxRAM"g -jar /opt/anaconda3/share/snpeff-4.3.1t-1/snpEff.jar "$snpEffDB" "$base_output"/variant_calling/All_variants_filtered_DP"$snp_coverage"QD"$snp_quality"GQ"$genotype_quality"_Biallelic.vcf > "$base_output"/variant_calling/All_variants_effects_DP"$snp_coverage"QD"$snp_quality"GQ"$genotype_quality"_Biallelic.vcf

