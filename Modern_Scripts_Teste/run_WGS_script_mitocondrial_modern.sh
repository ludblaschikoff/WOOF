#!/bin/bash
#Script to run entire analysis at once

##Execute as bash workflow.sh "reference" "threads" "pathtooutput" "basename" "snp_coverage" "snp_quality" "snpEffDB" "biosample" 

display_usage() {
echo '
1th argument is the path to the reference genome index
2th argument is the number of threads available to use. example "15"
3th argument is the complete path to the directory were results must be saved
4th argument is the sample name
5th argument is the minimum snp coverage. example "10.0"
6th argument is the minimum snp quality. example "20.0"
7th argument is the name of the snpEff database (ex: canis_mt or canfam3.1)
8th argument is the biosample name
9th argument is is the library name
10th argument is the minimum mapping quality. example "30"
11th argument is the minimum base coverage. example "10"'
}


#check if required arguments are there and display usage message
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ -z "$5" ] || [ -z "$6" ] || [ -z "$7" ] || [ -z "$8" ] [ -z "$9" ] [ -z "$10" ] [ -z "$11" ]; then
        printf "Please provide the arguments required for the script.\n\n"
        display_usage
        exit 1
fi

#variantes 
reference="$1"
threads="$2"
base_output="$3"
base_name="$4"
snp_coverage="$5"
snp_quality="$6"
snpEffDB="$7"
biosample="$8"
library_name="$9"
minmapqual="${10}"
base_coverage="${11}"

#########################
######## MAPPING ########
#########################

mkdir $base_output'/mapping'
cd "$base_output"/mapping

clean_fastq_file=../../preprocessing/cutadapt/"$base_name"_R2_trimmed.fq.gz

if [ -f "$clean_fastq_file" ]; then
    pair='true'
else
    pair='false'
fi

#bwa mem - mapping against dog reference

if [ $pair = 'true' ]; then
	bwa mem -t "$threads" "$reference" ../../preprocessing/cutadapt/"$base_name"_R1_trimmed.fq.gz ../../preprocessing/cutadapt/"$base_name"_R2_trimmed.fq.gz -o "$base_name".sam

else
	bwa mem -t "$threads" "$reference" ../../preprocessing/cutadapt/"$base_name"_trimmed.fq.gz -o "$base_name".sam
fi

#Create a file of mapped and unmapped reads (samtools view)

samtools view -bh "$base_name".sam -o "$base_name"_dog_mappedandunmapped.bam
samtools sort "$base_output"/mapping/"$base_name"_dog_mappedandunmapped.bam > "$base_output"/mapping/"$base_name"_dog_mappedandunmapped_sorted.bam
samtools index "$base_output"/mapping/"$base_name"_dog_mappedandunmapped_sorted.bam

rm $base_name'_dog_mappedandunmapped.bam'

#samtools (filtering by flag and qual)
samtools view -bh -f2 -F 4 -F 1024 -F 2048 -F 256 -q "$minmapqual" "$base_name".sam -o "$base_name".bam # -h guardar com head cabeçário e -b para guardar ficheiro .bam
samtools sort "$base_name".bam > "$base_name"_sorted.bam
samtools index "$base_name"_sorted.bam

#picartools (add read groups)
java -XX:ParallelGCThreads="$threads" -XX:ConcGCThreads="$threads" -jar $PICARD AddOrReplaceReadGroups VALIDATION_STRINGENCY="LENIENT" ID="$biosample" SM="$biosample" PU="PU" LB="$library_name" PL="illumina" I="$base_name"_sorted.bam O="$base_name"_RG.bam #add uma coluna de identificacao no .bam com identificador para depois pordermos remover os duplicados
samtools sort "$base_name"_RG.bam > "$base_name"_sorted_RG.bam
samtools index "$base_name"_sorted_RG.bam

#picartools (remove duplicates); reuse samtools index and samtools sort because we are not sure if the Picard desorganized the reads. 

java -XX:ParallelGCThreads="$threads" -XX:ConcGCThreads="$threads" -jar $PICARD MarkDuplicates VALIDATION_STRINGENCY="LENIENT" REMOVE_DUPLICATES="true" I="$base_name"_sorted_RG.bam O="$base_name"_no_dups.bam M="marked_dup_metrics.txt"
samtools sort "$base_name"_no_dups.bam > "$base_name"_sorted_no_dups.bam
samtools index "$base_name"_sorted_no_dups.bam

rm $base_name'.sam' $base_name'.bam' $base_name'_RG.bam' $base_name'_sorted_RG.bam' $base_name'_sorted_RG.bam.bai' $base_name'_no_dups.bam'   

qualimap bamqc -nt "$threads" -c -bam "$base_name"_sorted_no_dups.bam --java-mem-size=32G
 
#########################
#### VARIANT CALLING ####
#########################

mkdir $base_output'/variant_calling'


gatk HaplotypeCaller -R "$reference" --sample-ploidy 1 -I "$base_name"_sorted_no_dups.bam --bam-output "$base_name"_GATK_out.bam -O "$base_output"/variant_calling/"$base_name".vcf 	#--sample-ploidy 1 para dizer ao software que é um genoma haploide e quando houver posições heterozigoticas, ficar com o alelo mais frequente.

#vcf filter using GATK

cd $base_output'/variant_calling'

gatk VariantFiltration -R "$reference" -V "$base_output"/variant_calling/"$base_name".vcf --filter-name 'FAILED_qual' --filter-expression "QD < $snp_quality" --genotype-filter-name 'FAILED_DP' --genotype-filter-expression "DP < $snp_coverage" --set-filtered-genotype-to-no-call true -O "$base_output"/variant_calling/"$base_name"_filtered.vcf

#######

gatk SelectVariants -V "$base_output"/variant_calling/"$base_name"_filtered.vcf --select-type-to-exclude INDEL -O "$base_output"/variant_calling/"$base_name"_filtered_NO_INDELS.vcf	#para criar um ficheiro filtrado de SNPs e sem indels

bgzip -c "$base_output"/variant_calling/"$base_name"_filtered_NO_INDELS.vcf > "$base_output"/variant_calling/"$base_name"_filtered_NO_INDELS.vcf.gz

tabix -p vcf "$base_output"/variant_calling/"$base_name"_filtered_NO_INDELS.vcf.gz

#########

gatk SelectVariants -V "$base_output"/variant_calling/"$base_name"_filtered.vcf --exclude-filtered true --exclude-non-variants true --remove-unused-alternates true --restrict-alleles-to BIALLELIC -O "$base_output"/variant_calling/"$base_name"_PASS_ONLY_INDELS.vcf	#ficheiro SÓ das variantes que passaram, com INDELS

bgzip -c "$base_output"/variant_calling/"$base_name"_PASS_ONLY_INDELS.vcf > "$base_output"/variant_calling/"$base_name"_MT_PASS_ONLY_INDELS.vcf.gz 

tabix -p vcf "$base_output"/variant_calling/"$base_name"_MT_PASS_ONLY_INDELS.vcf.gz 

#########

gatk SelectVariants -V "$base_output"/variant_calling/"$base_name"_filtered.vcf --exclude-filtered true --select-type-to-exclude INDEL --exclude-non-variants true --remove-unused-alternates true --restrict-alleles-to BIALLELIC -O "$base_output"/variant_calling/"$base_name"_PASS_ONLY_noINDELS.vcf	#ficheiro SÓ das variantes que passaram, sem INDELS

bgzip -c "$base_output"/variant_calling/"$base_name"_PASS_ONLY_noINDELS.vcf > "$base_output"/variant_calling/"$base_name"_MT_PASS_ONLY_noINDELS.vcf.gz 

tabix -p vcf "$base_output"/variant_calling/"$base_name"_MT_PASS_ONLY_noINDELS.vcf.gz 

########

gatk SelectVariants -V "$base_output"/variant_calling/"$base_name"_filtered.vcf --exclude-filtered true --exclude-non-variants true --remove-unused-alternates true --select-type-to-include INDEL --restrict-alleles-to BIALLELIC -O "$base_output"/variant_calling/"$base_name"_INDELS.vcf #para criar um ficheiro SÓ de INDELS


######################################
####CONSENSUS SEQUENCE WITH INDELS####
######################################

mkdir $base_output'/consensus_sequence'
mkdir $base_output'/consensus_sequence/consensus_INDELS'
cd  $base_output'/consensus_sequence/consensus_INDELS'

bedtools genomecov -ibam "$base_output"/mapping/"$base_name"_sorted_no_dups.bam -bga > "$base_output"/consensus_sequence/"$base_name"_cov_regions.bed 

python /mnt/storage1/DATA/SCRIPTS/resolver_depth_nas_delecoes.py -v "$base_output"/variant_calling/"$base_name"_INDELS.vcf -b "$base_output"/consensus_sequence/"$base_name"_cov_regions.bed -o "$base_name"_cov_regions_DELok.bed  #para resolver o problema gerado entre a diferença de coverage gerada pelo bedtools na região das deleções e o valor de coverage anotado no vcf para essa deleção. Assim crio um novo ficheiro bed de cobertura corrigida

cat "$base_name"_cov_regions_DELok.bed | awk '$4 < "$base_coverage"' > "$base_name"_lowcov_regions.bed 

bedtools maskfasta -fi "$reference" -bed "$base_name"_lowcov_regions.bed -fo "$base_name"_MT_REFERENCE_with_Nscov.fa 


#################

zcat "$base_output"/variant_calling/"$base_name"_filtered_NO_INDELS.vcf.gz | grep "^#" > "$base_output"/variant_calling/"$base_name"_FAILED_ONLY.vcf # pega em todas as linhas que começam com "#", ou seja, o cabeçalho do ficheiro vcf

zcat "$base_output"/variant_calling/"$base_name"_filtered_NO_INDELS.vcf.gz | grep -v "^#" | grep "FAILED" >> "$base_output"/variant_calling/"$base_name"_FAILED_ONLY.vcf  #pega em todas as que não começam com "#", depois pega só nas que têm a palavra "FAILED" (ou seja, não passaram nos filtros) e adiciona ao cabeçalho que tinha guardado no comando anterior. Assim crio um ficheiro só de SNPs que não passaram na filtragem

bedtools maskfasta -fi "$base_name"_MT_REFERENCE_with_Nscov.fa -bed "$base_output"/variant_calling/"$base_name"_FAILED_ONLY.vcf -fo "$base_name"_MT_REFERENCE_with_Ns.fa #aqui pego no fasta com Ns anterior (baixa cobertura) e substituo posições FAILED por N.

#################

bcftools consensus -f "$base_name"_MT_REFERENCE_with_Ns.fa -o "$base_name"_MT_REFERENCE_with_Ns_SNPs_INDELS.fa "$base_output"/variant_calling/"$base_name"_MT_PASS_ONLY.vcf.gz 

######################################
#####CONSENSUS SEQUENCE noINDELS######
######################################

mkdir $base_output'/consensus_sequence/consensus_NOINDELS'
cd  $base_output'/consensus_sequence/consensus_NOINDELS'

cat "$base_output"/consensus_sequence/"$base_name"_cov_regions.bed  | awk '$4 < "$base_coverage"' > "$base_name"_lowcov_regions.bed

bedtools maskfasta -fi "$reference" -bed "$base_name"_lowcov_regions.bed -fo "$base_name"_MT_REFERENCE_with_Nscov.fa 

#####

zcat "$base_output"/variant_calling/"$base_name"_filtered_NO_INDELS.vcf.gz | grep "^#" > "$base_output"/variant_calling/"$base_name"_FAILED_ONLY.vcf # pega em todas as linhas que começam com "#", ou seja, o cabeçalho do ficheiro vcf

zcat "$base_output"/variant_calling/"$base_name"_filtered_NO_INDELS.vcf.gz | grep -v "^#" | grep "FAILED" >> "$base_output"/variant_calling/"$base_name"_FAILED_ONLY.vcf  #pega em todas as que não começam com "#", depois pega só nas que têm a palavra "FAILED" (ou seja, não passaram nos filtros) e adiciona ao cabeçalho que tinha guardado no comando anterior. Assim crio um ficheiro só de SNPs que não passaram na filtragem

bedtools maskfasta -fi "$base_name"_MT_REFERENCE_with_Nscov.fa -bed "$base_output"/variant_calling/"$base_name"_FAILED_ONLY.vcf -fo "$base_name"_MT_REFERENCE_with_Ns.fa #aqui pego no fasta com Ns anterior (baixa cobertura) e substituo posições FAILED por N.

######

bcftools consensus -f "$base_name"_MT_REFERENCE_with_Ns.fa -o "$base_name"_MT_REFERENCE_with_Ns_SNPs.fa "$base_output"/variant_calling/"$base_name"_MT_PASS_ONLY_noINDELS.vcf.gz 

