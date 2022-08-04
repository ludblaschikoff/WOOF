#!/bin/bash
#Script to run entire aDNA analysis at once

display_usage() {
echo '
1th argument is the path to the reference genome index
2rd argument is the minimum read length. example "30"
3nd argument is the mininum base quality for trimming. example "30"
4th argument is the minimum mapping quality. example "20"
5th argument is the path to the human genome index
6th argument is the number of threads available to use. example "15"
7th argument is the complete path to the directory were results must be saved
8th argument is the sample name
9th argument is the minimum snp coverage. example "3.0"
10th argument is the minimum snp quality. example "20.0"
11th argument is the name of the snpEff database (ex: canis_mt or CanFam3.1.95)
12th argument is the maximum read size (ex: 150bp)
13th argument is the path to the pig genome index
14th argument is the path to the chicken genome index
15th argument is the path to the cow genome index
16th argument is the path to the sheep genome index
17th argument is the path to the goat genome index
18th argument is the biosample name'
}


#check if required arguments are there and display usage message
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ -z "$5" ] || [ -z "$6" ] || [ -z "$7" ] || [ -z "$8" ] || [ -z "$9" ] || [ -z "$10" ] || [ -z "$11" ] || [ -z "$12" ] || [ -z "$13" ] || [ -z "$14" ] || [ -z "$15" ] || [ -z "$16" ] || [ -z "$17" ] || [ -z "$18" ]; then
        printf "Please provide the arguments required for the script.\n\n"
        display_usage
        exit 1
fi


#variantes 
reference="$1"
minreadlenght="$2"
trimqual="$3"
minmapqual="$4"
reference_human="$5"
threads="$6"
base_output="$7"
base_name="$8"
snp_coverage="$9"
snp_quality="${10}"
snpEffDB="${11}"
readsize="${12}"
reference_pig="${13}"
reference_chicken="${14}"
reference_cow="${15}"
reference_sheep="${16}"
reference_goat="${17}"
biosample="${18}"

#########################
######## MAPPING ########
#########################

cd "$base_output"../preprocessing

clean_fastq_file=cutadapt/"$base_name"_R2_trimmed.fq.gz

if [ -f "$clean_fastq_file" ]; then
    pair='true'
else
    pair='false'
fi


#Merge clean reads
if [ $pair = 'true' ]; then
        mkdir $base_output'../preprocessing/merge_reads'

        AdapterRemoval --file1 "$base_output"../preprocessing/cutadapt/"$base_name"_R1_trimmed.fq.gz --file2 "$base_output"../preprocessing/cutadapt/"$base_name"_R2_trimmed.fq.gz --basename "$base_output"../preprocessing/merge_reads/"$base_name"_output_paired --collapse --minlength "$minreadlenght" --minalignmentlength 11 --mm 1 --minquality "$trimqual" --trimqualities --trimns --maxns 10 --gzip   #adapter1 and adapter 2 are universal adapters sequences for Illumina platforms
fi


mkdir $base_output'/mapping'
cd "$base_output"/mapping

#bwa aln - mapping against dog reference

if [ $pair = 'true' ]; then
	bwa aln -t "$threads" -l 1024 -o 2 -n 0.01 "$reference" "$base_output"../preprocessing/merge_reads/"$base_name"_output_paired.collapsed.gz > "$base_output"/mapping/"$base_name"_dog.sai
	bwa samse "$reference" "$base_output"/mapping/"$base_name"_dog.sai "$base_output"../preprocessing/merge_reads/"$base_name"_output_paired.collapsed.gz > "$base_output"/mapping/"$base_name"_dog.sam
else
	bwa aln -t "$threads" -l 1024 -o 2 -n 0.01 "$reference" "$base_output"../preprocessing/cutadapt/"$base_name"_trimmed.fq.gz > "$base_output"/mapping/"$base_name"_dog.sai
	bwa samse "$reference" "$base_output"/mapping/"$base_name"_dog.sai "$base_output"../preprocessing/cutadapt/"$base_name"_trimmed.fq.gz > "$base_output"/mapping/"$base_name"_dog.sam
fi


#Create a file of mapped and unmapped reads (samtools view)

samtools view -bh "$base_output"/mapping/"$base_name"_dog.sam -o "$base_output"/mapping/"$base_name"_dog_mappedandunmapped.bam
samtools sort "$base_output"/mapping/"$base_name"_dog_mappedandunmapped.bam > "$base_output"/mapping/"$base_name"_dog_mappedandunmapped_sorted.bam
samtools index "$base_output"/mapping/"$base_name"_dog_mappedandunmapped_sorted.bam

#Filter only mapped reads (samtools view)
samtools view -bh -q "$minmapqual" "$base_output"/mapping/"$base_name"_dog.sam -o "$base_output"/mapping/"$base_name"_dog.bam 
samtools bam2fq "$base_output"/mapping/"$base_name"_dog.bam > "$base_output"/mapping/"$base_name"_dog.fq


#bwa aln - mapping against human reference

if [ $pair = 'true' ]; then
	bwa aln -t "$threads" -l 1024 -o 2 -n 0.01 "$reference_human" "$base_output"/mapping/"$base_name"_dog.fq > "$base_output"/mapping/"$base_name"_human.sai
	bwa samse "$reference_human" "$base_output"/mapping/"$base_name"_human.sai "$base_output"/mapping/"$base_name"_dog.fq > "$base_output"/mapping/"$base_name"_human.sam
else
	bwa aln -t "$threads" -l 1024 -o 2 -n 0.01 "$reference_human" "$base_output"/mapping/"$base_name"_dog.fq > "$base_output"/mapping/"$base_name"_human.sai
	bwa samse "$reference_human" "$base_output"/mapping/"$base_name"_human.sai "$base_output"/mapping/"$base_name"_dog.fq > "$base_output"/mapping/"$base_name"_human.sam
fi

#Filter only unmapped reads from human alignment

samtools view -bh -q 20 "$base_output"/mapping/"$base_name"_human.sam -o "$base_output"/mapping/"$base_name"_to_eliminate.bam -U "$base_output"/mapping/"$base_name"_no_human_contamination.bam 
samtools bam2fq "$base_output"/mapping/"$base_name"_no_human_contamination.bam > "$base_output"/mapping/"$base_name"_no_human_contamination.fq


#bwa aln - mapping against pig reference

if [ $pair = 'true' ]; then
	bwa aln -t "$threads" -l 1024 -o 2 -n 0.01 "$reference_pig" "$base_output"/mapping/"$base_name"_no_human_contamination.fq > "$base_output"/mapping/"$base_name"_pig.sai
	bwa samse "$reference_pig" "$base_output"/mapping/"$base_name"_pig.sai "$base_output"/mapping/"$base_name"_no_human_contamination.fq > "$base_output"/mapping/"$base_name"_pig.sam
else
	bwa aln -t "$threads" -l 1024 -o 2 -n 0.01 "$reference_pig" "$base_output"/mapping/"$base_name"_no_human_contamination.fq > "$base_output"/mapping/"$base_name"_pig.sai
	bwa samse "$reference_pig" "$base_output"/mapping/"$base_name"_pig.sai "$base_output"/mapping/"$base_name"_no_human_contamination.fq > "$base_output"/mapping/"$base_name"_pig.sam
fi

#Filter only unmapped reads from pig alignment

samtools view -bh -q 20 "$base_output"/mapping/"$base_name"_pig.sam -o "$base_output"/mapping/"$base_name"_to_eliminate.bam -U "$base_output"/mapping/"$base_name"_no_pig_human_contamination.bam 
samtools bam2fq "$base_output"/mapping/"$base_name"_no_pig_human_contamination.bam > "$base_output"/mapping/"$base_name"_no_pig_human_contamination.fq


#bwa aln - mapping against chicken reference

if [ $pair = 'true' ]; then
	bwa aln -t "$threads" -l 1024 -o 2 -n 0.01 "$reference_chicken" "$base_output"/mapping/"$base_name"_no_pig_human_contamination.fq > "$base_output"/mapping/"$base_name"_chicken.sai
	bwa samse "$reference_chicken" "$base_output"/mapping/"$base_name"_chicken.sai "$base_output"/mapping/"$base_name"_no_pig_human_contamination.fq > "$base_output"/mapping/"$base_name"_chicken.sam
else
	bwa aln -t "$threads" -l 1024 -o 2 -n 0.01 "$reference_chicken" "$base_output"/mapping/"$base_name"_no_pig_human_contamination.fq > "$base_output"/mapping/"$base_name"_chicken.sai
	bwa samse "$reference_chicken" "$base_output"/mapping/"$base_name"_chicken.sai "$base_output"/mapping/"$base_name"_no_pig_human_contamination.fq > "$base_output"/mapping/"$base_name"_chicken.sam
fi

#Filter only unmapped reads from chicken alignment

samtools view -bh -q 20 "$base_output"/mapping/"$base_name"_chicken.sam -o "$base_output"/mapping/"$base_name"_to_eliminate.bam -U "$base_output"/mapping/"$base_name"_no_chicken_pig_human_contamination.bam 
samtools bam2fq "$base_output"/mapping/"$base_name"_no_chicken_pig_human_contamination.bam > "$base_output"/mapping/"$base_name"_no_chicken_pig_human_contamination.fq


#bwa aln - mapping against cow reference

if [ $pair = 'true' ]; then
	bwa aln -t "$threads" -l 1024 -o 2 -n 0.01 "$reference_cow" "$base_output"/mapping/"$base_name"_no_chicken_pig_human_contamination.fq > "$base_output"/mapping/"$base_name"_cow.sai
	bwa samse "$reference_cow" "$base_output"/mapping/"$base_name"_cow.sai "$base_output"/mapping/"$base_name"_no_chicken_pig_human_contamination.fq > "$base_output"/mapping/"$base_name"_cow.sam
else
	bwa aln -t "$threads" -l 1024 -o 2 -n 0.01 "$reference_cow" "$base_output"/mapping/"$base_name"_no_chicken_pig_human_contamination.fq > "$base_output"/mapping/"$base_name"_cow.sai
	bwa samse "$reference_cow" "$base_output"/mapping/"$base_name"_cow.sai "$base_output"/mapping/"$base_name"_no_chicken_pig_human_contamination.fq > "$base_output"/mapping/"$base_name"_cow.sam
fi

#Filter only unmapped reads from cow alignment

samtools view -bh -q 20 "$base_output"/mapping/"$base_name"_cow.sam -o "$base_output"/mapping/"$base_name"_to_eliminate.bam -U "$base_output"/mapping/"$base_name"_no_cow_chicken_pig_human_contamination.bam 
samtools bam2fq "$base_output"/mapping/"$base_name"_no_cow_chicken_pig_human_contamination.bam > "$base_output"/mapping/"$base_name"_no_cow_chicken_pig_human_contamination.fq



#bwa aln - mapping against sheep reference

if [ $pair = 'true' ]; then
	bwa aln -t "$threads" -l 1024 -o 2 -n 0.01 "$reference_sheep" "$base_output"/mapping/"$base_name"_no_cow_chicken_pig_human_contamination.fq > "$base_output"/mapping/"$base_name"_sheep.sai
	bwa samse "$reference_sheep" "$base_output"/mapping/"$base_name"_sheep.sai "$base_output"/mapping/"$base_name"_no_cow_chicken_pig_human_contamination.fq > "$base_output"/mapping/"$base_name"_sheep.sam
else
	bwa aln -t "$threads" -l 1024 -o 2 -n 0.01 "$reference_sheep" "$base_output"/mapping/"$base_name"_no_cow_chicken_pig_human_contamination.fq > "$base_output"/mapping/"$base_name"_sheep.sai
	bwa samse "$reference_sheep" "$base_output"/mapping/"$base_name"_sheep.sai "$base_output"/mapping/"$base_name"_no_cow_chicken_pig_human_contamination.fq > "$base_output"/mapping/"$base_name"_sheep.sam
fi

#Filter only unmapped reads from sheep alignment

samtools view -bh -q 20 "$base_output"/mapping/"$base_name"_sheep.sam -o "$base_output"/mapping/"$base_name"_to_eliminate.bam -U "$base_output"/mapping/"$base_name"_no_sheep_cow_chicken_pig_human_contamination.bam 
samtools bam2fq "$base_output"/mapping/"$base_name"_no_sheep_cow_chicken_pig_human_contamination.bam > "$base_output"/mapping/"$base_name"_no_sheep_cow_chicken_pig_human_contamination.fq


#bwa aln - mapping against goat reference

if [ $pair = 'true' ]; then
	bwa aln -t "$threads" -l 1024 -o 2 -n 0.01 "$reference_goat" "$base_output"/mapping/"$base_name"_no_sheep_cow_chicken_pig_human_contamination.fq > "$base_output"/mapping/"$base_name"_goat.sai
	bwa samse "$reference_goat" "$base_output"/mapping/"$base_name"_goat.sai "$base_output"/mapping/"$base_name"_no_sheep_cow_chicken_pig_human_contamination.fq > "$base_output"/mapping/"$base_name"_goat.sam
else
	bwa aln -t "$threads" -l 1024 -o 2 -n 0.01 "$reference_goat" "$base_output"/mapping/"$base_name"_no_sheep_cow_chicken_pig_human_contamination.fq > "$base_output"/mapping/"$base_name"_goat.sai
	bwa samse "$reference_goat" "$base_output"/mapping/"$base_name"_goat.sai "$base_output"/mapping/"$base_name"_no_sheep_cow_chicken_pig_human_contamination.fq > "$base_output"/mapping/"$base_name"_goat.sam
fi

#Filter only unmapped reads from goat alignment

samtools view -bh -q 20 "$base_name"_goat.sam -o "$base_output"/mapping/"$base_name"_to_eliminate.bam -U "$base_name"_no_goat_sheep_cow_chicken_pig_human_contamination.bam # -h guardar com head cabeçário e -b para guardar ficheiro .bam
samtools bam2fq "$base_output"/mapping/"$base_name"_no_goat_sheep_cow_chicken_pig_human_contamination.bam > "$base_output"/mapping/"$base_name"_no_goat_sheep_cow_chicken_pig_human_contamination.fq


#bwa aln - mapping again against dog reference

if [ $pair = 'true' ]; then
	bwa aln -t "$threads" -l 1024 -o 2 -n 0.01 "$reference" "$base_output"/mapping/"$base_name"_no_goat_sheep_cow_chicken_pig_human_contamination.fq > "$base_output"/mapping/"$base_name".sai
	bwa samse "$reference" "$base_output"/mapping/"$base_name".sai "$base_output"/mapping/"$base_name"_no_goat_sheep_cow_chicken_pig_human_contamination.fq > "$base_output"/mapping/"$base_name".sam
else
	bwa aln -t "$threads" -l 1024 -o 2 -n 0.01 "$reference" "$base_output"/mapping/"$base_name"_no_goat_sheep_cow_chicken_pig_human_contamination.fq > "$base_output"/mapping/"$base_name".sai
	bwa samse "$reference" "$base_output"/mapping/"$base_name".sai "$base_output"/mapping/"$base_name"_no_goat_sheep_cow_chicken_pig_human_contamination.fq > "$base_output"/mapping/"$base_name".sam
fi


#Filter only mapped reads from dog alignment

samtools view -bh -q "$minmapqual" "$base_output"/mapping/"$base_name".sam -o "$base_output"/mapping/"$base_name".bam 
samtools sort "$base_name".bam > "$base_name"_sorted.bam
samtools index "$base_name"_sorted.bam

qualimap bamqc -nt "$threads" -c -bam "$base_name"_sorted.bam --java-mem-size=32G

gzip *.fq

rm $base_name'_to_eliminate.bam' $base_name'_dog.sai'  $base_name'_dog.sam' $base_name'_human.sai' $base_name'_human.sam' $base_name'_pig.sai' $base_name'_pig.sam' $base_name'_chicken.sai' $base_name'_chicken.sam' $base_name'_cow.sai' $base_name'_cow.sam' $base_name'_sheep.sam' $base_name'_sheep.sai' $base_name'_goat.sam' $base_name'_goat.sai'	$base_name'.sam' $base_name'.sai'	$base_name'.bam' 

#picartools (add read groups)
java -XX:ParallelGCThreads="$threads" -XX:ConcGCThreads="$threads" -jar $PICARD AddOrReplaceReadGroups ID="$biosample" SM="$biosample" PU="PU" LB="LB" PL="illumina" VALIDATION_STRINGENCY="LENIENT" I="$base_name"_sorted.bam O="$base_name"_RG.bam #add uma coluna de identificacao no .bam com identificador para depois pordermos remover os duplicados

#picartools (remove duplicates); reuse samtools index and samtools sort because we are not sure if the Picard desorganized the reads. 

samtools sort "$base_name"_RG.bam > "$base_name"_sorted_RG.bam
samtools index "$base_name"_sorted_RG.bam
java -XX:ParallelGCThreads="$threads" -XX:ConcGCThreads="$threads" -jar $PICARD MarkDuplicates VALIDATION_STRINGENCY="LENIENT" REMOVE_DUPLICATES="true" I="$base_name"_sorted_RG.bam O="$base_name"_no_dups.bam M="marked_dup_metrics.txt"

#samtools (filtering by flag and qual)

samtools sort "$base_name"_no_dups.bam > "$base_name"_sorted_no_dups.bam
samtools index "$base_name"_sorted_no_dups.bam
samtools view -F4 -q "$minmapqual" -bh "$base_name"_sorted_no_dups.bam -o "$base_name"_filtered.bam
samtools sort "$base_name"_filtered.bam > "$base_name"_filtered_sorted.bam
samtools index "$base_name"_filtered_sorted.bam


#########################
####### mapDamage #######
#########################

mkdir $base_output'/mapDamage'

mapDamage -l "$readsize" -d $base_output'/mapDamage' --rescale --rescale-out "$base_name"_filtered_sorted_mapdamage.bam -i "$base_name"_filtered_sorted.bam -r "$reference"

samtools index "$base_name"_filtered_sorted_mapdamage.bam

rm  $base_name'_sorted.bam.bai' $base_name'_RG.bam' $base_name'_sorted_RG.bam' $base_name'_sorted_RG.bam.bai' $base_name'_no_dups.bam' $base_name'_sorted_no_dups.bam' $base_name'_sorted_no_dups.bam.bai' 

qualimap bamqc -nt "$threads" -c -bam "$base_name"_filtered_sorted_mapdamage.bam --java-mem-size=32G

#Trim Bam file (trim 10 base pairs at the extremities)
bam trimBam "$base_name"_filtered_sorted_mapdamage.bam "$base_name"_trimBam.bam 10
samtools sort "$base_name"_trimBam.bam > "$base_name"_trimBam_sorted.bam
samtools index "$base_name"_trimBam_sorted.bam

#########################
#### VARIANT CALLING ####
#########################

mkdir $base_output'/variant_calling'

gatk HaplotypeCaller -R "$reference" -I "$base_name"_trimBam_sorted.bam --bam-output "$base_name"_GATK_out.bam --pcr-indel-model CONSERVATIVE --dont-use-soft-clipped-bases true --active-probability-threshold 0.002 --force-active true  -O "$base_output"/variant_calling/"$base_name".vcf

#vcf filter using GATK

cd $base_output'/variant_calling'
gatk VariantFiltration -R "$reference" -V "$base_output"/variant_calling/"$base_name".vcf --filter-name 'quality' --filter-expression "QD < $snp_quality" --filter-name 'readpos' --filter-expression "ReadPosRankSum < -1.0 && ReadPosRankSum > 1.0" --filter-name 'baserank' --filter-expression "AS_BaseQRankSum < -1.0 && AS_BaseQRankSum > 1.0" --genotype-filter-name 'DP' --genotype-filter-expression "DP < $snp_coverage" --set-filtered-genotype-to-no-call true -O "$base_output"/variant_calling/"$base_name"_filtered.vcf

gatk SelectVariants -V "$base_output"/variant_calling/"$base_name"_filtered.vcf --exclude-filtered true --exclude-non-variants true --remove-unused-alternates true --restrict-alleles-to BIALLELIC -O "$base_output"/variant_calling/"$base_name"_PASS_ONLY.vcf


########################
###### SNP effects #####
########################

java -jar /opt/anaconda3/share/snpeff-4.3.1t-1/snpEff.jar "$snpEffDB" "$base_output"/variant_calling/"$base_name"_PASS_ONLY.vcf > "$base_output"/variant_calling/"$base_name"_effects.vcf
