#!/bin/bash
#Script to run entire analysis at once

display_usage() {
echo '
1st argument is the path to the reference genome index
2nd argument is the number of threads available to use. example "15"
3rd argument is the complete path to the directory were results must be saved
4th argument is the sample name
5th argument is the biosample name
6th argument is the minimum mapping quality. example "20"
7th argument is the library name'
}


#check if required arguments are there and display usage message
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ -z "$5" ] || [ -z "$6" ] || [ -z "$7" ]; then
        printf "Please provide the arguments required for the script.\n\n"
        display_usage
        exit 1
fi


#variantes 
reference="$1"
threads="$2"
base_output="$3"
base_name="$4"
biosample="$5"
minmapqual="$6"
library_name="$7"


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

rm $base_name'.sam' $base_name'_dog_mappedandunmapped.bam' $base_name'.bam' $base_name'_RG.bam' $base_name'_sorted_RG.bam' $base_name'_sorted_RG.bam.bai'  $base_name'_no_dups.bam'

qualimap bamqc -nt "$threads" -c -bam "$base_name"_sorted_no_dups_.bam --java-mem-size=32G
 