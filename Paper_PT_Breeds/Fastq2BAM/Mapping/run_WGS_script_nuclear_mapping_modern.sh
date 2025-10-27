#!/bin/bash
#Script to run mapping against reference genome

#USAGE: bash path_to_script/run_WGS_script_nuclear_mapping_modern.sh /path/to/reference.fa threads path/to/output samplename biosamplename

display_usage() {
echo '
1st argument is the path to the reference genome index
2nd argument is the number of threads available to use. example "15"
3rd argument is the complete path to the directory were results must be saved
4th argument is the sample name
5th argument is the biosample'
}


#check if required arguments are there and display usage message
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ -z "$5" ]; then
        printf "Please provide the arguments required for the script.\n\n"
        display_usage
        exit 1
fi


#########################
####### VARIABLES #######
#########################

reference="$1"
threads="$2"
base_output="$3"
base_name="$4"
biosample="$5"

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
 

#samtools (filtering by flag and qual)
samtools view -bh -f2 -q 30 "$base_name".sam -o "$base_name".bam  # -f 2 - mapped in proper pair - cannot be set when read is not paired
samtools sort "$base_name".bam > "$base_name"_sorted.bam
samtools index "$base_name"_sorted.bam


#picartools (add read groups)
java -XX:ParallelGCThreads="$threads" -XX:ConcGCThreads="$threads" -jar $PICARD AddOrReplaceReadGroups VALIDATION_STRINGENCY="LENIENT" ID="$biosample" SM="$biosample" PU="PU" LB="LB" PL="illumina" I="$base_name"_sorted.bam O="$base_name"_RG.bam #add uma coluna de identificacao no .bam com identificador para depois pordermos remover os duplicados

#picartools (remove duplicates) 

java -XX:ParallelGCThreads="$threads" -XX:ConcGCThreads="$threads" -jar $PICARD MarkDuplicates VALIDATION_STRINGENCY="LENIENT" REMOVE_DUPLICATES="true" I="$base_name"_RG.bam O="$base_name"_no_dups.bam M="marked_dup_metrics.txt"
samtools index "$base_name"_no_dups.bam

rm $base_name'_sorted.bam' $base_name'_sorted.bam.bai' $base_name'_RG.bam' $base_name'.sam'

#Stats
qualimap bamqc -nt "$threads" -c -bam "$base_name"_no_dups.bam --java-mem-size=32G

echo "Done!"
