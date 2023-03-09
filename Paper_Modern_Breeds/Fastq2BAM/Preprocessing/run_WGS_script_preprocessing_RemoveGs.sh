#!/bin/bash
#Script to run reads quality control, filtering (pre-processing) and remove overrepresented polyG sequences.

#USAGE: bash path_to_script/run_WGS_script_preprocessing_RemoveGs.sh path/to/read1:path/to/read2 trimqual minreadlength threads path/to/output samplename

display_usage() {
echo '1st argument must be the path to the sample fastq file read1 [and read2 if paired] in the following format: "path/to/read1:" OR "path/to/read1:path/to/read2"
2nd argument is the mininum base quality for trimming
3rd argument is the minimum read length after adapter removal
4th argument is the number of threads available to use. example "15"
5th argument is the complete path to the directory were results must be saved
6th argument is the sample name)'
}


#check if required arguments are there and display usage message
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ -z "$5" ] || [ -z "$6" ] ; then
        printf "Please provide the arguments required for the script.\n\n"
        display_usage
        exit 1
fi

#########################
####### VARIABLES #######
#########################

#variables: read1 will be whatever I have to the left of the ":" ; read2 will be all I have to the right of the ":". ex: $1=P9306_1026_R1...fq:P9306_1026_R2.
read1=$(echo $1 | cut -d ":" -f1)
read2=$(echo $1 | cut -d ":" -f2)

pair='true'

if [ ${#read2} -eq 0 ]; then #[ ${#read2} -eq 0 ] means that the number of elements of the variable read 2 is equal to zero.
	pair='false'
fi

trimqual="$2"
minreadlength="$3"
threads="$4"    #use according to computational resources available
base_output="$5"
base_name="$6"

#########################
##### PREPROCESSING #####
#########################

#RUN THESE SEPARATE COMMANDS AND EVALUATE THE QUALITY OF THE READS BEFORE RUNNING THE SCRIPT

#mkdir $base_output
#mkdir $base_output'/preprocessing'
#mkdir $base_output/'preprocessing/fastqc_raw_reads'
#fastqc -t "$threads" -o "$base_output/preprocessing/fastqc_raw_reads" "$read1"
#fastqc -t "$threads" -o "$base_output/preprocessing/fastqc_raw_reads" "$read2"


#cutadapt 

mkdir $base_output'/preprocessing/cutadapt'

if [ $pair = 'true' ]; then
	cutadapt --pair-filter=any -b "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG" -B "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT" --nextseq-trim="$trimqual" -q "$trimqual","$trimqual" -j "$threads" -m "$minreadlength" --max-n 0 -o "$base_output"/preprocessing/cutadapt/"$base_name"_R1_trimmed.fq.gz -p "$base_output"/preprocessing/cutadapt/"$base_name"_R2_trimmed.fq.gz "$read1" "$read2"
else
	cutadapt --pair-filter=any -b "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG" --nextseq-trim="$trimqual" -q "$trimqual","$trimqual" -j "$threads" -m "$minreadlength" --max-n 0 -o "$base_output"/preprocessing/cutadapt/"$base_name"_trimmed.fq.gz "$read1"
fi


#FASTQC clean reads

mkdir $base_output'/preprocessing/fastqc_clean_reads'

if [ $pair = 'true' ]; then
	fastqc -t "$threads" -o "$base_output"/preprocessing/fastqc_clean_reads "$base_output"/preprocessing/cutadapt/"$base_name"_R1_trimmed.fq.gz
	fastqc -t "$threads" -o "$base_output"/preprocessing/fastqc_clean_reads "$base_output"/preprocessing/cutadapt/"$base_name"_R2_trimmed.fq.gz
else
	fastqc -t "$threads" -o "$base_output"/preprocessing/fastqc_clean_reads "$base_output"/preprocessing/cutadapt/"$base_name"_trimmed.fq.gz
fi

echo "Done!"
