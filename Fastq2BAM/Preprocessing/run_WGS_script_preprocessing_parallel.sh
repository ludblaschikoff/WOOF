#!/bin/bash
#Script to run reads quality control and filtering (pre-processing) in parallel

#USAGE: path_to_script/run_WGS_script_preprocessing_DOGPT12_46.sh path/to/read1:path/to/read2 trimqual minreadlength threads path/to/output samplename

display_usage() {
echo "
1nd argument is the complete path for the directory where fastq files are located.
2nd argument is the complete path to the directory where results must be saved
3rd argument is the mininum base quality for trimming. example "30"
4th argument is the minimum read length after adapter removal. example "80"
"
}


#check if required arguments are there and display usage message
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ]  ; then
        printf "Please provide the arguments required for the script.\n\n"
        display_usage
        exit 1
fi

#########################
####### VARIABLES #######
#########################

fastq_files_dir="$1"
base_output="$2"
trimqual="$3"
minreadlength="$4"		

threads=20    # xargs -n 1 -P $threads -> This means that you will activate 20 threads to run 20 different commands (each thread will run a sample you are giving). Each of those 20 individual commands can only use 1 thread, so the total is 20, so in fastqc we use -t 1. So if you have to process 50 samples, you give 50 threads when you run the script. The 50 samples run in parallel, taking the same time to run as just 1 sample, because they don't wait for each other.

#RUN THESE SEPARATE COMMANDS AND EVALUATE THE QUALITY OF THE READS BEFORE RUNNING THE SCRIPT

#mkdir $base_output
#mkdir $base_output'/preprocessing'
#mkdir $base_output/'preprocessing/fastqc_raw_reads'
#fastqc -t 50 "$fastq_files_dir"/"$base_name"_1.fastq.gz "$fastq_files_dir"/"$base_name"_2.fastq.gz -o "$base_output"/preprocessing/fastqc_raw_reads/ 


##################################
####### LOOP PREPROCESSING #######
##################################

ls $fastq_files_dir/*_1.fastq.gz | xargs -n 1 -P "$threads" -I {} sh -c 'base_name="$(basename {} | cut -d_ -f1)"; echo "$base_name";  mkdir "'$base_output'"/"$base_name"/preprocessing/fastqc_clean_reads/ ; mkdir "'$base_output'"/"$base_name"/preprocessing/cutadapt/ ; cutadapt -b "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG" -B "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT" -q '$trimqual','$trimqual' -m '$minreadlength' --max-n 0 -o "'$base_output'"/"$base_name"/preprocessing/cutadapt/"$base_name"_R1_trimmed.fq.gz -p "'$base_output'"/"$base_name"/preprocessing/cutadapt/"$base_name"_R2_trimmed.fq.gz "'$fastq_files_dir'"/"$base_name"_1.fastq.gz "'$fastq_files_dir'"/"$base_name"_2.fastq.gz > "'$base_output'"/"$base_name"/preprocessing/"$base_name".log; fastqc -t 1 "'$base_output'"/"$base_name"/preprocessing/cutadapt/"$base_name"_R1_trimmed.fq.gz "'$base_output'"/"$base_name"/preprocessing/cutadapt/"$base_name"_R2_trimmed.fq.gz -o "'$base_output'"/"$base_name"/preprocessing/fastqc_clean_reads >> "'$base_output'"/"$base_name"/preprocessing/"$base_name".log ' 		  #xargs (short for "eXtended ARGuments") is for running commands in parallel. "n-1" you are saying "take one file at a time".  "-I {}" this argument serves to replace the input (in this case the input is the file name) by a string. "sh –c" run the command as a shell, that way each command runs in a separate shell.
