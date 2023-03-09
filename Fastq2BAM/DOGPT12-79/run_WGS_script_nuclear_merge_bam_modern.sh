#!/bin/bash
#Script to merge bam files for samples with two sequencing runs.

#USAGE: bash run_WGS_script_nuclear_merge_bam_modern.sh <input1.bam> <input2.bam> threads path/output samplename

display_usage() {
echo '
1th argument is the path to the first bam file to merge
2th argument is the path to the second bam file to merge
3th argument is the number of threads available to use. example "15"
4th argument is the complete path to the directory were results must be saved
5th argument is the sample name'
}


#check if required arguments are there and display usage message
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ -z "$5" ] ; then
        printf "Please provide the arguments required for the script.\n\n"
        display_usage
        exit 1
fi

#############################
######### VARIABLES #########
#############################

bam1="$1"
bam2="$2"
threads="$3"
base_output="$4"
base_name="$5"


#########################
######### MERGE #########
#########################

cd $base_output

samtools merge --threads "$threads" "$base_name".bam "$bam1" "$bam2"

samtools sort "$base_name".bam > "$base_name"_sorted.bam
samtools index "$base_name"_sorted.bam

rm "$base_name".bam

qualimap bamqc -nt "$threads" -c -bam "$base_name"_sorted.bam --java-mem-size=32G

echo "Done!"
