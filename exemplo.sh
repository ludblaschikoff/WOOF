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
