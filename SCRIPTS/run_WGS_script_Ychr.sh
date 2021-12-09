#!/bin/bash
#Script to run entire aDNA analysis at once

display_usage() {
echo '
1th argument is the path to the reference chromosome Y  genome index
2th argument is the complete path to the directory were results must be saved
3th argument is the sample name'
}


#check if required arguments are there and display usage message
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] ; then
        printf "Please provide the arguments required for the script.\n\n"
        display_usage
        exit 1
fi

#variantes 
reference_chrY="$1"
base_output="$2"
base_name="$3"


#Fica com as reads que não mapearam contra o cromossoma nuclear

cd $base_output

samtools view -bh -f4 ../mapping/"$base_name".bam | samtools bam2fq - > "$base_name"_teste.fq

# Mapear com o Y cromossoma.

bwa mem -t40 "$reference_chrY" -p "$base_name"_teste.fq -o "$base_name".sam

# Ficar apenas com as que mapearam contra o cromossoma Y do cão.

samtools view -q 60 -F4 -bh "$base_name".sam -o "$base_name"_chrY.bam


#Para contar o número de reads que mapearam para cada sequência da referência:

samtools view "$base_name"_chrY.bam | cut -f 3 | sort | uniq -c > "$base_name"_chrY.txt

gzip "$base_name"_teste.fq 

rm "$base_name".sam
