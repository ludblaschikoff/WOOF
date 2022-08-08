#!/bin/bash
#Script to run entire analysis at once

display_usage() {
echo '
1st argument is the path to the reference genome index
2nd argument is the complete path for the directory where bam files are located.
3rd argument is the complete path to the directory where results must be saved
4th argument is the max RAM available for this job, in Gb (ex: 200)
'
}


#check if required arguments are there and display usage message
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] ; then
        printf "Please provide the arguments required for the script.\n\n"
        display_usage
        exit 1
fi


#variantes
reference="$1"
bam_files_dir="$2"
base_output="$3"
maxRAM="$4"				


#CALCULAMOS A RAM DISPONIVEL PARA CADA RUN

threads=30
ind_ram=$(( $maxRAM/$threads ))		#run é o numero de amostras a serem processadas. é o mesmo que threads.

ind_ram=$ind_ram'g'					##atualizo o valor da variavel $ind_ram para $ind_ram + g para o nome da variável levar a letra g.

#### BaseRecalibrator ####

ls $bam_files_dir/*.bam | xargs -n 1 -P $threads -I {} sh -c 'base_name="$(basename {} | cut -d_ -f1)"; echo "$base_name"; gatk --java-options "-Xmx'$ind_ram'" BaseRecalibrator -R '$reference' -I {} --known-sites /mnt/storage1/DATA/canis/SNP_array/variants_known_sites/All.vcf.gz -O "'$base_output'"/"$base_name"/NUCLEAR/mapping/recal_data_"$base_name".table ;  gatk --java-options -Xmx4G ApplyBQSR -R '$reference' -I {} --bqsr-recal-file "'$base_output'"/"$base_name"/NUCLEAR/mapping/recal_data_"$base_name".table -O "'$base_output'"/"$base_name"/NUCLEAR/mapping/BQSR_"$base_name"_no_dups_sorted.bam ' 		 # -P $runs é o numero de processadores que distribuo para cada amostra; #xargs (short for "eXtended ARGuments") serve para correr comandos em paralelo. is a command on Unix and most Unix-like operating systems used to build and execute commands from standard input. It converts input from standard input (ficheiros) into arguments to a command.

	
