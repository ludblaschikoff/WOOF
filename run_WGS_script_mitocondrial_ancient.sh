#!/bin/bash
#Script to run entire aDNA analysis at once

display_usage() {
echo '
1th argument is the path to the reference genome index
2th argument is the number of threads available to use. example "15"
3nd argument is the mininum base quality for trimming. example "30"
4rd argument is the minimum read length. example "30"
5th argument is the path to the reference_composite genome index
6th argument is the minimum mapping quality. example "20"
7th argument is the complete path to the directory were results must be saved
8th argument is the sample name
9th argument is the maximum read size (ex: 150bp)
10th argument is the biosample name
11th argument is the minimum SNP quality. example "20"'
}


#check if required arguments are there and display usage message
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ -z "$5" ] || [ -z "$6" ] || [ -z "$7" ] || [ -z "$8" ] || [ -z "$9" ] || [ -z "$10" ] || [ -z "$11" ] ; then
        printf "Please provide the arguments required for the script.\n\n"
        display_usage
        exit 1
fi

#variantes 
reference="$1"
threads="$2"
trimqual="$3"
minreadlenght="$4"
reference_composite="$5"
minmapqual="$6"
base_output="$7"
base_name="$8"
readsize="$9"
biosample="${10}"
minqual="${11}"



#########################
######## MAPPING ########
#########################

cd "$base_output"../preprocessing/

clean_fastq_file=cutadapt/"$base_name"_R2_trimmed.fq.gz				

if [ -f "$clean_fastq_file" ]; then
    pair='true'
else
    pair='false'
fi

#Merge clean reads
if [ $pair = 'true' ]; then
        mkdir $base_output'../preprocessing/merge_reads'

        AdapterRemoval --file1 "$base_output"../preprocessing/cutadapt/"$base_name"_R1_trimmed.fq.gz --file2 "$base_output"../preprocessing/cutadapt/"$base_name"_R2_trimmed.fq.gz --basename "$base_output"../preprocessing/merge_reads/"$base_name"_output_paired --collapse --minlength "$minreadlenght" --minalignmentlength 11 --mm 1 --minquality "$trimqual" --trimqualities --trimns --maxns 10 --gzip
fi


mkdir $base_output'/mapping'
cd "$base_output"/mapping

#bwa aln - mapping against dog reference

if [ $pair = 'true' ]; then
	bwa aln -t "$threads" -l 1024 -o 2 -n 0.01 "$reference" ../../preprocessing/merge_reads/"$base_name"_output_paired.collapsed.gz > "$base_output"/mapping/"$base_name"_dog.sai
	bwa samse "$reference" "$base_output"/mapping/"$base_name"_dog.sai ../../preprocessing/merge_reads/"$base_name"_output_paired.collapsed.gz > "$base_output"/mapping/"$base_name"_dog.sam
else
	bwa aln -t "$threads" -l 1024 -o 2 -n 0.01 "$reference" ../../preprocessing/cutadapt/"$base_name"_trimmed.fq.gz > "$base_output"/mapping/"$base_name"_dog.sai
	bwa samse "$reference" "$base_output"/mapping/"$base_name"_dog.sai ../../preprocessing/cutadapt/"$base_name"_trimmed.fq.gz > "$base_output"/mapping/"$base_name"_dog.sam
fi

#Create a file of mapped and unmapped reads (samtools view)

samtools view -bh "$base_output"/mapping/"$base_name"_dog.sam -o "$base_output"/mapping/"$base_name"_dog_mappedandunmapped.bam
samtools sort "$base_output"/mapping/"$base_name"_dog_mappedandunmapped.bam > "$base_output"/mapping/"$base_name"_dog_mappedandunmapped_sorted.bam
samtools index "$base_output"/mapping/"$base_name"_dog_mappedandunmapped_sorted.bam

#Filter only mapped reads (samtools view)
samtools view -bh -q "$minmapqual" "$base_output"/mapping/"$base_name"_dog.sam -o "$base_output"/mapping/"$base_name"_dog.bam 
samtools bam2fq "$base_output"/mapping/"$base_name"_dog.bam > "$base_output"/mapping/"$base_name"_dog.fq

#bwa aln - mapping against composite_genome

if [ $pair = 'true' ]; then
	bwa aln -t "$threads" -l 1024 -o 2 -n 0.01 "$reference_composite" "$base_output"/mapping/"$base_name"_dog.fq > "$base_output"/mapping/"$base_name"_no_contamination.sai
	bwa samse "$reference_composite" "$base_output"/mapping/"$base_name"_no_contamination.sai "$base_output"/mapping/"$base_name"_dog.fq > "$base_output"/mapping/"$base_name"_no_contamination.sam
else
	bwa aln -t "$threads" -l 1024 -o 2 -n 0.01 "$reference_composite" "$base_output"/mapping/"$base_name"_dog.fq > "$base_output"/mapping/"$base_name"_no_contamination.sai
	bwa samse "$reference_composite" "$base_output"/mapping/"$base_name"_no_contamination.sai "$base_output"/mapping/"$base_name"_dog.fq > "$base_output"/mapping/"$base_name"_no_contamination.sam
fi


#Filter only unmapped reads (samtools view)

samtools view -bh -q 20 "$base_output"/mapping/"$base_name"_no_contamination.sam -o "$base_output"/mapping/"$base_name"contamination.bam -U "$base_output"/mapping/"$base_name"_no_contamination.bam
samtools bam2fq "$base_output"/mapping/"$base_name"_no_contamination.bam > "$base_output"/mapping/"$base_name"_no_contamination.fq


#bwa aln - mapping again against dog reference

if [ $pair = 'true' ]; then
	bwa aln -t "$threads" -l 1024 -o 2 -n 0.01 "$reference" "$base_output"/mapping/"$base_name"_no_contamination.fq > "$base_output"/mapping/"$base_name".sai
	bwa samse "$reference" "$base_output"/mapping/"$base_name".sai "$base_output"/mapping/"$base_name"_no_contamination.fq > "$base_output"/mapping/"$base_name".sam
else
	bwa aln -t "$threads" -l 1024 -o 2 -n 0.01 "$reference" "$base_output"/mapping/"$base_name"_no_contamination.fq > "$base_output"/mapping/"$base_name".sai
	bwa samse "$reference" "$base_output"/mapping/"$base_name"_dog.sai "$base_output"/mapping/"$base_name"_no_contamination.fq > "$base_output"/mapping/"$base_name".sam
fi

#Number of TOTAL reads on sam file

samtools view -c "$base_output"/mapping/"$base_name".sam &> "$base_output"/mapping/"$base_name"_total_number_reads.txt

#Filter only mapped reads (samtools view)

samtools view -bh -q "$minmapqual" "$base_output"/mapping/"$base_name".sam -o "$base_output"/mapping/"$base_name".bam  	# -h guardar com head cabeçário e -b para guardar ficheiro .bam
samtools sort "$base_name".bam > "$base_name"_sorted.bam
samtools index "$base_name"_sorted.bam


#Creating summary reports
qualimap bamqc -nt "$threads" -c -bam "$base_name"_sorted.bam --java-mem-size=32G

gzip *.fq

rm $base_name'_dog_mappedandunmapped.bam' $base_name'_dog.sai' $base_name'_dog.sam' $base_name'_no_contamination.sai' $base_name'_no_contamination.sam'	$base_name'.sai'	$base_name'.sam' $base_name'.bam' 

#picartools (add read groups)

java -XX:ParallelGCThreads="$threads" -XX:ConcGCThreads="$threads" -jar $PICARD AddOrReplaceReadGroups VALIDATION_STRINGENCY="LENIENT" ID="$biosample" SM="$biosample" PU="PU" LB="LB" PL="illumina" I="$base_name"_sorted.bam O="$base_name"_RG.bam #add uma coluna de identificacao no .bam com identificador para depois pordermos remover os duplicados
samtools sort "$base_name"_RG.bam > "$base_name"_sorted_RG.bam  #reuse samtools index and samtools sort because we are not sure if the Picard desorganized the reads.
samtools index "$base_name"_sorted_RG.bam

#picartools (remove duplicates) 

java -XX:ParallelGCThreads="$threads" -XX:ConcGCThreads="$threads" -jar $PICARD MarkDuplicates VALIDATION_STRINGENCY="LENIENT" REMOVE_DUPLICATES="true" I="$base_name"_sorted_RG.bam O="$base_name"_no_dups.bam M="marked_dup_metrics.txt" &> logFile_removeduplicates.log
samtools sort "$base_name"_no_dups.bam > "$base_name"_sorted_no_dups.bam
samtools index "$base_name"_sorted_no_dups.bam


#samtools (filtering by flag and qual)

samtools view -q "$minmapqual" -bh "$base_name"_sorted_no_dups.bam -o "$base_name"_filtered.bam
samtools sort "$base_name"_filtered.bam > "$base_name"_filtered_sorted.bam
samtools index "$base_name"_filtered_sorted.bam



##############################
###Reads local realignment####
##############################

source /opt/anaconda3/etc/profile.d/conda.sh

conda activate GATK3 #comando para entrar no ambiente virtual

gatk3 -T RealignerTargetCreator -R $reference -I "$base_name"_filtered_sorted.bam -o targets.intervals
gatk3 -T IndelRealigner -R $reference -I "$base_name"_filtered_sorted.bam -targetIntervals targets.intervals -o "$base_name".final.bam --filter_bases_not_stored &> logFile_realignment.log

conda deactivate #comando para sair do ambiente virtual

samtools sort "$base_name".final.bam -o "$base_name".final.sort.bam
samtools index "$base_name".final.sort.bam

#########################
####### mapDamage #######
#########################

mkdir $base_output'/mapDamage'

mapDamage -l "$readsize" -d $base_output'/mapDamage' --rescale --rescale-out "$base_name"_mapdamage.bam -i "$base_name".final.sort.bam -r "$reference"

samtools index "$base_name"_mapdamage.bam

rm $base_name'_RG.bam' $base_name'_sorted_RG.bam' $base_name'_sorted_RG.bam.bai' $base_name'_no_dups.bam' $base_name'_sorted_no_dups.bam' $base_name'_sorted_no_dups.bam.bai' 

qualimap bamqc -nt "$threads" -c -bam "$base_name"_mapdamage.bam --java-mem-size=32G -outformat PDF:HTML
 

##########################
####CONSENSUS SEQUENCE####
##########################

mkdir $base_output'/consensus_sequence'
cd $base_output'/consensus_sequence'

angsd -doFasta 2  -out prefix.out -explode 1 -doCounts 1 -minQ $minqual -minMapQ $minmapqual -i "$base_output"/mapping/"$base_name"_mapdamage.bam -setMinDepth 3 -nThreads 20 -doDepth 1

