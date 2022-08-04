# Modern_Processing

In this directory you can find scripts and pipelines for the modern dogs analysis, starting from fastq files until BQSR.bam files.

* Please do twice the run_WGS_script_preprocessing.sh script if raw reads keep with polyG's sequences overrepresented in reverse read.

* OR DO

cutadapt -b "G{100}" -B "G{100}" -q "$trimqual","$trimqual" -j "$threads" -m "$minreadlength" --max-n 0 -o "$base_output"/preprocessing/cutadapt/"$base_name"_R1_trimmed.fq.gz -p "$base_output"/preprocessing/cutadapt/"$base_name"_R2_trimmed.fq.gz "$read1" "$read2"			#Cutadapt allows you to include N number of repeats by using the command '-b "G{100}". It simply trims away up to whatever you specify in the curly brackets.

* run_WGS_script_BQSR.sh is the last script do be run.

* Some samples have more than one run accession. So I decided to mapping them separately, merge bam files and do BQSR.
