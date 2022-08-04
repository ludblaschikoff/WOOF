No preprocessing dessas amostras modernas ao correr o fastqc das raw reads me deparei com polyG's sequencias overrepresented na read reverse.

Removi utilizando o script run_WGS_script_preprocessing_PointerProject, onde no cutadapt fiz: 

cutadapt --pair-filter=any -b "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG" -B "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG" -q "$trimqual","$trimqual" -j "$threads" -m "$minreadlength" --max-n 0 -o "$base_output"/preprocessing/cutadapt/"$base_name"_R1_trimmed.fq.gz -p "$base_output"/preprocessing/cutadapt/"$base_name"_R2_trimmed.fq.gz "$read1" "$read2"

Porém, algumas amostras continuaram com poly G´s. Então repeti o mesmo script.

Porém, a parte testei o seguinte comando cutadapt numa amostra sem nenhum processamento:

cutadapt -b "G{100}" -B "G{100}" -q "$trimqual","$trimqual" -j "$threads" -m "$minreadlength" --max-n 0 -o "$base_output"/preprocessing/cutadapt/"$base_name"_R1_trimmed.fq.gz -p "$base_output"/preprocessing/cutadapt/"$base_name"_R2_trimmed.fq.gz "$read1" "$read2"			#Cutadapt allows you to include N number of repeats by using the command '-b "G{100}". It simply trims away up to whatever you specify in the curly brackets.

E o resultado foi logo amostras sem polyG´s sem precisar repetir o cutadapt.

Não refiz a análise porque o resultado acaba por ser o mesmo ou muito semelhante. Mas da próxima vez que tiver esse problema, convém fazer o segundo comando do cutadapt.
