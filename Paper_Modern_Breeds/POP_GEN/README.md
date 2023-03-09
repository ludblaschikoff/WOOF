# SNP calling

*variant_calling_loop_modern_nuclear.sh*

From Base Quality Recalibrated BAM files I perform variant calling for modern samples dataset. Here I use GATK (v4.1.7.0) HaplotypeCaller to call variants and generate individual gVCF files from each BAM file. Next all gVCFs are combined into a single file with GATK's GenomicsDBImport. I run joint genotyping with GATK's GenotypeGVCFs. Finally, Gatk SelectVariants is used to call variants based on a set of filters.

*stats.sh*

This script outputs several files containing basic VCF statistics, including allele frequencies, mean depth per individual and per site, site quality, proportion of missing data, etc. These files will be used to generate plots in a RMarkdown report that can help you choose the settings to filter your VCF file.

*VCF_stats.rmd*

Plotting of several VCF statistics before and after filtering. Run stats.sh to generate the necessary input files. To render the RMD script run:
R --> rmarkdown::render("VCF_stats.rmd", clean=TRUE, output_format="html_document")
Command-line --> Rscript -e "rmarkdown::render('VCF_stats.rmd', clean=TRUE)"
