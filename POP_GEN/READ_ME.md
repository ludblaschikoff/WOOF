# SNP calling

*stats.sh*

This script outputs several files containing basic VCF statistics, including allele frequencies, mean depth per individual and per site, site quality, proportion of missing data, etc. These files will be used to generate plots in a RMarkdown report that can help you choose the settings to filter your VCF file.

*VCF_stats.rmd*

Plotting of several VCF statistics before and after filtering. Run stats.sh to generate the necessary input files. To render the RMD script run:
R --> rmarkdown::render("VCF_stats.rmd", clean=TRUE, output_format="html_document")
Command-line --> Rscript -e "rmarkdown::render('VCF_stats.rmd', clean=TRUE)"
