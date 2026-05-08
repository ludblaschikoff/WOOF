##https://cran.r-project.org/web/packages/StAMPP/StAMPP.pdf
#by Luke Pembleton 
install.packages("StAMPP")
install.packages("pinfsc50")
install.packages("vcfR")
install.packages("adegenet")

library(pinfsc50)
library(vcfR)
library(adegenet)
library(StAMPP)


##Read vcf file into into R 
vcf<-read.vcfR("dog.kingship_clean_recode_samplesize_NewFID_vcf.vcf", verbose = FALSE) 
##Convert vcf into genlight object 
x <- vcfR2genlight(vcf) 

### convert genlight object to a matrix in stampp format 
x2 <- as.matrix(x) #convert genlight object to matrix 
sample <- row.names(x2) #sample names 

#pop.names <- pop(x) #extract ploidy info from genlight object (however this is not available when imported via a vcf file) 
#pop.names <- #provide a vector here of the same length of the number of samples, with a corresponding population name/id
pop.names <- c("CLWD","PTWD","PTWD","EMD","PTWD","PTWD","SMCD","SMCD","SMCD","SMCD","CLWD","CLWD","CLWD","CLWD","CLWD","PWHLwire","PWHMwire","PWHMwire","PWHSwire","PWHSwire","PWHSshort","ALTM","ALTM","ALTM","ALTM","ALTM","EMD","EMD","EMD","EMD","PTSD","PTSD","PTSD","PTSD","PPOI","PPOI","PPOI","ALGB","ALGB","ALGB","ALGB", "PPOI")
pop.names <- c("CLWD","PTWD","PTWD","EMD","PTWD","PTWD","SMCD","SMCD","SMCD","SMCD","CLWD","CLWD","CLWD","CLWD","CLWD","PWH","PWH","PWH","PWH","PWH","PWH","ALTM","ALTM","ALTM","ALTM","ALTM","EMD","EMD","EMD","EMD","PTSD","PTSD","PTSD","PTSD","PPOI","PPOI","PPOI","ALGB","ALGB","ALGB","ALGB", "PPOI")
#e.g. c("popA", "popA", "popB", "popB", "popB", "popC", "popC") 

ploidy <- ploidy(x) #extract ploidy info from genlight object 
x2 = x2 * (1/ploidy) #convert allele counts to frequency 
x2[is.na(x2)] = NaN 
format <- vector(length = length(sample))
#format id for the genotype data
format[1:length(format)] = "freq"  

x.stampp <- as.data.frame(cbind(sample, pop.names, ploidy, format, x2)) #convert to basic r data.frame suitable to stamppConvert 

geno <- stamppConvert(x.stampp, 'r') 

#you should now be able to run all the stampp commands with 'geno' as the input object 
# e.g. fst <- stamppFst(geno)

#Calculate Nei's Genetic Distance between population

distance.mat <- stamppNeisD(geno, pop = TRUE)

##Export to Phylip Format

stamppPhylip(distance.mat, file = "nei_distances_StAmpp2")
