argv <- commandArgs(T)
INPUT<- argv[1]
OUTPUT<- argv[2]
GENOME<- argv[3]

#source("01_scripts/Rscripts/fix_sniffles_delly.R")
source("/storage/homefs/kw23y068/software/scripts/R_Delly/fix_sniffles_delly.R")

fix_sniffles(input_vcf=INPUT, output_vcf=OUTPUT, refgenome = GENOME)
