# Extract SNP locations from large VCF, and save them to an R object for easy access.

# Run this script like this from the terminal:
# Rscript snp_positions_from_vcf.R --args xxxxxxx.vcf xxxxxx.rds

args <- commandArgs(trailingOnly = TRUE)

filein <- args[2]
fileout <- args[3]

library(VariantAnnotation)

myvcf <- readVcf(filein, param = ScanVcfParam(geno = NA))

saveRDS(myvcf, file = fileout)
