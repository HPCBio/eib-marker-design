# Functions for getting various marker statistics.
# The vcf argument is a VCF object imported into R using VariantAnnotation.
# The markers argument is a character vector of names of markers to analyze.
# The refgenome argument is a FaFile object pointing to the FASTA of the
# reference genome.
# The flanking_bp argument indicates how many basepairs to either side of the
# SNP should be considered.
# rr is a GRanges object indicating SNP locations.

getKaspRange <- function(rr, flanking_bp = 50){
  if(!all(BiocGenerics::width(rr) == 1)){
    stop("Ranges are too wide to be SNPs.")
  }
  return(GenomicRanges::promoters(rr, upstream = flanking_bp,
                                  downstream = flanking_bp + 1))
}

gcContent <- function(vcf, markers, refgenome, flanking_bp = 50){
  rr <- SummarizedExperiment::rowRanges(vcf)[markers]
  rr2 <- getKaspRange(rr, flanking_bp = flanking_bp)
  seq <- Rsamtools::scanFa(refgenome, rr2, as = "DNAStringSet")
  tally <- Biostrings::letterFrequency(seq, c("ATW", "GCS"))
  return(tally[,2] / rowSums(tally))
}
