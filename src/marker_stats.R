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
  gc <- tally[,2] / rowSums(tally)
  names(gc) <- markers
  return(gc)
}

nFlankingSNPs <- function(vcf, markers, flanking_bp = 50){
  rr <- SummarizedExperiment::rowRanges(vcf)
  rr2 <- getKaspRange(rr[markers], flanking_bp = flanking_bp)
  fo <- GenomicRanges::findOverlaps(rr2, rr, type = "any")
  nhits <- S4Vectors::countLnodeHits(fo) - 1
  names(nhits) <- markers
  return(nhits)
}

MAF <- function(vcf, markers){
  
}

LD <- function(){}

# Function to format SNPs for KASP assay, annotating flanking SNPs.
formatKasp <- function(vcf, markers, refgenome, flanking_bp = 50){
  rr <- SummarizedExperiment::rowRanges(vcf)
  rr2 <- getKaspRange(rr[markers], flanking_bp = flanking_bp)
  seq <- Rsamtools::scanFa(refgenome, rr2, as = "DNAStringSet")
  fo <- GenomicRanges::findOverlaps(rr2, rr, type = "any")
  for(i in seq_along(markers)){
    toshift <- BiocGenerics::start(rr2)[i] - 1
    theseSNPs <- rr[S4Vectors::subjectHits(fo)[S4Vectors::queryHits(fo) == i]]
    pos <- BiocGenerics::start(theseSNPs) - toshift
    ambig <- Biostrings::mergeIUPACLetters(paste0(theseSNPs$REF, unlist(theseSNPs$ALT)))
    seq[[i]] <- Biostrings::replaceLetterAt(seq[[i]], pos, ambig)
  }
  outstrings <- paste0(XVector::subseq(seq, start = 1, width = flanking_bp),
                       "[",
                       XVector::subseq(seq, start = flanking_bp + 1, width = 1),
                       "]",
                       XVector::subseq(seq, start = flanking_bp + 2, width = flanking_bp))
  return(data.frame(SNP_ID = markers,
                    Sequence = outstrings,
                    stringsAsFactors = FALSE))
}
