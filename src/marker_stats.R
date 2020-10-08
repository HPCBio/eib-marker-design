# Functions for getting various marker statistics.
# The vcf argument is a VCF object imported into R using VariantAnnotation.
# The markers argument is a character vector of names of markers to analyze.
# The refgenome argument is a FaFile object pointing to the FASTA of the
# reference genome.
# The flanking_bp argument indicates how many basepairs to either side of the
# SNP should be considered.
# rr is a GRanges object indicating SNP locations.

# Packages are specified for Bioconductor functions here just to make everything
# a little more robust to the user having loaded other packages.

# Function to get a GRanges object representing the flanking regions for a set
# of SNPs, given a GRanges object representing the SNPs.
getKaspRange <- function(rr, flanking_bp = 50){
  if(!all(BiocGenerics::width(rr) == 1)){
    stop("Ranges are too wide to be SNPs.")
  }
  return(GenomicRanges::promoters(rr, upstream = flanking_bp,
                                  downstream = flanking_bp + 1))
}

# Function to correct the chromosome names from a TASSEL VCF to match the
# reference genome.  x is a character vector, and a character vector of the same
# length is returned.
# !! You can modify this function if chromosome names need to be modified in
# some different way.  Alternatively, you can build a similar function to pass
# to the fixfn argument of the below functions.
fixTasselNames <- function(x){
  gsub("^OM", "chrom", x)
}

# Function to get a GRanges object from a VCF, but with the names corrected.
rowRanges_correctedSeqnames <- function(vcf, markers = rownames(vcf),
                                        fixfn = fixTasselNames){
  rr <- SummarizedExperiment::rowRanges(vcf)[markers]
  rr1 <- GenomicRanges::GRanges(fixfn(GenomeInfoDb::seqnames(rr)),
                                IRanges::ranges(rr))
  S4Vectors::mcols(rr1) <- S4Vectors::mcols(rr)
  return(rr1)
}

# Function to get the proportion GC content for the flanking regions for a set
# of SNPs.
gcContent <- function(vcf, markers, refgenome, flanking_bp = 50,
                      fixfn = fixTasselNames){
  rr <- rowRanges_correctedSeqnames(vcf, markers, fixfn)
  rr2 <- getKaspRange(rr, flanking_bp = flanking_bp)
  seq <- Rsamtools::scanFa(refgenome, rr2, as = "DNAStringSet")
  tally <- Biostrings::letterFrequency(seq, c("ATW", "GCS"))
  gc <- tally[,2] / rowSums(tally)
  names(gc) <- markers
  return(gc)
}

# Function to count the number of SNPs in the flanking region for each SNP.
nFlankingSNPs <- function(vcf, markers, flanking_bp = 50){
  rr <- SummarizedExperiment::rowRanges(vcf)
  rr2 <- getKaspRange(rr[markers], flanking_bp = flanking_bp)
  fo <- GenomicRanges::findOverlaps(rr2, rr, type = "any")
  nhits <- S4Vectors::countLnodeHits(fo) - 1
  names(nhits) <- markers
  return(nhits)
}

# MAF <- function(vcf, markers){
#   
# }
# 
# LD <- function(){}

# Function to format SNPs for KASP assay, annotating flanking SNPs.
formatKasp <- function(vcf, markers, refgenome, flanking_bp = 50,
                       fixfn = fixTasselNames){
  rr <- rowRanges_correctedSeqnames(vcf, markers = rownames(vcf), fixfn = fixfn)
  rr2 <- getKaspRange(rr[markers], flanking_bp = flanking_bp)
  seq <- Rsamtools::scanFa(refgenome, rr2, as = "DNAStringSet")
  fo <- GenomicRanges::findOverlaps(rr2, rr, type = "any")
  kaspstart <- BiocGenerics::start(rr2)
  for(i in seq_along(markers)){
    toshift <- kaspstart[i] - 1
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
