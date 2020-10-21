# Evaluate sets of markers near a QTL.

# numgen: numeric genotype matrix with SNPs in rows and individuals in columns.
# vcf: a VCF object containing SNP metadata, with rows corresponding to those in numgen.

LD <- function(numgen, vcf){
  rr <- SummarizedExperiment::rowRanges(vcf)
  if(is.null(rr$paramRangeID)){
    stop("Need VCF with paramRangeID column from using GRanges to specify regions.")
  }
  if(!identical(rownames(numgen), names(rr))){
    stop("Rows in numgen and vcf need to match.")
  }
  targets <- levels(rr$paramRangeID)
  ldlist <- lapply(targets,
                   function(snp){
                     theserows <- which(rr$paramRangeID == snp)
                     snprow <- theserows[names(rr[theserows]) == sub("_[ACGT]+$", "", snp)]
                     sapply(names(rr[theserows]),
                            function(x){
                              if(sd(numgen[x,]) == 0){
                                return(0)
                              } else {
                                return(cor(numgen[snprow,], numgen[x,],
                                           use = "pairwise.complete.obs", method = "pearson") ^ 2)
                              }
                            })
                   })
  ldlist <- as(ldlist, "NumericList")
  names(ldlist) <- targets
  return(ldlist)
}
