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

# Get Pearson's correlation coefficient of phenotype with markers
# markers: a character vector naming markers of interest
# traits: a character vector indicating which trait goes with each marker
# traittab: a dataframe, with rows corresponding to columns of numgen,
# containing trait values.  Column names should match traits.
phenoCorr <- function(numgen, vcf, markers, traits, traittab,
                      samplecol = "DRS"){
  rr <- SummarizedExperiment::rowRanges(vcf)
  if(is.null(rr$paramRangeID)){
    stop("Need VCF with paramRangeID column from using GRanges to specify regions.")
  }
  if(!identical(rownames(numgen), names(rr))){
    stop("Rows in numgen and vcf need to match.")
  }
  if(!all(traits %in% colnames(traittab))){
    stop("Traits need to match column names")
  }
  if(length(markers) != length(traits)){
    stop("markers and traits must be same length (one trait per marker).")
  }
  if(!all(markers %in% rr$paramRangeID)){
    stop("Marker names not found in paramRangeID.")
  }
  if(is.null(traittab[[samplecol]]) ||
     !all(traittab[[samplecol]] %in% colnames(numgen))){
    stop("Make sure samplecol points to column with sample names, and that these match column names of numgen.")
  }
  numgen <- numgen[, traittab[[samplecol]]]
  
  corlist <- lapply(seq_along(markers),
                    function(i){
                      theserows <- which(rr$paramRangeID == markers[i])
                      sapply(names(rr[theserows]),
                      function(x){
                        cor(traittab[[traits[i]]],
                            numgen[x,],
                            use = "pairwise.complete.obs", method = "pearson")
                      })
                    })
  corlist <- as(corlist, "NumericList")
  names(corlist) <- markers
  return(corlist)
}

# Determine which allele is positively associated with the trait at each SNP
# corrval is a vector of Pearson's correlation coefficient, matching markers.
whichAllele <- function(markers, vcf, corrval){
  if(length(markers) != length(corrval)){
    stop("markers and corrval must be same length")
  }
  rr <- SummarizedExperiment::rowRanges(vcf)
  out <- rr[markers,]$REF
  altfill <- rr[markers[which(corrval > 0)],]$ALT
  if(all(lengths(altfill) == 1)){
    altfill <- unlist(altfill)
  } else {
    altfill <- sapply(altfill, function(x) x[1])
  }
  out[which(corrval > 0)] <- altfill
  names(out) <- markers
  return(out)
}
whichAlleleList <- function(numgen, vcf, markers, traits, traittab,
                        samplecol = "DRS"){
  cr <- phenoCorr(numgen, vcf, markers, traits, traittab, samplecol)
  al_list <- lapply(cr,
                     function(cors){
                       whichAllele(names(cors), vcf, cors)
                     })
  al_list <- DNAStringSetList(al_list)
  return(al_list)
}
