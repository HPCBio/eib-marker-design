# Function to get a numeric SNP matrix from a VCF.
# Biallelic SNPs assumed.

getNumGeno <- function(file){
  gtmat <- VariantAnnotation::readGT(file)
  ct <- stringr::str_count(gtmat, "1")
  ms <- stringr::str_detect(gtmat, "\\.")
  ct[ms] <- NA_integer_
  out <- matrix(ct, nrow = nrow(gtmat), ncol = ncol(gtmat),
                dimnames = dimnames(gtmat))
  return(out)
}

# Function to make a distance matrix for identifying clones.
# Downsample SNPs to speed up processing.
interIndividualDist <- function(numgen, minMAF = 0.05, nsnp = 1e4, ploidy = 2){
  freq <- rowMeans(numgen, na.rm = TRUE) / ploidy
  numgen <- numgen[freq >= minMAF,]
  randsnp <- sample(nrow(numgen), nsnp)
  ind_dist <- dist(t(numgen[randsnp,]))
  return(ind_dist)
}

# remove clonal individuals from dataset
removeClones <- function(numgen, ind_dist, threshold = 10){
  toremove <- c()
  clone_index <- which(as.matrix(ind_dist) < threshold, arr.ind = TRUE)
  clone_index <- clone_index[clone_index[,1] != clone_index[,2],]
  while(nrow(clone_index) > 0){
    thisind <- max(clone_index[,2])
    toremove <- c(toremove, thisind)
    clone_index <-
      clone_index[clone_index[,1] != thisind & clone_index[,2] != thisind, , drop = FALSE]
  }
  if(length(toremove) > 0){
    return(numgen[, -toremove])
  } else {
    return(numgen)
  }
}
