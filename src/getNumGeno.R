# Function to get a numeric SNP matrix from a VCF.
# Biallelic SNPs assumed.

getNumGeno <- function(file, ranges = NULL){
  if(is.null(ranges)){
    param <- VariantAnnotation::ScanVcfParam()
  } else {
    param <- VariantAnnotation::ScanVcfParam(which = ranges)
  }
  gtmat <- VariantAnnotation::readGT(file, param = param)
  ct <- stringr::str_count(gtmat, "1")
  ms <- stringr::str_detect(gtmat, "\\.")
  ct[ms] <- NA_integer_
  out <- matrix(ct, nrow = nrow(gtmat), ncol = ncol(gtmat),
                dimnames = dimnames(gtmat))
  
  # Add sample names if missed
  if(is.null(colnames(gtmat))){
    hdr <- VariantAnnotation::scanVcfHeader(file)
    sam <- VariantAnnotation::samples(hdr)
    sam <- sub("^.*/", "", sam)
    sam <- sub("\\..*", "", sam)
    if(!anyDuplicated(sam)){
      colnames(out) <- sam
    }
  }
  
  return(out)
}

# Function to make a distance matrix for identifying clones.
# Downsample SNPs to speed up processing.
interIndividualIBS <- function(numgen, minMAF = 0.05, nsnp = 1e4, ploidy = 2){
  freq <- rowMeans(numgen, na.rm = TRUE) / ploidy
  numgen <- numgen[freq >= minMAF,]
  if(nsnp > nrow(numgen)) nsnp <- nrow(numgen)
  randsnp <- sample(nrow(numgen), nsnp)
  ind_dist <- dist(t(numgen[randsnp,]), method = "manhattan")
  ncp <- nsnp * ploidy
  ibs <- (ncp - ind_dist) / ncp
  return(ibs)
}

# remove clonal individuals from dataset
removeClones <- function(numgen, ibs, threshold = 0.97){
  toremove <- c()
  clone_index <- which(as.matrix(ibs) > threshold, arr.ind = TRUE)
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
