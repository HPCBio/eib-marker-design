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
