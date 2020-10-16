# Functions to assess a set of markers for how well it can differentiate accessions

# freqmat is a matrix of allele frequencies with markers in rows and populations in columns.

# Function to get expected heterozygosity, from a vector of allele frequencies,
# assuming biallelic markers.
Expected_het <- function(freq){
  return(1 - freq ^ 2 - (1 - freq) ^ 2)
}

# Get Jost's D from a matrix of allele frequencies.
JostD <- function(freqmat){
  meanfreq <- rowMeans(freqmat, na.rm = TRUE)
  Ht <- Expected_het(meanfreq)
  Hs_mat <- apply(freqmat, 2, Expected_het)
  Hs <- rowMeans(Hs_mat, na.rm = TRUE)
  return((Ht - Hs) / (1 - Hs))
}

# The geometric mean of pairwise Jost's D could be a good way to see if we have markers that distinguish all pops
DiffScore <- function(freqmat){
  npop <- ncol(freqmat)
  if(npop < 2) stop("Can't do population differentiation with less than two populations.")
  npairs <- sum(seq_len(npop - 1))
  p <- 1
  dmat <- matrix(NA_real_, nrow = nrow(freqmat), ncol = npairs,
                 dimnames = list(rownames(freqmat), NULL))
  for(i in 1:(npop - 1)){
    for(j in (i + 1):npop){
      dmat[,p] <- JostD(freqmat[,c(i, j)])
      p <- p + 1
    }
  }
  
  score <- exp(mean(log(colMeans(dmat, na.rm = TRUE))))
  return(score)
}
