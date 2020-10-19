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

# Geometric mean of expected heterozygosity across populations
DivScore <- function(freqmat){
  hmat <- Expected_het(freqmat)
  score <- exp(mean(log(colMeans(hmat, na.rm = TRUE))))
  return(score)
}

# Simulated annealing algorithm to find an ideal marker set
findMarkerSet <- function(freqmat, nSNP = 20,
                          within_pop_weight = 1, between_pop_weight = 1,
                          T0 = 0.5, reps = 100,
                          maxrounds = 100){
  totSNP <- nrow(freqmat)
  thisset <- sample(totSNP, nSNP)
  thisscore <- DiffScore(freqmat[thisset,]) ^ between_pop_weight *
    DivScore(freqmat[thisset,]) ^ within_pop_weight
  bestset <- thisset
  bestscore <- thisscore
  Temp <- T0
  temps_used <- numeric(maxrounds)
  scores_by_round <- numeric(maxrounds)
  step <- T0 / maxrounds
  
  for(i in seq_len(maxrounds)){
    message(paste("Simulated annealing round", i))
    switched <- FALSE
    for(j in seq_len(reps)){
      newset <- sample(thisset, nSNP - 1)
      newset <- c(newset,
                  sample(seq_len(totSNP)[-newset], 1))
      newscore <- DiffScore(freqmat[newset,]) ^ between_pop_weight *
        DivScore(freqmat[newset,]) ^ within_pop_weight
      if(newscore > thisscore ||
         runif(1) < Temp - (thisscore - newscore)){
        thisset <- newset
        thisscore <- newscore
        switched <- TRUE
        if(thisscore > bestscore){
          bestset <- thisset
          bestscore <- thisscore
        }
      }
    }
    temps_used[i] <- Temp
    scores_by_round[i] <- bestscore
    Temp <- Temp - step
    if(!switched) break
  }
  
  return(list(Set = rownames(freqmat)[bestset],
              Temperatures = temps_used[1:i],
              Scores = scores_by_round[1:i]))
}
