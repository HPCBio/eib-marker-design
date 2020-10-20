# Functions to assess a set of markers for how well it can differentiate accessions

# freqmat is a matrix of allele frequencies with markers in rows and populations in columns.

# Function to get expected heterozygosity, from a vector of allele frequencies,
# assuming biallelic markers.
Expected_het <- function(freq){
  return(1 - freq ^ 2 - (1 - freq) ^ 2)
}

# Get Jost's D from a matrix of allele frequencies.
# Global Jost's D is calculated.  To get pairwise Jost's D, do two columns at a time.
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

# Simulated annealing algorithm to find an ideal marker set.
# nSNP: The number of SNPs to have in the set.
# within_pop_weight and between_pop_weight indicate how important it is to have
# variation within and between populations, respectively.
# T0: Starting temperature.  The higher this is, the more the algorithm will
# randomly explore the solution space.
# reps: Number of sets to test at each temperature.  More reps means a longer
# computational time but more thorough exploration of solutions.
# rho: The factor by which temperature decreases with each round.  The higher
# this is, the more the algorithm will randomly explore the solution space.
# maxrounds: The maximum number of rounds of simulated annealing to perform
# before stopping.  The algorithm will also stop if no new solutions are found
# within a round.
# min_dist_bp: The minimum distance, in basepairs, between markers in the set.
findMarkerSet <- function(freqmat, nSNP = 20,
                          within_pop_weight = 1, between_pop_weight = 0.5,
                          T0 = 0.2, reps = 500, rho = 0.95,
                          maxrounds = 100,
                          min_dist_bp = 1e5){
  # setup
  totSNP <- nrow(freqmat)
  Temp <- T0
  temps_used <- numeric(maxrounds)
  scores_by_round <- numeric(maxrounds)
  
  # function to generate a score for a set of SNPs
  scorefn <- function(snps){
    num <- log(DiffScore(freqmat[snps,])) * between_pop_weight +
      log(DivScore(freqmat[snps,])) * within_pop_weight
    den <- between_pop_weight + within_pop_weight
    return(exp(num / den))
  }
  
  # function to determine if a SNP is too physically close to any already in the set
  chromosomes <- sub("_[[:digit:]]+$", "", rownames(freqmat))
  positions <- as.integer(sub("^.+_", "", rownames(freqmat)))
  
  snp_pos_ok <- function(snp, set){
    samechr <- set[chromosomes[set] == chromosomes[snp]]
    return(all(abs(positions[samechr] - positions[snp]) >= min_dist_bp))
  }
  
  # initial set of SNPs
  thisset <- integer(nSNP)
  for(s in seq_len(nSNP)){
    snp <- sample(totSNP, 1)
    while(!snp_pos_ok(snp, thisset[seq_len(s - 1)])){
      snp <- sample(totSNP, 1)
    }
    thisset[s] <- snp
  }
  thisscore <- scorefn(thisset)
  bestset <- thisset
  bestscore <- thisscore
  
  # simulated annealing
  for(i in seq_len(maxrounds)){
    message(paste("Simulated annealing round", i))
    switched <- FALSE
    for(j in seq_len(reps)){
      # remove one SNP
      newset <- sample(thisset, nSNP - 1)
      # find a new SNP and confirm it is not too close to any in the set
      repeat{
        newsnp <- sample(seq_len(totSNP)[-newset], 1)
        if(snp_pos_ok(newsnp, newset)) break
      }
      newset <- c(newset, newsnp)
      # score the set of SNPs by diversity and differentiation
      newscore <- scorefn(newset)
      # determine whether to move to this new solution
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
    Temp <- Temp * rho
    if(!switched) break
  }
  
  return(list(Set = rownames(freqmat)[bestset],
              Temperatures = temps_used[1:i],
              Scores = scores_by_round[1:i]))
}
