jc <- function(set1, set2) {
  length(intersect(set1, set2)) / length(union(set1, set2))
} 


#' Function to calculate the bootstrap estimate of the Jaccard Index of two sets
#' 
#' @param set1 vector with all elements of first set
#' @param set2 vectpr with all elements of second set
#' @param num_bs number of bootstrap samples, default to 500
#' @param seed seed for random number generator, default to 123
#' 
#' @return a list with 3 entries
#'         1) jc_mean: the mean Jaccard index from the bootstrap samples
#'         2) jc_se: the bootstrap standard error for the Jaccard Index
#'         3) jc_samples: a vector with Jaccard Indices calculated from all bootstrap samples
jarccard_bs <- function(set1, set2, num_bs = 500, seed = 123) {
  
  set.seed(seed)
  
  union_set <- union(set1, set2)
  intersect_set <- intersect(set1, set2)
  
  samples <- lapply(1:num_bs, function(i) table(sample(union_set, length(union_set), replace = T)))
  samples <- sapply(samples, "[", union_set)
  rownames(samples) <- union_set
  samples[is.na(samples)] <- 0 
 
  jaccard <- colSums(samples[intersect_set, ]) / colSums(samples)
  
  res <- list(
    jc_mean = mean(jaccard),
    jc_se = sd(jaccard) / sqrt(length(jaccard)),
    jc_samples = jaccard
  )
  
  res
  
}


#' Function to calculate the bootstrap estimate of the Jaccard Index of two sets
#' 
#' @param set1 vector with all elements of first set
#' @param set2 vector with all elements of second set
#' @param n number of items to choose from set2
#' @param nsamples number of sampling runs, default to 500
#' @param seed seed for random number generator, default to 123
#' 
#' @return a list with 3 entries
#'         1) jc_mean: the mean Jaccard index from the bootstrap samples
#'         2) jc_se: the bootstrap standard error for the Jaccard Index
#'         3) jc_samples: a vector with Jaccard Indices calculated from all bootstrap samples
jarccard_resample <- function(set1, set2, n, nsampling = 500, seed = 123) {
  
  set.seed(seed)
  
  if (n >= length(set2)) {
    stop("Number of items to sample n has to be smaller than number of items in set 2!")
  } 
 
  union_set <- union(set1, set2)
  intersect_set <- intersect(set1, set2)
  set1_ex <- setdiff(set1, set2)
  
  samples <- lapply(1:nsampling, function(i) c(set1_ex, sample(set2, n, replace = F)))
  samples <- sapply(samples, function(xx) union_set %in% xx)
  rownames(samples) <- union_set

  jaccard <- colSums(samples[intersect_set, ]) / colSums(samples)
  
  res <- list(
    jc_mean = mean(jaccard),
    jc_sd = sd(jaccard),
    jc_se = sd(jaccard) / sqrt(length(jaccard)),
    jc_med = median(jaccard),
    jc_mad = mad(jaccard),
    jc_samples = jaccard
  )
  
  res
  
}


# Test with simulated data

## define sets and and calculate Jaccard index
set1 <- LETTERS[1:6]
set2_orig <- LETTERS[5:11]
jc_orig <- jc(set1, set2_orig)

## calculate resampled Jaccard index based on random samples from all letters
jc_resampling <- jarccard_resample(set1, set2 = LETTERS, 
                                    n = length(set2_orig), 
                                    nsampling = 10000, 
                                    seed = 123)

## get ECDF from samples JC indices
## and plot it
jc_cdf <- ecdf(jc_resampling$jc_samples)

plot(jc_cdf, main = "", xlab = 'sampled Jaccard indices', ylab = "empirical CDF")
lims <- par("usr")
points(jc_orig, jc_cdf(jc_orig), pch = 20, col = 'red')
segments(lims[1], jc_cdf(jc_orig), jc_orig, jc_cdf(jc_orig), col = 'red', lty = 2)
segments(jc_orig, jc_cdf(jc_orig), jc_orig, lims[3], col = 'red', lty = 2)

## calculate one-sided p-value from ECDF 
p <- 1-jc_cdf(jc_orig)
message("Resampling p-value: ", p)
