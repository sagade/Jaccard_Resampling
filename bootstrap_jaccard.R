#' Helper function to calculate the Jaccard index from two sets
#' 
#' @param set1 vector with elements for set 1
#' @param set2 vector with elements for set2
#' 
#' @return the Jaccard index for the two sets
jaccard <- function(set1, set2) {
  length(intersect(set1, set2)) / length(union(set1, set2))
} 


#' Helper function to calculate a p-value based on a statistic and a sampling
#' based null distribution
#' 
#' @param x the statistic to be tested: only one value
#' @param y vector with the sampling based null distribution
#' @param alternative a character string specifying the alternative hypothesis
#' @param plot boolean indicating whether the empirical distriubtion function
#'               used to calculate the p-value should be plotted
#' @param ... additional argument to plot for the cdf 
#' 
#' @return the p-value from the test
resample_test <- function(x, y, 
                          alternative = c("two.sides", "less", "greater"), 
                          plot = F,
                          ...) {
  
  alternative = match.arg(alternative)
  x <- x[1]
  
  cdf <- ecdf(y)
  
  if (plot) {
    plot(cdf, ...)
    lims <- par("usr")
    points(x, cdf(x), pch = 20, col = 'red')
    segments(lims[1], cdf(x), x, cdf(x), col = 'red', lty = 2)
    segments(x, cdf(x), x, lims[3], col = 'red', lty = 2)
  }
  
  ## calculate p-value
  p <- switch(alternative,
              two.sided = 2 * min(cdf(x), 1 - cdf(x)),
              less = cdf(x),
              greater = 1 - cdf(x))
  
  p
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
jaccard_bs <- function(set1, set2, num_bs = 500, seed = 123) {
  
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
#'         1) jc_mean: the mean Jaccard index from all re-sampling runs
#'         2) jc_sd: the accard index sd from all re-sampling runs
#'         3) jc_se: the resampling standard error for the Jaccard Index
#'         4) jc_med: the median Jaccard index from all re-sampling runs
#'         5) jc_mad: the index MAD from all re-sampling runs
#'         6) jc_samples: a vector with Jaccard Indices calculated from all bootstrap samples
jaccard_resample <- function(set1, set2, n, nsampling = 500, seed = 123) {
  
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
jc_orig <- jaccard(set1, set2_orig)

## calculate resampled Jaccard index based on random samples from all letters
jc_resampl<- jaccard_resample(set1, set2 = LETTERS, 
                              n = length(set2_orig), 
                              nsampling = 10000, 
                              seed = 353)

## test the Jaccard index aggainst the sampled distribution of random Jaccard
## indices. Test if the original index is higher than random sampling
p <- resample_test(x = jc_orig, y = jc_resampl$jc_samples, 
                   alternative = "greater",
                   plot = T, 
                   main = "Jaccard index ECDF", 
                   xlab = "random sample indices", 
                   ylab = "empirical distribution function")
message("p-value from re-sampling based test: ", signif (p, 3))
