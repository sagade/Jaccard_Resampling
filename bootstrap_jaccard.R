#' @name bootsttrap_jaccard.R
#' @version 0.1.1
#' @author Stephan Gade


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
#' @param set2 vector with all elements of second set
#' @param n_samples number of bootstrap samples, default to 500
#' @param seed seed for random number generator, default to 123
#' 
#' @return a list with 3 entries
#'         1) jc_mean: the mean Jaccard index from the bootstrap samples
#'         2) jc_se: the bootstrap standard error for the Jaccard Index
#'         3) jc_samples: a vector with Jaccard Indices calculated from all bootstrap samples
jaccard_bs <- function(set1, set2, n_samples = 500, seed = 123) {
  
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


#' Function to calculate a Null distribution for the Jaccard index based on
#' random sampling from a Null set. 
#' 
#' @details The function can be used to calculate the Jaccard index between 
#' a fixed set1 and random samplings from a set of items (the Null set) to
#' calculate a Null distribution of random Jaccard indices. This Null
#' distribution can then be used to test a Jaccard index of set1 against a real
#' set2.
#'  
#' @param set1 vector with all elements of first set
#' @param null_set vector with all elements of the null set
#' @param n_items number of items to choose from the null set
#' @param n_samples number of sampling runs, default to 500
#' @param seed seed for random number generator, default to 123
#' 
#' @return a list with following entries
#'         1) jc_mean: the mean Jaccard index from all re-sampling runs
#'         2) jc_sd: the accard index sd from all re-sampling runs
#'         3) jc_se: the resampling standard error for the Jaccard Index
#'         4) jc_med: the median Jaccard index from all re-sampling runs
#'         5) jc_mad: the index MAD from all re-sampling runs
#'         6) jc_samples: a vector with Jaccard Indices calculated from all bootstrap samples
jaccard_resample <- function(set1, null_set, n_items, n_samples = 500, seed = 123) {
  
  set.seed(seed)
  
  if (n_items >= length(null_set)) {
    stop("Number of items to sample 'n_items' has to be smaller than number of items in null set!")
  } 
 
  # calculate the union and intersect of the set of interest (set1) 
  # and the Null set
  union_set <- union(set1, null_set)
  intersect_set <- intersect(set1, null_set)
  
  # Draw random samples from the Null set, number of samples to draw is 
  # specified by n_samples. We will draw without replacement.
  sampling_list <- lapply(1:n_samples, function(i) sample(null_set, n_items, replace = F))
   
  # Calculate two matrices: 
  #   - s_union: matrix of 0 and 1 with rows the items in the union set and
  #              columns the sampling runs containing a 1 if an item ie either
  #              part of set1 or of the corresponding random sample from the 
  #              Null set, 0 otherwise
  #   - s_intersect: matrix of 0 and 1 with rows corresponding to the items
  #                  in the intersect between set1 and the Null set and columns
  #                  to the sampling runs. Here, 1 indicates an item from the
  #                  intersect set was present in the random sample from the Null
  #                  set and 0 otherwise.
  s_union     <- sapply(sampling_list, function(xx) union_set %in% c(set1, xx))
  s_intersect <- sapply(sampling_list, function(xx) intersect_set %in% xx)
 
  # Building up the two matrices from the list of random samplings from the Null
  # set, the Jaccard index between set1 and each random sample can now simply be
  # calculated by the element-wise division of two column sums. The first column
  # sum is calculated for the matrix holding only items from the intersect set
  # returning the intersection of set1 and all sampling runs from the Null set. 
  # Consequently the second column sum returns the union of set1 and all samples
  # from the Null set. So, the element-wise divsion of these tow column sums 
  # gives the vectors of Jaccard indices from all sampling runs.
  jaccard <- colSums(s_intersect) / colSums(s_union)
  
  # Assemble the result list with the samples Jaccard indices and some
  # summary statistics.
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
# set1 <- LETTERS[1:6]
set1 <- c("1", "2", "3", "A", "B", "E", "F", "R")
set2 <- LETTERS[5:11]
jc_orig <- jaccard(set1, set2)

## calculate resampled Jaccard index based on random samples from all letters
jc_resampl<- jaccard_resample(set1 = set1, 
                              null_set = LETTERS, 
                              n_items = length(set2), 
                              n_samples = 10000, 
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



# Run time test
library(microbenchmark)
n_samples <- 1000
microbenchmark(
  jc_standard = ({
    set.seed(123)
    jc <- sapply(1:n_samples, function(i) jaccard(set1, sample(LETTERS, length(set2), replace = F)))
  }),
  jc_resampl = jaccard_resample(set1, null_set = LETTERS, n_items = length(set2), n_samples = n_samples, seed = 123)
)
