
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
