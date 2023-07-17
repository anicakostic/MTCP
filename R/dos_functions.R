#' DOS sequence calculation
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param p_seq  A sequence of p-values
#' @param alpha  DOS sequence power parameter
#' @param exc  Number of excluded values from the beginning
#' @return A list containing
#' * `dos_seq` the DOS sequence
#'
#' * `cp_loc` the location of the maximum of the DOS sequence
#' @export
dos_fun <- function(p_seq, alpha = 1, exc = 0) {
  p_seq <- sort(p_seq)
  n <- length(p_seq)
  t <- (p_seq[2 * (1:(n/2))] - 2 * p_seq[1:(n/2)])/((1:(n/2))/n)^alpha

  ret <- NULL
  ret$dos_seq <- t
  ret$cp_loc <- ifelse(max(t) > 0, which.max(t), which.max(t[1:(n/2)]))
  return(ret)
}

#' Storey's estimator function
#'
#' This function calculates the Storey's estimator
#' of the false null proportion, given the tuning parameter values
#' and the sequence of p-values
#'
#' @param lambda  A vector of tuning parameter values
#' @param p_seq  A sequence of p-values
#' @param mod Whether a modified Storey's estimator is to be calculated (see Blanchard, G.,
#'  & Roquain, É. (2009). Adaptive false discovery rate control under independence and dependence.
#'   Journal of Machine Learning Research, 10, 2837–2871.)
#'
#' @return A list containing
#' * `lambda` Input vector lambda
#'
#'* `est` The corresponding vector of estimates
#' @export
storey_pi1est <- function(lambda, p_seq, mod = FALSE) {
  p_seq <- sort(p_seq)
  n <- length(p_seq)
  N <- length(lambda)
  est <- rep(0, N)
  for (i in 1:N) {
    est[i] <- (sum(p_seq <= lambda[i])/n - lambda[i] - mod/n)/(1 - lambda[i])
    est[i] <- max(0, est[i])
  }
  ret <- NULL
  ret$lambda <- lambda
  ret$est <- est

  return(ret)

}




#' Generating p-values from the Gaussian model
#'
#' This function generates one-sided p-values from the Gaussian mean testing
#' where the p-values are calculated assuming that the statistics
#' have standard normal distribution under the null
#' and a mean greater than zero under the alternative.
#'
#' @param n  The length of the output vector of p-values
#' @param null_mean  A vector of mean values under the null, or a single
#' value if it is the same for all variables
#' @param alt_mean A vector of mean values under the alternative, or a single
#' value if it is the same for all variables
#' @param pi1 A proportion of false null p-values
#' @param rho A correlation parameter between 0 and 1, such that the correlation between
#' each pair of p-values is equal to rho
#'
#' @return A vector of p-values
#'
#' @export
sim_pval <- function(n, null_mean, alt_mean, pi1, rho) {
  if (length(null_mean) == 1)
    null_mean <- rep(null_mean, n - floor(n * pi1))
  if (length(alt_mean) == 1)
    alt_mean <- rep(alt_mean, floor(n * pi1))
  if (rho > 0)
    x <- sqrt(rho) * rnorm(1) + sqrt(1 - rho) * rnorm(n, c(alt_mean, null_mean)) else x <- rnorm(n, c(alt_mean, null_mean))

    return(1 - pnorm(x))
}


