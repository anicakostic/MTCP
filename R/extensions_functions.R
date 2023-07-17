#' Hoang & Dickhaus randomised p-values calculation
#'
#' This function produces a sequence of randomised p-values
#' using the method proposed in
#' Hoang, A.-T., & Dickhaus, T. (2020).
#' On the usage of randomized p-values in the Schweder–Spjøtvoll estimator.
#' Annals of the Institute of Statistical Mathematics, 74, 289–319.
#'
#' @param p_seq  A sequence of p-values
#' @param lambda Parameter for Storey's estimator and for randomisation
#' @param c Threshold parameter for calculating the randomised p-values.
#' If NULL, the threshold is calculated by optimising a function defined in Eq. 8 in the paper.
#'
#' @return A list containing
#' * `random_p_seq` the sequence of randomised p-values
#'
#' * `pi1_est` the corresponding Storey's estimator of the false null proportion based on the randomised sequence and for given lamdba
#'
#' @export
hd_est=function(p_seq,lambda=0.5,c=NULL)
{
  p_seq=sort(p_seq)
  if(is.null(c))
  {
    eval_points=sort(c(p_seq,p_seq[p_seq<lambda]/lambda))
    g_seq=rep(0,length(eval_points))
    for(i in 1:length(eval_points))
    {
      g_seq[i]=lambda*sum(p_seq>=eval_points[i])+ sum(p_seq/lambda<eval_points[i])
    }
    c=eval_points[which.max(g_seq)]
  }
  random_p_seq=runif(length(p_seq))*(p_seq>=c) + p_seq/c*(p_seq<c)
  res=NULL
  res$random_p_seq=random_p_seq
  res$pi1_est=storey_pi1est(lambda,random_p_seq)$est
  return(res)

}


#' Adaptive Benjamini-Hochberg procedure
#'
#' This function calculates the FDR and the power of the adaptive BH procedure
#' for a given value of the proportion estimator
#' (Benjamini, Y., & Hochberg, Y. (2000). On the Adaptive Control of the False
#' Discovery Rate in Multiple Testing with Independent Statistics.
#'  Journal of Educational and Behavioral Statistics, 25(1), 60.)
#'
#'
#' @param p_seq  A sorted sequence of p-values
#' @param pi0 Estimated true null proportion
#' @param alpha FDR control level
#' @param sig A vector matching p-values, zero if the corresponding p-value
#' is true null, and 1 if it is false null
#'
#' @return A list containing
#' * `fdr` the realised proportion of false discoveries among all discoveries
#'
#' * `power` the proportion of true discoveries among all false null p-values
#'
#' @export
bh_adaptive <- function(p_seq, pi0, alpha, sig) {
  n <- length(p_seq)
  res <- NULL
  bh_est <- ifelse(sum(p_seq < (1:n) / n * alpha / pi0) == 0, 0,
                   max(which(p_seq <= (1:n) / n * alpha / pi0)))
  res$est <- bh_est
  if (bh_est == 0) {
    res$fdr <- 0
    res$power <- 0
  } else {
    res$fdr <- sum(1 - sig[1:bh_est]) / bh_est
    res$power <- sum(sig[1:bh_est]) / sum(sig)
  }
  return(res)
}


