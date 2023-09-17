#' Jiang & Doerge estimator calculation
#'
#' This function calculates the true null proportion estimates
#' using the method proposed in Jiang, H., & Doerge, R. W. (2008)
#'
#' @param p_seq  A sequence of p-values
#' @param N  The number of bootstrap samples
#'
#' @return Estimated true null proportion
#' @references Jiang, H., & Doerge, R. W. (2008).
#' Estimating the Proportion of True Null Hypotheses
#' for Multiple Comparisons. Cancer Informatics, 6, 25--32
#' https://doi.org/10.1177/117693510800600001
#' @examples
#' p_seq <- sim_pval(100, 0, 3, 0.1, 0)
#' JD_est(p_seq, 1000)
#' @export
JD_est <- function(p_seq, N) {
  n <- length(p_seq)
  B <- c(5, 10, 20, 50, 100)
  pi0_estB <- rep(0, length(B))
  mse <- rep(0, length(B))
  for (j in 1:length(B)) {

    t <- seq(0, 1, by = 1/B[j])
    ns <- nb <- rep(0, length(t) - 1)
    for (i in 1:(length(t) - 1)) {
      ns[i] <- sum(p_seq >= t[i] & p_seq < t[i + 1])
      nb[i] <- sum(p_seq >= t[i])
    }
    i_st <- min(which(ns <= nb/(B[j] - (1:length(ns)) + 1)))
    pi0_estB[j] <- mean((nb/(1 - t[-length(t)])/length(p_seq))[(i_st -
                                                                  1):(B[j])])

    pi0star_b <- rep(0, N)
    for (b in 1:N) {
      bootst_p <- sample(p_seq, n, replace = TRUE)
      for (i in 1:(length(t) - 1)) {
        ns[i] <- sum(bootst_p >= t[i] & bootst_p < t[i + 1])
        nb[i] <- sum(bootst_p >= t[i])
      }
      i_st <- min(which(ns <= nb/(B[j] - (1:length(ns)) + 1)))

      pi0star_b[b] <- mean((nb/(1 - t[-length(t)])/length(bootst_p))[(i_st -
                                                                        1):(B[j])])

    }

    mse[j] <- mean((pi0star_b - pi0_estB[j])^2)

  }

  whichB <- which.min(mse)
  return(min(pi0_estB[whichB], 1))


}

#' Meinshausen & Rice estimator calculation
#'
#' This function calculates the false null proportion estimates
#' using the method proposed in Meinshausen, N., & Rice, J. (2006)
#'
#'
#'
#' @param p_seq  A sequence of p-values
#' @param alpha  For the confidence interval of level 1-alpha
#'
#' @return Estimated false null proportion
#' @references "Meinshausen, N., & Rice, J. (2006). Estimating the proportion of false null
#'  hypotheses among a large number of independently tested hypotheses.
#'  Annals of Statistics, 34(1), 373–393" https://doi.org/10.1214/009053605000000741
#' @examples
#' p_seq <- sim_pval(100, 0, 3, 0.1, 0)
#' MR_est(p_seq)
#' @export
MR_est <- function(p_seq, alpha = 0.05) {
  n <- length(p_seq)
  eps <- 0.0001
  a <- sqrt(2 * n * log(log(n)))
  b <- 2 * log(log(n)) + log(log(log(n))) / 2 - log(4 * pi) / 2
  beta <- (-log(-log(1 - alpha)) + b) / a

  t <- seq(1 / n, 1 - 1 / n, by = eps)

  ecd <- ecdf(p_seq)
  seq <- na.exclude((ecd(t) - t - beta * sqrt(t * (1 - t))) / (1 - t))
  return(max(floor(max(seq) * n), 0) / n)
}

#' Langaas, Lindqvist & Ferkingstad estimator calculation
#'
#' Estimator of the true null proportion at the longest constant interval in the Grenander estimator,
#' as proposed in Langaas, M., Lindqvist, B. H., & Ferkingstad, E. (2005).
#' Estimating the Proportion of True Null Hypotheses, with Application to DNA Microarray Data.
#'  Journal of the Royal Statistical Society. Series B (Statistical Methodology), 67(4), 555–572.
#'
#' @param p_seq  A sequence of p-values
#' @references Langaas, M., Lindqvist, B. H., & Ferkingstad, E. (2005).
#' Estimating the Proportion of True Null Hypotheses, with Application to DNA Microarray Data.
#'  Journal of the Royal Statistical Society. Series B (Statistical Methodology), 67(4), 555–572.
#'  https://doi.org/10.1111/j.1467-9868.2005.00515.x
#' @return Estimated true null proportion
#' @examples
#' p_seq <- sim_pval(100, 0, 3, 0.1, 0)
#' LLF_Grenander_est(p_seq)
#'
#' @export
LLF_Grenander_est <- function(p_seq) {
  res <- fdrtool::grenander(ecdf(p_seq))
  f_le_1 <- which(res$f.knots <= 1)
  max_idx <- which.max(diff(res$x.knots[f_le_1]))

  if (length(max_idx) > 0) {
    return(min(1, res$f.knots[f_le_1[max_idx]]))
  }

  return(NULL)
}


#' Broberg's MGF estimator calculation
#'
#' This function calculates the estimates of the true null proportion using the method 'MGF'
#' proposed in Broberg, P. (2005). A comparative review of estimates of the proportion
#' unchanged genes and the false discovery rate. BMC Bioinformatics, 6, 199. The code belongs to Per Broberg and
#' it was a part of the now deprecated package 'SAGx' in BioConductor. 'SAGx' can be downloaded from
#' https://www.bioconductor.org/packages//2.10/bioc/html/SAGx.html
#'
#'
#' @param p_seq  A sequence of p-values
#' @return Estimated true null proportion
#' @references  Broberg, P. (2005). A comparative review of estimates of the proportion
#' unchanged genes and the false discovery rate. BMC Bioinformatics, 6, 199.
#' https://doi.org/10.1186/1471-2105-6-199
#' @examples
#' p_seq <- sim_pval(100, 0, 3, 0.1, 0)
#' MGF_est(p_seq)
#'
#' @export
MGF_est <- function(p_seq){

  R.s <- function(p_seq, s){
    sum(exp(s*p_seq))/length(p_seq)
  }
  ss <- seq(0.01, 0.1, by = 0.001)
  R.s.s <- vector(length = length(ss) )
  R.s.s <- apply(as.data.frame(ss), 1, function(x) R.s(p_seq, x))
  g <- function(x) (exp(x) - 1)/x
  r.n <- vector(length = length(ss))
  g.n <- g(ss)

  r.n.est <- function(x){
    temp.vect <- vector(length =length(ss))
    temp.vect[1] <- x
    for(i in 2:length(ss)) temp.vect[i] <- (R.s.s[i]*g.n[i-1] - R.s.s[i-1]*g.n[i] + temp.vect[i-1]*(g.n[i] - R.s.s[i]))/(g.n[i-1] - R.s.s[i-1])
    temp.vect
  }

  target <- function(x) sd((R.s.s - r.n.est(x))/(g.n - r.n.est(x)))/mean((R.s.s - r.n.est(x))/(g.n - r.n.est(x)))

  optimal.value <- optimize(target, c(1,min(R.s.s)))$minimum

  r.n <- r.n.est(optimal.value)

  p0.0 <- min(mean((R.s.s-r.n)/(g.n-r.n)),1)

  return(p0.0)

}

#' Broberg's PRE estimator calculation
#'
#' This function calculates the estimates of the true null proportion using the method 'MGF'
#' proposed in Broberg, P. (2005). A comparative review of estimates of the proportion
#' unchanged genes and the false discovery rate. BMC Bioinformatics, 6, 199. The code belongs to Per Broberg and
#' it was a part of the now deprecated package 'SAGx' in BioConductor. 'SAGx' can be downloaded from
#' https://www.bioconductor.org/packages//2.10/bioc/html/SAGx.html
#'
#'
#' @param p_seq  A sequence of p-values
#'
#' @return Estimated true null proportion
#' @references Broberg, P. (2005). A comparative review of estimates of the proportion
#' unchanged genes and the false discovery rate. BMC Bioinformatics, 6, 199.
#' https://doi.org/10.1186/1471-2105-6-199
#' @examples
#' p_seq <- sim_pval(100, 0, 3, 0.1, 0)
#' PRE_est(p_seq)
#' @export
PRE_est <- function(p_seq)
{

  R.s <- function(p_seq, s){
    sum(exp(s*p_seq))/length(p_seq)
  }
  ss <- seq(0.01, 0.1, by = 0.001)
  R.s.s <- vector(length = length(ss) )
  R.s.s <- apply(as.data.frame(ss), 1, function(x) R.s(p_seq, x))
  g <- function(x) (exp(x) - 1)/x
  r.n <- vector(length = length(ss))
  g.n <- g(ss)

  counts <- hist(p_seq, plot = FALSE, breaks = 'FD')

  K <- length(counts$mids)
  M <- matrix(ncol = K, nrow = K)
  lambda <- 0.1
  for(i in 1:K) for(j in 1:K) M[i,j] <- dnorm((counts$mids[i]-counts$mids[j])/lambda)/lambda
  row.sum <- matrix(rep(1, ncol(M))%*%M)
  M <- sweep(M, 1, row.sum, FUN = "/")
  mu <- M%*%counts$counts
  x1 <- counts$mids;x2  <- x1^2;x3 <- x1^3;x4 <- x1^4
  glm.test <- glm(counts$counts ~ x1 + x2 + x3 + x4, family  = poisson(), offset = log(mu))
  spl <- smooth.spline(counts$mids, glm.test$fitted, df = 4)
  x1 <- p_seq;x2  <- x1^2;x3 <- x1^3;x4 <- x1^4
  fs <- length(counts$counts)*predict(spl, p_seq)$y/sum(counts$counts)


  p0.1 <- min(fs)

  p0.1 <- min(1, quantile(fs, probs = c(0.25)))

  if(p0.1 < 1){
    r.n[1] <- (R.s.s[1] - p0.1 * g.n[1])/(1 - p0.1)
    for(i in 2:length(ss)) r.n[i] <- (R.s.s[i]*g.n[i-1] - R.s.s[i-1]*g.n[i] + r.n[i-1]*(g.n[i] - R.s.s[i]))/(g.n[i-1] - R.s.s[i-1])
  }

  return(p0.1)
}


#' Hwang et al. Slope Difference estimator calculation
#'
#' This function calculates the estimates of the true null proportion using the
#' Slope Difference estimator proposed in Hwang, Y. T. et al. (2014).
#'
#' @param p_seq  A sequence of p-values
#'
#' @return Estimated true null proportion
#' @references Hwang, Y. T., Kuo, H. C., Wang,
#'  C. C., & Lee, M. F. (2014). Estimating the number of true null hypotheses
#'   in multiple hypothesis testing. Statistics and Computing, 24(3), 399–416.
#'   https://doi.org/10.1007/s11222-013-9377-5
#' @examples
#' p_seq <- sim_pval(100, 0, 3, 0.1, 0)
#' Hwang_est(p_seq)
#' @export
Hwang_est <- function(p_seq) {
  n <- length(p_seq)
  p_seq <- sort(p_seq)
  k <- max(which(p_seq < 1/2))
  sd_seq <- (1 - p_seq) / (n + 1 - 1:n) - p_seq / (1:n)
  J <- which.max(sd_seq[1:k])
  return(min((n + 1 - J) / (1 - p_seq[J]) - 1, n) / n)
}



#' Turkheimer et al. Change-point based estimator
#'
#' This function calculates the estimates of the true null proportion using the
#' Change-point based estimator proposed in
#'  Turkheimer, F. E. et al. (2001)
#'
#' @param p_seq  A sequence of p-values
#' @param alpha Significance level for iterative uniformity testing
#'
#' @return Estimated true null proportion
#' @references Turkheimer, F. E., Smith, C. B., & Schmidt, K. (2001).
#'  Estimation of the number of “true” null hypotheses in multivariate
#'   analysis of neuroimaging data. NeuroImage, 13(5), 920–930.
#'   https://doi.org/10.1006/nimg.2001.0764
#' @examples
#' p_seq <- sim_pval(100, 0, 3, 0.1, 0)
#' turkheimer_est(p_seq, 0.05)
#' @export
turkheimer_est <- function(p_seq, alpha) {
  p_seq <- sort(p_seq)
  oldlen <- length(p_seq)
  alphas <- c(0.15, 0.1, 0.05, 0.025, 0.01)
  quantiles <- c(0.973, 1.073, 1.224, 1.358, 1.518)
  thr <- quantiles[which(alphas == alpha)]
  k <- 0

  repeat {
    n <- length(p_seq)
    c_min <- max((1:n) / (n + 1) - p_seq)
    c_min_mod <- (c_min + 0.4 / n) * (sqrt(n) + 0.2 + 0.68 / sqrt(n))

    if (c_min_mod < thr) break

    p_seq <- p_seq[-1]
    k <- k + 1
  }

  newlen <- length(p_seq)
  br <- 1:newlen

  vars <- br * (newlen - br + 1) / (newlen + 1)^2 / (newlen + 2)

  data_pval <- data.frame("pvals" = p_seq, "ind" = (1 + k):oldlen/ oldlen)
  wls_model <- lm(pvals ~ ind, data = data_pval, weights = vars)

  return(1 / wls_model$coefficients[2])
}
