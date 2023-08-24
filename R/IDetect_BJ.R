#' Berk-Jones statistic on subsets
#'
#' Given a sequence of sorted p-values and the start and endpoint indices \eqn{s} and \eqn{e},
#' check if the values \eqn{p_{(s)},..., p_{(e)}} form a sample from uniform \eqn{U[p_{(s-1)}, p_{(e+1)}]}
#' distribution.
#'
#' @param ind A vector containing the starting and the ending point of the interval
#' @param p_seq The sorted sequence of p-values
#' @param mbj Boolean, if `TRUE`, the modified Berk-Jones statistic is used,
#' otherwise the standard form
#'
#' @references Li, J., & Siegmund, D. (2015). Higher criticism:
#'  $p$-values and criticism. The Annals of Statistics, 43(3), 1323–1350.
#'   https://doi.org/10.1214/15-aos1312
#'
#' @return Returns a list of 2 elements
#' \item{val}{the value of the Berk-Jones statistic}
#' \item{loc}{the location in the corresponding sequence where the maximum is achieved}
#'
max_BJ <- function(ind = c(1, length(p_seq)), p_seq, mbj = FALSE) {
  s = ind[1]
  e = ind[2]

  if (e - s < 3) {
    return(c(0, 0))
  }

  if (e == length(p_seq)) {
    p_seq <- c(p_seq, 1)
  }

  if (s == 1) {
    p_seq <- c(0, p_seq)
    s = s + 1

    if (ind[2] == length(p_seq))
      e = e + 1
  }

  i <- (s + 1):e

  scaled_p <- (p_seq[i] - p_seq[s - 1]) / (p_seq[e + 1] - p_seq[s - 1])

  if (mbj == TRUE) {
    BJ_seq <- (scaled_p < i / (e - s)) * ((i - s) * log((i - s) / (e - s) / scaled_p) - ((i - s) / (e - s) - scaled_p))
  } else {
    BJ_seq <- (scaled_p < i / (e - s)) * ((i - s) * log((i - s) / (e - s) / scaled_p) + (e - i) * log((e - i) / (e - s) / (1 - scaled_p)))
  }

  ip.max <- which.max(BJ_seq)

  return(c(s + ip.max - 1, BJ_seq[ip.max]))
}


#' IDetect method with the Berk-Jones statistic
#'
#' Given a sequence of sorted p-values, perform segmentation into uniform segments using
#' the IDetect method with Berk-Jones statistic. The code provided is a slight
#'  modification of the code for the `IDetect` procedure from the `breakfast` package
#'
#' @param p_seq The sorted sequence of p-values
#' @param thr The threshold for the Berk-Jones statistic. If `NULL` the threshold is calculated through simulations
#' @param points The number of points added to the interval at each step (the minimum number of points)
#' @param mode Set to 'left' or 'right' to indicate the type of expanding intervals to consider:
#' 'left' for left-expanding intervals and 'right' for right-expanding intervals
#' @param mbj Boolean, if `TRUE`, the modified Berk-Jones statistic is used,
#' otherwise the standard form
#' @param alpha The significance level for calculating the threshold
#' @param N The number of repetitions for calculating the threshold
#'
#' @references Li, J., & Siegmund, D. (2015). Higher criticism:
#'  $p$-values and criticism. The Annals of Statistics, 43(3), 1323–1350.
#'   https://doi.org/10.1214/15-aos1312,
#'
#' Anastasiou, A., & Fryzlewicz, P. (2022).
#' Detecting multiple generalized change-points by isolating single ones.
#' Metrika, 85(2), 141–174. https://doi.org/10.1007/s00184-021-00821-6
#'
#' @return Returns a list of 2 elements
#' \item{changepoint}{The locations of the estimated change-points}
#' \item{full_information}{A dataframe with the information about the Berk-Jones statistic on the subintervals
#' at each step of the algorithm.
#' The endpoints of the interval considered are in the first two columns, the location of the maximum in the third,
#' and the value of the Berk-Jones statistic at the fourth column}
#' @export
#' @example
#' set.seed(1)
#' p_seq = c(runif(100,max = 0.1),runif(300, max = 0.3), runif(600))
#' plot(sort(p_seq))
#' model = idetect.th_BJ(sort(p_seq), mode = "left")
#' model$changepoint
#' model$full_information
idetect.th_BJ <- function(p_seq, thr = NULL , points = 5, mode = "left", mbj = FALSE, alpha = 0.05, N = 1000)
{
  n <- length(p_seq)

  if (is.null(thr)) {
    thr = BJ_threshold_calc(n, points, N, alpha, mbj = mbj)
  }

  if(mode == "left")
  {
    res = idetect.th_BJ_le(p_seq, thr, points = points, s = 1, e = length(p_seq),  k_l = 1, mbj = mbj)
  }

  if(mode == "right")
  {
    res = idetect.th_BJ_re(p_seq, thr, points = points, s = 1, e = length(p_seq),  k_r = 1, mbj = mbj)
  }

  return(res)


}


#' IDetect method with right-expanding intervals and the Berk-Jones statistic
#'
#' Given a sequence of sorted p-values, perform segmentation into uniform segments using
#' the IDetect method with Berk-Jones statistic. We consider right-expanding intervals.
#'
#' @param p_seq The sorted sequence of p-values
#' @param thr The threshold for the Berk-Jones statistic
#' @param points The number of points added to the interval at each step (the minimum number of points)
#' @param s The starting point of the current interval observed for change (for the recursion)
#' @param e The endpoint of the current interval observed for change (for the recursion)
#' @param k_r The current step in expanding the right edge of the subintervals within the (s,e) range.
#' @param mbj Boolean, if `TRUE`, the modified Berk-Jones statistic is used,
#' otherwise the standard form
#'
#' @return Returns a list of 2 elements
#' \item{changepoint}{The locations of the estimated change-points}
#' \item{full_information}{A dataframe with the information about the Berk-Jones statistic on the subintervals
#' at each step of the algorithm.
#' The endpoints of the interval considered are in the first two columns, the location of the maximum in the third,
#' and the value of the Berk-Jones statistic at the fourth column}
#' @example
#' set.seed(1)
#' p_seq = c(runif(100,max = 0.1),runif(300, max = 0.3), runif(600))
#' plot(sort(p_seq))
#' model_r = idetect.th_BJ_re(sort(p_seq), thr = 11)
#' model_r$changepoint
#' model_r$full_information
#' @keywords internal
idetect.th_BJ_re <- function(p_seq, thr, points = 5, s = 1, e = length(p_seq),  k_r = 1, mbj = FALSE) {

  l <- length(p_seq)
  Res <- matrix(0, 1, 4)
  points <- as.integer(points)
  r_e_points <- seq(points, l, points)
  chp <- 0
  if (e - s < points) {
      Res_fin <- matrix(0, 1, 4)
      cpt <- integer(0)
    } else {
      pos_r <- numeric()
      BJ_r <- numeric()
      right_points=c(r_e_points[r_e_points>s & r_e_points<e],e)
      rur <- length(right_points)
      while ( (chp == 0) & (k_r < rur)) {
        ind <- c(s, right_points[k_r])
        tmp <- max_BJ(ind, p_seq, mbj = mbj)
        pos_r[k_r] <- tmp[1]
        BJ_r[k_r] <- tmp[2]
        Res <- rbind(Res, c(s,right_points[k_r],pos_r[k_r],BJ_r[k_r]))
        if (BJ_r[k_r] > thr) {
          chp <- pos_r[k_r]
        } else {
          k_r <- k_r + 1
        }
      }
      if (chp == 0) {
        while ( (chp == 0) &  (k_r <= rur)) {
          ind <- c(s, right_points[k_r])
          tmp <- max_BJ(ind, p_seq, mbj = mbj)
          pos_r[k_r] <- tmp[1]
          BJ_r[k_r] <- tmp[2]
          Res <- rbind(Res, c(s,right_points[k_r],pos_r[k_r],BJ_r[k_r]))
          if (BJ_r[k_r] > thr) {
            chp <- pos_r[k_r]
          }
          else {k_r <- k_r+1}
        }
      }
      if (chp != 0) {

        r <- idetect.th_BJ_re(p_seq, thr, points = points, s = chp + 1, e = e, k_r = 1)
        cpt <- c(chp, r[[1]])
        Res_fin <- rbind(Res, r[[2]])
      } else {
        cpt <- chp
        Res_fin <- Res
      }
    }
    cpt <- cpt[cpt != 0]
    Res_fin <- Res_fin[which(Res_fin[,3] != 0), , drop = FALSE]

  return(list(changepoints = sort(cpt), full_information = Res_fin))
}

#' IDetect method with left-expanding intervals and the Berk-Jones statistic
#'
#' Given a sequence of sorted p-values, perform segmentation into uniform segments using
#' the IDetect method with Berk-Jones statistic. We consider left-expanding intervals.
#'
#' @param p_seq The sorted sequence of p-values
#' @param thr The threshold for the Berk-Jones statistic
#' @param points The number of points added to the interval at each step (the minimum number of points)
#' @param s The starting point of the current interval observed for change (for the recursion)
#' @param e The endpoint of the current interval observed for change (for the recursion)
#' @param k_l The current step in expanding the left edge of the subintervals within the (s,e) range.
#' @param mbj Boolean, if `TRUE`, the modified Berk-Jones statistic is used,
#' otherwise the standard form
#' @return Returns a list of 2 elements
#' \item{changepoint}{The locations of the estimated change-points}
#' \item{full_information}{A dataframe with the information about the Berk-Jones statistic on the subintervals
#' at each step of the algorithm.
#' The endpoints of the interval considered are in the first two columns, the location of the maximum in the third,
#' and the value of the Berk-Jones statistic at the fourth column}
#' @example
#' set.seed(1)
#' p_seq = c(runif(100,max = 0.1),runif(300, max = 0.3), runif(600))
#' plot(sort(p_seq))
#' model_l = idetect.th_BJ_le(sort(p_seq), thr = 11)
#' model_l$changepoint
#' model_l$full_information[ model_l$full_information[,4] > thr, ]
#' @keywords internal
idetect.th_BJ_le <- function(p_seq, thr, points = 5, s = 1, e = length(p_seq), k_l = 1, mbj = FALSE) {
  l <- length(p_seq)
  Res <- matrix(0, 1, 4)
  points <- as.integer(points)
  l_e_points <- seq(l, points, -points)
  chp <- 0

  if (e - s < points) {
    Res_fin <- matrix(0, 1, 4)
    cpt <- integer(0)
  } else {
    pos_l <- numeric()
    BJ_l <- numeric()
    left_points <- c(l_e_points[l_e_points > s & l_e_points < e], s)
    lur <- length(left_points)

    while ((chp == 0) & (k_l < lur)) {
      ind <- c(left_points[k_l], e)
      tmp <- max_BJ(ind, p_seq, mbj = mbj)
      pos_l[k_l] <- tmp[1]
      BJ_l[k_l] <- tmp[2]
      Res <- rbind(Res, c(left_points[k_l], e, pos_l[k_l], BJ_l[k_l]))

      if (BJ_l[k_l] > thr) {
        chp <- pos_l[k_l]
      } else {
        k_l <- k_l + 1
      }
    }

    if (chp == 0) {
      while ((chp == 0) & (k_l <= lur)) {
        ind <- c(left_points[k_l], e)
        tmp <- max_BJ(ind, p_seq, mbj = mbj)
        pos_l[k_l] <- tmp[1]
        BJ_l[k_l] <- tmp[2]
        Res <- rbind(Res, c(left_points[k_l], e, pos_l[k_l], BJ_l[k_l]))

        if (BJ_l[k_l] > thr) {
          chp <- pos_l[k_l]
        } else {
          k_l <- k_l + 1
        }
      }
    }

    if (chp != 0) {
      r <- idetect.th_BJ_le(p_seq, s = s, e = chp, points = points, k_l = 1, thr = thr)
      cpt <- c(chp, r[[1]])
      Res_fin <- rbind(Res, r[[2]])
    } else {
      cpt <- chp
      Res_fin <- Res
    }
  }

  cpt <- cpt[cpt != 0]
  Res_fin <- Res_fin[which(Res_fin[, 3] != 0), , drop = FALSE]

  return(list(changepoints = sort(cpt), full_information = Res_fin))
}



#' Simulating the threshold for the IDetect Berk-Jones procedure
#'
#' Calculate the threshold to be used by `idetect.th_BJ_re` or `idetect.th_BJ_le`.
#' The threshold is calculated as an upper quantile of the sequence of test statistics,
#'  such that the procedure under the null (for uniform p-values) does not detect any change-points.
#'
#' @param n The length of the p-value sequence
#' @param points The length of the smallest interval considered (step parameter)
#' @param N The number of repetitions
#' @param alpha The upper quantile to be used
#' @param mbj Boolean, if `TRUE`, the modified Berk-Jones statistic is used,
#' otherwise the standard form
#' @references  Anastasiou, A., & Fryzlewicz, P. (2022).
#' Detecting multiple generalized change-points by isolating single ones.
#' Metrika, 85(2), 141–174. https://doi.org/10.1007/s00184-021-00821-6
#'
#' @return Returns the estimated threshold value
#' @export
BJ_threshold_calc <- function(n, points, N, alpha, mbj = FALSE) {
  max_seq <- rep(0, N)

  for (i in 1:N) {
    r_e_points <- seq(points, n, points)
    values <- rep(0, length(r_e_points))
    x <- sort(runif(n))

    for (j in 1:length(r_e_points)) {
      values[j] <- max_BJ(ind = c(1, r_e_points[j]), x, mbj = mbj)[2]
    }

    max_seq[i] <- max(values)
  }

  return(quantile(max_seq, probs = 1 - alpha))
}


#' Local FDR estimator using IDetect Berk-Jones procedure
#'
#' Calculate the estimate of the local FDR using the IDetect Berk-Jones procedure
#' with right or left expanding intervals
#'
#' @param p_seq The sequence of sorted p-values
#' @param thr The threshold for the procedure. If `NULL`, then calculate using the
#' `BJ_threshold_calc` function.
#' @param points The length of the smallest interval considered (step parameter)
#' @param mode Set to 'left' for the procedure with left expanding intervals, or 'right'
#' for the procedure with right expanding intervals
#' @param mbj Boolean, if `TRUE`, the modified Berk-Jones statistic is used,
#' otherwise the standard form
#' @param plot Boolean, if `TRUE` a plot of the estimated local FDR will be produced
#' @param alpha The significance level for calculating the threshold
#' @param N The number of repetitions for calculating the threshold
#'
#' @return Returns a list of 2 elements
#' \item{lfdr}{The estimated lfdr values corresponding to the p-values sequence}
#' \item{groups}{A dataframe with the information about the resulting grouping of the p-values
#' based on their significance}
#'
#' @references  Anastasiou, A., & Fryzlewicz, P. (2022).
#' Detecting multiple generalized change-points by isolating single ones.
#' Metrika, 85(2), 141–174. https://doi.org/10.1007/s00184-021-00821-6
#'
#' @export
#' @example
#' set.seed(6)
#' p_seq = 1-pnorm(c(rnorm(50,3),rnorm(100, 2), rnorm(200)))
#' p_seq = sort(p_seq)
#' thr = BJ_threshold_calc(length(p_seq),points=5,N=1000,alpha = 0.05)
#' lfdr_right = ID_BJ_group(p_seq, thr = 10, mod = 'right', plot = TRUE)
#' lfdr_left =  ID_BJ_group(p_seq, thr = 10, mod = 'left', plot = TRUE)
#' lfdr_left
ID_BJ_group <- function(p_seq, thr = NULL, points = 5, mode = 'left', mbj = FALSE,  plot = FALSE, alpha = 0.05, N = 1000) {
  n = length(p_seq)

  cps = idetect.th_BJ(sort(p_seq), thr = thr, points = points, mode = mode, mbj = mbj)$changepoints

  hist = hist(p_seq, breaks = c(0, p_seq[cps], 1), plot = FALSE)
  lfdr_values = 1 / hist$density * hist$density[length(hist$density)]

  lfdr_est = c()
  div_cps = c(0, cps, n)

  for (i in 1:(length(div_cps) - 1)) {
    lfdr_est = c(lfdr_est, rep(lfdr_values[i], diff(div_cps)[i]))
  }

  if (plot == TRUE) {
    plot(lfdr_est, col = "red", lwd = 2, type = "l", ylim = c(0, 1), ylab = "", xlab = '', main = "local FDR estimate")
  }

  res = list()
  res$lfdr = lfdr_est
  res$groups = data.frame(start = div_cps[-length(div_cps)] + 1, end = div_cps[-1])
  return(res)
}


