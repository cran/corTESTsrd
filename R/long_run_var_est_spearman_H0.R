#' Estimator of long-run variance for Spearman's Rho
#'
#' Estimates the long-run variance of the estimator of Spearman's Rho between short-range dependent observations
#' of independent random variables. The expression comes from the asymptotic normal distribution of the estimator, see Corollary 2 in Lun et al. (2022).
#'
#' @usage long_run_var_est_spearman_H0(x, y,
#'            kernelf=function(z) {return(ifelse(abs(z) <= 1,(1 - z^2)^2, 0))},
#'            bwf=function(nnn){3*nnn^(1/4)})
#'
#' @param x numeric input vector.
#' @param y numeric input vector.
#' @param kernelf kernel-function that should be used in the estimation procedure.
#' @param bwf function for choosing the bandwidth based on the sample size that is used in the estimation procedure.
#'
#' @references D. Lun, S. Fischer, A. Viglione, and G. Blöschl, Significance testing of rank cross-correlations between autocorrelated time series with short-range dependence, submitted to Journal of Applied Statistics, 2022.
#'
#' @return Estimate of long-run variance of estimator.
#' @keywords internal
#' @export
#'
#' @examples
#' long_run_var_est_spearman_H0(x=rnorm(50),y=rnorm(50))
long_run_var_est_spearman_H0 <- function(x, y,
                                          kernelf=function(z) {return(ifelse(abs(z) <= 1,(1 - z^2)^2, 0))},
                                          bwf=function(nnn){3*nnn^(1/4)}) {

  if(length(x)!=length(y)) {
    stop("x and y should be vectors with the same length")
  }

  n = length(x)

  spearman_acf_x = rankacf(x,lag.max = n-2)
  spearman_acf_y = rankacf(y,lag.max = n-2)
  bn = bwf(n)

  sigma2_est = 1+2*sum(kernelf(c(1:(n-2))/bn)*spearman_acf_x[c(1:(n-2))]*spearman_acf_y[c(1:(n-2))])

  return(sigma2_est)
}
