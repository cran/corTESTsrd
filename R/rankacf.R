#' Estimator of Spearman autocorrelations
#'
#' Computes Spearman-autocorrelations when the observations contain \code{NA}s.
#' If the observations do not contain \code{NA}s, the function is consistent with the function acf applied to ranks.
#'
#' @usage rankacf(x,lag.max=length(x)-2)
#'
#' @param x numeric input vector.
#' @param lag.max maximum lag for which acf is estimated, automatically limited to one less than number of observations, should be positive integer.
#'
#' @return Estimated spearman autocorrelation up to lag.max.
#' @keywords internal
#' @export
#' @import stats
#'
#' @examples
#' x = rnorm(10)
#' rankacf(x)
#' acf(rank(x),plot=FALSE,lag.max = length(x)-2)$acf[-1,,1]
rankacf = function(x,lag.max=length(x)-2) {
  if(lag.max>length(x)) {lag.max = length(x)-1}
  r = rank(x,na.last="keep")
  acfvals = rep(NA,lag.max)
  nobs = sum(!is.na(x))
  n = length(x)
  m = mean(r,na.rm=T)
  gamma_0 = (1/nobs)*sum((r-m)^2,na.rm=T)
  for(j in c(1:lag.max)) {
    acfvals[j] = (1/nobs)*sum((r[1:(n-j)]-m)*(r[(1+j):(n)]-m),na.rm=T)/gamma_0
  }
  names(acfvals) = c(1:lag.max)

  return(acfvals)
}
