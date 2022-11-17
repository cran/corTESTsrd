#' Significance test of rank cross-correlations
#'
#' Significance test of Spearman's Rho or Kendall's Tau between time series of short-range dependent random variables.
#' The test is based on the asymptotic normal distributions of the estimators.
#
#' @usage corTESTsrd(x, y,
#'            iid=TRUE, method="spearman",
#'            alternative="two.sided",
#'            kernelf=function(z){return(ifelse(abs(z) <= 1, (1 - z^2)^2, 0))},
#'            bwf=function(n){3*n^(1/4)})
#'
#' @param x numeric input vector.
#' @param y numeric input vector.
#' @param iid logical, if TRUE, observations are assumed to be iid, if FALSE observations are assumed to be short-range dependent and the long-run variance of the estimator is estimated from the observations.
#' @param method a character string, indicating which correlation coefficient should be used for the test. One of "spearman" or "kendall", cannot be abbreviated.
#' @param alternative a character string indicating the alternative hypothesis. Must be one of "two.sided", "greater" or "less", cannot be abbreviated.
#' @param kernelf a function that is used in the estimation procedure. The default kernel-function is a quartic kernel. Should be a vectorized function.
#' @param bwf a function for choosing the bandwidth, based on the sample size \eqn{n}, that should be used in the estimation procedure. Default is \eqn{3n^{1/4}}{3*n^(1/4)}, \eqn{b_n=o(n^{1/2})}{bn=o(n^(1/2))} must hold.
#'
#' @details Calculates an estimate of the rank correlation coefficient between the inputs x and y, which are assumed to be evenly spaced time series with equal time-increments,
#' and performs a significance test for the rank correlation coefficient with \eqn{\mathcal{H}_0: \rho_S/\tau=0}{H0: \rho/\tau=0} against an alternative specified by the user. The function returns the estimate of the rank correlation coefficient
#' and a p-value. Missing observations (\code{NA}) are allowed, but will prompt a warning. Ties are not allowed. \cr \cr
#' The test statistic and the corresponding p-value are based on the distribution of the respective estimator under the assumption of independence between the inputs x and y, and an additional assumption regarding
#' the dependence structure of the inputs on their own past. The distribution of the test statistic is modelled as a normal distribution. \cr \cr
#' If the option iid is TRUE, the inputs are assumed to be realizations of independent and identically distributed random variables. In this case the asymptotic variance of the test statistic is given by \eqn{\frac{1}{n-1}}{(1)/(n-1)}
#' for Spearman's Rho and as \eqn{\frac{2(2n+5)}{9n(n-1)}}{(2(2n+5))/(9n(n-1))} for Kendall's Tau, see Gibbons and Chakraborti (2003), equations 3.13 and 2.29 in chapter 11, respectively. \cr \cr
#' If the option iid is FALSE, the inputs are assumed to be realizations of short-range dependent random variables (see Corollary 1 in Lun et al., 2022). The asymptotic variance of the test statistic is modelled as
#' \eqn{\frac{1}{n}(1+2\sum_{h=1}^{\infty} \rho_S^X(h) \rho_S^Y(h))}{(1/n)(1+2\sum \rho_X(h) \rho_Y(h))} for Spearman's Rho and as \eqn{\frac{4}{9n}(1+2\sum_{h=1}^{\infty} \rho_S^X(h) \rho_S^Y(h))}{(4/(9n))(1+2\sum \rho_X(h) \rho_Y(h))} for Kendall's Tau.
#' Here \eqn{\rho_S^X(h)}{\rho_X(h)} refers to the Spearman autocorrelation of the first input x for lag \eqn{h}, and the analogue applies to \eqn{\rho_S^X(h)}{\rho_Y(h)}. In this case the asymptotic variance of the test statistic is estimated
#' (see Corollary 2 in Lun et al., 2022). For this estimation procedure a kernel-function together with a bandwidth is used, which can be specified by the user.
#'
#' @references J. D. Gibbons, and S. Chakraborti, Nonparametric statistical inference (4th Edition). CRC press, 2003. \cr \cr
#' D. Lun, S. Fischer, A. Viglione, and G. Blöschl, Significance testing of rank cross-correlations between autocorrelated time series with short-range dependence, Journal of Applied Statistics, 2022, 1-17. \doi{10.1080/02664763.2022.2137115}.
#'
#' @return Estimate of rank correlation coefficient and p-value of corresponding hypothesis test.
#' @import stats
#' @export
#'
#' @examples
#' #Demonstration
#' sam_size = 50
#' nsim = 1000
#'
#' pval_iid <- rep(NA, nsim)
#' pval_srd <- rep(NA, nsim)
#' #iid-simulation: if we have iid observations the modified test
#' #is able to maintain the desired significance level
#' for(j in c(1:nsim)) {
#'   x <- rnorm(n=sam_size)
#'   y <- rnorm(n=sam_size)
#'   pval_iid[j] <- corTESTsrd(x, y, iid=TRUE, method="spearman")[2]
#'   pval_srd[j] <- corTESTsrd(x, y, iid=FALSE, method="spearman")[2]
#' }
#' sum(pval_iid <= 0.05)/nsim
#' sum(pval_srd <= 0.05)/nsim
#'
#' #ar(1)-simulation: if we have srd-observations the modified test
#' #counteracts the inflation of type-I-errors
#' for(j in c(1:nsim)) {
#'   x <- as.numeric(arima.sim(model=list(ar=c(0.8)), n=sam_size))
#'   y <- as.numeric(arima.sim(model=list(ar=c(0.8)), n=sam_size))
#'   pval_iid[j] <- corTESTsrd(x, y, iid=TRUE, method="spearman")[2]
#'   pval_srd[j] <- corTESTsrd(x, y, iid=FALSE, method="spearman")[2]
#' }
#' sum(pval_iid <= 0.05)/nsim
#' sum(pval_srd <= 0.05)/nsim
#'
#' #the test can be made more conservative be choosing a bigger bandwidth,
#' #but this decreases the power
#' bwfbig <- function(n) {10*(n^(1/4))}
#' #ar(1)-simulation: if we have srd-observations the modified test
#' #counteracts the inflation of type-I-errors
#' for(j in c(1:nsim)) {
#'   x <- as.numeric(arima.sim(model=list(ar=c(0.8)), n=sam_size))
#'   y <- as.numeric(arima.sim(model=list(ar=c(0.8)), n=sam_size))
#'   pval_iid[j] <- corTESTsrd(x, y, iid=TRUE, method="spearman")[2]
#'   pval_srd[j] <- corTESTsrd(x, y, iid=FALSE, method="spearman", bwf=bwfbig)[2]
#' }
#' sum(pval_iid<=0.05)/nsim
#' sum(pval_srd<=0.05)/nsim
#'
corTESTsrd <- function(x, y, iid=TRUE, method="spearman", alternative="two.sided",
                       kernelf=function(z){return(ifelse(abs(z) <= 1, (1 - z^2)^2, 0))},
                       bwf=function(n){3*n^(1/4)}) {

  if(!is.vector(x)|!is.vector(y)) {
    stop("Please only provide numeric vectors as inputs for x and y")
  }
  if(length(x)!=length(y)) {
    stop("x and y should be vectors with the same length")
  }
  if(any(duplicated(x[!is.na(x)]))|any(duplicated(y[!is.na(y)]))) {
    stop("Ties are not permitted in the input vectors")
  }

  if(!is.logical(iid)) {
    stop("The paramter iid should be either TRUE or FALSE, indicating if observations are assumed to be independent and identically distributed or short-range dependent")
  }
  if(!method %in% c("spearman","kendall")) {
    stop("method should either be 'spearman' or 'kendall'")
  }
  if(!alternative %in% c("two.sided", "less", "greater")) {
    stop("alternative should either be 'two.sided', 'less' or 'greater'")
  }
  if(!is.function(kernelf)) {
    stop("kernelf should be a vectorized function")
  }
  if(!is.function(bwf)) {
    stop("bwf should be a function")
  }

  if(any(is.na(x)) | any(is.na(y))) {
    warning("your data contains NAs")
  }

  n = sum(!is.na(x)&!is.na(y))
  if( (length(x)==0)|(length(y)==0)|(n<4) ) {
    stop("Less than 4 observations available for calculation")
  }

  corr_est = cor(x,y,method=method,use="pairwise.complete.obs")
  if(iid==TRUE) {
    var_est = ifelse(method=="spearman",spearman_var_independent(n),kendall_var_independent(n))
  }
  else if(iid==FALSE) {
    var_est = ifelse(method=="spearman",
                     long_run_var_est_spearman_H0(x=x,y=y,kernelf = kernelf,bwf=bwf)*(1/n),
                     long_run_var_est_kendall_H0(x=x,y=y,kernelf = kernelf,bwf=bwf)*(1/n))
  }

  if(alternative=="less") {
    p_val = pnorm(corr_est,mean=0,sd=sqrt(var_est) )
  }
  else if(alternative=="greater") {
    p_val = 1-pnorm(corr_est,mean=0,sd=sqrt(var_est) )
  }
  else if(alternative=="two.sided") {
    p_val = 2*(1-pnorm(abs(corr_est),mean=0,sd=sqrt(var_est) ) )
  }

  return(c(rho=corr_est,pval=p_val))
}
