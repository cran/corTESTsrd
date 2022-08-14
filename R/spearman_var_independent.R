#' Asymptotic variance of estimator of Spearman's Rho between independent iid random variables
#'
#' Calculates the asymptotic variance of the estimator of Spearman's Rho between iid observations
#' of independent random variables. The expression comes from the asymptotic normal distribution of the estimator,
#' see e.g. equation 3.13 in chapter 11 of Gibbons and Chakraborti (2003).
#'
#' @usage spearman_var_independent(n)
#'
#' @param n number of observations, should be integer bigger than 0.
#'
#' @references J. D. Gibbons, and S. Chakraborti, Nonparametric statistical inference (4th Edition). CRC press, 2003.
#'
#' @return Asymptotic variance of estimator.
#' @keywords internal
#' @export
#'
#' @examples
#' spearman_var_independent(n=50)
spearman_var_independent <- function(n) {1/(n - 1)}
