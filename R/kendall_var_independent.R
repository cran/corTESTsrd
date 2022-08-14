#' Asymptotic variance of estimator of Kendall's Tau between independent iid random variables
#'
#' Calculates the asymptotic variance of the estimator of Kendall's Tau between iid observations
#' of independent random variables. The expression comes from the asymptotic normal distribution of the estimator,
#' see e.g. equation 2.30 in chapter 11 of Gibbons and Chakraborti (2003).
#'
#' @usage kendall_var_independent(n)
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
#' kendall_var_independent(n=50)
kendall_var_independent <- function(n) {(2*(2*n + 5))/(9*n*(n - 1))}
