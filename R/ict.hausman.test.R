

#' Hausman Specification Test for Two List Experiment Regression Fit Objects
#'
#' @param ml Maximum likelihood model fit, for example from ictreg(method = "ml")
#' @param nls NLS model fit, for example from ictreg(method = "nls")
#' @param abs Flag to override error when Hausman statistic is negative, which may indicate misspecification. Set to \code{FALSE} to override.
#' @param psd Flag to override error when variance-covariance difference is non-positive semidefinite, suggesting misspecification.  Set to \code{TRUE} to override.
#'
#' @return List containing Hausman statistic, p-value, and the degrees of freedom of the test.
ict.hausman.test <- function(ml, nls, abs = FALSE, psd = FALSE){
  if(length(coef(ml)) != length(coef(nls))){
    stop("Please provide two ictreg fit objects that used the same number of covariates.")
  }
  psd.check <- min(eigen(vcov(nls) - vcov(ml))$values) < 0
  stat <-
    as.numeric((coef(ml) - coef(nls)) %*%
                 ginv((vcov(nls) - vcov(ml))) %*%
                 (coef(ml) - coef(nls)))
  # we found that the following check for positive semi-definiteness led to 
  # significant overrejection of the null
  if (psd.check & psd) stop("The variance-covariance difference is non-positive semidefinite, suggesting misspecification.  Set psd = FALSE to compute anyway.")  
  if (stat < 0 & !abs) stop("Hausman test statistic is negative, suggesting possible misspecification.  Set abs = TRUE to compute anyway.")  
  df <- length(coef(ml))
  p <- pchisq(q = abs(stat), df = df, lower.tail = FALSE)
  return(structure(list(stat = stat, df = df, p = p), class = "ict.hausman.test"))
}

print.ict.hausman.test <- function(x, ...){
  
  cat("\nTest for List Experiment Model Misspecification\n\n")
  
  cat("Hausman statistic = ", x$stat, ", df = ", x$df, ", p-value = ", x$p, ".", sep = "")
  
  cat("\n")
  
}
