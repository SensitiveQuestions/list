

#' Hausman Specification Test for Two List Experiment Regression Fit Objects
#'
#' @param ml 
#' @param nls 
#'
#' @return
#' @export
#'
#' @examples
ict.hausman.test <- function(ml, nls, ginv = FALSE){
  if(length(coef(ml)) != length(coef(nls))){
    stop("Please provide two ictreg fit objects that used the same number of covariates.")
  }
  if (ginv) {
    require(MASS)
    stat <-
      as.numeric((coef(ml) - coef(nls)) %*%
                   ginv((vcov(nls) - vcov(ml))) %*%
                   (coef(ml) - coef(nls))) 
  } else {
    stat <-
      as.numeric((coef(ml) - coef(nls)) %*%
                   solve((vcov(nls) - vcov(ml))) %*%
                   (coef(ml) - coef(nls)))
    if (stat < 0) stop("Hausman test statistic is negative, suggesting misspecification.  Set ginv = TRUE to compute anyway.")    
  }
  df <- length(coef(ml)) - 1
  p <- pchisq(q = stat, df = df, lower.tail = FALSE)
  return(structure(list(stat = stat, df = df, p = p), class = "ict.hausman.test"))
}

#' @export
print.ict.hausman.test <- function(x, ...){
  
  cat("\nTest for List Experiment Model Misspecification\n\n")
  
  cat("Hausman statistic = ", x$stat, ", df = ", x$df, ", p-value = ", x$p, ".", sep = "")
  
  cat("\n")
  
}
