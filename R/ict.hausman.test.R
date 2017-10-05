

#' Hausman Specification Test for Two List Experiment Regression Fit Objects
#'
#' @param fit1 
#' @param fit2 
#'
#' @return
#' @export
#'
#' @examples
ict.hausman.test <- function(fit1, fit2){
  if(length(coef(fit1)) != length(coef(fit2))){
    stop("Please provide two ictreg fit objects that used the same number of covariates.")
  }
  stat <-
    as.numeric((coef(fit1) - coef(fit2)) %*%
                 solve((vcov(fit1) - vcov(fit2))) %*%
                 (coef(fit1) - coef(fit2)))
  df <- length(coef(fit1))
  p <- pchisq(q = abs(stat), df = df, lower.tail = FALSE)
  return(structure(list(stat = stat, df = df, p = p), class = "ict.hausman.test"))
}

#' @export
print.ict.hausman.test <- function(x, ...){
  
  cat("\nTest for List Experiment Model Misspecification\n\n")
  
  cat("Hausman statistic = ", x$stat, ", df = ", x$df, ", p-value = ", x$p, ".", sep = "")
  
  cat("\n")
  
}
