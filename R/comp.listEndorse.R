#' Comparing List and Endorsement Experiment Data
#' 
#' Function to conduct a statistical test with the null hypothesis that there
#' is no difference between the correlation coefficients between list
#' experiment and endorsement experiment data.
#' 
#' This function allows the user to calculate the correlation between list and
#' endorsement experiment data within the control group and the treatment
#' group, and to conduct a statistical test with the null hypothesis of no
#' difference between the two correlation coefficients.
#' 
#' @param y.endorse A numerical matrix containing the response data for the
#' endorsement experiment.
#' @param y.list A numerical vector containing the response data for a list
#' experiment.
#' @param treat A numerical vector containing the binary treatment status for
#' the experiments. The treatment assignment must be the same for both
#' experiments to compare across experiments.
#' @param n.draws Number of Monte Carlo draws.
#' @param alpha Confidence level for the statistical test.
#' @param endorse.mean A logical value indicating whether the mean endorsement
#' experiment response is taken across questions.
#' @param method The method for calculating the correlation, either Pearson's
#' rho or Kendall's tau.
#' @return \code{comp.listEndorse} returns a list with four elements: the
#' correlation statistic (rho or tau) for the treatment group as
#' \code{cor.treat}, the correlation statistic for the control group as
#' \code{cor.control}, the p.value for the statistical test comparing the two
#' correlation statistics as \code{p.value}, and the bootstrapped confidence
#' interval of the difference as \code{ci}.
#' @author Graeme Blair, UCLA, \email{graeme.blair@ucla.edu}
#' and Kosuke Imai, Princeton University, \email{kimai@princeton.edu}
#' @references Blair, Graeme, Jason Lyall and Kosuke Imai. (2014) ``Comparing
#' and Combining List and Experiments: Evidence from Afghanistan."  American
#' Journal of Political Science. available at
#' \url{http://imai.princeton.edu/research/comp.html}
#' @keywords models regression
#' @export comp.listEndorse
comp.listEndorse <- function(y.endorse, y.list, treat, n.draws = 10000, alpha = .05, endorse.mean = FALSE,
                             method = "pearson") {

  if (endorse.mean == FALSE) {
    y.endorse.mean <- apply(y.endorse, 1, function(x) mean(x, na.omit = T))
  } else {
    y.endorse.mean <- y.endorse
  }
  
  calc.corr <- function(y.endorse.mean, y.list, treat) {
    
    no.na.vec <- !(is.na(y.endorse.mean) | is.na(y.list))
    
    cor.treat <- cor(y.endorse.mean[treat == 1 & no.na.vec == 1],
                     y.list[treat==1 & no.na.vec == 1], method = method)
    cor.control <- cor(y.endorse.mean[treat == 0 & no.na.vec == 1],
                       y.list[treat==0 & no.na.vec == 1], method = method)
    return(list(cor.treat = cor.treat, cor.control = cor.control))
  }

  calc.d <- function(cor.treat, cor.control) {
    
    d <- cor.treat - cor.control
    
    return(d)
  }

  cor <- calc.corr(y.endorse.mean, y.list, treat)

  d.bootstrap <- rep(NA, n.draws)
  for(i in 1:n.draws) {
    sample <- sample(1:length(treat), length(treat), replace = T)
    cor.bootstrap <- calc.corr(y.endorse.mean[sample], y.list[sample], treat[sample])
    d.bootstrap[i] <- calc.d(cor.treat = cor.bootstrap$cor.treat,
                               cor.control = cor.bootstrap$cor.control)
  }

  return(list(cor.treat = cor$cor.treat, cor.control = cor$cor.control,
              p.value = mean(d.bootstrap <= 0),
              ci = c(quantile(d.bootstrap, alpha/2), quantile(d.bootstrap, (1-alpha)/2))))
  
}


