\name{ictregBayes} 
\alias{ictregBayes}
\title{Item Count Technique} 
\description{ 
Function to conduct multivariate regression analyses of survey data with the item count technique, also known as the list experiment and the unmatched count technique.
} 
\usage{
ictregBayes(formula, data = parent.frame(), treat = "treat", J, 
            constrained.single = "full", 
	    constrained.multi = TRUE, fit.start = "lm", 
	    n.draws = 10000, burnin = 5000, thin = 0, delta.start, 
	    psi.start, Sigma.start, Phi.start, delta.mu0, psi.mu0, 
	    delta.A0, psi.A0, Sigma.df, Sigma.scale, Phi.df, Phi.scale, 
            delta.tune, psi.tune, gamma.tune, zeta.tune, 
	    formula.mixed, group.mixed, 
	    verbose = TRUE, sensitive.model = "logit", 
	    df = 5, endorse.options, ...)
} 

\arguments{ 
\item{formula}{An object of class "formula": a symbolic description of the model to be fitted.} 
\item{data}{A data frame containing the variables in the model}
\item{treat}{Name of treatment indicator as a string. For single sensitive item models, this refers to a binary indicator, and for multiple sensitive item models it refers to a multi-valued variable with zero representing the control condition. This can be an integer (with 0 for the control group) or a factor (with "control" for the control group).}
\item{J}{Number of non-sensitive (control) survey items. This will be set automatically to the maximum value of the outcome variable in the treatment group if no input is sent by the user.}
\item{constrained.single}{A string indicating whether the control group parameters are constrained to be equal in the single sensitive item design, either setting all parameters to be equal (\code{full}) or only the intercept (\code{intercept}). If neither, set to \code{none}.}
\item{constrained.multi}{A logical value indicating whether the non-sensitive item count is included as a predictor in the sensitive item fits for the multiple sensitive item design.}
\item{fit.start}{Fit method for starting values. The options are \code{lm}, \code{glm}, \code{nls}, and \code{ml}, which use OLS, logistic regression, non-linear least squares, and maximum likelihood estimation to generate starting values, respectively. The default is \code{lm}.}
\item{n.draws}{Number of MCMC iterations after the burnin.}
\item{burnin}{The number of initial MCMC iterations that are discarded.}
\item{thin}{The interval of thinning, in which every other (\code{thin} = 1) or more iterations are discarded in the output object}
\item{delta.start}{Optional starting values for the sensitive item fit. This should be a vector with the length of the number of covariates for the single sensitive item design, and either a vector or a list with a vector of starting values for each of the sensitive items. The default runs an \code{ictreg} fit with the method set by the \code{fit.start} option.}
\item{psi.start}{Optional starting values for the control items fit. This should be a vector of length the number of covariates for the constrained models. The default runs an \code{ictreg} fit with the method set by the \code{fit.start} option.}
\item{Sigma.start}{Optional starting values for Sigma parameter for mixed effects models for sensitive item.}
\item{Phi.start}{Optional starting values for the Phi parameter for mixed effects models for control item.}
\item{delta.mu0}{Optional vector of prior means for the sensitive item fit parameters, a vector of length the number of covariates.}
\item{psi.mu0}{Optional vector of prior means for the control item fit parameters, a vector of length the number of covariates.}
\item{delta.A0}{Optional matrix of prior precisions for the sensitive item fit parameters, a matrix of dimension the number of covariates.}
\item{psi.A0}{Optional matrix of prior precisions for the control items fit parameters, a matrix of dimension the number of covariates.}
\item{Sigma.df}{Optional prior degrees of freedom parameter for mixed effects models for sensitive item.}
\item{Sigma.scale}{Optional prior scale parameter for mixed effects models for sensitive item.}
\item{Phi.df}{Optional prior degress of freedom parameter for mixed effects models for control item.}
\item{Phi.scale}{Optional prior scale parameter for mixed effects models for control item.}
\item{delta.tune}{A required vector of tuning parameters for the Metropolis algorithm for the sensitive item fit. This must be set and refined by the user until the acceptance ratios are approximately .4 (reported in the output).}
\item{psi.tune}{A required vector of tuning parameters for the Metropolis algorithm for the control item fit. This must be set and refined by the user until the acceptance ratios are approximately .4 (reported in the output).}
\item{gamma.tune}{An optional vector of tuning parameters for the Metropolis algorithm for the control item fit for the random effects. This can be set and refined by the user until the acceptance ratios are approximately .4 (reported in the output).}
\item{zeta.tune}{An optional vector of tuning parameters for the Metropolis algorithm for the sensitive item fit for the random effects. This can be set and refined by the user until the acceptance ratios are approximately .4 (reported in the output).}
\item{formula.mixed}{To specify a mixed effects model, include this formula object for the group-level fit. ~1 allows intercepts to vary, and including covariates in the formula allows the slopes to vary also.}
\item{group.mixed}{A numerical group indicator specifying which group each individual belongs to for a mixed effects model.}
\item{verbose}{A logical value indicating whether model diagnostics are printed out during fitting.}
\item{sensitive.model}{A logical value indicating which model is used for the sensitive item fit,  logistic regression (\code{logit}, default), robit regression (\code{robit}), or probit regression (\code{probit}).}
\item{df}{Degrees of freedom for the robit model for the sensitive item fit, only used if \code{robit} is set to \code{TRUE}.}
\item{endorse.options}{A list of inputs and options for running the combined list experiment and endorsement experiment model. Options documented more fully in \code{endorse} package.}
\item{...}{further arguments to be passed to NLS regression commands.}
}

\details{This function allows the user to perform regression analysis on data from the item count technique, also known as the list experiment and the unmatched count technique using a Bayesian MCMC algorithm.

Unlike the maximum likelihood and least squares estimators in the \code{ictreg} function, the Metropolis algorithm for the Bayesian MCMC estimators in this function must be tuned to work correctly. The \code{delta.tune} and \code{psi.tune} are required, and the values, one for each estimated parameter, will need to be manipulated. The output of the \code{ictregBayes} function, and of the \code{summary} function run on an \code{ictregBayes} object display the acceptance ratios from the Metropolis algorithm. If these values are far from 0.4, the tuning parameters should be changed until the ratios approach 0.4.

For the single sensitive item design, the model can constrain all control parameters to be equal (\code{constrained = "full"}), or just the intercept (\code{constrained = "intercept"}) or all the control fit parameters can be allowed to vary across the potential sensitive item values (\code{constrained = "none"}).

For the multiple sensitive item design, the model can include the estimated number of affirmative responses to the control items as a covariate in the sensitive item model fit (\code{constrained} set to \code{TRUE}) or exclude it (\code{FALSE}).

The function also allows the user to perform combined list experiment and endorsement experiment regression. Setting \code{endorse.options} to a list with the options from the \code{endorse} package for endorsement experiment regression, the function will return the combined model in which the relationship between covariates and the sensitive item in the list experiment model is set to be identical to the relationship between covariates and support for endorsers in the endorsement experiment model. For more details on endorsement experiment regression, see the help for the \code{endorse} package.

Convergence is at times difficult to achieve, so we recommend running multiple chains from overdispersed starting values by, for example, running an MLE or linear model using the ictreg() function, and then generating a set of overdispersed starting values using those estimates and their estimated variance-covariance matrix. An example is provided below for each of the possible designs. Running \code{summary()} after such a procedure will output the Gelman-Rubin convergence statistics in addition to the estimates. If the G-R statistics are all below 1.1, the model is said to have converged.

}

\value{\code{ictregBayes} returns an object of class "ictregBayes".  The function \code{summary} is used to obtain a table of the results, using the \code{coda} package. Two attributes are also included, the data ("x"), the call ("call"), which can be extracted using the command, e.g., attr(ictregBayes.object, "x").

\item{mcmc}{an object of class "mcmc" that can be analyzed using the \code{coda} package.}
\item{x}{the design matrix}
\item{multi}{a logical value indicating whether the data included multiple sensitive items.}
\item{constrained}{a logical or character value indicating whether the control group parameters are constrained to be equal in the single sensitive item design, and whether the non-sensitive item count is included as a predictor in the sensitive item fits for the multiple sensitive item design.}
\item{delta.start}{Optional starting values for the sensitive item fit. This should be a vector with the length of the number of covariates. The default runs an \code{ictreg} fit with the method set by the \code{fit.start} option.}
\item{psi.start}{Optional starting values for the control items fit. This should be a vector of length the number of covariates. The default runs an \code{ictreg} fit with the method set by the \code{fit.start} option.}
\item{delta.mu0}{Optional vector of prior means for the sensitive item fit parameters, a vector of length the number of covariates.}
\item{psi.mu0}{Optional vector of prior means for the control item fit parameters, a vector of length the number of covariates.}
\item{delta.A0}{Optional matrix of prior precisions for the sensitive item fit parameters, a matrix of dimension the number of covariates.}
\item{psi.A0}{Optional matrix of prior precisions for the control items fit parameters, a matrix of dimension the number of covariates.}
\item{delta.tune}{A required vector of tuning parameters for the Metropolis algorithm for the sensitive item fit. This must be set and refined by the user until the acceptance ratios are approximately .4 (reported in the output).}
\item{psi.tune}{A required vector of tuning parameters for the Metropolis algorithm for the control item fit. This must be set and refined by the user until the acceptance ratios are approximately .4 (reported in the output).}
\item{J}{Number of non-sensitive (control) survey items set by the user or detected.}
\item{treat.labels}{a vector of the names used by the \code{treat} vector for the sensitive item or items. This is the names from the \code{treat} indicator if it is a factor, or the number of the item if it is numeric.} 
\item{control.label}{a vector of the names used by the \code{treat} vector for the control items. This is the names from the \code{treat} indicator if it is a factor, or the number of the item if it is numeric.}
\item{call}{the matched call}

If the data includes multiple sensitive items, the following object is also included:
\item{treat.values}{a vector of the values used in the \code{treat} vector for the sensitive items, either character or numeric depending on the class of \code{treat}. Does not include the value for the control status}
}

\references{
Blair, Graeme and Kosuke Imai. (2012) ``Statistical Analysis of List Experiments." 
  Political Analysis, Vol. 20, No 1 (Winter). available at 
  \url{http://imai.princeton.edu/research/listP.html}

Imai, Kosuke. (2011) ``Multivariate Regression Analysis for the Item
  Count Technique.''  Journal of the American Statistical Association, 
  Vol. 106, No. 494 (June), pp. 407-416. available at
  \url{http://imai.princeton.edu/research/list.html}

Blair, Graeme, Jason Lyall and Kosuke Imai. (2013) ``Comparing and Combining
  List and Experiments: Evidence from Afghanistan." 
  Working paper. available at 
  \url{http://imai.princeton.edu/research/comp.html}}

\author{Graeme Blair, Princeton University,
\email{gblair@princeton.edu} and Kosuke Imai, Princeton University,
\email{kimai@princeton.edu} 
}

\examples{

data(race)

\dontrun{

## Multiple chain MCMC list experiment regression
## starts with overdispersed MLE starting values

## Standard single sensitive-item design

## Control item parameters fully constrained

mle.estimates <- ictreg(y ~ male + college + age + south, data = race)

draws <- mvrnorm(n = 3, mu = coef(mle.estimates), 
  Sigma = vcov(mle.estimates) * 9)

bayesDraws.1 <- ictregBayes(y ~ male + college + age + south, data = race, 
  delta.start = draws[1, 1:5], psi.start = draws[1, 6:10], burnin = 10000, 
  n.draws = 100000, delta.tune = diag(.002, 5), psi.tune = diag(.00025, 5),
  constrained.single = "full")

bayesDraws.2 <- ictregBayes(y ~ male + college + age + south, data = race, 
  delta.start = draws[2, 1:5], psi.start = draws[2, 6:10], burnin = 10000, 
  n.draws = 100000, delta.tune = diag(.002, 5), psi.tune = diag(.00025, 5),
  constrained.single = "full")

bayesDraws.3 <- ictregBayes(y ~ male + college + age + south, data = race, 
  delta.start = draws[3, 1:5], psi.start = draws[3, 6:10], burnin = 10000, 
  n.draws = 100000, delta.tune = diag(.002, 5), psi.tune = diag(.00025, 5),
  constrained.single = "full")

bayesSingleConstrained <- as.list(bayesDraws.1, bayesDraws.2, bayesDraws.3)

summary(bayesSingleConstrained)

## Control item parameters unconstrained

mle.estimates <- ictreg(y ~ male + college + age + south, data = race,
  constrained = FALSE)

draws <- mvrnorm(n = 3, mu = coef(mle.estimates), 
  Sigma = vcov(mle.estimates) * 9)

bayesDraws.1 <- ictregBayes(y ~ male + college + age + south, data = race, 
  delta.start = draws[1, 1:5], psi.start = list(psi0 = draws[1, 6:10], 
  psi1 = draws[1, 11:15]), burnin = 10000, n.draws = 100000, 
  delta.tune = diag(.002, 5), 
  psi.tune = list(psi0 = diag(.0017, 5), psi1 = diag(.00005, 5)),
  constrained.single = "none")

bayesDraws.2 <- ictregBayes(y ~ male + college + age + south, data = race, 
  delta.start = draws[2, 1:5], psi.start = list(psi0 = draws[2, 6:10], 
  psi1 = draws[2, 11:15]), burnin = 10000, n.draws = 100000, 
  delta.tune = diag(.002, 5), 
  psi.tune = list(psi0 = diag(.0017, 5), psi1 = diag(.00005, 5)),
  constrained.single = "none")

bayesDraws.3 <- ictregBayes(y ~ male + college + age + south, data = race, 
  delta.start = draws[3, 1:5], psi.start = list(psi0 = draws[3, 6:10], 
  psi1 = draws[3, 11:15]), burnin = 10000, n.draws = 100000, 
  delta.tune = diag(.002, 5), 
  psi.tune = list(psi0 = diag(.0017, 5), psi1 = diag(.00005, 5)),
  constrained.single = "none")

bayesSingleUnconstrained <- as.list(bayesDraws.1, bayesDraws.2, bayesDraws.3)

summary(bayesSingleUnconstrained)

## Control item parameters constrained except intercept

mle.estimates <- ictreg(y ~ male + college + age + south, data = race,
  constrained = TRUE)

draws <- mvrnorm(n = 3, mu = coef(mle.estimates), 
  Sigma = vcov(mle.estimates) * 9)

bayesDraws.1 <-  ictregBayes(y ~ male + college + age + south, data = race, 
  delta.start = draws[1, 1:5], psi.start = c(draws[1, 6:10],0),
  burnin = 10000, n.draws = 100000, delta.tune = diag(.002, 5),
  psi.tune = diag(.0004, 6), constrained.single = "intercept")

bayesDraws.2 <-  ictregBayes(y ~ male + college + age + south, data = race, 
  delta.start = draws[2, 1:5], psi.start = c(draws[2, 6:10],0),
  burnin = 10000, n.draws = 100000, delta.tune = diag(.002, 5),
  psi.tune = diag(.0004, 6), constrained.single = "intercept")

bayesDraws.3 <-  ictregBayes(y ~ male + college + age + south, data = race, 
  delta.start = draws[3, 1:5], psi.start = c(draws[3, 6:10],0),
  burnin = 10000, n.draws = 100000, delta.tune = diag(.002, 5),
  psi.tune = diag(.0004, 6), constrained.single = "intercept")

bayesSingleInterceptOnly <- as.list(bayesDraws.1, bayesDraws.2, bayesDraws.3)

summary(bayesSingleInterceptOnly)

## Multiple sensitive item design

## Constrained (estimated control item count not included in sensitive fit)

mle.estimates.multi <- ictreg(y ~ male + college + age + south, data = multi,
  constrained = TRUE)

draws <- mvrnorm(n = 3, mu = coef(mle.estimates.multi), 
  Sigma = vcov(mle.estimates.multi) * 9)

bayesMultiDraws.1 <- ictregBayes(y ~ male + college + age + south, 
  data = multi, delta.start = list(draws[1, 6:10], draws[1, 11:15]), 
  psi.start = draws[1, 1:5], burnin = 10000, n.draws = 100000, 
  delta.tune = diag(.002, 5), psi.tune = diag(.001, 5), 
  constrained.multi = TRUE)

bayesMultiDraws.2 <- ictregBayes(y ~ male + college + age + south, 
  data = multi, delta.start = list(draws[2, 6:10], draws[2, 11:15]), 
  psi.start = draws[2, 1:5], burnin = 10000, n.draws = 100000, 
  delta.tune = diag(.002, 5), psi.tune = diag(.001, 5), 
  constrained.multi = TRUE)

bayesMultiDraws.3 <- ictregBayes(y ~ male + college + age + south, 
  data = multi, delta.start = list(draws[3, 6:10], draws[3, 11:15]), 
  psi.start = draws[3, 1:5], burnin = 10000, n.draws = 100000, 
  delta.tune = diag(.002, 5), psi.tune = diag(.001, 5), 
  constrained.multi = TRUE)

bayesMultiConstrained <- as.list(bayesMultiDraws.1, bayesMultiDraws.2, 
  bayesMultiDraws.3)

summary(bayesMultiConstrained)

## Unconstrained (estimated control item count is included in sensitive fit)

mle.estimates.multi <- ictreg(y ~ male + college + age + south, data = multi,
  constrained = FALSE)

draws <- mvrnorm(n = 3, mu = coef(mle.estimates.multi), 
  Sigma = vcov(mle.estimates.multi) * 9)

bayesMultiDraws.1 <- ictregBayes(y ~ male + college + age + south, 
  data = multi, delta.start = list(draws[1, 6:10], draws[1, 11:15]), 
  psi.start = draws[1, 1:5], burnin = 50000, n.draws = 300000, 
  delta.tune = diag(.0085, 6), psi.tune = diag(.00025, 5), 
  constrained.multi = FALSE)

bayesMultiDraws.2 <- ictregBayes(y ~ male + college + age + south, 
  data = multi, delta.start = list(draws[2, 6:10], draws[2, 11:15]), 
  psi.start = draws[2, 1:5], burnin = 50000, n.draws = 300000, 
  delta.tune = diag(.0085, 6), psi.tune = diag(.00025, 5), 
  constrained.multi = FALSE)

bayesMultiDraws.3 <- ictregBayes(y ~ male + college + age + south, 
  data = multi, delta.start = list(draws[3, 6:10], draws[3, 11:15]), 
  psi.start = draws[3, 1:5], burnin = 50000, n.draws = 300000, 
  delta.tune = diag(.0085, 6), psi.tune = diag(.00025, 5), 
  constrained.multi = FALSE)

bayesMultiUnconstrained <- as.list(bayesMultiDraws.1, bayesMultiDraws.2, 
  bayesMultiDraws.3)

summary(bayesMultiUnconstrained)

## Mixed effects models

## Varying intercepts

mle.estimates <- ictreg(y ~ male + college + age + south, data = race)

draws <- mvrnorm(n = 3, mu = coef(mle.estimates), 
  Sigma = vcov(mle.estimates) * 9)

bayesDraws.1 <- ictregBayes(y ~ male + college + age + south, data = race, 
  delta.start = draws[1, 1:5], psi.start = draws[1, 6:10], burnin = 100, 
  n.draws = 1000, delta.tune = diag(.002, 5), psi.tune = diag(.00025, 5),
  constrained.single = "full", group.mixed = "state", formula.mixed = ~ 1)

bayesDraws.2 <- ictregBayes(y ~ male + college + age + south, data = race, 
  delta.start = draws[2, 1:5], psi.start = draws[2, 6:10], burnin = 10000, 
  n.draws = 100000, delta.tune = diag(.002, 5), psi.tune = diag(.00025, 5),
  constrained.single = "full", group.mixed = "state", formula.mixed = ~ 1)

bayesDraws.3 <- ictregBayes(y ~ male + college + age + south, data = race, 
  delta.start = draws[3, 1:5], psi.start = draws[3, 6:10], burnin = 10000, 
  n.draws = 100000, delta.tune = diag(.002, 5), psi.tune = diag(.00025, 5),
  constrained.single = "full", group.mixed = "state", formula.mixed = ~ 1)

bayesMixed <- as.list(bayesDraws.1, bayesDraws.2, bayesDraws.3)

summary(bayesMixed)

}
}
 
\seealso{
\code{\link{predict.ictreg}} for fitted values
}

\keyword{models}
\keyword{regression}