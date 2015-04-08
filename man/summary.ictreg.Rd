\name{summary.ictreg} 
\alias{summary.ictreg} 
\title{Summary Method for the Item Count Technique} 
\description{ 
Function to summarize results from list experiment regression based on the ictreg() function, and to produce proportions of liars estimates.
} 
\usage{

  \method{summary}{ictreg}(object, boundary.proportions = FALSE, n.draws = 10000, ...) 

} 

\arguments{ 
\item{object}{Object of class inheriting from "ictreg"}
\item{boundary.proportions}{A switch indicating whether, for models with ceiling effects, floor effects, or both (indicated by the \code{floor = TRUE}, \code{ceiling = TRUE} options in ictreg), the conditional probability of lying and the population proportions of liars are calculated.}
\item{n.draws}{For quasi-Bayesian approximation based predictions, specify the number of Monte Carlo draws.}
\item{...}{further arguments to be passed to or from other methods.}
}

\details{
\code{predict.ictreg} produces a summary of the results from an \code{ictreg} object. It displays the coefficients, standard errors, and fit statistics for any model from \code{ictreg}.

\code{predict.ictreg} also produces estimates of the conditional probability of lying and of the population proportion of liars for boundary models from \code{ictreg()} if \code{ceiling = TRUE} or \code{floor = TRUE}. 

The conditional probability of lying for the ceiling model is the probability that a respondent with true affirmative views of all the sensitive and non-sensitive items lies and responds negatively to the sensitive item. The conditional probability for the floor model is the probability that a respondent lies to conceal her true affirmative views of the sensitive item when she also holds true negative views of all the non-sensitive items. In both cases, the respondent may believe her privacy is not protected, so may conceal her true affirmative views of the sensitive item. 
}

\seealso{\code{\link{ictreg}} for model fitting}

\references{
Blair, Graeme and Kosuke Imai. (2012) ``Statistical Analysis of List Experiments." 
  Political Analysis, Vol. 20, No 1 (Winter). available at 
  \url{http://imai.princeton.edu/research/listP.html}

Imai, Kosuke. (2011) ``Multivariate Regression Analysis for the Item
  Count Technique.''  Journal of the American Statistical Association, 
  Vol. 106, No. 494 (June), pp. 407-416. available at
  \url{http://imai.princeton.edu/research/list.html}}

\author{Graeme Blair, Princeton University,
\email{gblair@princeton.edu} and Kosuke Imai, Princeton University,
\email{kimai@princeton.edu} 
}

\examples{
data(race)
\dontrun{
# Fit standard design ML model with ceiling effects
# Replicates Table 7 Columns 3-4 in Blair and Imai (2012)

ceiling.results <- ictreg(y ~ age + college + male + south, treat = "treat", 
		   	  J = 3, data = affirm, method = "ml", fit.start = "nls",
			  ceiling = TRUE, ceiling.fit = "bayesglm",
			  ceiling.formula = ~ age + college + male + south)

# Summarize fit object and generate conditional probability 
# of ceiling liars the population proportion of ceiling liars,
# both with standard errors.
# Replicates Table 7 Columns 3-4 last row in Blair and Imai (2012)

summary(ceiling.results, boundary.proportions = TRUE)
}

}
 
\keyword{models}
\keyword{regression}