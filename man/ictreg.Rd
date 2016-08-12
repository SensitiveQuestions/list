\name{ictreg} 
\alias{ictreg} 
\alias{ict}
\alias{list}
\title{Item Count Technique} 
\description{ 
Function to conduct multivariate regression analyses of survey data with the item count technique, also known as the list experiment and the unmatched count technique.
} 
\usage{
ictreg(formula, data = parent.frame(), treat = "treat", J, method = "ml", 
       weights, h = NULL, group = NULL, matrixMethod = "efficient",
       overdispersed = FALSE, constrained = TRUE, floor = FALSE, 
       ceiling = FALSE, ceiling.fit = "glm", floor.fit = "glm", 
       ceiling.formula = ~ 1, floor.formula = ~ 1, fit.start = "lm", 
       fit.nonsensitive = "nls", multi.condition = "none", 
       maxIter = 5000, verbose = FALSE, ...)
} 

\arguments{ 
\item{formula}{An object of class "formula": a symbolic description of the model to be fitted.} 
\item{data}{A data frame containing the variables in the model}
\item{treat}{Name of treatment indicator as a string. For single sensitive item models, this refers to a binary indicator, and for multiple sensitive item models it refers to a multi-valued variable with zero representing the control condition. This can be an integer (with 0 for the control group) or a factor (with "control" for the control group).}
\item{J}{Number of non-sensitive (control) survey items.}
\item{method}{ Method for regression, either \code{ml} for the Maximum Likelihood (ML) estimation with the Expectation-Maximization algorithm; \code{lm} for linear model estimation; or \code{nls} for the Non-linear Least Squares (NLS) estimation with the two-step procedure.}
\item{weights}{Name of the weights variable as a string, if weighted regression is desired. Not implemented for the ceiling/floor models, multiple sensitive item design, or for the modified design.}
\item{h}{Auxiliary data functionality. Optional named numeric vector with length equal to number of groups. Names correspond to group labels and values correspond to auxiliary moments.}
\item{group}{Auxiliary data functionality. Optional character vector of group labels with length equal to number of observations.}
\item{matrixMethod}{Auxiliary data functionality. Procedure for estimating optimal weighting matrix for generalized method of moments. One of "efficient" for two-step feasible and "cue" for continuously updating. Default is "efficient". Only relevant if \code{h} and \code{group} are specified.}
\item{overdispersed}{Indicator for the presence of overdispersion. If \code{TRUE}, the beta-binomial model is used in the EM algorithm, if \code{FALSE} the binomial model is used. Not relevant for the \code{NLS} or \code{lm} methods.}
\item{constrained}{A logical value indicating whether the control group parameters are constrained to be equal.  Not relevant for the \code{NLS} or \code{lm} methods}
\item{floor}{A logical value indicating whether the floor liar model should be used to adjust for the possible presence of respondents dishonestly reporting a negative preference for the sensitive item among those who hold negative views of all the non-sensitive items.}
\item{ceiling}{A logical value indicating whether the ceiling liar model should be used to adjust for the possible presence of respondents dishonestly reporting a negative preference for the sensitive item among those who hold affirmative views of all the non-sensitive items.}
\item{ceiling.fit}{Fit method for the M step in the EM algorithm used to fit the ceiling liar model. \code{glm} uses standard logistic regression, while \code{bayesglm} uses logistic regression with a weakly informative prior over the parameters.}
\item{floor.fit}{Fit method for the M step in the EM algorithm used to fit the floor liar model. \code{glm} uses standard logistic regression, while \code{bayesglm} uses logistic regression with a weakly informative prior over the parameters.}
\item{ceiling.formula}{Covariates to include in ceiling liar model. These must be a subset of the covariates used in \code{formula}.}
\item{floor.formula}{Covariates to include in floor liar model. These must be a subset of the covariates used in \code{formula}.}
\item{fit.start}{Fit method for starting values for standard design \code{ml} model. The options are \code{lm}, \code{glm}, and \code{nls}, which use OLS, logistic regression, and non-linear least squares to generate starting values, respectively. The default is \code{nls}.}
\item{fit.nonsensitive}{Fit method for the non-sensitive item fit for the \code{nls} method and the starting values for the \code{ml} method for the \code{modified} design. Options are \code{glm} and \code{nls}, and the default is \code{nls}.}
\item{multi.condition}{For the multiple sensitive item design, covariates representing the estimated count of affirmative responses for each respondent can be included directly as a level variable by choosing \code{level}, or as indicator variables for each value but one by choosing \code{indicators}. The default is \code{none}.}
\item{maxIter}{Maximum number of iterations for the Expectation-Maximization algorithm of the ML estimation.  The default is 5000.} 
\item{verbose}{ a logical value indicating whether model diagnostics are printed out during fitting.}
\item{...}{further arguments to be passed to NLS regression commands.}
}

\details{This function allows the user to perform regression analysis on data from the item count technique, also known as the list experiment and the unmatched count technique.

Three list experiment designs are accepted by this function: the standard design; the multiple sensitive item standard design; and the modified design proposed by Corstange (2009). 

For the standard design, three methods are implemented in this function: the linear model; the Maximum Likelihood (ML) estimation for the Expectation-Maximization (EM) algorithm; the nonlinear least squares (NLS) estimation with the two-step procedure both proposed in Imai (2010); and the Maximum Likelihood (ML) estimator in the presence of two types of dishonest responses, "ceiling" and "floor" liars. The ceiling model, floor model, or both, as described in Blair and Imai (2010) can be activated by using the \code{ceiling} and \code{floor} options. The constrained and unconstrained ML models presented in Imai (2010) are available through the \code{constrained} option, and the user can specify if overdispersion is present in the data for the no liars models using the \code{overdispersed} option to control whether a beta-binomial or binomial model is used in the EM algorithm to model the item counts.

The modified design and the multiple sensitive item design are automatically detected by the function, and only the binomial model without overdispersion is available. }

\value{\code{ictreg} returns an object of class "ictreg".  The function \code{summary} is used to obtain a table of the results.  The object \code{ictreg} is a list that contains the following components.  Some of these elements are not available depending on which method is used (\code{lm}, \code{nls} or \code{ml}), which design is used (\code{standard}, \code{modified}), whether multiple sensitive items are include (\code{multi}), and whether the constrained model is used (\code{constrained = TRUE}).

\item{par.treat}{point estimate for effect of covariate on item count fitted on treatment group}
\item{se.treat}{standard error for estimate of effect of covariate on item count fitted on treatment group}
\item{par.control}{point estimate for effect of covariate on item count fitted on control group}
\item{se.control}{standard error for estimate of effect of covariate on item count fitted on control group}
\item{coef.names}{variable names as defined in the data frame}
\item{design}{call indicating whether the \code{standard} design as proposed in Imai (2010) or thee \code{modified} design as proposed in Corstange (2009) is used}
\item{method}{call of the method used}
\item{overdispersed}{call indicating whether data is overdispersed}
\item{constrained}{call indicating whether the constrained model is used}
\item{boundary}{call indicating whether the floor/ceiling boundary models are used}
\item{multi}{indicator for whether multiple sensitive items were included in the data frame}
\item{call}{the matched call}
\item{data}{the \code{data} argument}
\item{x}{the design matrix}
\item{y}{the response vector}
\item{treat}{the vector indicating treatment status}
\item{J}{Number of non-sensitive (control) survey items set by the user or detected.}
\item{treat.labels}{a vector of the names used by the \code{treat} vector for the sensitive item or items. This is the names from the \code{treat} indicator if it is a factor, or the number of the item if it is numeric.} 
\item{control.label}{a vector of the names used by the \code{treat} vector for the control items. This is the names from the \code{treat} indicator if it is a factor, or the number of the item if it is numeric.}

For the maximum likelihood models, an additional output object is included:
\item{pred.post}{posterior predicted probability of answering "yes" to the sensitive item. The weights from the E-M algorithm.}

For the floor/ceiling models, several additional output objects are included:
\item{ceiling}{call indicating whether the assumption of no ceiling liars is relaxed, and ceiling parameters are estimated}
\item{par.ceiling}{point estimate for effect of covariate on whether respondents who answered affirmatively to all non-sensitive items and hold a true affirmative opinion toward the sensitive item lied and reported a negative response to the sensitive item }
\item{se.ceiling}{standard error for estimate for effect of covariate on whether respondents who answered affirmatively to all non-sensitive items and hold a true affirmative opinion toward the sensitive item lied and reported a negative response to the sensitive item}
\item{floor}{call indicating whether the assumption of no floor liars is relaxed, and floor parameters are estimated}
\item{par.ceiling}{point estimate for effect of covariate on whether respondents who answered negatively to all non-sensitive items and hold a true affirmative opinion toward the sensitive item lied and reported a negative response to the sensitive item }
\item{se.ceiling}{standard error for estimate for effect of covariate on whether respondents who answered negatively to all non-sensitive items and hold a true affirmative opinion toward the sensitive item lied and reported a negative response to the sensitive item}
\item{coef.names.ceiling}{variable names from the ceiling liar model fit, if applicable}
\item{coef.names.floor}{variable names from the floor liar model fit, if applicable}

For the multiple sensitive item design, the \code{par.treat} and \code{se.treat} vectors are returned as lists of vectors, one for each sensitive item.

For the unconstrained model, the \code{par.control} and \code{se.control} output is replaced by:
\item{par.control.phi0}{point estimate for effect of covariate on item count fitted on treatment group}
\item{se.control.phi0}{standard error for estimate of effect of covariate on item count fitted on treatment group}
\item{par.control.phi1}{point estimate for effect of covariate on item count fitted on treatment group}
\item{se.control.phi1}{standard error for estimate of effect of covariate on item count fitted on treatment group}

Depending upon the estimator requested by the user, model fit statistics are also included:
\item{llik}{the log likelihood of the model, if \code{ml} is used}
\item{resid.se}{the residual standard error, if \code{nls} or \code{lm} are used. This will be a scalar if the standard design was used, and a vector if the multiple sensitive item design was used}
\item{resid.df}{the residual degrees of freedom, if \code{nls} or \code{lm} are used. This will be a scalar if the standard design was used, and a vector if the multiple sensitive item design was used}

When using the auxiliary data functionality, the following objects are included:
\item{aux}{logical value indicating whether estimation incorporates auxiliary moments}
\item{nh}{integer count of the number of auxiliary moments}
\item{wm}{procedure used to estimate the optimal weight matrix}
\item{J.stat}{numeric value of the Sargan Hansen overidentifying restriction test statistic}
\item{overid.p}{corresponding p-value for the Sargan Hansen test}

}

\references{
Blair, Graeme and Kosuke Imai. (2012) ``Statistical Analysis of List Experiments." 
  Political Analysis. Forthcoming. available at 
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

set.seed(1)

# Calculate list experiment difference in means

diff.in.means.results <- ictreg(y ~ 1, data = race, 
	       	      treat = "treat", J=3, method = "lm")

summary(diff.in.means.results)

# Fit linear regression
# Replicates Table 1 Columns 1-2 Imai (2011); note that age is divided by 10

lm.results <- ictreg(y ~ south + age + male + college, data = race, 
	       	      treat = "treat", J=3, method = "lm")

summary(lm.results)

# Fit two-step non-linear least squares regression
# Replicates Table 1 Columns 3-4 Imai (2011); note that age is divided by 10

nls.results <- ictreg(y ~ south + age + male + college, data = race, 
	       	      treat = "treat", J=3, method = "nls")

summary(nls.results)

\dontrun{

# Fit EM algorithm ML model with constraint
# Replicates Table 1 Columns 5-6, Imai (2011); note that age is divided by 10

ml.constrained.results <- ictreg(y ~ south + age + male + college, data = race, 
		       	  	 treat = "treat", J=3, method = "ml", 
				 overdispersed = FALSE, constrained = TRUE)

summary(ml.constrained.results)

# Fit EM algorithm ML model with no constraint
# Replicates Table 1 Columns 7-10, Imai (2011); note that age is divided by 10

ml.unconstrained.results <- ictreg(y ~ south + age + male + college, data = race, 
			    	   treat = "treat", J=3, method = "ml", 
				   overdispersed = FALSE, constrained = FALSE)

summary(ml.unconstrained.results)

# Fit EM algorithm ML model for multiple sensitive items
# Replicates Table 3 in Blair and Imai (2010)

multi.results <- ictreg(y ~ male + college + age + south + south:age, treat = "treat", 
	      	 	J = 3, data = multi, method = "ml", 
			multi.condition = "level")

summary(multi.results)

# Fit standard design ML model
# Replicates Table 7 Columns 1-2 in Blair and Imai (2010)

noboundary.results <- ictreg(y ~ age + college + male + south, treat = "treat",
		      	     J = 3, data = affirm, method = "ml", 
			     overdispersed = FALSE)

summary(noboundary.results)

# Fit standard design ML model with ceiling effects alone
# Replicates Table 7 Columns 3-4 in Blair and Imai (2010)

ceiling.results <- ictreg(y ~ age + college + male + south, treat = "treat", 
		   	  J = 3, data = affirm, method = "ml", fit.start = "nls",
			  ceiling = TRUE, ceiling.fit = "bayesglm",
			  ceiling.formula = ~ age + college + male + south)

summary(ceiling.results)

# Fit standard design ML model with floor effects alone
# Replicates Table 7 Columns 5-6 in Blair and Imai (2010)

floor.results <- ictreg(y ~ age + college + male + south, treat = "treat", 
	      	 	J = 3, data = affirm, method = "ml", fit.start = "glm", 
			floor = TRUE, floor.fit = "bayesglm",
			floor.formula = ~ age + college + male + south)

summary(floor.results)

# Fit standard design ML model with floor and ceiling effects
# Replicates Table 7 Columns 7-8 in Blair and Imai (2010)

both.results <- ictreg(y ~ age + college + male + south, treat = "treat", 
	     	       J = 3, data = affirm, method = "ml", 
		       floor = TRUE, ceiling = TRUE, 
		       floor.fit = "bayesglm", ceiling.fit = "bayesglm",
		       floor.formula = ~ age + college + male + south,
		       ceiling.formula = ~ age + college + male + south)

summary(both.results)

}
}
 
\seealso{
\code{\link{predict.ictreg}} for fitted values
}

\keyword{models}
\keyword{regression}