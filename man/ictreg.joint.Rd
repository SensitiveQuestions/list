\name{ictreg.joint} 
\alias{ictreg.joint}
\title{Item Count Technique: Outcome Models}
\description{
Function to conduct multivariate regression analyses of survey data with the item count technique, also known as the list experiment, using predicted responses from list experiments as predictors in outcome regression models.
}
\usage{
ictreg.joint(formula, data = parent.frame(), treat = "treat", J, 
                    outcome = "outcome", outcome.reg = "logistic",
                    constrained = FALSE,
                    maxIter = 1000)
}
\arguments{
\item{formula}{An object of class "formula": a symbolic description of the model to be fitted.}
\item{data}{A data frame containing the variables in the model}
\item{treat}{Name of treatment indicator as a string. For single sensitive item models, this refers to a binary indicator, and for multiple sensitive item models it refers to a multi-valued variable with zero representing the control condition. This can be an integer (with 0 for the control group) or a factor (with "control" for the control group).}
\item{J}{Number of non-sensitive (control) survey items.}
\item{outcome}{Name of outcome indicator as a string.}
\item{outcome.reg}{Model for outcome regression. Options are "logistic" or "linear;" default is "logistic".}
\item{constrained}{A logical value indicating whether the control group parameters are constrained to be equal. Default is FALSE.}
\item{maxIter}{Maximum number of iterations for the Expectation-Maximization algorithm of the ML estimation.  The default is 1000.}
}

\details{This function allows the user to perform regression analysis on survey data with the item count technique, also known as the list experiment, using predicted responses from list experiments as predictors in outcome regression models.
}

\value{\code{ictreg.joint} returns an object of class "ictreg.joint". The function \code{summary} is used to obtain a table of the results. The object \code{ictreg.joint} is a list that contains the following components. 
\item{par.treat}{point estimate for effect of covariate on item count fitted on treatment group}
\item{se.treat}{standard error for estimate of effect of covariate on item count fitted on treatment group}
\item{par.control}{point estimate for effect of covariate on item count fitted on control group}
\item{se.control}{standard error for estimate of effect of covariate on item count fitted on control group}
\item{par.outcome}{point estimate for effect of covariate and sensitive item on outcome}
\item{se.outcome}{standard error for estimate of effect of covariate and sensitive item on outcome}
\item{coef.names}{variable names as defined in the data frame}
\item{constrained}{call indicating whether the constrained model is used}
\item{call}{the matched call}
\item{data}{the \code{data} argument}
\item{outcome.reg}{the \code{outcome.reg} argument}
\item{x}{the design matrix}
\item{y}{the response vector}
\item{treat}{the vector indicating treatment status}
\item{J}{Number of non-sensitive (control) survey items set by the user or detected.}
\item{treat.labels}{a vector of the names used by the \code{treat} vector for the sensitive item or items. This is the names from the \code{treat} indicator if it is a factor, or the number of the item if it is numeric.} 
\item{control.label}{a vector of the names used by the \code{treat} vector for the control items. This is the names from the \code{treat} indicator if it is a factor, or the number of the item if it is numeric.}
}

\references{
Imai, Kosuke, Bethany Park, and Kenneth F. Greene. (2014) ``Using the Predicted Responses 
from List Experiments as Explanatory Variables in Regression Models.'' available 
at \url{http://imai.princeton.edu/research/files/listExp.pdf}
}


\examples{
\dontrun{
data(mexico)
loyal <- mexico[mexico$mex.loyal == 1,]
notloyal <- mexico[mexico$mex.loyal == 0,]

## Logistic outcome regression
## (effect of vote-selling on turnout)
## This replicates Table 4 in Imai et al. 2014

loyalreg <- ictreg.joint(formula = mex.y.all ~ mex.male + mex.age + mex.age2 + mex.education +  
                         mex.interest + mex.married +
                         mex.wealth + mex.urban + mex.havepropoganda + mex.concurrent, data = loyal,
                         treat = "mex.t", outcome = "mex.votecard", J = 3, constrained = TRUE,
                         outcome.reg = "logistic", maxIter = 1000)
summary(loyalreg)

## Linear outcome regression
## (effect of vote-selling on candidate approval)
## This replicates Table 5 in Imai et al. 2014

approvalreg <- ictreg.joint(formula = mex.y.all ~ mex.male + mex.age + mex.age2 +
                            mex.education +
                            mex.interest + mex.married +
                            mex.urban + 
                            mex.cleanelections + mex.cleanelectionsmiss +
                            mex.havepropoganda +
                            mex.wealth + mex.northregion +
                            mex.centralregion + mex.metro + mex.pidpriw2 + 
			    mex.pidpanw2 + mex.pidprdw2,
                            data = mexico, treat = "mex.t", outcome = "mex.epnapprove",
                            J = 3, constrained = TRUE,
                            outcome.reg = "linear", maxIter = 1000)


summary(approvalreg)
}
}
