\name{predict.ictreg.joint} 
\alias{predict.ictreg.joint} 
\title{Predict Method for Item Count Technique, Outcome Regressions} 
\description{ 
Function to calculate predictions and uncertainties of predictions from estimates from multivariate regression analysis of survey data with the item count technique, using predicted responses to list experiments as predictors in outcome regressions.
} 
\usage{

\method{predict}{ictreg.joint}(object, newdata, newdata.diff, se.fit = FALSE, 
   	interval = c("none","confidence"), level = .95, avg = FALSE, 
   	sensitive.value = c("0", "1", "both"), sensitive.diff = FALSE,
	return.draws = FALSE, predict.sensitive = FALSE, ...) 

} 

\arguments{ 
\item{object}{Object of class inheriting from "ictreg.joint"}
\item{newdata}{An optional data frame containing data that will be used to make predictions from. If omitted, the data used to fit the regression are used.}
\item{newdata.diff}{An optional data frame used to compare predictions with predictions from the data in the provided newdata data frame.}
\item{se.fit}{A switch indicating if standard errors are required.}
\item{interval}{Type of interval calculation.}
\item{level}{Significance level for confidence intervals.}
\item{avg}{A switch indicating if the mean prediction and associated statistics across all obserations in the dataframe will be returned instead of predictions for each observation.}
\item{sensitive.value}{User-specified value for the sensitive item.}
\item{sensitive.diff}{A switch indicating if the difference in predictions when the sensitive item = 1 and when the sensitive item = 0 is calculated.}
\item{return.draws}{A switch indicating if the draws from the simulations used to generate predictions will be returned.}
\item{predict.sensitive}{A switch indicating whether predictions from the sensitive item model are generated.}
\item{...}{further arguments to be passed to or from other methods.}
}

\details{ \code{predict.ictreg.joint} produces predicted values, obtained by evaluating the regression function in the
frame newdata (which defaults to model.frame(object)). By using sensitive.value, users must set
the value of z -- the latent response to the sensitive item -- to be either zero or
one, depending on the prediction that the user requires.

Two additional types of mean prediction are also available. The first, if a newdata.diff data frame
is provided by the user, calculates the mean predicted values across two datasets, as well as the
mean difference in predicted value. Standard errors and confidence intervals are also added. For 
newdata.diff predictions, sensitive.value must be set to 1 or 0, not "both" (and sensitive.diff
must also be set to FALSE). 
Users may also set the logical sensitive.diff to TRUE and sensitive.value to "both", which will output 
the mean predicted values across all observations for z = 0 as well as z = 1, in addition to the mean 
difference in predicted value. Standard errors and confidence intervals are also added. For
difference predictions (sensitive.diff and newdata.diff), the option avg must be set to TRUE.

Users can also use the predict.sensitive = TRUE option to generate predictions of responses
to the sensitive item, with standard errors and confidence intervals.

NOTE: In order to generate predictions from user-provided data frames (newdata and newdata.diff),
users MUST run models using \code{ictreg.joint} on data that does not contain any missingness. 
Further, the data frames provided to \code{predict.ictreg.joint} must also not contain any missingness.
}

\value{\code{predict.ictreg.joint} produces a vector of predictions or a matrix of predictions and bounds with column names fit, lwr, and upr if interval is set. If sensitive.value = "both", \code{predict.ictreg.joint} will produce a list, where the first element corresponds to when the sensitive item = 0 and the second element corresponds to when the sensitive item = 1. If sensitive.diff = TRUE, the third element in the list corresponds to the difference (sensitive = 0 subtracted from sensitive = 1). If se.fit is TRUE, a list with the following components is returned:

\item{fit}{vector or matrix as above.}
\item{se.fit}{standard error of prediction(s)}

If return.draws is TRUE, the list includes

\item{draws.predict}{A matrix of draws from a multivariate normal distribution with mean equal to the vector of estimated coefficients from the outcome regression model, and sigma equal to the variance-covariance matrix of the outcome regression model. Rows are observations; colums are 10,000 draws. If sensitive.value = both, will be a list of two elements where each element is a matrix as described; the first matrix will be for when the sensitive item = 0, the second matrix will be for when the sensitive item = 1. If newdata.diff is provided, \code{draws.predict} will be a list of two elements where each element is a matrix as described; the first matrix will correspond to the newdata data frame; the second matrix will correspond to the newdata.diff data frame.}
\item{draws.mean}{The \code{draws.predict} matrix averaged over all observations; a vector of 10,000 draws. If sensitive.value = both, will be a list of two elements where each element is a vector as described; the first matrix will be for when the sensitive item = 0, the second matrix will be for when the sensitive item = 1. If newdata.diff is provided, \code{draws.mean} will be a list of two elements where each element is a matrix as described; the first matrix will correspond to the newdata data frame; the second matrix will correspond to the newdata.diff data frame.}
\item{sens.diff}{If sensitive.diff = TRUE, a vector of 10,000 draws generated from subtracting the first item in \code{draws.mean} from the second item. A vector of 10,000 draws.}

If predict.sensitive = TRUE, the list also includes

\item{fitsens}{a vector of predictions and bounds with column names fit, lwr, and upr if interval is set, generated from the sensitive item model.}
\item{draws.predict.sens}{A matrix of draws from a multivariate normal distribution with mean equal to the vector of estimated coefficients from the sensitive item model, and sigma equal to the variance-covariance matrix of the sensitive item model. Rows are observations; colums are 10,000 draws (only returned if return.draws is TRUE). If newdata.diff is provided, this will be a list of two matrices as described. The first will correspond to newdata, and the second to newdata.diff.}
\item{draws.mean.sens}{The \code{draws.predict.sens} matrix averaged over all observations; a vector of 10,000 draws (only returned if return.draws is TRUE). If newdata.diff is provided, this will be a list of two matrices as described. The first will correspond to newdata, and the second to newdata.diff.}
}

\references{
Imai, Kosuke, Bethany Park, and Kenneth F. Greene. (2014) ``Using the Predicted Responses 
from List Experiments as Explanatory Variables in Regression Models.'' available 
at \url{http://imai.princeton.edu/research/files/listExp.pdf}
}


\examples{

data(mexico)
loyal <- mexico[mexico$mex.loyal == 1,]
notloyal <- mexico[mexico$mex.loyal == 0,]

\dontrun{

## Logistic outcome regression
## (effect of vote-selling on turnout)
## This replicates Table 4 in Imai et al. 2014

loyalreg <- ictreg.joint(formula = mex.y.all ~ mex.male + mex.age + mex.age2 + mex.education +  
                         mex.interest + mex.married +
                         mex.wealth + mex.urban + mex.havepropoganda + mex.concurrent, data = loyal,
                         treat = "mex.t", outcome = "mex.votecard", J = 3, constrained = TRUE,
                         outcome.reg = "logistic", maxIter = 1000)


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

## Generate predicted probability of turnout, averaged over the whole sample,
## for vote sellers (z = 1), non-vote sellers (z = 0), and the difference
## between vote sellers and non-vote sellers, in the sample of party supporters.
## This replicates the results in the righthand panel of Figure 2 in Imai et al. 2014

loyalpred <- predict.ictreg.joint(loyalreg, se.fit = TRUE, interval = "confidence", 
					level = 0.95, avg = TRUE, 
					sensitive.value = "both", 
					sensitive.diff = TRUE, return.draws = TRUE,
    					predict.sensitive = TRUE)

loyalpred$fit

## View predicted probability of vote selling, in the sample of party supporters.
## This replicates the results in the lefthand panel of Figure 2 in Imai et al. 2014

loyalpred$fitsens


}
}
