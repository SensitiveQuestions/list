\name{plot.predict.ictreg} 
\alias{plot.predict.ictreg} 
\title{Plot Method for the Item Count Technique} 
\description{ 
Function to plot predictions and confidence intervals of predictions from estimates from multivariate regression analysis of survey data with the item count technique.
} 
\usage{

  \method{plot}{predict.ictreg}(x, labels = NA, axes.ict = TRUE,
                                xlim = NULL, ylim = NULL, xlab = NULL, 
				ylab = "Estimated Proportion",
                                axes = F, pch = 19, xvec = NULL, ...) 

} 

\arguments{ 
\item{x}{object or set of objects of class inheriting from "predict.ictreg". Either a single object from an \code{ictreg()} model fit or multiple \code{predict} objects combined with the c() function. }
\item{labels}{a vector of labels for each prediction, plotted at the x axis.}
\item{axes.ict}{a switch indicating if custom plot axes are to be used with the user-provided estimate \code{labels}.}
\item{xlim}{a title for the y axis.}
\item{ylim}{a title for the y axis.}
\item{xlab}{a title for the x axis.}
\item{ylab}{a title for the y axis.}
\item{axes}{an indicator for whether default plot axes are included.}
\item{pch}{either an integer specifying a symbol or a single character to be used as the default in plotting points.}
\item{xvec}{a vector of x values at which the proportions will be printed.}
\item{...}{Other graphical parameters to be passed to the \code{plot()} command are accepted.}
}

\details{\code{plot.predict.ictreg} produces plots with estimated population proportions of respondents answering the sensitive item in a list experiment in the affirmative, with confidence intervals. 

The function accepts a set of \code{predict.ictreg} objects calculated in the following manner:

\code{predict(ictreg.object, avg = TRUE, interval = "confidence")}

For each average prediction, a point estimate and its confidence interval is plotted at equally spaced intervals. The x location of the points can be manipulated with the \code{xvec} option.

Either a single predict object can be plotted, or a group of them combined with \code{c(predict.object1, predict.object2)}. Predict objects with the \code{newdata.diff} option, which calculates the mean difference in probability between two datasets, and the \code{direct.glm} option, which calculates the mean difference between the mean predicted support for the sensitive item in the list experiment and in a direct survey item, can also be plotted in the same way as other \code{predict} objects.
}

\seealso{\code{\link{ictreg}} for model fitting and \code{\link{predict.ictreg}} for predictions based on the model fits.}


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
race.south <- race.nonsouth <- race
race.south[, "south"] <- 1
race.nonsouth[, "south"] <- 0

\dontrun{

# Fit EM algorithm ML model with constraint
ml.constrained.results <- ictreg(y ~ south + age + male + college, 
   data = race, treat = "treat", J=3, method = "ml", 
   overdispersed = FALSE, constrained = TRUE)

# Calculate average predictions for respondents in the South 
# and the the North of the US for the MLE model, replicating the 
# estimates presented in Figure 1, Imai (2011)
avg.pred.south.mle <- predict(ml.constrained.results, 
   newdata = race.south, avg = TRUE, interval = "confidence")
avg.pred.nonsouth.mle <- predict(ml.constrained.results, 
   newdata = race.nonsouth, avg = TRUE, interval = "confidence")

# A plot of a single estimate and its confidence interval
plot(avg.pred.south.mle, labels = "South")

# A  plot of the two estimates and their confidence intervals
# use c() to combine more than one predict object for plotting
plot(c(avg.pred.south.mle, avg.pred.nonsouth.mle), labels = c("South", "Non-South"))

# The difference option can also be used to simultaneously
# calculate separate estimates of the two sub-groups
# and the estimated difference. This can also be plotted.

avg.pred.diff.mle <- predict(ml.constrained.results, 
   newdata = race.south, newdata.diff = race.nonsouth,
   se.fit = TRUE, avg = TRUE)

plot(avg.pred.diff.mle)

# Social desirability bias plots

# Estimate logit for direct sensitive question

data(mis)

mis.list <- subset(mis, list.data == 1)

mis.sens <- subset(mis, sens.data == 1)

# Fit EM algorithm ML model

fit.list <- ictreg(y ~ age + college + male + south,
   J = 4, data = mis.list, method = "ml")

# Fit logistic regression with directly-asked sensitive question

fit.sens <- glm(sensitive ~ age + college + male + south, 
   data = mis.sens, family = binomial("logit"))

# Predict difference between response to sensitive item
# under the direct and indirect questions (the list experiment).
# This is an estimate of the revealed social desirability bias
# of respondents. See Blair and Imai (2010).

avg.pred.social.desirability <- predict(fit.list, 
   direct.glm = fit.sens, se.fit = TRUE)

plot(avg.pred.social.desirability)

}
}
 
\keyword{models}
\keyword{regression}