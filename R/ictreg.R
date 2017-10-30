#' Item Count Technique
#' 
#' Function to conduct multivariate regression analyses of survey data with the
#' item count technique, also known as the list experiment and the unmatched
#' count technique.
#' 
#' This function allows the user to perform regression analysis on data from
#' the item count technique, also known as the list experiment and the
#' unmatched count technique.
#' 
#' Three list experiment designs are accepted by this function: the standard
#' design; the multiple sensitive item standard design; and the modified design
#' proposed by Corstange (2009).
#' 
#' For the standard design, three methods are implemented in this function: the
#' linear model; the Maximum Likelihood (ML) estimation for the
#' Expectation-Maximization (EM) algorithm; the nonlinear least squares (NLS)
#' estimation with the two-step procedure both proposed in Imai (2010); and the
#' Maximum Likelihood (ML) estimator in the presence of two types of dishonest
#' responses, "ceiling" and "floor" liars. The ceiling model, floor model, or
#' both, as described in Blair and Imai (2010) can be activated by using the
#' \code{ceiling} and \code{floor} options. The constrained and unconstrained
#' ML models presented in Imai (2010) are available through the
#' \code{constrained} option, and the user can specify if overdispersion is
#' present in the data for the no liars models using the \code{overdispersed}
#' option to control whether a beta-binomial or binomial model is used in the
#' EM algorithm to model the item counts.
#' 
#' The modified design and the multiple sensitive item design are automatically
#' detected by the function, and only the binomial model without overdispersion
#' is available.
#' 
#' @aliases ictreg ict list
#' @param formula An object of class "formula": a symbolic description of the
#' model to be fitted.
#' @param data A data frame containing the variables in the model
#' @param treat Name of treatment indicator as a string. For single sensitive
#' item models, this refers to a binary indicator, and for multiple sensitive
#' item models it refers to a multi-valued variable with zero representing the
#' control condition. This can be an integer (with 0 for the control group) or
#' a factor (with "control" for the control group).
#' @param J Number of non-sensitive (control) survey items.
#' @param method Method for regression, either \code{ml} for the Maximum
#' Likelihood (ML) estimation with the Expectation-Maximization algorithm;
#' \code{lm} for linear model estimation; or \code{nls} for the Non-linear
#' Least Squares (NLS) estimation with the two-step procedure.
#' @param weights Name of the weights variable as a string, if weighted
#' regression is desired. Not implemented for the ceiling/floor models,
#' multiple sensitive item design, or for the modified design.
#' @param h Auxiliary data functionality. Optional named numeric vector with
#' length equal to number of groups. Names correspond to group labels and
#' values correspond to auxiliary moments.
#' @param group Auxiliary data functionality. Optional character vector of
#' group labels with length equal to number of observations.
#' @param matrixMethod Auxiliary data functionality. Procedure for estimating
#' optimal weighting matrix for generalized method of moments. One of
#' "efficient" for two-step feasible and "cue" for continuously updating.
#' Default is "efficient". Only relevant if \code{h} and \code{group} are
#' specified.
#' @param overdispersed Indicator for the presence of overdispersion. If
#' \code{TRUE}, the beta-binomial model is used in the EM algorithm, if
#' \code{FALSE} the binomial model is used. Not relevant for the \code{NLS} or
#' \code{lm} methods.
#' @param constrained A logical value indicating whether the control group
#' parameters are constrained to be equal.  Not relevant for the \code{NLS} or
#' \code{lm} methods
#' @param floor A logical value indicating whether the floor liar model should
#' be used to adjust for the possible presence of respondents dishonestly
#' reporting a negative preference for the sensitive item among those who hold
#' negative views of all the non-sensitive items.
#' @param ceiling A logical value indicating whether the ceiling liar model
#' should be used to adjust for the possible presence of respondents
#' dishonestly reporting a negative preference for the sensitive item among
#' those who hold affirmative views of all the non-sensitive items.
#' @param ceiling.fit Fit method for the M step in the EM algorithm used to fit
#' the ceiling liar model. \code{glm} uses standard logistic regression, while
#' \code{bayesglm} uses logistic regression with a weakly informative prior
#' over the parameters.
#' @param floor.fit Fit method for the M step in the EM algorithm used to fit
#' the floor liar model. \code{glm} uses standard logistic regression, while
#' \code{bayesglm} uses logistic regression with a weakly informative prior
#' over the parameters.
#' @param ceiling.formula Covariates to include in ceiling liar model. These
#' must be a subset of the covariates used in \code{formula}.
#' @param floor.formula Covariates to include in floor liar model. These must
#' be a subset of the covariates used in \code{formula}.
#' @param fit.start Fit method for starting values for standard design
#' \code{ml} model. The options are \code{lm}, \code{glm}, and \code{nls},
#' which use OLS, logistic regression, and non-linear least squares to generate
#' starting values, respectively. The default is \code{nls}.
#' @param fit.nonsensitive Fit method for the non-sensitive item fit for the
#' \code{nls} method and the starting values for the \code{ml} method for the
#' \code{modified} design. Options are \code{glm} and \code{nls}, and the
#' default is \code{nls}.
#' @param multi.condition For the multiple sensitive item design, covariates
#' representing the estimated count of affirmative responses for each
#' respondent can be included directly as a level variable by choosing
#' \code{level}, or as indicator variables for each value but one by choosing
#' \code{indicators}. The default is \code{none}.
#' @param maxIter Maximum number of iterations for the Expectation-Maximization
#' algorithm of the ML estimation.  The default is 5000.
#' @param verbose a logical value indicating whether model diagnostics are
#' printed out during fitting.
#' @param ... further arguments to be passed to NLS regression commands.
#' @return \code{ictreg} returns an object of class "ictreg".  The function
#' \code{summary} is used to obtain a table of the results.  The object
#' \code{ictreg} is a list that contains the following components.  Some of
#' these elements are not available depending on which method is used
#' (\code{lm}, \code{nls} or \code{ml}), which design is used (\code{standard},
#' \code{modified}), whether multiple sensitive items are include
#' (\code{multi}), and whether the constrained model is used (\code{constrained
#' = TRUE}).
#' 
#' \item{par.treat}{point estimate for effect of covariate on item count fitted
#' on treatment group} \item{se.treat}{standard error for estimate of effect of
#' covariate on item count fitted on treatment group} \item{par.control}{point
#' estimate for effect of covariate on item count fitted on control group}
#' \item{se.control}{standard error for estimate of effect of covariate on item
#' count fitted on control group} \item{coef.names}{variable names as defined
#' in the data frame} \item{design}{call indicating whether the \code{standard}
#' design as proposed in Imai (2010) or thee \code{modified} design as proposed
#' in Corstange (2009) is used} \item{method}{call of the method used}
#' \item{overdispersed}{call indicating whether data is overdispersed}
#' \item{constrained}{call indicating whether the constrained model is used}
#' \item{boundary}{call indicating whether the floor/ceiling boundary models
#' are used} \item{multi}{indicator for whether multiple sensitive items were
#' included in the data frame} \item{call}{the matched call} \item{data}{the
#' \code{data} argument} \item{x}{the design matrix} \item{y}{the response
#' vector} \item{treat}{the vector indicating treatment status} \item{J}{Number
#' of non-sensitive (control) survey items set by the user or detected.}
#' \item{treat.labels}{a vector of the names used by the \code{treat} vector
#' for the sensitive item or items. This is the names from the \code{treat}
#' indicator if it is a factor, or the number of the item if it is numeric.}
#' \item{control.label}{a vector of the names used by the \code{treat} vector
#' for the control items. This is the names from the \code{treat} indicator if
#' it is a factor, or the number of the item if it is numeric.}
#' 
#' For the maximum likelihood models, an additional output object is included:
#' \item{pred.post}{posterior predicted probability of answering "yes" to the
#' sensitive item. The weights from the E-M algorithm.}
#' 
#' For the floor/ceiling models, several additional output objects are
#' included: \item{ceiling}{call indicating whether the assumption of no
#' ceiling liars is relaxed, and ceiling parameters are estimated}
#' \item{par.ceiling}{point estimate for effect of covariate on whether
#' respondents who answered affirmatively to all non-sensitive items and hold a
#' true affirmative opinion toward the sensitive item lied and reported a
#' negative response to the sensitive item } \item{se.ceiling}{standard error
#' for estimate for effect of covariate on whether respondents who answered
#' affirmatively to all non-sensitive items and hold a true affirmative opinion
#' toward the sensitive item lied and reported a negative response to the
#' sensitive item} \item{floor}{call indicating whether the assumption of no
#' floor liars is relaxed, and floor parameters are estimated}
#' \item{par.ceiling}{point estimate for effect of covariate on whether
#' respondents who answered negatively to all non-sensitive items and hold a
#' true affirmative opinion toward the sensitive item lied and reported a
#' negative response to the sensitive item } \item{se.ceiling}{standard error
#' for estimate for effect of covariate on whether respondents who answered
#' negatively to all non-sensitive items and hold a true affirmative opinion
#' toward the sensitive item lied and reported a negative response to the
#' sensitive item} \item{coef.names.ceiling}{variable names from the ceiling
#' liar model fit, if applicable} \item{coef.names.floor}{variable names from
#' the floor liar model fit, if applicable}
#' 
#' For the multiple sensitive item design, the \code{par.treat} and
#' \code{se.treat} vectors are returned as lists of vectors, one for each
#' sensitive item.
#' 
#' For the unconstrained model, the \code{par.control} and \code{se.control}
#' output is replaced by: \item{par.control.phi0}{point estimate for effect of
#' covariate on item count fitted on treatment group}
#' \item{se.control.phi0}{standard error for estimate of effect of covariate on
#' item count fitted on treatment group} \item{par.control.phi1}{point estimate
#' for effect of covariate on item count fitted on treatment group}
#' \item{se.control.phi1}{standard error for estimate of effect of covariate on
#' item count fitted on treatment group}
#' 
#' Depending upon the estimator requested by the user, model fit statistics are
#' also included: \item{llik}{the log likelihood of the model, if \code{ml} is
#' used} \item{resid.se}{the residual standard error, if \code{nls} or
#' \code{lm} are used. This will be a scalar if the standard design was used,
#' and a vector if the multiple sensitive item design was used}
#' \item{resid.df}{the residual degrees of freedom, if \code{nls} or \code{lm}
#' are used. This will be a scalar if the standard design was used, and a
#' vector if the multiple sensitive item design was used}
#' 
#' When using the auxiliary data functionality, the following objects are
#' included: \item{aux}{logical value indicating whether estimation
#' incorporates auxiliary moments} \item{nh}{integer count of the number of
#' auxiliary moments} \item{wm}{procedure used to estimate the optimal weight
#' matrix} \item{J.stat}{numeric value of the Sargan Hansen overidentifying
#' restriction test statistic} \item{overid.p}{corresponding p-value for the
#' Sargan Hansen test}
#' @author Graeme Blair, UCLA, \email{graeme.blair@ucla.edu}
#' and Kosuke Imai, Princeton University, \email{kimai@princeton.edu}
#' @seealso \code{\link{predict.ictreg}} for fitted values
#' @references Blair, Graeme and Kosuke Imai. (2012) ``Statistical Analysis of
#' List Experiments."  Political Analysis. Forthcoming. available at
#' \url{http://imai.princeton.edu/research/listP.html}
#' 
#' Imai, Kosuke. (2011) ``Multivariate Regression Analysis for the Item Count
#' Technique.'' Journal of the American Statistical Association, Vol. 106, No.
#' 494 (June), pp. 407-416. available at
#' \url{http://imai.princeton.edu/research/list.html}
#' @keywords models regression
#' @examples
#' 
#' 
#' data(race)
#' 
#' set.seed(1)
#' 
#' # Calculate list experiment difference in means
#' 
#' diff.in.means.results <- ictreg(y ~ 1, data = race, 
#' 	       	      treat = "treat", J=3, method = "lm")
#' 
#' summary(diff.in.means.results)
#' 
#' # Fit linear regression
#' # Replicates Table 1 Columns 1-2 Imai (2011); note that age is divided by 10
#' 
#' lm.results <- ictreg(y ~ south + age + male + college, data = race, 
#' 	       	      treat = "treat", J=3, method = "lm")
#' 
#' summary(lm.results)
#' 
#' # Fit two-step non-linear least squares regression
#' # Replicates Table 1 Columns 3-4 Imai (2011); note that age is divided by 10
#' 
#' nls.results <- ictreg(y ~ south + age + male + college, data = race, 
#' 	       	      treat = "treat", J=3, method = "nls")
#' 
#' summary(nls.results)
#' 
#' \dontrun{
#' 
#' # Fit EM algorithm ML model with constraint
#' # Replicates Table 1 Columns 5-6, Imai (2011); note that age is divided by 10
#' 
#' ml.constrained.results <- ictreg(y ~ south + age + male + college, data = race, 
#' 		       	  	 treat = "treat", J=3, method = "ml", 
#' 				 overdispersed = FALSE, constrained = TRUE)
#' 
#' summary(ml.constrained.results)
#' 
#' # Fit EM algorithm ML model with no constraint
#' # Replicates Table 1 Columns 7-10, Imai (2011); note that age is divided by 10
#' 
#' ml.unconstrained.results <- ictreg(y ~ south + age + male + college, data = race, 
#' 			    	   treat = "treat", J=3, method = "ml", 
#' 				   overdispersed = FALSE, constrained = FALSE)
#' 
#' summary(ml.unconstrained.results)
#' 
#' # Fit EM algorithm ML model for multiple sensitive items
#' # Replicates Table 3 in Blair and Imai (2010)
#' 
#' multi.results <- ictreg(y ~ male + college + age + south + south:age, treat = "treat", 
#' 	      	 	J = 3, data = multi, method = "ml", 
#' 			multi.condition = "level")
#' 
#' summary(multi.results)
#' 
#' # Fit standard design ML model
#' # Replicates Table 7 Columns 1-2 in Blair and Imai (2010)
#' 
#' noboundary.results <- ictreg(y ~ age + college + male + south, treat = "treat",
#' 		      	     J = 3, data = affirm, method = "ml", 
#' 			     overdispersed = FALSE)
#' 
#' summary(noboundary.results)
#' 
#' # Fit standard design ML model with ceiling effects alone
#' # Replicates Table 7 Columns 3-4 in Blair and Imai (2010)
#' 
#' ceiling.results <- ictreg(y ~ age + college + male + south, treat = "treat", 
#' 		   	  J = 3, data = affirm, method = "ml", fit.start = "nls",
#' 			  ceiling = TRUE, ceiling.fit = "bayesglm",
#' 			  ceiling.formula = ~ age + college + male + south)
#' 
#' summary(ceiling.results)
#' 
#' # Fit standard design ML model with floor effects alone
#' # Replicates Table 7 Columns 5-6 in Blair and Imai (2010)
#' 
#' floor.results <- ictreg(y ~ age + college + male + south, treat = "treat", 
#' 	      	 	J = 3, data = affirm, method = "ml", fit.start = "glm", 
#' 			floor = TRUE, floor.fit = "bayesglm",
#' 			floor.formula = ~ age + college + male + south)
#' 
#' summary(floor.results)
#' 
#' # Fit standard design ML model with floor and ceiling effects
#' # Replicates Table 7 Columns 7-8 in Blair and Imai (2010)
#' 
#' both.results <- ictreg(y ~ age + college + male + south, treat = "treat", 
#' 	     	       J = 3, data = affirm, method = "ml", 
#' 		       floor = TRUE, ceiling = TRUE, 
#' 		       floor.fit = "bayesglm", ceiling.fit = "bayesglm",
#' 		       floor.formula = ~ age + college + male + south,
#' 		       ceiling.formula = ~ age + college + male + south)
#' 
#' summary(both.results)
#' 
#' }
#' 
ictreg <- function(formula, data = parent.frame(), treat = "treat", J, method = "ml", weights, 
                   h = NULL, group = NULL, matrixMethod = "efficient", robust = FALSE, error = "none", 
                   overdispersed = FALSE, constrained = TRUE, floor = FALSE, ceiling = FALSE, 
                   ceiling.fit = "glm", floor.fit = "glm", ceiling.formula = ~ 1, floor.formula = ~ 1, 
                   fit.start = "lm", fit.nonsensitive = "nls", multi.condition = "none", maxIter = 5000, verbose = FALSE, ...){

  ictreg.call <- match.call()

  # set up data frame, with support for standard and modified responses
  mf <- match.call(expand.dots = FALSE)
  
  # make all other call elements null in mf <- NULL in next line
  mf$method <- mf$maxIter <- mf$verbose <- mf$fit.start <- mf$J <- mf$design <- 
    mf$treat <- mf$weights <- mf$constrained <- mf$overdispersed <- mf$floor <- 
    mf$ceiling <- mf$ceiling.fit <- mf$fit.nonsensitive <- mf$floor.fit <- 
    mf$multi.condition <- mf$floor.formula <- mf$ceiling.formula <- mf$h <-
    mf$group <- mf$matrixMethod <- mf$robust <- mf$error <- NULL
  mf[[1]] <- as.name("model.frame")
  mf$na.action <- 'na.pass'
  mf <- eval.parent(mf)
  
  if(floor == TRUE | ceiling == TRUE) {
    boundary <- TRUE
  } else {
    boundary <- FALSE
  }
  
  # define design, response data frames
  x.all <- model.matrix.default(attr(mf, "terms"), mf)
  y.all <- model.response(mf)
  
  if(class(y.all)=="matrix") {
    design <- "modified"
  } else {
    design <- "standard"
  }
  
  if(missing("weights") == FALSE) {
      weighted <- TRUE
      w.all <- data[, paste(weights)]
  } else {
      weighted <- FALSE
      w.all <- rep(1, length(y.all))
  }

  if (method == "nls") fit.start <- "nls"
  
  # extract number of non-sensitive items from the response matrix
  if(design=="modified") J <- ncol(y.all) - 1
  
  # set -777 code for treatment status
#  if(design=="modified"){
#    for(j in 1:J) y.all[,j] <- ifelse(!is.na(y.all[,J+1]), -777, y.all[,j])
#    y.all[,J+1] <- ifelse(is.na(y.all[,J+1]), -777, y.all[,J+1])
#  }
  
  # list-wise missing deletion
  na.x <- apply(is.na(x.all), 1, sum)
  if(design=="standard") na.y <- is.na(y.all)
  if(design=="modified") {
    na.y <- apply(is.na(y.all[,1:J]), 1, sum)
    na.y[is.na(y.all[,J+1]) == FALSE] <- 0
  }
  na.w <- is.na(w.all)

  ## treatment indicator for subsetting the dataframe
  if(design=="standard") t <- data[na.x==0 & na.y==0 & na.w==0, paste(treat)]
  if(design=="modified") t <- as.numeric(is.na(y.all[na.x == 0 & na.y == 0 & na.w==0, J+1]) == FALSE)

  if(class(t) == "factor") {

    levels(t) <- tolower(levels(t))
    
    if (length(which(levels(t) == "control")) == 1) {
      t <- relevel(t, ref = "control")
    } else {
      warning("Note: using the first level as the control condition, but it is not labeled 'control'.")
    }
        
    condition.labels <- levels(t)
    t <- as.numeric(t) - 1
    treatment.labels <- condition.labels[2:length(condition.labels)]
    control.label <- condition.labels[1]
    
  } else {
    ##condition.labels <- as.character(paste("Sensitive Item",sort(unique(t))))
    ##condition.labels[which(condition.labels == "Sensitive Item 0")] <- "control"
    ##treatment.labels <- condition.labels[condition.labels != "control"]

    condition.labels <- sort(unique(t))
    treatment.labels <- condition.labels[condition.labels != 0]
    control.label <- 0
  }
  
  # list wise delete
  if(design=="standard") y.all <- y.all[na.x==0 & na.y==0 & na.w==0]
  if(design=="modified") y.all <- y.all[na.x==0 & na.y==0 & na.w==0,]
  x.all <- x.all[na.x==0 & na.y==0, , drop = FALSE]
  w.all <- w.all[na.x==0 & na.y==0 & na.w==0]

  ## so that the output data has the same dimension as x.all and y.all
  data <- data[na.x==0 & na.y==0 & na.w==0, , drop = FALSE]
  
  # set up data objects for y and x for each group from user input
  x.treatment <- x.all[t != 0, , drop = FALSE]
  y.treatment <- subset(y.all, t != 0)
  w.treatment <- subset(w.all, t != 0)
  x.control <- x.all[t == 0 , , drop = FALSE]
  y.control <- subset(y.all, t==0)
  w.control <- subset(w.all, t==0)

  if(design=="standard" & missing("J")) {
    J <- max(y.treatment) - 1
    cat(paste("Note: number of control items set to", J, "\n"))
  }

  condition.values <- sort(unique(t))
  treatment.values <- 1:length(condition.values[condition.values!=0])
  
  if(length(treatment.values) > 1) {
    multi <- TRUE
  } else {
    multi <- FALSE
  }

  if(weighted == TRUE) {
    if(design == "modified")
      stop("Weighted list experiment regression is not supported (yet) for the modified design.")
    if(boundary == TRUE)
      stop("Weighted list experiment regression is not supported (yet) for the ceiling/floor design.")
    if(multi == TRUE)
      stop("Weighted list experiment regression is not supported (yet) for the multi-item design.")
  }

  # robust ml functionality -- check conditions
  if (robust) {
    if (multi.condition != "none") 
      stop("The robust ML functionality is not yet supported for multiple sensitive item designs.")
    if (method != "ml")
      stop("You must specify method as 'ml' to use the robust ML functionality.")
  }

  # auxiliary data functionality -- check conditions
  aux.check <- !(is.null(h) & is.null(group))

  if (aux.check) {
    if (multi.condition != "none") 
      stop("The auxiliary data functionality is not yet supported for multiple sensitive item designs.")
    if (method != "nls")
      stop("The auxiliary data functionality is currently supported only for the nonlinear least squares method.")
    if (robust)
      stop("The auxiliary data functionality is not yet compatible with the robust ML functionality.")
  }

  n <- nrow(x.treatment) + nrow(x.control)

  coef.names <- colnames(x.all)

  nPar <- ncol(x.all)

  intercept.only <- ncol(x.all)==1 & sum(x.all[,1]==1) == n

  if (intercept.only == TRUE) {
    x.vars <- "1"
  } else {
    x.vars <- coef.names[-1]
  }
  
  logistic <- function(x) exp(x)/(1+exp(x))
  logit <- function(x) return(log(x)-log(1-x))

  if (design == "standard" & method == "lm") {
    
    treat <- matrix(NA, nrow = n, ncol = length(treatment.values))
    
    for (m in 1:length(treatment.values)) 
      treat[ , m] <- as.numeric(t == m)
      
    x.all.noint <- x.all[, -1, drop = FALSE]

    x.all.lm <- x.all
    for (m in 1:length(treatment.values))
      x.all.lm <- cbind(x.all.lm, x.all * treat[, m])

    if (intercept.only == TRUE) {
      fit.lm <- lm(y.all ~ treat, weights = w.all)
    } else {
      ## fit.lm <- lm(y.all ~ x.all.noint * treat)
      fit.lm <- lm(y.all ~ x.all.lm - 1, weights = w.all)
    }
    
    vcov <- vcovHC(fit.lm)

    par.control <- coef(fit.lm)[1:length(coef.names)]
    se.control <- sqrt(diag(vcov))[1:length(coef.names)]

    names(par.control) <- names(se.control) <- coef.names
    
    par.treat <- se.treat <- list()
    
    for (m in 1:length(treatment.values)) {
      if (intercept.only == TRUE) {
        par.treat[[m]] <- coef(fit.lm)[nPar + m]
        se.treat[[m]] <- sqrt(diag(vcov))[nPar + m]
        names(par.treat[[m]]) <- names(se.treat[[m]]) <- "(Intercept)"
      } else {
        
        par.treat[[m]] <- coef(fit.lm)[(nPar + (m-1) * nPar + 1) : (nPar + m * nPar)]
        se.treat[[m]] <- sqrt(diag(vcov))[(nPar + (m-1) * nPar + 1) : (nPar + m * nPar)]
        
        names(par.treat[[m]]) <- names(se.treat[[m]]) <- coef.names
      }
      
    }

    if(multi == FALSE) {
      par.treat <- par.treat[[1]]
      se.treat <- se.treat[[1]]
    }
    
    sum.fit.lm <- summary(fit.lm)
      
    resid.se <- sum.fit.lm$sigma
    resid.df <- sum.fit.lm$df[2]
    
    if (multi == TRUE) {
      return.object <- list(par.treat=par.treat, se.treat=se.treat,
                            par.control=par.control, se.control=se.control,
                            vcov=vcov, treat.values = treatment.values, treat.labels = treatment.labels,
                            J=J, coef.names=coef.names,
                            resid.se = resid.se, resid.df = resid.df,  design = design, method = method,
                            multi = multi, boundary = boundary, call = match.call(),
                            data = data, x = x.all, y = y.all, treat = t)
      if(weighted == TRUE)
          return.object$weights <- w.all
    } else {
      return.object <- list(par.treat=par.treat, se.treat=se.treat,
                            par.control=par.control, se.control=se.control,
                            vcov=vcov, J=J,  coef.names=coef.names, resid.se = resid.se,
                            resid.df = resid.df,  design = design, method = method,
                            multi = multi, boundary = boundary, call = match.call(),
                            data = data, x = x.all, y = y.all, treat = t)
      if(weighted == TRUE)
          return.object$weights <- w.all
    }
    
  } else if (design=="standard" & method != "lm") {
    
    if (error == "none") {
      
      ##
      # Run two-step estimator
      
      vcov.twostep.std <- function(treat.fit, control.fit, J, y, treat, x) {
        
        n <- length(y)
        y1 <- y[treat == 1]
        y0 <- y[treat == 0]
        x1 <- x[treat == 1, , drop = FALSE]
        x0 <- x[treat == 0, , drop = FALSE]
        delta <- coef(treat.fit)
        gamma <- coef(control.fit)
        
        m1 <- c((y1 - J*logistic(x1 %*% gamma) - logistic(x1 %*% delta)) *
                  logistic(x1 %*% delta)/(1+exp(x1 %*% delta))) * x1
        m0 <- c((y0 - J*logistic(x0 %*% gamma)) * J *
                  logistic(x0 %*% gamma)/(1+exp(x0 %*% gamma))) * x0
        Em1 <- t(m1) %*% m1 / n
        Em0 <- t(m0) %*% m0 / n
        F <- adiag(Em1, Em0)
        Gtmp <- c(logistic(x1 %*% delta)/(1 + exp(x1 %*% delta))) * x1
        G1 <- -t(Gtmp) %*% Gtmp / n  
        Gtmp <- c(sqrt(J*logistic(x1 %*% delta)*logistic(x1 %*% gamma)/
                         ((1 + exp(x1 %*% delta))*(1 + exp(x1 %*% gamma))))) * x1
        G2 <- -t(Gtmp) %*% Gtmp / n
        Gtmp <- c(J*logistic(x0 %*% gamma)/(1 + exp(x0 %*% gamma))) * x0
        G3 <- -t(Gtmp) %*% Gtmp / n
        invG1 <- solve(G1, tol = 1e-20)
        invG3 <- solve(G3, tol = 1e-20)
        invG <- rbind(cbind(invG1, - invG1 %*% G2 %*% invG3),
                      cbind(matrix(0, ncol = ncol(G1), nrow = nrow(G3)), invG3))
        
        return(invG %*% F %*% t(invG) / n)
        
      }
      
      ## fit to the control group
      fit.glm.control <- glm(cbind(y.control, J-y.control) ~ x.control - 1,
                             family = binomial(logit), weights = w.control)
      
      coef.glm.control <- coef(fit.glm.control)
      names(coef.glm.control) <- paste("beta", 1:length(coef.glm.control), sep = "")
      
      if(fit.start == "nls") {
        if(is.null(ictreg.call$control)){
          fit.control <- nls( as.formula(paste("I(y.control/J) ~ logistic(x.control %*% c(",
                                               paste(paste("beta", 1:length(coef.glm.control),
                                                           sep=""), collapse= ","), "))")) ,
                              start = coef.glm.control, weights = w.control, control =
                                nls.control(maxiter=maxIter, warnOnly=TRUE), ... = ...)
        } else {
          fit.control <- nls( as.formula(paste("I(y.control/J) ~ logistic(x.control %*% c(",
                                               paste(paste("beta", 1:length(coef.glm.control),
                                                           sep=""), collapse= ","), "))")) ,
                              weights = w.control,
                              start = coef.glm.control, ... = ...)
        }
        
        fit.control.coef <- summary(fit.control)$parameters[,1]
        
      } else if (fit.start == "glm") {
        fit.control <- fit.glm.control
        fit.control.coef <- coef(fit.glm.control)
      } else if (fit.start == "lm") {
        fit.control <- lm(y.control ~ x.control - 1, weights = w.control)
        fit.control.coef <- coef(fit.glm.control)
      }
      
      fit.treat <- vcov.nls <- se.twostep <- par.treat.nls.std <- list()
      
      for(m in 1:length(treatment.values)){
        
        curr.treat <- t[t!=0] == treatment.values[m]
        
        x.treatment.curr <- x.treatment[curr.treat, , drop = F]
        w.treatment.curr <- w.treatment[curr.treat]
        
        ## calculate the adjusted outcome    
        y.treatment.pred <- y.treatment[curr.treat] - logistic(x.treatment.curr
                                                               %*% fit.control.coef)*J
        
        ## fit to the treated		
        y.treatment.pred.temp <- ifelse(y.treatment.pred>1, 1, y.treatment.pred)
        y.treatment.pred.temp <- ifelse(y.treatment.pred.temp<0, 0, y.treatment.pred.temp)
        
        alpha <- mean(y.treatment.pred.temp)
        
        y.treatment.start <- ifelse(y.treatment.pred.temp >
                                      quantile(y.treatment.pred.temp, alpha), 1, 0)
        
        try(fit.glm.treat <- glm(cbind(y.treatment.start, 1 - y.treatment.start) ~ x.treatment.curr - 1,
                                 family = binomial(logit), weights = w.treatment.curr), silent = F)
        try(coef.glm.treat <- coef(fit.glm.treat), silent = T)
        try(names(coef.glm.treat) <- paste("beta", 1:length(coef.glm.treat), sep = ""), silent = T)
        
        if(fit.start == "nls") {
          if(exists("coef.glm.treat")) {
            if (is.null(ictreg.call$control)) {
              fit.treat[[m]] <- nls( as.formula(paste("y.treatment.pred ~ logistic(x.treatment.curr %*% c(",
                                                      paste(paste("beta", 1:length(coef.glm.treat),
                                                                  sep=""), collapse= ","), "))")) ,
                                     start = coef.glm.treat, weights = w.treatment.curr, control =
                                       nls.control(maxiter = maxIter, warnOnly = TRUE), ... = ...)
            } else {
              fit.treat[[m]] <- nls( as.formula(paste("y.treatment.pred ~ logistic(x.treatment.curr %*% c(",
                                                      paste(paste("beta", 1:length(coef.glm.treat),
                                                                  sep=""), collapse= ","), "))")) ,
                                     start = coef.glm.treat, weights = w.treatment.curr, ... = ...)
            }
          } else {
            if (is.null(ictreg.call$control)) {
              fit.treat[[m]] <- nls( as.formula(paste("y.treatment.pred ~ logistic(x.treatment.curr %*% c(",
                                                      paste(paste("beta", 1:length(coef.glm.treat),
                                                                  sep=""), collapse= ","), "))")),
                                     weights = w.treatment.curr,
                                     control = nls.control(maxiter = maxIter, warnOnly = TRUE),
                                     ... = ...)
            } else {
              fit.treat[[m]] <- nls( as.formula(paste("y.treatment.pred ~ logistic(x.treatment.curr %*% c(",
                                                      paste(paste("beta", 1:length(coef.glm.treat),
                                                                  sep=""), collapse= ","), "))")) ,
                                     weights = w.treatment.curr,
                                     ... = ...)
            }
          }
          
        } else if (fit.start == "glm") {
          fit.treat[[m]] <- fit.glm.treat
        } else if (fit.start == "lm") {
          fit.treat[[m]] <- lm(y.treatment.pred ~ x.treatment.curr - 1, weights = w.treatment.curr)
        }
        
        sample.curr <- t==0 | t==treatment.values[m]
        
        vcov.nls[[m]] <- vcov.twostep.std(fit.treat[[m]], fit.control, J,
                                          y.all[sample.curr], t[sample.curr]>0, x.all[sample.curr, , drop = FALSE])
        
        se.twostep[[m]] <- sqrt(diag(vcov.nls[[m]]))
        par.treat.nls.std[[m]] <- coef(fit.treat[[m]])
        
      }
      
      if(multi == FALSE) {
        fit.treat <- fit.treat[[1]]
        vcov.nls <- vcov.nls[[1]]
        se.twostep <- se.twostep[[1]]
        par.treat.nls.std <- par.treat.nls.std[[m]]
      }
      
      par.control.nls.std <- coef(fit.control)
      
    } else if (error == "topcode") {
      
      ## NLS top-coded error model
      
      #\l[Y_i - pJ - T_i\{ p + (1-p) \E(Z_i \mid X_i) \}- (1-p) \E(Y_i^\ast \mid X_i) \r]^2
      
      sse_nls_topcoded <- function(par, J, y, treat, x){
        pstar <- par[1]
        p <- exp(pstar) / (1 + exp(pstar))
        
        k <- ncol(x)
        coef.h <- par[2:(k + 1)]
        coef.g <- par[(k + 2):(2 * k + 1)]
        gX <- logistic(x %*% coef.g)
        hX <- logistic(x %*% coef.h)

        sse <- sum((y - (p * J + treat * (p + (1 - p) * gX) + (1 - p) * J * hX)) ^ 2)
        
        return(sse)
      }
      
      k <- ncol(x.all)
      
      start <- runif(k*2 + 1, min = -.5, max = .5)
      
      NLSfit <- optim(par = start, 
                      fn = sse_nls_topcoded, J = J, y = y.all,
                      treat = t, x = x.all, hessian = TRUE, control = list(maxit = maxIter))
      
      vcov.nls <- solve(-NLSfit$hessian, tol = 1e-20)
      se.nls <- sqrt(diag(vcov.nls))
      
      par.treat <- NLSfit$par[(k + 2):(2 * k + 1)]
      par.control <- NLSfit$par[2:(k + 1)]
      se.treat <- se.nls[(k + 2):(2 * k + 1)]
      se.control <- se.nls[2:(k + 1)]
      
      p.est <- logistic(NLSfit$par[1])
      p.ci <- c("lwr" = logistic(NLSfit$par[1] - qnorm(.975)*se.nls[1]),
                "upr" = logistic(NLSfit$par[1] + qnorm(.975)*se.nls[1]))
      
    } else if (error == "uniform") {
      
      test <- 5
      
      ## NLS uniform error model
      
      #\frac{p_0 (1-T_i) J}{2} + T_i \l\{\frac{p_1 (J+1)}{2} +
      #    (1-p_1)\E(Z_i \mid X_i)\r\}  \nonumber \\
      # + \{(1-T_i)(1 - p_0) + T_i(1-p_1)\} \E(Y_i^\ast \mid X_i)
      
      # For the parameter p, we use the logit transformation by defining a new parameter p* = log {p / (1-p)}.  Then, substitute p = exp(p*) / { 1 + exp(p*)} into equation 11.
      
      sse_nls_uniform <- function(par, J, y, treat, x){
        p0star <- par[1]
        p1star <- par[2]
        p0 <- exp(p0star) / (1 + exp(p0star))
        p1 <- exp(p1star) / (1 + exp(p1star))
        
        k <- ncol(x)
        coef.h <- par[3:(k + 2)]
        coef.g <- par[(k + 3):(2 * k + 2)]
        gX <- logistic(x %*% coef.g)
        hX <- logistic(x %*% coef.h)
        
        # sse <- sum((y - ((p0 * (1 - treat) * J) / 2 + 
        #                    treat * ((p1 * (J + 1)) / 2 + 
        #                               (1 - p1) * mean(gX)) +
        #               ((1 - treat) * (1 - p0) + 
        #                  treat * (1 - p1)) * J * mean(hX))) ^ 2)
        
        sse <- sum((y - ((p0 * (1 - treat) * J) / 2 + 
                    treat * ((p1 * (J + 1)) / 2 + 
                               (1 - p1) * gX) +
                    ((1 - treat) * (1 - p0) + 
                       treat * (1 - p1)) * J * hX)) ^ 2)
        
        return(sse)
      }
      
      k <- ncol(x.all)
      
      start <- c(0, 0, 0, 1, 0, 1) #runif(k*2 + 2)
      
      NLSfit <- optim(par = start, 
                      fn = sse_nls_uniform, J = J, y = y.all,
                      treat = t, x = x.all, hessian = TRUE, control = list(maxit = maxIter))
      
      vcov.nls <- solve(-NLSfit$hessian, tol = 1e-20)
      se.nls <- sqrt(diag(vcov.nls))
      
      par.treat <- NLSfit$par[(k + 3):(2 * k + 2)]
      par.control <- NLSfit$par[3:(k + 2)]
      se.treat <- se.nls[(k + 3):(2 * k + 2)]
      se.control <- se.nls[3:(k + 2)]
      
      p0.est <- logistic(NLSfit$par[1])
      p0.ci <- c("lwr" = logistic(NLSfit$par[1] - qnorm(.975)*se.nls[1]),
                "upr" = logistic(NLSfit$par[1] + qnorm(.975)*se.nls[1]))
      
      p1.est <- logistic(NLSfit$par[2])
      p1.ci <- c("lwr" = logistic(NLSfit$par[2] - qnorm(.975)*se.nls[2]),
                "upr" = logistic(NLSfit$par[2] + qnorm(.975)*se.nls[2]))
      
    }
      
    if(method=="ml" & robust == FALSE) {
      
      if(boundary == FALSE) {
      
        if (multi == FALSE) {

          coef.control.start <- par.control.nls.std
          coef.treat.start <- par.treat.nls.std
          
          ## ##
          ## Run EM estimator if requested
          
          ##
          ## observed-data log-likelihood for beta binomial
          ##
          obs.llik.std <- function(par, J, y, treat, x, wt, const = FALSE) {
            
            k <- ncol(x)
            if (const) {
              coef.h0 <- coef.h1 <- par[1:k]
              rho.h0 <- rho.h1 <- logistic(par[(k+1)])
              coef.g <- par[(k+2):(2*k+1)]
            } else {
              coef.h0 <- par[1:k]
              rho.h0 <- logistic(par[(k+1)])
              coef.h1 <- par[(k+2):(2*k+1)]
              rho.h1 <- logistic(par[(2*k+2)])
              coef.g <- par[(2*k+3):(3*k+2)]
            }
            
            h0X <- logistic(x %*% coef.h0)
            h1X <- logistic(x %*% coef.h1)
            gX <- logistic(x %*% coef.g)
            
            ind10 <- ((treat == 1) & (y == 0))
            ind1J1 <- ((treat == 1) & (y == (J+1)))
            ind1y <- ((treat == 1) & (y > 0) & (y < (J+1)))
            
            if (sum(ind10) > 0) {
              p10 <- sum(wt[ind10] * (log(1-gX[ind10]) + dBB(x = 0, bd = J, mu = h0X[ind10],
                                                sigma = rho.h0/(1-rho.h0), log = TRUE)))
            } else {
              p10 <- 0
            }
            
            if (sum(ind1J1) > 0) {
              p1J1 <- sum(wt[ind1J1] * (log(gX[ind1J1]) + dBB(J, bd = J, mu = h1X[ind1J1],
                                                sigma = rho.h1/(1-rho.h1), log = TRUE)))
            } else {
              p1J1 <- 0
            }
            
            if (sum(ind1y) > 0) {
              p1y <- sum(wt[ind1y] * (log(gX[ind1y]*dBB(y[ind1y]-1, bd = J, mu = h1X[ind1y], sigma =
                                           rho.h1/(1-rho.h1), log = FALSE) + (1-gX[ind1y])*
                             dBB(y[ind1y], bd = J, mu = h0X[ind1y],
                                 sigma = rho.h0/(1-rho.h0), log = FALSE))))
            } else {
              p1y <- 0
            }
            
            if (sum(treat == 0) > 0) {
              p0y <- sum(wt[!treat] * (log(gX[!treat]*dBB(y[!treat], bd = J, mu = h1X[!treat],
                                            sigma = rho.h1/(1-rho.h1), log = FALSE) +
                             (1-gX[!treat])*dBB(y[!treat], bd = J, mu = h0X[!treat],
                                                sigma = rho.h0/(1-rho.h0), log = FALSE))))
            } else {
              p0y <- 0
            }
            
            return(p10+p1J1+p1y+p0y)
            
          }
          
          ##
          ##  Observed data log-likelihood for binomial
          ##
          obs.llik.binom.std <- function(par, J, y, treat, x, wt, const = FALSE) {
            
            k <- ncol(x)
            if (const) {
              coef.h0 <- coef.h1 <- par[1:k]
              coef.g <- par[(k+1):(2*k)]
            } else {
              coef.h0 <- par[1:k]
              coef.h1 <- par[(k+1):(2*k)]
              coef.g <- par[(2*k+1):(3*k)]
            }
            
            h0X <- logistic(x %*% coef.h0)
            h1X <- logistic(x %*% coef.h1)
            gX <- logistic(x %*% coef.g)
            
            ind10 <- ((treat == 1) & (y == 0))
            ind1J1 <- ((treat == 1) & (y == (J+1)))
            ind1y <- ((treat == 1) & (y > 0) & (y < (J+1)))
            
            if (sum(ind10) > 0) {
              p10 <- sum(wt[ind10] * (log(1-gX[ind10]) + dbinom(x = 0, size = J, prob = h0X[ind10], log = TRUE)))
            } else {
              p10 <- 0
            }
            if (sum(ind1J1) > 0) {
              p1J1 <- sum(wt[ind1J1] * (log(gX[ind1J1]) + dbinom(J, size = J, prob = h1X[ind1J1], log = TRUE)))
            } else {
              p1J1 <- 0
            }
            if (sum(ind1y) > 0) {
              p1y <- sum(wt[ind1y] * (log(gX[ind1y]*dbinom(y[ind1y]-1, size = J, prob = h1X[ind1y],
                                              log = FALSE) + (1-gX[ind1y])*
                             dbinom(y[ind1y], size = J, prob = h0X[ind1y], log = FALSE))))
            } else {
              p1y <- 0
            }
            
            if (sum(treat == 0) > 0) {
              p0y <- sum(wt[!treat] * (log(gX[!treat]*dbinom(y[!treat], size = J, prob = h1X[!treat], log = FALSE) +
                             (1-gX[!treat])*dbinom(y[!treat], size = J, prob = h0X[!treat], log = FALSE))))
            } else {
              p0y <- 0
            }
            
            return(p10+p1J1+p1y+p0y)
          }
          
          ##
          ## EM algorithm
          ##
          
          ## Estep for beta-binomial
          Estep.std <- function(par, J, y, treat, x, const = FALSE) {
            
            k <- ncol(x)
            if (const) {
              coef.h0 <- coef.h1 <- par[1:k]
              rho.h0 <- rho.h1 <- logistic(par[(k+1)])
              coef.g <- par[(k+2):(2*k+1)]
            } else {
              coef.h0 <- par[1:k]
              rho.h0 <- logistic(par[(k+1)])
              coef.h1 <- par[(k+2):(2*k+1)]
              rho.h1 <- logistic(par[(2*k+2)])
              coef.g <- par[(2*k+3):(3*k+2)]
            }
            
            h0X <- logistic(x %*% coef.h0)
            h1X <- logistic(x %*% coef.h1)
            gX <- logistic(x %*% coef.g)
            
            ind <- !((treat == 1) & ((y == 0) | (y == (J+1))))
            w <- rep(NA, length(y))
            w[ind] <- gX[ind]*dBB((y-treat)[ind], bd = J, mu = h1X[ind], sigma = rho.h1/(1-rho.h1), log = FALSE) /
              (gX[ind]*dBB((y-treat)[ind], bd = J, mu = h1X[ind], sigma = rho.h1/(1-rho.h1), log = FALSE) +
               (1-gX[ind])*dBB(y[ind], bd = J, mu = h0X[ind], sigma = rho.h0/(1-rho.h0), log = FALSE))
            
            w[(treat == 1) & (y == 0)] <- 0
            w[(treat == 1) & (y == (J+1))] <- 1
            
            return(w)
          }
          
          
          ## Estep for binomial
          Estep.binom.std <- function(par, J, y, treat, x, const = FALSE) {
            
            k <- ncol(x)
            if (const) {
              coef.h0 <- coef.h1 <- par[1:k]
              coef.g <- par[(k+1):(2*k)]
            } else {
              coef.h0 <- par[1:k]
              coef.h1 <- par[(k+1):(2*k)]
              coef.g <- par[(2*k+1):(3*k)]
            }
            
            h0X <- logistic(x %*% coef.h0)
            h1X <- logistic(x %*% coef.h1)
            gX <- logistic(x %*% coef.g)
            
            ind <- !((treat == 1) & ((y == 0) | (y == (J+1))))
            w <- rep(NA, length(y))
            w[ind] <- gX[ind]*dbinom((y-treat)[ind], size = J, prob = h1X[ind], log = FALSE) /
              (gX[ind]*dbinom((y-treat)[ind], size = J, prob = h1X[ind], log = FALSE) +
               (1-gX[ind])*dbinom(y[ind], size = J, prob = h0X[ind], log = FALSE))
            w[(treat == 1) & (y == 0)] <- 0
            w[(treat == 1) & (y == (J+1))] <- 1
            
            return(w)
          }
          
          ## Mstep 1: weighted MLE for logistic regression
          wlogit.fit.std <- function(y, treat, x, w, par = NULL, wt) {
              ## wt is survey weights
              
            yrep <- rep(c(1,0), each = length(y))
            xrep <- rbind(x, x)
            wrep <- c(w, 1-w)
            wtrep <- c(wt, wt)
            return(glm(cbind(yrep, 1-yrep) ~ xrep - 1, weights = wrep * wtrep, family = binomial(logit), start = par))
            
          }
          
          ## Mstep 2: weighted MLE for beta-binomial regression
          wbb.fit.std <- function(J, y, treat, x, w, par0 = NULL, par1 = NULL) {
            
            Not0 <- ((treat == 1) & (y == (J+1)))
            y0 <- y[!Not0]
            x0 <- x[!Not0,]
            w0 <- 1-w[!Not0]
            fit0 <- vglm(cbind(y0, J-y0) ~ x0, betabinomial, weights = w0, coefstart = par0)
            
            Not1 <- ((treat == 1) & (y == 0))
            y1 <- y
            y1[treat == 1] <- y1[treat == 1] - 1
            y1 <- y1[!Not1]
            x1 <- x[!Not1]
            w1 <- w[!Not1]
            fit1 <- vglm(cbind(y1, J-y1) ~ x1, betabinomial, weights = w1, coefstart = par1)
            
            return(list(fit1 = fit1, fit0 = fit0))
          }
          
          ## Mstep2: weighted MLE for binomial regression
          wbinom.fit.std <- function(J, y, treat, x, w, par0 = NULL, par1 = NULL) {
            
            Not0 <- ((treat == 1) & (y == (J+1)))
            y0 <- y[!Not0]
            x0 <- x[!Not0,]
            w0 <- 1-w[!Not0]
            fit0 <- glm(cbind(y0, J-y0) ~ x0, family = binomial(logit), weights = w0, start = par0)
            
            Not1 <- ((treat == 1) & (y == 0))
            y1 <- y
            y1[treat == 1] <- y1[treat == 1] - 1
            y1 <- y1[!Not1]
            x1 <- x[!Not1,]
            w1 <- w[!Not1]
            fit1 <- glm(cbind(y1, J-y1) ~ x1, family = binomial(logit), weights = w1, start = par1)
            
            return(list(fit1 = fit1, fit0 = fit0))
          }
          
          ##
          ## Running the EM algorithm
          ##
                    
          if(constrained==F) {
            
            ## Run unconstrained model
            
            if (overdispersed==T) {
              par <- c(rep(c(coef.control.start, 0.5), 2), coef.treat.start)
            } else {
              par <- c(rep(coef.control.start, 2), coef.treat.start)
            }
            
            ## max number of iterations
            pllik <- -Inf
            
            if (overdispersed==T) {
              llik <- obs.llik.std(par, J = J, y = y.all, treat = t, x = x.all, wt = w.all)
            } else {
              llik <- obs.llik.binom.std(par, J = J, y = y.all, treat = t, x = x.all, wt = w.all)
            }
            
            Not0 <- (t & (y.all == (J+1)))
            Not1 <- (t & (y.all == 0))
            
            counter <- 0
            while (((llik - pllik) > 10^(-8)) & (counter < maxIter)) {
              
              if(overdispersed==T) {
                
                w <- Estep.std(par, J, y.all, t, x.all)
                
                lfit <- wlogit.fit.std(y.all, t, x.all, w, par = par[(nPar*2+3):length(par)], wt = w.all)
              
                y1 <- y.all
                
                y1[t==1] <- y1[t==1]-1

                if (intercept.only == TRUE) {
                  fit0 <- vglm(cbind(y.all[!Not0], J-y.all[!Not0]) ~ 1,
                               betabinomial, weights = (1-w[!Not0]) * w.all[!Not0], coefstart = par[1:2])
                  fit1 <- vglm(cbind(y1[!Not1], J-y1[!Not1]) ~ 1,
                               betabinomial, weights = w[!Not1] * w.all[!Not1],
                               coefstart = par[3:4])
                  par <- c(coef(fit0), coef(fit1), coef(lfit))
                  
                } else {
                  fit0 <- vglm(cbind(y.all[!Not0], J-y.all[!Not0]) ~ x.all[!Not0, -1, drop = FALSE],
                               betabinomial, weights = (1-w[!Not0]) * w.all[!Not0],
                               coefstart = par[c(1,(nPar+1),2:(nPar))])
                  fit1 <- vglm(cbind(y1[!Not1], J-y1[!Not1]) ~ x.all[!Not1, -1, drop = FALSE] ,
                               betabinomial, weights = w[!Not1] * w.all[!Not1],
                               coefstart = par[c((nPar+2),(2*nPar+2),(nPar+3):(2*nPar+1))])
                  par <- c(coef(fit0)[c(1,3:(nPar+1),2)], coef(fit1)[c(1,3:(nPar+1),2)], coef(lfit))
                }
                                
              } else {
                
                w <- Estep.binom.std(par, J, y.all, t, x.all)

                lfit <- wlogit.fit.std(y.all, t, x.all, w, par = par[(nPar*2+1):length(par)], wt = w.all)
                
                fit0 <- glm(cbind(y.all[!Not0], J-y.all[!Not0]) ~ x.all[!Not0,] - 1,
                            family = binomial(logit), weights = (1-w[!Not0]) * w.all[!Not0], start = par[1:(nPar)])
                
                y1 <- y.all
                
                y1[t==1] <- y1[t==1]-1
                
                fit1 <- glm(cbind(y1[!Not1], J-y1[!Not1]) ~ x.all[!Not1,] - 1,
                            family = binomial(logit), weights = w[!Not1] * w.all[!Not1],
                            start = par[(nPar+1):(2*nPar)])
                
                par <- c(coef(fit0), coef(fit1), coef(lfit))
                
              }
              
              pllik <- llik
              
              if(verbose==T)
                cat(paste(counter, round(llik, 4), "\n"))
              
              if (overdispersed==T) {
                llik <- obs.llik.std(par, J = J, y = y.all, treat = t, x = x.all, wt = w.all)
              } else {
                llik <- obs.llik.binom.std(par, J = J, y = y.all, treat = t, x = x.all, wt = w.all)
              }
              
              counter <- counter + 1
              
              if (llik < pllik) 
                warning("log-likelihood is not monotonically increasing.")
              ## NOTE CHANGED FROM STOP()
              
              if(counter == (maxIter-1))
                warning("number of iterations exceeded maximum in ML.")
              
            }
            
            ## getting  standard errors
            
            if (overdispersed==T) {
              MLEfit <- optim(par, obs.llik.std, method = "BFGS", J = J, y = y.all,
                              treat = t, x = x.all, wt = w.all, hessian = TRUE, control = list(maxit = 0))
              vcov.mle <- solve(-MLEfit$hessian, tol = 1e-20)
              se.mle <- sqrt(diag(vcov.mle))
            } else {
              MLEfit <- optim(par, obs.llik.binom.std, method = "BFGS", J = J, y = y.all,
                              treat = t, x = x.all, wt = w.all, hessian = TRUE, control = list(maxit = 0))
              vcov.mle <- solve(-MLEfit$hessian, tol = 1e-20)
              se.mle <- sqrt(diag(vcov.mle))
            }
            
          } else { # end of unconstrained model
            
            ## constrained model

            if (overdispersed==T) {
              par <- c(coef.control.start, 0.5, coef.treat.start)
            } else {
              par <- c(coef.control.start, coef.treat.start)
            }
            
            pllik.const <- -Inf
            
            if (overdispersed==T) {
              llik.const <- obs.llik.std(par, J = J, y = y.all, treat = t,
                                         x = x.all, wt = w.all, const = TRUE)
            } else {
              llik.const <- obs.llik.binom.std(par, J = J, y = y.all, treat = t,
                                               x = x.all, wt = w.all, const = TRUE)
            }
            
            counter <- 0
            while (((llik.const - pllik.const) > 10^(-8)) & (counter < maxIter)) {
              
              if (overdispersed==T) {
                
                w <- Estep.std(par, J, y.all, t, x.all, const = TRUE)
                lfit <- wlogit.fit.std(y.all, t, x.all, w, par = par[(nPar+2):(nPar*2+1)], wt = w.all)
                
              } else {
                
                w <- Estep.binom.std(par, J, y.all, t, x.all, const = TRUE)
                lfit <- wlogit.fit.std(y.all, t, x.all, w, par = par[(nPar+1):(nPar*2)], wt = w.all)
                
              }
              
              y.var <- as.character(formula)[[2]]
 
              data.all <- as.data.frame(cbind(y.all, x.all))
              names(data.all)[1] <- y.var
              
              dtmp <- rbind(cbind(data.all, w, t, w.all), cbind(data.all, w, t, w.all)[t==1, ])
              
              dtmp[((dtmp$t == 1) & ((1:nrow(dtmp)) <= n)), paste(y.var)] <-
                dtmp[((dtmp$t == 1) & ((1:nrow(dtmp)) <= n)), paste(y.var)] - 1
              
              dtmp$w[((dtmp$t == 1) & ((1:nrow(dtmp)) > n))] <-
                1 - dtmp$w[((dtmp$t == 1) & ((1:nrow(dtmp)) > n))]
              
              dtmp$w[dtmp$t == 0] <- 1

              dtmp$w <- dtmp$w * dtmp$w.all
              
              dtmp <- dtmp[dtmp$w > 0, ]
              
              if (overdispersed==T) {

                if (intercept.only == TRUE) {

                  fit <- vglm(as.formula(paste("cbind(", y.var, ", J-", y.var, ") ~ 1")),
                              betabinomial, weights = dtmp$w,
                              coefstart = par[1:2], data = dtmp)
                  
                  par <- c(coef(fit), coef(lfit))

                } else {
                  
                  fit <- vglm(as.formula(paste("cbind(", y.var, ", J-", y.var, ") ~ ",
                                               paste(x.vars, collapse=" + "))),
                              betabinomial, weights = dtmp$w,
                              coefstart = par[c(1,(nPar+1),2:(nPar))], data = dtmp)
                  
                  par <- c(coef(fit)[c(1,3:(nPar+1),2)], coef(lfit))

                }
                  
              } else {

                ## GB: This is where the error comes in for survey weights.

                fit <- glm(as.formula(paste("cbind(", y.var, ", J-", y.var, ") ~ ",
                                            paste(x.vars, collapse=" + "))),
                           family = binomial(logit), weights = dtmp$w,
                           start = par[1:(nPar)], data = dtmp)
                
                par <- c(coef(fit), coef(lfit))
                
              }
              
              pllik.const <- llik.const
              
              if(verbose==T)
                cat(paste(counter, round(llik.const, 4), "\n"))
              
              if (overdispersed==T) {
                llik.const <- obs.llik.std(par, J = J, y = y.all, treat = t,
                                           x = x.all, wt = w.all, const = TRUE)
              } else {
                llik.const <- obs.llik.binom.std(par, J = J, y = y.all, treat = t,
                                                 x = x.all, wt = w.all, const = TRUE)
              }
              
              counter <- counter + 1
              if (llik.const < pllik.const)
                warning("log-likelihood is not monotonically increasing.")
              ## NOTE CHANGED FROM STOP()
              
              if(counter == (maxIter-1))
                warning("number of iterations exceeded maximum in ML")
              
            }
            
            if (overdispersed==T) {
              
              MLEfit <- optim(par, obs.llik.std, method = "BFGS", J = J,
                              y = y.all, treat = t, x = x.all, wt = w.all, const = TRUE,
                              hessian = TRUE, control = list(maxit = 0))
              
              vcov.mle <- solve(-MLEfit$hessian, tol = 1e-20)
              se.mle <- sqrt(diag(vcov.mle))
              
            } else {
              
              MLEfit <- optim(par, obs.llik.binom.std, method = "BFGS", J = J,
                              y = y.all, treat = t,  x = x.all, wt = w.all, const = TRUE,
                              hessian = TRUE, control = list(maxit = 0))
              
              vcov.mle <- solve(-MLEfit$hessian, tol = 1e-20)
              se.mle <- sqrt(diag(vcov.mle))
              
            }

          } # end of constrained model

        } else if (multi == TRUE) {

          ##
          ## Start multi-treatment standard design model
          
          ## test start values
          ## par <- rep(0, 3*k + 2*(J+1))
            
          obs.llik.multi <- function(par, J, y, treat, x, multi = "none") {

            treat.values <- sort(unique(treat[treat!=0]))
            k <- ncol(x)

            ## pull out all g coefs, including alpha_yt and beta_t
            par.g <- par[(k+1):length(par)]

            ## pull out h coefs and construct h vector
            coef.h <- par[1:k]
            hX <- logistic(x %*% coef.h)

            ## separate coefs for alpha_yt and beta_t and construct
            ## g_t vector for each treatment condition by filling in a
            ## matrix with columns representing possible values of Y_i(0),
            ## rows individuals, and it is in a list where g_t = gX[[t]]
            coef.g.b <- coef.g.a <- gX <- list()
            for (m in 1:length(treat.values)) {

              if (multi == "indicators") {
                coef.g.b[[m]] <- par.g[((m-1)*(k+J)+1):((m-1)*(k+J)+k)]
                coef.g.a[[m]] <- c(par.g[((m-1)*(k+J)+k+1):((m-1)*(k+J)+k+J)],0)
              } else if (multi == "level") {
                coef.g.b[[m]] <- par.g[((m-1)*(k+1)+1):((m-1)*(k+1)+k)]
                coef.g.a[[m]] <- par.g[((m-1)*(k+1)+k+1)]
              } else if (multi == "none") {
                coef.g.b[[m]] <- par.g[((m-1)*k+1):((m-1)*k+k)]
              }

              gX[[m]] <- matrix(NA, nrow = length(y), ncol = J+1)

              if (multi == "indicators") {
                for(y.iter in 0:J) {
                  gX[[m]][,y.iter + 1] <- logistic(coef.g.a[[m]][y.iter + 1] + x %*% coef.g.b[[m]])
                }
              } else if (multi == "level") {
                for(y.iter in 0:J) {
                  gX[[m]][,y.iter + 1] <- logistic(y.iter * coef.g.a[[m]] + x %*% coef.g.b[[m]])
                }
              } else if (multi == "none") {
                for(y.iter in 0:J) {
                  gX[[m]][,y.iter + 1] <- logistic(x %*% coef.g.b[[m]])
                }
              }
            }

            ## contribution from control group
            
            ind0y <- (treat==0)
            p0y <- sum(dbinom(x = y[ind0y], size = J, prob = hX[ind0y], log = TRUE))

            ## contribution from treatment groups
            
            pm0 <- pmJ1 <- rep(NA, length(treat.values))
            pmy <- matrix(NA, nrow = length(treat.values), ncol = J)

            for (m in 1:length(treat.values)) {
              
              indm0 <- ((treat == m) & (y == 0))
              pm0[m] <- sum(dbinom(x = 0, size = J, prob = hX[indm0], log = TRUE) +
                            log(1 - gX[[m]][indm0, 0 + 1]))

              indmJ1 <- ((treat == m) & (y == (J+1)))
              pmJ1[m] <- sum(dbinom(x = J, size = J, prob = hX[indmJ1], log = TRUE) +
                             log(gX[[m]][indmJ1, J + 1]))

              for (y.iter in 1:J) {
                indmy <- ((treat == m) & (y == y.iter))
                pmy[m,y.iter] <- sum(log(gX[[m]][indmy, (y.iter - 1) + 1] *
                                         dbinom(x = y.iter - 1, size = J,
                                                prob = hX[indmy], log = FALSE) +
                 (1-gX[[m]][indmy, y.iter + 1]) * dbinom(x = y.iter, size = J,
                                                         prob = hX[indmy], log = FALSE)))
              }
                         
            }

            return(p0y + sum(pm0) + sum(pmJ1) + sum(pmy))
            
          }

          Estep.multi <- function(par, J, y, treat, x, m, multi = "none") {

            treat.values <- sort(unique(treat[treat!=0]))
            k <- ncol(x)

            ## pull out all g coefs, including alpha_yt and beta_t
            par.g <- par[(k+1):length(par)]

            ## pull out h coefs and construct h vector
            coef.h <- par[1:k]
            hX <- logistic(x %*% coef.h)

            if (multi == "indicators") {
              coef.g.b <- par.g[((m-1)*(k+J)+1):((m-1)*(k+J)+k)]
              coef.g.a <- c(par.g[((m-1)*(k+J)+k+1):((m-1)*(k+J)+k+J)],0)
            } else if (multi == "level") {
              coef.g.b <- par.g[((m-1)*(k+1)+1):((m-1)*(k+1)+k)]
              coef.g.a <- par.g[((m-1)*(k+1)+k+1)]
            } else if (multi == "none") {
              coef.g.b <- par.g[((m-1)*k+1):((m-1)*k+k)]
            }

            ## construct g_t matrix, where columns represent values of y_i(0),
            ## since g_t varies in y_i(0) through alpha_ty. Rows are indivs.

            ## in constrained model, g_t is constant across y_i(0) values
            
            gX <- matrix(NA, nrow = length(y), ncol = J+1)

            if (multi == "indicators") {
              for(y.iter in 0:J) {
                gX[,y.iter + 1] <- logistic(coef.g.a[y.iter + 1] + x %*% coef.g.b)
              }
            } else if (multi == "level") {
              for(y.iter in 0:J) {
                gX[,y.iter + 1] <- logistic(y.iter * coef.g.a + x %*% coef.g.b)
              }
            } else if (multi == "none") {
              for(y.iter in 0:J) {
                gX[,y.iter + 1] <- logistic(x %*% coef.g.b)
              }
            }

            w <- rep(NA, length(y))

            ind0 <- ((y==0) & (treat==m))
            w[ind0] <- 0

            indJ1 <- ((y==(J+1)) & (treat==m))
            w[indJ1] <- 1

            ## for values of Y_i(0) not 0 or (J+1), construct the weight by
            ## looping through all possible values of Y_i(0) and adding weights
            ## for individuals with that Y_i(0)
            
            for (y.iter in 1:J) {
              
              indother <- ((ind0==FALSE) & (indJ1==FALSE) & (y==y.iter) & (treat==m))
              w[indother] <- (gX[indother, (y.iter - 1) + 1] *
                              dbinom(x = y.iter - 1, size = J, prob = hX[indother],
                                     log = FALSE)) /
                                (gX[indother, (y.iter - 1) + 1] *
                                 dbinom(x = y.iter - 1, size = J, prob = hX[indother],
                                        log = FALSE) +
                                 (1 - gX[indother, y.iter + 1]) *
                                 dbinom(x = y.iter, size = J, prob = hX[indother],
                                        log = FALSE) )
              
            }
            
            return(w)
            
          }

          ##
          ## start EM algorithm

          ## get starting values from NLS
          
          par <- par.control.nls.std
          if (multi.condition == "none") {
            par <- c(par, do.call(c, par.treat.nls.std))
          } else if (multi.condition == "indicators") {
            for (m in 1:length(treatment.values)){
              par <- c(par, par.treat.nls.std[[m]], rep(0, J))
            }
          } else if (multi.condition == "level") {
            for (m in 1:length(treatment.values)){
              par <- c(par, par.treat.nls.std[[m]], 0)
            }
          }
           
          pllik <- -Inf
          
          llik <- obs.llik.multi(par, J = J, y = y.all, treat = t, x = x.all,
                                 multi = multi.condition)

          ##
          ## create temporary datasets for M step

          ## start fit data for H including control group data
          
          ytmpH <- y.all[t==0]
          xtmpH <- x.all[t==0, , drop = FALSE]
          
          ## now loop over treatment conditions to construct data for g_t fit
          ## and add treatment condition data (m times) to data for h fit

          ytmpG <- xtmpG <- dvtmpG <- list()

          for (m in 1:length(treatment.values)){

            ##
            ## set up temporary data for g fit
            
            ytmpG[[m]] <- c(y.all[t==m]-1, y.all[t==m])
            xtmpG[[m]] <- rbind(x.all[t==m, , drop = FALSE], x.all[t==m, , drop = FALSE])
            dvtmpG[[m]] <- c(rep(1, sum(t==m)), rep(0, sum(t==m)))

            if (multi.condition == "indicators") {
              
              ## create matrix of indicator variables for non-sensitive y_i(0) value
              ## omits indicator for value == J              
              y.ind <- matrix(NA, nrow = sum(t==m) * 2, ncol = J)
              for (y.iter in 0:(J-1))
                y.ind[, y.iter + 1] <- ytmpG[[m]] == y.iter
              
              ## x vector for g_t fit includes covariates and y_i(0) indicator matrix
              ## the constrained model does not include these indicators in the g_t M fits              
              xtmpG[[m]] <- cbind(xtmpG[[m]], y.ind)
              
            } else if (multi.condition == "level") {

              xtmpG[[m]] <- cbind(xtmpG[[m]], ytmpG[[m]])
              
            }
            
            ##
            ## set up temporary data for h fit
            
            ytmpH <- c(ytmpH, y.all[t==m]-1, y.all[t==m])
            xtmpH <- rbind(xtmpH, x.all[t==m, , drop = FALSE], x.all[t==m, , drop = FALSE])
            
          }

          counter <- 0
          while (((llik - pllik) > 10^(-8)) & (counter < maxIter)) {
            
            wtmpH <- rep(1, sum(t==0))
            
            lfit <- list()
            
            for (m in 1:length(treatment.values)){

              ## get E step values for this treatment condition
              w <- Estep.multi(par, J, y.all, t, x.all, m, multi = multi.condition)

              ## create weight vectors for temporary data for fits

              wtmpG <- c(w[t==m], 1-w[t==m])

              ## for h, we are adding weights to the vector of weights = 1 created earlier
              
              wtmpH <- c(wtmpH, w[t==m], 1-w[t==m])

              if (multi.condition == "none")
                par.start <- par[(nPar + (m-1) * nPar + 1):(nPar + (m-1) * nPar + nPar)]
              else if (multi.condition == "indicators")
                par.start <- par[(nPar + (m - 1) * (nPar + J) + 1):
                                 (nPar + (m - 1) * (nPar + J) + nPar + J)]
              else if (multi.condition == "level")
                par.start <- par[(nPar + (m - 1) * (nPar + 1) + 1):
                                 (nPar + (m - 1) * (nPar + 1) + nPar + 1)]

              ## run g_t fits
              lfit[[m]] <- glm(cbind(dvtmpG[[m]][wtmpG>0], 1 - dvtmpG[[m]][wtmpG>0]) ~
                               xtmpG[[m]][wtmpG>0, , drop = FALSE] - 1, weights = wtmpG[wtmpG>0],
                               family = binomial(logit), start = par.start,
                               control = glm.control(maxit = maxIter))

            }

            ## run h fit
            fit <- glm(cbind(ytmpH[wtmpH>0], J - ytmpH[wtmpH>0]) ~ xtmpH[wtmpH>0, , drop = FALSE] - 1,
                       family = binomial(logit), weights = wtmpH[wtmpH>0], start = par[1:(nPar)])

            ## put coefficients back together
            par <- coef(fit)
            par.treat <- list()
            for (m in 1:length(treatment.values)) {
               par <- c(par, coef(lfit[[m]]))
               par.treat[[m]] <- coef(lfit[[m]])
            }

            pllik <- llik
            
            if(verbose==T)
              cat(paste(counter, round(llik, 4), "\n"))
            
            llik <- obs.llik.multi(par, J = J, y = y.all, treat = t, x = x.all,
                                   multi = multi.condition)
            
            counter <- counter + 1
            if (llik < pllik)
              warning("log-likelihood is not monotonically increasing.")
            ## NOTE CHANGED FROM STOP()
            
            if(counter == (maxIter-1))
              warning("number of iterations exceeded maximum in ML")
            
          }

          MLEfit <- optim(par, obs.llik.multi, method = "BFGS", J = J, y = y.all, treat = t,
                          x = x.all, hessian = TRUE, control = list(maxit = 0),
                          multi = multi.condition)
          
          vcov.mle <- solve(-MLEfit$hessian, tol = 1e-20)
          
          se.mle <- sqrt(diag(vcov.mle))
  
        } # end of multi treat regression
        
      } else if (boundary == TRUE) {
                
        coef.control.start <- par.control.nls.std
        coef.treat.start <- par.treat.nls.std

        ##
        ##  Observed data log-likelihood for binomial
        ##
        
        obs.llik.binom.boundary <- function(par, J, y, treat, x, x.floor, x.ceiling, ceiling.fit,
                                            floor.fit, ceiling, floor) {
          
          k <- ncol(x)
          k.ceiling <- ncol(x.ceiling)
          k.floor <- ncol(x.floor)
          
          coef.h <- par[1:k]
          coef.g <- par[(k+1):(2*k)]
          coef.ql <- par[(2*k+1):(2*k + k.floor)] 
          coef.qu <- par[(2*k+k.floor+1):(2*k+k.floor+k.ceiling)] 
          
          hX <- logistic(x %*% coef.h)
          gX <- logistic(x %*% coef.g)
          quX <- logistic(x.ceiling %*% coef.qu)
          qlX <- logistic(x.floor %*% coef.ql)
          
          ind0y <- ((treat == 0) & (y < (J+1)))
          ind10 <- ((treat == 1) & (y == 0))
          ind11 <- ((treat == 1) & (y == 1))
          ind1y <- ((treat == 1) & (y > 1) & (y < J))
          ind1J <- ((treat == 1) & (y == J))
          ind1J1 <- ((treat == 1) & (y == (J+1)))
          
          p0y <- p10 <- p11 <- p1y <- p1J <- p1J1 <- 0
          
          if(sum(ind0y) > 0){
            p0y <- sum(dbinom(x = y[ind0y], size = J, prob = hX[ind0y], log = TRUE))
          }
          
          if(sum(ind10) > 0){
            p10 <- sum(log((1-gX[ind10]) + gX[ind10] * qlX[ind10])
                       + dbinom(x = 0, size = J, prob = hX[ind10], log = TRUE))
          }
          
          if (sum(ind11) > 0) {
            p11 <- sum(log( gX[ind11] * (1-qlX[ind11]) *
                           dbinom(x = 0, size = J, prob = hX[ind11], log = FALSE) +
                           (1-gX[ind11]) * dbinom(x = 1, size = J, prob = hX[ind11], log = FALSE) ))
          }
          
          if (sum(ind1y) > 0) {
            p1y <- sum(log( gX[ind1y] * dbinom(x = y[ind1y]-1, size = J, prob = hX[ind1y], log = FALSE) +
                           (1-gX[ind1y]) * dbinom(x = y[ind1y], size = J, prob = hX[ind1y], log = FALSE) ))
          }
          
          if (sum(ind1J)) {
            p1J <- sum(log( gX[ind1J] * (dbinom(x = J-1, size = J, prob = hX[ind1J], log = FALSE) +
                                         quX[ind1J] * dbinom(x = J, size = J, prob = hX[ind1J], log = FALSE)) +
                           (1-gX[ind1J]) * dbinom(x = J, size = J, prob = hX[ind1J], log = FALSE) ))
          }
          
          if (sum(ind1J1) > 0) {
            p1J1 <- sum(log(1-quX[ind1J1]) + log(gX[ind1J1]) + dbinom(x = J, size = J, prob = hX[ind1J1], log = TRUE))
          }
          
          if (ceiling.fit=="bayesglm" & ceiling == TRUE) {
            p.prior.ceiling <- sum(dcauchy(x = coef.qu,
                                           scale = rep(2.5, length(coef.qu)), log = TRUE))
          } else {
            p.prior.ceiling <- 0
          }
          
          if (floor.fit=="bayesglm" & floor == TRUE) {
            p.prior.floor <- sum(dcauchy(x = coef.ql,
                                         scale = rep(2.5, length(coef.ql)), log = TRUE))
          } else {
            p.prior.floor <- 0
          }
					
          return(p0y  + p10  + p11  + p1y  + p1J  + p1J1 + p.prior.ceiling + p.prior.floor)
        }
				
        ##
        ##  Observed data log-likelihood for binomial for optim -- shortened par
        ##
        obs.llik.binom.optim.boundary <- function(par, J, y, treat, x, x.floor, x.ceiling,
                                                  ceiling.fit, floor.fit, floor, ceiling) {
          
          k <- ncol(x)
          k.ceiling <- ncol(x.ceiling)
          k.floor <- ncol(x.floor)
          coef.h <- par[1:k]
          coef.g <- par[(k+1):(2*k)]
          if(floor==TRUE) {
            coef.ql <- par[(2*k+1):(2*k + k.floor)]
            if(ceiling==TRUE){
              coef.qu <- par[(2*k+k.floor+1):(2*k+k.floor+k.ceiling)] 
            } else {
              coef.qu <- c(-Inf, rep(0, (k.ceiling-1)))
            }
          } else {
            coef.ql <- c(-Inf, rep(0, (k.floor-1)))
            if(ceiling==TRUE){
              coef.qu <- par[(2*k+1):(2*k + k.ceiling)] 
            } else {
              coef.qu <- c(-Inf, rep(0, (k.ceiling-1)))
            }
          }
          
          hX <- logistic(x %*% coef.h)
          gX <- logistic(x %*% coef.g)
          quX <- logistic(x.ceiling %*% coef.qu)
          qlX <- logistic(x.floor %*% coef.ql)
          
          ind0y <- ((treat == 0) & (y < (J+1)))
          ind10 <- ((treat == 1) & (y == 0))
          ind11 <- ((treat == 1) & (y == 1))
          ind1y <- ((treat == 1) & (y > 1) & (y < J))
          ind1J <- ((treat == 1) & (y == J))
          ind1J1 <- ((treat == 1) & (y == (J+1)))
          
          p0y <- p10 <- p11 <- p1y <- p1J <- p1J1 <- 0
          
          if(sum(ind0y) > 0){
            p0y <- sum(dbinom(x = y[ind0y], size = J, prob = hX[ind0y], log = TRUE))
          }
          
          if(sum(ind10) > 0){
            p10 <- sum(log((1-gX[ind10]) + gX[ind10] * qlX[ind10]) + dbinom(x = 0, size = J, prob = hX[ind10], log = TRUE))
          }
          
          if(sum(ind11) > 0){
            p11 <- sum(log( gX[ind11] * (1-qlX[ind11]) * dbinom(x = 0, size = J, prob = hX[ind11], log = FALSE) + (1-gX[ind11]) * dbinom(x = 1, size = J, prob = hX[ind11], log = FALSE) ))
          }
          
          if(sum(ind1y) > 0){
            p1y <- sum(log( gX[ind1y] * dbinom(x = y[ind1y]-1, size = J, prob = hX[ind1y], log = FALSE) + (1-gX[ind1y]) * dbinom(x = y[ind1y], size = J, prob = hX[ind1y], log = FALSE) ))
          }
          
          if(sum(ind1J)){
            p1J <- sum(log( gX[ind1J] * (dbinom(x = J-1, size = J, prob = hX[ind1J], log = FALSE) + quX[ind1J] * dbinom(x = J, size = J, prob = hX[ind1J], log = FALSE)) + (1-gX[ind1J]) * dbinom(x = J, size = J, prob = hX[ind1J], log = FALSE) ))
          }
          
          if(sum(ind1J1) > 0){
            p1J1 <- sum(log(1-quX[ind1J1]) + log(gX[ind1J1]) + dbinom(x = J, size = J, prob = hX[ind1J1], log = TRUE))
          }
          
          if(ceiling.fit=="bayesglm" & ceiling==TRUE) {
            p.prior.ceiling <- sum(dcauchy(x = coef.qu,
                                         scale = rep(2.5, length(coef.qu)), log = TRUE))
          } else {
            p.prior.ceiling <- 0
          }
          
          if(floor.fit=="bayesglm" & floor==TRUE) {
            p.prior.floor <- sum(dcauchy(x = coef.ql,
                                         scale = rep(2.5, length(coef.ql)), log = TRUE))
          } else {
            p.prior.floor <- 0
          }
          
          return(p0y  + p10  + p11  + p1y  + p1J  + p1J1 + p.prior.ceiling + p.prior.floor)
          
        }
        
        ##
        ## EM algorithm
        ##
        
        ## Estep for beta-binomial
        Estep.boundary <- function(par, J, y, treat, x, product) {

          ## not currently used
          ## needs to be updated for various changes, including ceiling/floor formulaes
          
          k <- ncol(x)
          coef.h <- par[1:k]
          rho.h <- logistic(par[(k+1)])
          coef.g <- par[(k+2):(2*k+1)]
          coef.ql <- par[(2*k+2):(3*k+1)] 
          coef.qu <- par[(3*k+2):(4*k+1)] 
          
          hX <- logistic(x %*% coef.h)
          gX <- logistic(x %*% coef.g)
          quX <- logistic(x %*% coef.qu)
          qlX <- logistic(x %*% coef.ql)
          
          ind0 <- (y==0)
          ind1 <-  (y==1)
          indJ <- (y==J)
          indJ1 <- (y==(J+1))
          indother <- (y > 1 & y < J)
          
          w <- rep(NA, length(y))
          
          if(product == FALSE){
            
            w[ind0] <- dBB(0, bd = J, mu = hX[ind0], sigma = rho.h/(1-rho.h), log = FALSE)*gX[ind0]*qlX[ind0] /
              (dBB(0, bd = J, mu = hX[ind0], sigma = rho.h/(1-rho.h), log = FALSE) * (gX[ind0]*qlX[ind0] + (1-gX[ind0]) ) )
            
            w[ind1] <- dBB(0, bd = J, mu = hX[ind1], sigma = rho.h/(1-rho.h), log = FALSE)*gX[ind1]*(1-qlX[ind1]) / ( dBB(0, bd = J, mu = hX[ind1], sigma = rho.h/(1-rho.h), log = FALSE) * gX[ind1] * (1-qlX[ind1]) + dBB(1, bd = J, mu = hX[ind1], sigma = rho.h/(1-rho.h), log = FALSE)*(1-gX[ind1]) )
            
            w[indother] <- gX[indother]*dBB((y-1)[indother], bd = J, mu = hX[indother], sigma = rho.h/(1-rho.h), log = FALSE) /
              (gX[indother]*dBB((y-1)[indother], bd = J, mu = hX[indother], sigma = rho.h/(1-rho.h), log = FALSE) +
               (1-gX[indother])*dBB(y[indother], bd = J, mu = hX[indother], sigma = rho.h/(1-rho.h), log = FALSE))
            
            w[indJ] <- gX[indJ]*( dBB(J-1, bd = J, mu = hX[indJ], sigma = rho.h/(1-rho.h), log = FALSE) + quX[indJ]*dBB(J, bd = J, mu = hX[indJ], sigma = rho.h/(1-rho.h), log = FALSE) ) / ( gX[indJ] * ( dBB(J-1, bd = J, mu = hX[indJ], sigma = rho.h/(1-rho.h), log = FALSE) + quX[indJ] * dBB(J, bd = J, mu = hX[indJ], sigma = rho.h/(1-rho.h), log = FALSE)) + (1-gX[indJ])*dBB(J, bd = J, mu = hX[indJ], sigma = rho.h/(1-rho.h), log = FALSE) )
            
            w[indJ1] <- 1
            
          }
          
          if(product == TRUE){
            
            w[ind0] <- 0
            
            w[ind1] <- dBB(0, bd = J, mu = hX[ind1], sigma = rho.h/(1-rho.h), log = FALSE)*gX[ind1]*(1-qlX[ind1]) / ( dBB(0, bd = J, mu = hX[ind1], sigma = rho.h/(1-rho.h), log = FALSE) * gX[ind1] * (1-qlX[ind1]) + dBB(1, bd = J, mu = hX[ind1], sigma = rho.h/(1-rho.h), log = FALSE)*(1-gX[ind1]) )
            
            w[indother] <- gX[indother]*dBB((y-1)[indother], bd = J, mu = hX[indother], sigma = rho.h/(1-rho.h), log = FALSE) / (gX[indother]*dBB((y-1)[indother], bd = J, mu = hX[indother], sigma = rho.h/(1-rho.h), log = FALSE) + (1-gX[indother])*dBB(y[indother], bd = J, mu = hX[indother], sigma = rho.h/(1-rho.h), log = FALSE))
            
            w[indJ] <- gX[indJ]*dBB(J - 1, bd = J, mu = hX[indJ], sigma = rho.h/(1-rho.h), log = FALSE) / (gX[indJ]*(dBB(J - 1, bd = J, mu = hX[indJ], sigma = rho.h/(1-rho.h), log = FALSE) + quX[indJ]*dBB(J, bd = J, mu = hX[indJ], sigma = rho.h/(1-rho.h), log = FALSE)) + (1-gX[indJ])*dBB(J, bd = J, mu = hX[indJ], sigma = rho.h/(1-rho.h), log = FALSE))
            
            w[indJ1] <- 1
            
          }
          
          return(w)
          
        }
        
        ## Estep for binomial
        Estep.binom.boundary <- function(par, J, y, treat, x, x.floor, x.ceiling, product) {
          
          k <- ncol(x)
          k.ceiling <- ncol(x.ceiling)
          k.floor <- ncol(x.floor)
          coef.h <- par[1:k]
          coef.g <- par[(k+1):(2*k)]
          coef.ql <- par[(2*k+1):(2*k+k.floor)] 
          coef.qu <- par[(2*k+k.floor+1):(2*k+k.floor+k.ceiling)] 
          
          hX <- logistic(x %*% coef.h)
          gX <- logistic(x %*% coef.g)
          quX <- logistic(x.ceiling %*% coef.qu)
          qlX <- logistic(x.floor %*% coef.ql)
          
          ind0 <- (y==0)
          ind1 <-  (y==1)
          indJ <- (y==J)
          indJ1 <- (y==(J+1))
          indother <- (y > 1 & y < J)
          
          w <- rep(NA, length(y))
          
          if(product == FALSE){
            
            w[ind0] <- (dbinom(x = 0, size = J, prob = hX[ind0], log = FALSE) * gX[ind0] * qlX[ind0]) / (dbinom(x = 0, size = J, prob = hX[ind0], log = FALSE) * (gX[ind0] * qlX[ind0] + (1-gX[ind0])))
            
            w[ind1] <- (dbinom(x = 0, size = J, prob = hX[ind1], log = FALSE) * gX[ind1] * (1-qlX[ind1])) / (dbinom(x = 0, size = J, prob = hX[ind1], log = FALSE) * gX[ind1] * (1-qlX[ind1]) + dbinom(x = 1, size = J, prob = hX[ind1], log = FALSE) * (1-gX[ind1]))
            
            w[indother] <- (gX[indother] * dbinom(x = y[indother]-1, size = J, prob = hX[indother], log = FALSE)) / (gX[indother] * dbinom(x = y[indother]-1, size = J, prob = hX[indother], log = FALSE) + (1-gX[indother]) * dbinom(x = y[indother], size = J, prob = hX[indother], log = FALSE) )
            
            w[indJ] <- (gX[indJ] * (dbinom(x = J-1, size = J, prob = hX[indJ], log = FALSE) + quX[indJ] * dbinom(x = J, size = J, prob = hX[indJ], log = FALSE))) / (gX[indJ] * (dbinom(x = J-1, size = J, prob = hX[indJ], log = FALSE) + quX[indJ] * dbinom(x = J, size = J, prob = hX[indJ], log = FALSE)) + (1-gX[indJ]) * dbinom(x = J, size = J, prob = hX[indJ], log = FALSE))
            
            w[indJ1] <- 1
            
          }
          
          if(product == TRUE){
            
            w[ind0] <- 0
            
            w[ind1] <- (dbinom(x = 0, size = J, prob = hX[ind1], log = FALSE) * gX[ind1] * (1-qlX[ind1])) / (dbinom(x = 0, size = J, prob = hX[ind1], log = FALSE) * gX[ind1] * (1-qlX[ind1]) + dbinom(x = 1, size = J, prob = hX[ind1], log = FALSE) * (1-gX[ind1]))
            
            w[indother] <- (gX[indother] * dbinom(x = y[indother]-1, size = J, prob = hX[indother], log = FALSE)) / (gX[indother] * dbinom(x = y[indother]-1, size = J, prob = hX[indother], log = FALSE) + (1-gX[indother]) * dbinom(x = y[indother], size = J, prob = hX[indother], log = FALSE) )
            
            w[indJ] <- (gX[indJ] * dbinom(x = J-1, size = J, prob = hX[indJ], log = FALSE)) / (gX[indJ] * (dbinom(x = J-1, size = J, prob = hX[indJ], log = FALSE) + quX[indJ] * dbinom(x = J, size = J, prob = hX[indJ], log = FALSE)) + (1-gX[indJ]) * dbinom(x = J, size = J, prob = hX[indJ], log = FALSE))
        
            w[indJ1] <- 1
            
          }
          
          return(w)
          
        }
        
        ## Mstep 1: weighted MLE for logistic regression
        wlogit.fit.boundary <- function(y, x, w, par = NULL, maxIter = 50000) {
          yrep <- rep(c(1,0), each = length(y))
          xrep <- rbind(x, x)
          wrep <- c(w, 1-w)
          return(glm(cbind(yrep, 1-yrep) ~ xrep - 1, weights = wrep, family = binomial(logit), start = par, control = glm.control(maxit = maxIter)))
          
        }
        
        ## Mstep2: weighted MLE for binomial regression
        wbinom.fit.boundary <- function(J, y, treat, x, w, par0 = NULL, par1 = NULL) {
          
          Not0 <- ((treat == 1) & (y == (J+1)))
          y0 <- y[!Not0]
          x0 <- x[!Not0,]
          w0 <- 1-w[!Not0]
          fit0 <- glm(cbind(y0, J-y0) ~ x0, family = binomial(logit), weights = w0, start = par0)
          
          Not1 <- ((treat == 1) & (y == 0))
          y1 <- y
          y1[treat == 1] <- y1[treat == 1] - 1
          y1 <- y1[!Not1]
          x1 <- x[!Not1,]
          w1 <- w[!Not1]
          fit1 <- glm(cbind(y1, J-y1) ~ x1, family = binomial(logit), weights = w1, start = par1)
          
          return(list(fit1 = fit1, fit0 = fit0))
        }
        
        ##
        ## Running the EM algorithm
        ##

        if (overdispersed == T)
          stop("The ceiling-floor liar model does not support overdispersion. Set overdispersed = F.")
                
        y.var <- as.character(formula)[[2]]

        x.ceiling.all <- model.matrix(as.formula(ceiling.formula), as.data.frame(x.all))
        x.ceiling.treatment <- model.matrix(as.formula(ceiling.formula), as.data.frame(x.treatment))
        x.ceiling.control <- model.matrix(as.formula(ceiling.formula), as.data.frame(x.control))
        
        nPar.ceiling <- ncol(x.ceiling.all)
        
        intercept.only.ceiling <- ncol(x.ceiling.all)==1 & sum(x.ceiling.all[,1]==1) == n
        
        coef.names.ceiling <- colnames(x.ceiling.all)
        
        if (intercept.only.ceiling == TRUE) {
          x.vars.ceiling <- "1"
        } else {
          x.vars.ceiling <- coef.names.ceiling[-1]
        }
        
        x.floor.all <- model.matrix(as.formula(floor.formula), as.data.frame(x.all))
        x.floor.treatment <- model.matrix(as.formula(floor.formula), as.data.frame(x.treatment))
        x.floor.control <- model.matrix(as.formula(floor.formula), as.data.frame(x.control))
        nPar.floor <- ncol(x.floor.all)
        
        intercept.only.floor <- ncol(x.floor.all)==1 & sum(x.floor.all[,1]==1) == n
        
        coef.names.floor <- colnames(x.floor.all)
        
        if (intercept.only.floor == TRUE) {
          x.vars.floor <- "1"
        } else {
          x.vars.floor <- coef.names.floor[-1]
        }
        
        data.all <- as.data.frame(cbind(y.all, x.all))
        names(data.all)[1] <- y.var
        
        data.treatment <- as.data.frame(cbind(y.treatment, x.treatment))
        names(data.treatment)[1] <- y.var
        
        ## set start values
        
        data.ceiling.all <- as.data.frame(cbind(y.all, x.ceiling.all))
        names(data.ceiling.all)[1] <- y.var
       
        if (ceiling == TRUE){
    
          data.ceiling.treatment <- as.data.frame(cbind(y.treatment, x.ceiling.treatment))
          names(data.ceiling.treatment)[1] <- y.var
          
          dtmpS <- data.treatment[data.ceiling.treatment[,c(y.var)]==J |
                                  data.ceiling.treatment[,c(y.var)]==(J+1), ]
          
          dtmpS$dv <- dtmpS[,c(y.var)] == J
          
          coef.ceiling.start <- coef(glm(as.formula(paste("cbind(dv, 1-dv) ~ ",
                                                          paste(x.vars.ceiling, collapse=" + "))),
                                         family = binomial(logit), data = dtmpS,
                                         control = glm.control(maxit = maxIter)))
        } else {
          coef.ceiling.start <- c(-Inf, rep(0, (nPar.ceiling-1)))
        }

        data.floor.all <- as.data.frame(cbind(y.all, x.floor.all))
        names(data.floor.all)[1] <- y.var
        
        if (floor == TRUE){
          
          data.floor.treatment <- as.data.frame(cbind(y.treatment, x.floor.treatment))
          names(data.floor.treatment)[1] <- y.var
          
          dtmpS <- data.floor.treatment[data.floor.treatment[,c(y.var)]==0 |
                                        data.floor.treatment[,c(y.var)]==1, ]
          
          dtmpS$dv <- dtmpS[,c(y.var)] == 1
          
          coef.floor.start <- coef(glm(as.formula(paste("cbind(dv, 1-dv) ~ ",
                                                        paste(x.vars.floor, collapse=" + "))),
                                       family = binomial(logit), data = dtmpS,
                                       control = glm.control(maxit = maxIter)))
          
        } else {
          coef.floor.start <- c(-Inf, rep(0, (nPar.floor-1)))
        }
        
        par <- c(coef.control.start, coef.treat.start, coef.floor.start, coef.ceiling.start)
       
        ## calculate starting log likelihood

        pllik <- -Inf
        
        llik <- obs.llik.binom.boundary(par, J = J, y = y.all, treat = t, x = x.all,
                                        x.ceiling = x.ceiling.all, x.floor = x.floor.all,
                                        ceiling.fit = ceiling.fit, floor.fit = floor.fit,
                                        ceiling = ceiling, floor = floor)
        
        ## begin EM iterations
        
        counter <- 0
        while (((llik - pllik) > 10^(-8)) & (counter < maxIter)) {
          
          w <- Estep.binom.boundary(par, J, y.all, t, x.all,
                                    x.floor.all, x.ceiling.all, product = FALSE)
          w.prod <- Estep.binom.boundary(par, J, y.all, t, x.all,
                                         x.floor.all, x.ceiling.all, product = TRUE)
          
          lfit <- wlogit.fit.boundary(y.all[t==1], x.all[t==1, , drop = FALSE], w[t==1],
                             par = par[(nPar+1):(nPar*2)], maxIter = maxIter)
          
          dtmp <- rbind(data.all[t==0, ], data.all[t==1, ], data.all[t==1, ], data.all[t==1, ])
          
          dtmp$w <- c(rep(1, sum(t==0)), (w-w.prod)[t==1], (1-w)[t==1], w.prod[t==1])
          
          dtmp[(sum(c(t==0, t==1, t==1))+1):nrow(dtmp),paste(y.var)] <-
            dtmp[(sum(c(t==0, t==1, t==1))+1):nrow(dtmp),paste(y.var)] - 1
          
          dtmp <- dtmp[dtmp$w > 0, ]
          
          ## temporary data for ceiling logit
          dtmpC <- rbind(data.ceiling.all[t==1 & y.all==(J+1),], data.ceiling.all[t==1 & y.all==J, ])
          
          dtmpC[,paste(y.var)] <- c(rep(0,sum(t==1 & y.all==(J+1))), rep(1,sum(t==1 & y.all==J)))
          
          dtmpC$w <- c( w.prod[t==1 & y.all==(J+1)], (w - w.prod)[t==1 & y.all==J])
          
          dtmpC <- dtmpC[dtmpC$w > 0, ]
          
          ## temporary data for floor logit
          dtmpF <- rbind(data.floor.all[t==1 & y.all==1,], data.floor.all[t==1 & y.all==0, ])
          
          dtmpF[,paste(y.var)] <- c(rep(0,sum(t==1 & y.all==1)), rep(1,sum(t==1 & y.all==0)))
          
          dtmpF$w <- c(w.prod[t==1 & y.all==1], (w-w.prod)[t==1 & y.all==0])
          
          dtmpF <- dtmpF[dtmpF$w > 0, ]
                      
          fit <- glm(as.formula(paste("cbind(", y.var, ", J-", y.var, ") ~ ",
                                      paste(x.vars, collapse=" + "))),
                     family = binomial(logit), weights = dtmp$w, start = par[1:(nPar)],
                     data = dtmp)
          
          if (ceiling==TRUE) {
            ##coef.qufit.start <- par[(3*nPar+1):(4*nPar)]
            coef.qufit.start <- par[(2*nPar + nPar.floor + 1):(2*nPar + nPar.floor + nPar.ceiling)]
          } else {
            coef.qufit.start <- rep(0, nPar.ceiling)
          }

          if (ceiling == TRUE) {
            if (ceiling.fit=="glm") {
              qufit <- glm(as.formula(paste("cbind(", y.var, ", 1-", y.var, ") ~ ",
                                            paste(x.vars.ceiling, collapse=" + "))),
                           weights = dtmpC$w, family = binomial(logit),
                           start = coef.qufit.start, data = dtmpC,
                           control = glm.control(maxit = maxIter))
            } else if(ceiling.fit=="bayesglm") {
              
              if (intercept.only.ceiling == F) {
                qufit <- bayesglm.internal(as.formula(paste("cbind(", y.var, ", 1-", y.var, ") ~ ",
                                                   paste(x.vars.ceiling, collapse=" + "))),
                                  weights = dtmpC$w, family = binomial(logit),
                                  start = coef.qufit.start, data = dtmpC,
                                  control = glm.control(maxit = maxIter), scaled = F)
                
              } else {
                qufit <- bayesglm.internal(as.formula(paste("cbind(", y.var, ", 1-", y.var, ") ~ 1")),
                                  weights = dtmpC$w, family = binomial(logit),
                                  start = coef.qufit.start, data = dtmpC,
                                  control = glm.control(maxit = maxIter), scaled = F)
                
              }
            }
          }
          
          if(floor==TRUE) {
            ##coef.qlfit.start <- par[(2*nPar+1):(3*nPar)]
            coef.qlfit.start <- par[(2*nPar+1):(2*nPar + nPar.floor)]
          } else {
            coef.qlfit.start <- rep(0, nPar.floor)
          }

          if (floor == TRUE) {
            if(floor.fit=="glm") {
              qlfit <- glm(as.formula(paste("cbind(", y.var, ", 1-", y.var, ") ~ ",
                                            paste(x.vars.floor, collapse=" + "))), weights = dtmpF$w,
                           family = binomial(logit), start = coef.qlfit.start, data = dtmpF,
                           control = glm.control(maxit = maxIter))
            } else if(floor.fit=="bayesglm") {
              if (intercept.only.floor == F) {
                qlfit <- bayesglm.internal(as.formula(paste("cbind(", y.var, ", 1-", y.var, ") ~ ",
                                                   paste(x.vars.floor, collapse=" + "))),
                                  weights = dtmpF$w, family = binomial(logit),
                                  start = coef.qlfit.start, data = dtmpF,
                                  control = glm.control(maxit = maxIter), scaled = F)
              } else {
                qlfit <- bayesglm.internal(as.formula(paste("cbind(", y.var, ", 1-", y.var, ") ~ 1")),
                                  weights = dtmpF$w, family = binomial(logit),
                                  start = coef.qlfit.start, data = dtmpF,
                                  control = glm.control(maxit = maxIter), scaled = F)
              }
            }
          }
          
          ## reset parameters based on floor / ceiling options
          
          if(floor==TRUE & ceiling==TRUE) {
            par <- c(coef(fit), coef(lfit),coef(qlfit), coef(qufit))
          } else if(floor==FALSE & ceiling==TRUE) {
            par <- c(coef(fit), coef(lfit),  c(-Inf, rep(0, (nPar.floor-1))), coef(qufit))
          } else if(floor==TRUE & ceiling==FALSE) {
            par <- c(coef(fit), coef(lfit),coef(qlfit),  c(-Inf, rep(0, (nPar.ceiling-1))))
          }
          
          pllik <- llik
          
          if(verbose==T) cat(paste(counter, round(llik, 4), "\n"))
          
          llik <- obs.llik.binom.boundary(par, J = J, y = y.all, treat = t, x = x.all,
                                          x.ceiling = x.ceiling.all, x.floor = x.floor.all,
                                          ceiling.fit = ceiling.fit, floor.fit = floor.fit,
                                          ceiling = ceiling, floor = floor)
          			  
          counter <- counter + 1
          if (llik < pllik)
            warning("log-likelihood is not monotonically increasing.")
          ## NOTE CHANGED FROM STOP()
          
          if(counter == (maxIter-1))
            warning("number of iterations exceeded maximum in ML")
				  				  
        }
          
        if(floor==FALSE & ceiling==TRUE) {
          par <- par[c((1:(2*nPar)), (((2*nPar + nPar.floor)+1):(2*nPar + nPar.floor + nPar.ceiling)))]
        } else if(floor==TRUE & ceiling==FALSE) {
          par <- par[1:(2*nPar + nPar.floor)]
        }
        
        MLEfit <- optim(par, obs.llik.binom.optim.boundary, method = "BFGS", J = J, y = y.all,
                        treat = t,  x = x.all, x.ceiling = x.ceiling.all, x.floor = x.floor.all,
                        ceiling.fit = ceiling.fit,
                        floor.fit  = floor.fit, floor = floor, ceiling = ceiling,
                        hessian = TRUE, control = list(maxit = 0))
        
        if(ceiling.fit=="bayesglm" & ceiling==TRUE) {
          p.prior.ceiling <- sum(dcauchy(x = coef(qufit),
                                         scale = rep(2.5, length(coef(qufit))), log = TRUE))
        } else {
          p.prior.ceiling <- 0
        }
        
        if(floor.fit=="bayesglm" & floor==TRUE) {
          p.prior.floor <- sum(dcauchy(x = coef(qlfit),
                                       scale = rep(2.5, length(coef(qlfit))), log = TRUE))
        } else {
          p.prior.floor <- 0
        }
        
        llik <- llik - p.prior.floor - p.prior.ceiling
        
        vcov.mle <- solve(-MLEfit$hessian, tol = 1e-20)
        se.mle <- sqrt(diag(vcov.mle))
        
      } # end of boundary ml
      
    } # end of em algorithm
    
  } else if (design=="modified") {

    ## end of standard design
    
    ##
    ## Run two-step estimator for modified design
    
    ## fit to the control group
    
    fit.control <- list()
    par.control.list <- list()
    pi.pred <- matrix(NA, nrow(x.treatment), J)
    
    for (j in 1:J) {
      
      fit.glm.control <- glm(y.control[,j] ~ x.control - 1, family = binomial(logit), weights = w.control)
      coef.glm.control <- coef(fit.glm.control)
      names(coef.glm.control) <- paste("beta", 1:length(coef.glm.control), sep = "")
      
      if (fit.nonsensitive == "glm") {
        fit.control[[j]] <- fit.glm.control
      } else if (fit.nonsensitive == "nls") {

        tmpY <- y.control[, j]
        
        fit.control[[j]] <- nls( as.formula(paste("tmpY ~ logistic( x.control %*% c(",
                                                  paste(paste("beta", 1:length(coef.glm.control),
                                                              sep=""), collapse= ","), "))")) ,
                                start = coef.glm.control,
                                weights = w.control,
                                control = nls.control(maxiter=maxIter, warnOnly=TRUE))
      }
      
      par.control.list[[j]] <- coef(fit.control[[j]])
      
      pi.pred[,j] <- logistic(x.treatment %*% par.control.list[[j]])
      
    }
    
    y.treatment.pred <- y.treatment[,J+1] - apply(pi.pred, 1, sum)
    
    y.treatment.pred.temp <- ifelse(y.treatment.pred > 1, 1, y.treatment.pred)
    y.treatment.pred.temp <- ifelse(y.treatment.pred.temp < 0, 0, y.treatment.pred.temp)
    
    alpha <- mean(y.treatment.pred.temp)
    
    y.treatment.start <- ifelse(y.treatment.pred.temp >=
                                quantile(y.treatment.pred.temp, alpha), 1, 0)
    
    try(fit.glm.control <- glm(cbind(y.treatment.start, 1 - y.treatment.start) ~ x.treatment - 1, family = binomial(logit),
                               weights = w.treatment))
    try(coef.glm.control <- coef(fit.glm.control), silent = T)
    try(names(coef.glm.control) <- paste("beta", 1:length(coef.glm.control), sep = ""), silent = T)
    
    if(exists("coef.glm.control")) {
      try(fit.treat <- nls( as.formula(paste("y.treatment.pred ~ logistic( x.treatment %*% c(",
                                             paste(paste("beta", 1:length(coef.glm.control),sep=""),
                                                   collapse= ","), "))")) , start = coef.glm.control,
                           weights = w.treatment,
                           control = nls.control(maxiter=maxIter, warnOnly=TRUE)), silent = T)
    } else {
      try(fit.treat <- nls( as.formula(paste("y.treatment.pred ~ logistic( x.treatment %*% c(",
                                             paste(paste("beta", 1:length(coef.glm.control),sep=""),
                                                   collapse= ","), "))")) ,
                           weights = w.treatment,
                           control = nls.control(maxiter=maxIter, warnOnly=TRUE)), silent = T)
    }

    if(!exists("fit.treat"))
      fit.treat <- lm(y.treatment.pred ~ x.treatment - 1, weights = w.treatment)
      
    vcov.twostep.modified <- function(coef.treat, coef.control, J,
                                      x1, y1, x0, y0, fit.nonsensitive = "nls") {
      
      nPar <- length(coef.treat)
      
      n <- length(c(y1,y0))
      
      Htmp <- c(logistic(x1 %*% coef.treat) / (1 + exp(x1 %*% coef.treat))) * x1
      
      H <- - t(Htmp) %*% Htmp / n
      
      L <- matrix(NA, ncol=(J*nPar), nrow=nPar)
      
      for(j in 1:J){
        
        Ltmp <- c(sqrt(logistic(x1 %*% coef.control[[j]])*
                       logistic(x1 %*% coef.treat)/((1 + exp(x1 %*% coef.control[[j]]))*
                                                    (1 + exp(x1 %*% coef.treat))))) * x1
        
        L[,((j-1)*nPar+1):((j)*nPar)] <- -t(Ltmp) %*% Ltmp / n
        
      }
      
      M <- matrix(0, ncol = (J*nPar), nrow = (J*nPar))
      
      if (fit.nonsensitive=="glm") {
        
        for(j in 1:J){
        			
          Mtmp <- sqrt(c(exp(x0 %*% coef.control[[j]])/(1+exp(x0 %*% coef.control[[j]])) -
                         exp(2*x0 %*% coef.control[[j]])/(1+exp(x0 %*% coef.control[[j]]))^2))*
                           x0
          
          M[((j-1)*nPar+1):((j)*nPar),((j-1)*nPar+1):((j)*nPar)] <- -t(Mtmp) %*% Mtmp / n
          
        }
        
      } else if (fit.nonsensitive == "nls") {
        
        for(j in 1:J){
          
          Mtmp <- c(logistic(x0 %*% coef.control[[j]])/(1 + exp(x0 %*% coef.control[[j]]))) * x0
          
          M[((j-1)*nPar+1):((j)*nPar),((j-1)*nPar+1):((j)*nPar)] <- -t(Mtmp) %*% Mtmp / n
          
        }
        
      }
      
      invH <- solve(H, tol = 1e-20)
      invM <- solve(M, tol = 1e-20)
      
      invG <- rbind( cbind(invH, -invH %*% L %*% invM),
                    cbind(matrix(0, nrow(invM), ncol(invH)), invM)  ) 
      
      gtmp <- matrix(NA, nrow=nrow(x1), ncol=J)
      for(j in 1:J)
        gtmp[,j] <- logistic(x1 %*% coef.control[[j]])
      
      gtmpsum <- as.vector(apply(gtmp, 1, sum))
      
      g <- c((y1[,J+1] - gtmpsum - logistic(x1 %*% coef.treat))*
             logistic(x1 %*% coef.treat)/(1+exp(x1 %*% coef.treat))) * x1
      
      gg <- t(g) %*% g / n
      
      h <- list()
      
      if (fit.nonsensitive == "glm"){
        for (j in 1:J) {
          h[[j]] <- c((y0[,j] - logistic(x0 %*% coef.control[[j]]))) * x0
        }
      } else if (fit.nonsensitive == "nls") {
        for (j in 1:J) {
          h[[j]] <- c((y0[,j] - logistic(x0 %*% coef.control[[j]]))*
                      logistic(x0 %*% coef.control[[j]])/(1+exp(x0 %*% coef.control[[j]]))) * x0
        }	
      }
      
      Fh <- matrix(NA, nrow = J*nPar, ncol = J*nPar)
      for(j in 1:J) {
        for(k in 1:J){
          Fh[((j-1)*nPar+1):(j*nPar), ((k-1)*nPar+1):(k*nPar)] <- t(h[[j]]) %*% h[[k]] / n
        }
      }
      
      F <- adiag(gg, Fh)
    
      return(invG %*% F %*% t(invG) / n)
      
    }
    
    par.treat.nls.mod <- coef(fit.treat)
    par.control.nls.mod <- do.call(c, par.control.list)
    
    if(method == "nls") {
      vcov.twostep <- vcov.twostep.modified(par.treat.nls.mod, par.control.list,
                                            J, x.treatment, y.treatment, x.control, y.control,
                                            fit.nonsensitive)
      se.twostep <- sqrt(diag(vcov.twostep))
    }
    
    ##
    ## Run maximum likelihood method
    
    if(method == "ml"){

      coef.control.start <- par.control.nls.mod
      coef.treat.start <- par.treat.nls.mod
      
      R.poisbinomR <- function(k, p) {								
        colSum <- 0
        for(i in 1:k){
          if((k-i)>=1) colSum <- colSum + (-1)^(i+1) * sum( (p/(1-p))^i ) * R.poisbinomR(k - i, p)
          if((k-i)==0) colSum <- colSum + (-1)^(i+1) * sum( (p/(1-p))^i )
          if((k-i)<0) colSum <- colSum + 0						
        }
        if(k>0) return((1/k) * colSum)	
        else if(k==0) return(1)
      }
      
      ## Poisson-Binomial density function
      ## takes y scalar and p vector of Bernoulli success probabilities
      dpoisbinomR <- function(y, p) {
        
        return( R.poisbinomR(y, p) * prod(1 - p) )
        
      }
      
      dpoisbinomC <- function(y, p) {
        
        k <- length(p)
        return(.C("dpoisbinom", as.integer(y), as.integer(k), as.double(p),
                  res = double(1), PACKAGE = "list")$res)
        
      }
      
      R.poisbinomC <- function(k, p) {
        
        return(.C("RpoisbinomReturn", as.integer(k), as.double(p),
                  as.integer(length(p)), res = double(1), PACKAGE = "list")$res)
        
      }
      
      R.poisbinomE <- function(k, p) {
        
        maxk <- max(k)
        return(.C("RpoisbinomEff", as.integer(maxk), 
                  as.double(p), as.integer(length(p)), Rs = double(maxk+1),
                  PACKAGE = "list")$Rs[k+1])
        
      }
      
      dpoisbinomE <- function(y, p) {
        
        return(R.poisbinomE(y, p)*prod(1-p))
        
      }
      
      R.poisbinomM <- function(k, p) {
        
        m <- ncol(p)
        if (m != length(k))
          stop("the dimension of k match with the column dimension of p")
        return(.C("RpoisbinomEffMatrix", as.integer(k), as.integer(max(k)),
                  as.double(p), as.integer(nrow(p)), as.integer(m),
                  Rs = double(m), PACKAGE = "list")$Rs)
        
      }
      
      dpoisbinomM <- function(y, p) {
        
        return(R.poisbinomM(y, p)*apply(1-p, 2, prod))
        
      }
      
      obs.llik.modified <- function(par, y, X, treat, wt) {
        
        J <- ncol(y) - 1
        
        nPar <- length(par) / (J+1)
        
        n.treat <- nrow(y[treat==1,])
        
        x.treat <- X[treat==1,]
        y.treat <- y[treat==1,]
        
        pi <- matrix(NA, ncol = nrow(X), nrow = J+1)
        for(j in 1:(J+1))
          pi[j,] <- logistic(X %*% par[((j-1)*nPar+1):(j*nPar)])
        
        llik.treat <- sum(wt[t==1] * log(dpoisbinomM(y.treat[,J+1], pi[,t==1])))
        
        x.control <- X[treat==0,]
        y.control <- y[treat==0,]
        
        llik.control <- 0
        for(j in 1:J)
          llik.control <- llik.control + sum(wt[t==0] * (y.control[,j]*log(pi[j,t==0]) +
                                             (1-y.control[,j])*log(1-pi[j,t==0])))
        
        llik <- llik.treat + llik.control
        
        return(llik)
        
      }
      
      Estep.modified <- function(par, y, X, treat) {
        
        J <- ncol(y) - 1
        
        nPar <- length(par) / (J+1)
        
        x.treat <- X[treat==1, , drop = FALSE]
        y.treat <- y[treat==1, , drop = FALSE]
        
        n.treat <- nrow(x.treat)
        
        pi <- matrix(NA, ncol = n.treat, nrow = J+1)
        for(j in 1:(J+1))
          pi[j,] <- logistic(x.treat %*% par[((j-1)*nPar+1):(j*nPar)])
        
        all0.cond <- y.treat[,J+1]==0
        all1.cond <- y.treat[,J+1]==(J+1)
        either.cond <- all0.cond==TRUE | all1.cond==TRUE
        
        rpb <- matrix(NA, nrow = n.treat, ncol = J+1)
        rpb.inv <- matrix(NA, nrow = n.treat, ncol = J+1)
        
        rpb.inv.vector <- R.poisbinomM(y.treat[!either.cond, J+1], pi[,!either.cond])
        
        for(j in 1:(J+1)){
          rpb[!either.cond ,j] <- R.poisbinomM(y.treat[!either.cond,J+1]-1, pi[-j,!either.cond])
          rpb.inv[!either.cond,j] <- rpb.inv.vector
        }
        
        w <- t(pi) * rpb / ((1-t(pi)) * rpb.inv)
        w[all0.cond, ] <- matrix(0, ncol = ncol(w), nrow = sum(all0.cond))
        w[all1.cond, ] <- matrix(1, ncol = ncol(w), nrow = sum(all1.cond))
        
        return(w)
        
      }
      
      par <- c(coef.control.start, coef.treat.start)
      
      nPar <- length(par) / (J+1)
      
      pllik <- -Inf
      
      llik <- obs.llik.modified(par, y.all, x.all, t, wt = w.all)
      
      counter <- 0
      while (((llik - pllik) > 10^(-8)) & (counter < maxIter)) {
        
        w <- Estep.modified(par, y.all, x.all, t)
        
        y.var <- as.character(formula)[[2]]
                
        coef.list <- list()
        
        for(j in 1:(J+1)){
          
          dtmp <- as.data.frame(rbind(cbind(1, x.treatment, w[,j], 1),
                                      cbind(y.control[,j], x.control, 1, 0),
                                      cbind(0, x.treatment, 1-w[,j], 1)))

          dtmp.weights <- c(w.treatment, w.control, w.treatment)

          if (intercept.only == TRUE)
            names(dtmp) <- c("y","(Intercept)","w","treat")
          else
            names(dtmp) <- c("y","(Intercept)",x.vars,"w","treat")
          
          if(j==(J+1)) {
              dtmp <- dtmp[dtmp$treat==1,]
              dtmp.weights <- dtmp.weights[dtmp$treat==1]
          }
          
          if (j<(J+1)) {
            trials <- 1
          } else {
            trials <- j
          }
          
          fit <- glm(as.formula(paste("cbind(y, 1-y) ~ ", paste(x.vars, collapse=" + "))),
                     family = binomial(logit), weights = dtmp$w * dtmp.weights,
                     start = par[((j-1)*nPar+1):(j*nPar)], data = dtmp)
          
          coef.list[[j]] <- coef(fit)
          
        }
        
        par <- do.call(c, coef.list)
        
        pllik <- llik
        
        if(verbose==T)
          cat(paste(counter, round(llik, 4), "\n"))
        
        llik <- obs.llik.modified(par, y.all, x.all, t)
        
        counter <- counter + 1
        
        if (llik < pllik)
          warning("log-likelihood is not monotonically increasing.")
        ## NOTE CHANGED FROM STOP()
        
        if(counter == (maxIter-1))
          warning("number of iterations exceeded maximum in ML.")
        
      } # end ml while loop
      
      MLEfit <- optim(par, fn = obs.llik.modified, method = "BFGS", y = y.all, X = x.all,
                      treat = t, hessian = TRUE, control = list(maxit = 0))
      
      vcov.mle <- solve(-MLEfit$hessian, tol = 1e-20)
      se.mle <- sqrt(diag(vcov.mle))
      
    } # end ml for modified design
    
  } # end modified design
  
  # measurement error models
  if (error == "topcode" & method == "ml") {
    
    ####################
    # TOP CODING ERROR #
    ####################
    obs.llik.top <- function(all.pars, formula, data, treat, J) {

      mf <- model.frame(formula, data)
      n  <- nrow(mf)

      Y  <- model.response(mf)
      X  <- model.matrix(formula, data)
      Tr <- data[row.names(mf), treat]
      k  <- ncol(X)
      beta <- all.pars[1:k]
      gamma <- all.pars[(k+1):(2*k)]

      log.odds <- all.pars[2*k + 1]
      p0 <- logistic(log.odds)

      Xb <- X %*% beta
      Xg <- X %*% gamma

      lliks <- ifelse(Y == J + 1, log(p0 + (1 - p0) * logistic(Xb) * logistic(Xg)^J), 
        ifelse(Y == J & Tr == 0, log(p0 + (1 - p0) * logistic(Xg)^J),  
          ifelse(Y == 0 & Tr == 1, log((1 - p0) * (1 - logistic(Xb)) * (1 - logistic(Xg))^J), 
        ifelse(Y %in% 1:J & Tr == 1, log(
          (1 - p0) * (logistic(Xb) * choose(J, Y-1) * logistic(Xg)^(Y-1) * (1 - logistic(Xg))^(J-Y+1) + 
            (1 - logistic(Xb)) * choose(J, Y) * logistic(Xg)^Y * (1 - logistic(Xg))^(J-Y))), 
          log((1 - p0) * choose(J, Y) * logistic(Xg)^Y * (1 - logistic(Xg))^(J-Y))))))
     
     sum(lliks) 

    }

    topcodeE <- function(formula, data, treat, J, p0, params) {
      
      mf <- model.frame(formula, data)
      n  <- nrow(mf)

      Y  <- model.response(mf)
      X  <- model.matrix(formula, data)
      Tr <- data[row.names(mf), treat]
      k  <- ncol(X)
      beta <- params[1:k]
      gamma <- params[(k+1):(2*k)]

      Xb <- X %*% beta
      Xg <- X %*% gamma

      # probability of top coding error
      xi <- ifelse(Y == J + 1, 
        p0/(p0 * (1 - logistic(Xb) * logistic(Xg)^J) + logistic(Xb)*logistic(Xg)^J), 
          ifelse(Tr == 0 & Y == J, 
            p0/(p0 * (1 - logistic(Xg)^J) + logistic(Xg)^J), 
              0
            )
        )

      # probability of sensitive trait
      eta <- ifelse(Tr == 1 & Y == 0, 0,
        ifelse(Y == J + 1, 
          ((1-p0) * logistic(Xb) * logistic(Xg)^J + p0 * logistic(Xb))/
            (p0 + (1-p0)*logistic(Xb)*logistic(Xg)^J),
          logistic(Xb) * choose(J, Y-1) * logistic(Xg)^(Y-1) * (1 - logistic(Xg))^(J-Y+1)/ 
            (logistic(Xb) * choose(J, Y-1) * logistic(Xg)^(Y-1) * (1 - logistic(Xg))^(J-Y+1) + 
              (1 - logistic(Xb)) * choose(J, Y) * logistic(Xg)^Y * (1 - logistic(Xg))^(J - Y))
          )
        )

      # eta <- ifelse(Tr == 0, 0, eta)

      # latent number of control items
      yzeta0 <- ifelse(Tr == 0 & Y == J, 
        J * logistic(Xg)^J/(p0 * (1-logistic(Xg)^J) + logistic(Xg)^J), 0)

      yzeta0y <- rep(0, n)
      for (y in 1:(J - 1)) {
        tmp <- ifelse(Tr == 0 & Y == J, 
          y * p0 * choose(J, y) * logistic(Xg)^y * (1 - logistic(Xg))^(J - y)/
            (p0 * (1 - logistic(Xg)^J) + logistic(Xg)^J), 
              ifelse(Tr == 0 & Y == y, y, 0))
        yzeta0y <- yzeta0y + tmp
      }

      yzeta1J <- ifelse(Tr == 1 & Y == J + 1, 
        J * ((1-p0) * logistic(Xb) * logistic(Xg)^J + p0*logistic(Xg)^J)/(p0 + (1-p0)*logistic(Xb)*logistic(Xg)^J), 
          ifelse(Tr == 1 & Y == J, 
        J * (1 - logistic(Xb)) * logistic(Xg)^J/
          ((1 - logistic(Xb)) * logistic(Xg)^J + logistic(Xb) * J * logistic(Xg)^(J-1) * (1 - logistic(Xg))), 
            0))

      yzeta1y <- rep(0, n)
      for (y in 1:(J - 1)) {
        tmp <- ifelse(Tr == 1 & Y == y + 1, 
          y * logistic(Xb)*choose(J, Y-1)*logistic(Xg)^(Y-1)*(1-logistic(Xg))^(J-Y+1)/
          (logistic(Xb)*choose(J, Y-1)*logistic(Xg)^(Y-1)*(1-logistic(Xg))^(J-Y+1) + 
            (1-logistic(Xb))*choose(J, Y)*logistic(Xg)^Y*(1-logistic(Xg))^(J-Y)),
          ifelse(Tr == 1 & Y == y, 
            y * (1-logistic(Xb))*choose(J, Y)*logistic(Xg)^Y*(1-logistic(Xg))^(J-Y)/
          ((1-logistic(Xb))*choose(J, Y)*logistic(Xg)^Y*(1-logistic(Xg))^(J-Y) + 
            logistic(Xb)*choose(J, Y-1)*logistic(Xg)^(Y-1)*(1-logistic(Xg))^(J-Y+1)), 
          ifelse(Tr == 1 & Y == J + 1, 
            (y * p0 * choose(J, y) * logistic(Xg)^y * (1 - logistic(Xg))^(J - y))/
              (p0 + (1-p0)*logistic(Xb)*logistic(Xg)^J), 
                0)))
        yzeta1y <- yzeta1y + tmp

      }

      yzeta <- yzeta0 + yzeta0y + yzeta1J + yzeta1y

      return(list(xi = xi, eta = eta, yzeta = yzeta))
      
    }

    topcodeM <- function(formula, data, treat, J, xi, eta, yzeta) {

      mf <- model.frame(formula, data)
      n  <- nrow(mf)

      Y  <- model.response(mf)
      X  <- model.matrix(formula, data)
      Tr <- data[row.names(mf), treat]
      k  <- ncol(X)

      if (length(attr(attr(mf, "terms"), "term.labels")) == 0) {

        betaX <- rbind(X[Tr == 1, ], X[Tr == 1, ])
        betaY <- rep(0:1, each = sum(Tr))

        beta.fit <- glm(cbind(betaY, 1 - betaY) ~ 1, weights = c(1 - eta[Tr == 1], eta[Tr == 1]), 
          family = binomial("logit"), control = list(maxit = 5000))

        gammaX <- rbind(X, X)
        gammaY <- rep(0:1, each = n)

        gamma.fit <- glm(cbind(gammaY, 1 - gammaY) ~ 1, weights = c(J - yzeta, yzeta), 
          family = binomial("logit"), control = list(maxit = 5000))

      } else {

        betaX <- rbind(X[Tr == 1, ], X[Tr == 1, ])
        betaY <- rep(0:1, each = sum(Tr))

        beta.fit <- glm(cbind(betaY, 1 - betaY) ~ betaX - 1, weights = c(1 - eta[Tr == 1], eta[Tr == 1]), 
          family = binomial("logit"), control = list(maxit = 5000))

        gammaX <- rbind(X, X)
        gammaY <- rep(0:1, each = n)

        gamma.fit <- glm(cbind(gammaY, 1 - gammaY) ~ gammaX - 1, weights = c(J - yzeta, yzeta), 
          family = binomial("logit"), control = list(maxit = 5000))
      }

      return(list(beta.fit = beta.fit, gamma.fit = gamma.fit, 
        pars = c(coef(beta.fit), coef(gamma.fit)), p0 = mean(xi)))

    }

    topcodeEM <- function(formula, data, treat, J, p0 = 0.05, params = NULL, eps = 1e-08, maxIter) {

      mf <- model.frame(formula, data)
      n  <- nrow(mf)

      Y  <- model.response(mf)
      X  <- model.matrix(formula, data)
      Tr <- data[row.names(mf), treat]
      k  <- ncol(X)

      if (is.null(params)) {
        params <- c(coef.treat.start, coef.control.start)
      }

      beta  <- params[1:k]
      gamma <- params[(k+1):(2*k)]

      Xb <- X %*% beta
      Xg <- X %*% gamma

      Estep0 <- topcodeE(formula, data, treat, J, p0, params)
      Mstep0 <- topcodeM(formula, data, treat, J, 
        eta = Estep0$eta, yzeta = Estep0$yzeta, xi = Estep0$xi)

      Estep <- topcodeE(formula, data, treat, J, p0, params)
      Mstep <- topcodeM(formula, data, treat, J, 
        eta = Estep$eta, yzeta = Estep$yzeta, xi = Estep$xi)

      par.holder <- c(Mstep$pars, Mstep$p0)

      iteration <- 1

      while ((iteration == 1 | max(abs(par.holder - c(Mstep$par, Mstep$p0))) > eps) & 
        iteration <= maxIter) {

        par.holder <- c(Mstep$pars, Mstep$p0)

        obs.llik0 <- obs.llik.top(all.pars = c(Mstep$pars, log(Mstep$p0/(1-Mstep$p0))), 
          formula, data, treat, J)

        Estep <- topcodeE(formula, data, treat, J, p0 = Mstep$p0, params = Mstep$pars)
        Mstep <- topcodeM(formula, data, treat, J, 
          eta = Estep$eta, yzeta = Estep$yzeta, xi = Estep$xi)

        obs.llik1 <- obs.llik.top(all.pars = c(Mstep$pars, log(Mstep$p0/(1-Mstep$p0))), 
          formula, data, treat, J)

        if(signif(obs.llik1) < signif(obs.llik0)) warning("Observed-data likelihood is not monotonically increasing.")
        iteration <- iteration + 1

      }

      if (iteration == maxIter) warning("Maximum number of iterations reached.")

      optim.out <- optim(fn = obs.llik.top, hessian = TRUE, 
        par = c(Mstep$pars, log(Mstep$p0/(1 - Mstep$p0))), 
        formula = formula, data = data, treat = treat, J = J, 
        control = list(fnscale = -1))

      vcov <- vcov.flip <- solve(-optim.out$hessian)
      std.errors <- sqrt(diag(vcov))
      vcov.flip[(1:k), (1:k)] <- vcov[(k+1):(2*k), (k+1):(2*k)] 
      vcov.flip[(k+1):(2*k), (k+1):(2*k)] <- vcov[(1:k), (1:k)]

      return(
        list(
          par.treat   = Mstep$par[1:k], 
          par.control = Mstep$par[(k+1):(2*k)],
          se.treat    = std.errors[1:k],
          se.control  = std.errors[(k+1):(2*k)],
          vcov        = vcov.flip, 
          p.est       = mean(Estep$xi),
          p.ci        = c("lwr" = logistic(optim.out$par[2*k+1] - qnorm(0.975)*std.errors[2*k+1]),
            "upr" = logistic(optim.out$par[2*k+1] + qnorm(0.975)*std.errors[2*k+1])), 
          lp.est      = optim.out$par[2*k+1], 
          lp.se       = std.errors[2*k+1], 
          iterations  = iteration, 
          obs.llik    = optim.out$value
          )
        )

    }

    topcode.out <- topcodeEM(formula = formula, data = data, treat = treat, J = J, maxIter = maxIter)
    topcode.par.treat <- topcode.out$par.treat
    topcode.par.control <- topcode.out$par.control
    topcode.vcov <- topcode.out$vcov
    topcode.se.treat <- topcode.out$se.treat
    topcode.se.control <- topcode.out$se.control
    llik.topcode <- topcode.out$obs.llik   
    topcode.p <- topcode.out$p.est
    topcode.p.ci <- topcode.out$p.ci
    topcode.lp <- topcode.out$lp.est
    topcode.lp.se <- topcode.out$lp.se
    topcode.iter <- topcode.out$iterations


  } 

  if (error == "uniform" & method == "ml") {

    ########################
    # UNIFORM CODING ERROR #
    ########################
    obs.llik.uniform <- function(all.pars, formula, data, treat, J) {

      mf <- model.frame(formula, data)
      n  <- nrow(mf)

      Y  <- model.response(mf)
      X  <- model.matrix(formula, data)
      Tr <- data[row.names(mf), treat]
      k  <- ncol(X)

      beta <- all.pars[1:k]
      gamma <- all.pars[(k+1):(2*k)]

      Xb <- X %*% beta
      Xg <- X %*% gamma

      p0 <- logistic(all.pars[2*k + 1])
      p1 <- logistic(all.pars[2*k + 2])

      lliks <- ifelse(Y == J + 1, 
        log((1 - p1) * logistic(Xb) * logistic(Xg)^J + p1/(J + 2)), 
        ifelse(Y == 0 & Tr == 1, 
          log((1 - p1) * (1 - logistic(Xb)) * (1 - logistic(Xg))^J + p1/(J + 2)), 
        ifelse(Y %in% 1:J & Tr == 1, 
          log((1 - p1) * (logistic(Xb) * choose(J, Y-1) * logistic(Xg)^(Y-1) * (1 - logistic(Xg))^(J-Y+1) + 
            (1 - logistic(Xb)) * choose(J, Y) * logistic(Xg)^Y * (1 - logistic(Xg))^(J-Y)) + p1/(J + 2)), 
        log((1 - p0) * choose(J, Y) * logistic(Xg)^Y * (1 - logistic(Xg))^(J-Y) + p0/(J+1)))))

      sum(lliks)

    }

    uniformE <- function(formula, data, treat, J, p0, p1, params) {

      mf <- model.frame(formula, data)
      n  <- nrow(mf)

      Y  <- model.response(mf)
      X  <- model.matrix(formula, data)
      Tr <- data[row.names(mf), treat]
      k  <- ncol(X)

      beta  <- params[1:k]
      gamma <- params[(k+1):(2*k)]

      Xb <- X %*% beta
      Xg <- X %*% gamma

      # probability of error
      xi <- ifelse(Y == J + 1, p1/(J+2)/(p1/(J+2) + (1-p1) * logistic(Xb) * logistic(Xg)^J), 
        ifelse(Tr == 1 & Y == 0, p1/(J+2)/(p1/(J+2) + (1-p1) * (1 - logistic(Xb)) * (1 - logistic(Xg))^J), 
          ifelse(Tr == 1 & Y %in% 1:J, 
            p1/(J+2)/(p1/(J+2) + (1-p1) * (logistic(Xb) * choose(J, Y-1) * logistic(Xg)^(Y-1) * 
              (1 - logistic(Xg))^(J-Y+1) + (1 - logistic(Xb)) * choose(J, Y) * logistic(Xg)^Y * 
                (1 - logistic(Xg))^(J-Y))), 
            p0/(J+1)/(p0/(J+1) + (1-p0) * choose(J, Y) * logistic(Xg)^Y * (1 - logistic(Xg))^(J-Y)))))

      # probability of sensitive trait
      eta <- ifelse(Y == J + 1, 
         (p1/(J + 2) * logistic(Xb) + (1 - p1) * logistic(Xb) * logistic(Xg)^J)/
          ((1 - p1) * logistic(Xb) * logistic(Xg)^J + p1/(J + 2)), 
        ifelse(Tr == 1 & Y == 0, 
          p1/(J + 2) * logistic(Xb)/
           ((1 - p1) * (1 - logistic(Xb)) * (1 - logistic(Xg))^J + p1/(J + 2)), 
        ifelse(Tr == 1 & Y %in% 1:J, 
          (p1/(J + 2) + (1 - p1) * choose(J, Y - 1) * logistic(Xg)^(Y - 1) * (1 - logistic(Xg))^(J - Y + 1)) * logistic(Xb)/
            (p1/(J + 2) + (1 - p1) * (logistic(Xb) * choose(J, Y - 1) * logistic(Xg)^(Y - 1) * (1 - logistic(Xg))^(J - Y + 1) + 
              (1 - logistic(Xb)) * choose(J, Y) * logistic(Xg)^Y * (1 - logistic(Xg))^(J - Y))), 
            0)
          )
        )

      # expected values for control obs
      yzeta0 <- Y * (p0/(J + 1) + (1 - p0)) * choose(J, Y) * logistic(Xg)^Y * (1 - logistic(Xg))^(J - Y)/
        (p0/(J + 1) + (1 - p0) * choose(J, Y) * logistic(Xg)^Y * (1 - logistic(Xg))^(J - Y))
      for (y in 0:J) {
        yzeta0 <- yzeta0 + ifelse(Y == y, 0, y * 
          p0/(J + 1) * choose(J, y) * logistic(Xg)^y * (1 - logistic(Xg))^(J - y)/
          (p0/(J + 1) + (1 - p0) * choose(J, Y) * logistic(Xg)^Y * (1 - logistic(Xg))^(J - Y)))
      }

      # expected values for treated obs
      yzetaJ1 <- J * ifelse(Y == J + 1, ((1-p1) * logistic(Xb) * logistic(Xg)^J + p1/(J+2) * logistic(Xg)^J)/
        (p1/(J+2) + (1-p1) * logistic(Xb) * logistic(Xg)^J), 
          ifelse(Tr == 1 & Y == J, ((1-p1) * (1 - logistic(Xb)) * logistic(Xg)^J + p1/(J+2) * logistic(Xg)^J)/
            (p1/(J+2) + (1-p1) * (logistic(Xb) * J * logistic(Xg)^(J-1) * (1 - logistic(Xg)) + 
              (1 - logistic(Xb)) * logistic(Xg)^J)), 
          ifelse(Tr == 1 & Y == 0, (p1/(J+2) * logistic(Xg)^J)/ 
            (p1/(J+2) + (1-p1) * (1 - logistic(Xb)) * (1 - logistic(Xg))^J), 
            ifelse(Tr == 1 & Y %in% 1:(J-1), p1/(J+2) * logistic(Xg)^J/
        (p1/(J+2) + (1-p1) * ((1 - logistic(Xb)) * choose(J, Y) * logistic(Xg)^Y * (1 - logistic(Xg))^(J - Y) + 
          logistic(Xb) * choose(J, Y-1) * logistic(Xg)^(Y-1) * (1 - logistic(Xg))^(J-Y+1))), 
              0))))

      yzetay1 <- rep(0, n) 
      for (y in 1:(J - 1)) {
        tmp <- ifelse(Tr == 1 & Y == y + 1, 
          y * (p1/(J + 2) + (1 - p1) * logistic(Xb)) * choose(J, y) * logistic(Xg)^y * (1 - logistic(Xg))^(J - y)/
              (p1/(J + 2) + (1 - p1) * (logistic(Xb) * choose(J, y) * logistic(Xg)^y * (1 - logistic(Xg))^(J - y) + 
                (1 - logistic(Xb)) * choose(J, y + 1) * logistic(Xg)^(y + 1) * (1 - logistic(Xg))^(J - y - 1))), 
        ifelse(Tr == 1 & Y == y, 
          y * (p1/(J + 2) + (1 - p1) * (1 - logistic(Xb))) * choose(J, y) * logistic(Xg)^y * (1 - logistic(Xg))^(J - y)/
              (p1/(J + 2) + (1 - p1) * (logistic(Xb) * choose(J, y - 1) * logistic(Xg)^(y - 1) * (1 - logistic(Xg))^(J - y + 1) + 
                (1 - logistic(Xb)) * choose(J, y) * logistic(Xg)^y * (1 - logistic(Xg))^(J - y))), 
        ifelse(Tr == 1 & Y == 0, 
          y * p1/(J + 2) * choose(J, y) * logistic(Xg)^y * (1 - logistic(Xg))^(J - y)/
              (p1/(J + 2) + (1 - p1) * (1 - logistic(Xb)) * (1 - logistic(Xg))^J), 
        ifelse(Tr == 0, 0, 
          y * p1/(J + 2) * choose(J, y) * logistic(Xg)^y * (1 - logistic(Xg))^(J - y)/
              (p1/(J + 2) + (1 - p1) * (
                (1 - logistic(Xb)) * choose(J, Y) * logistic(Xg)^Y * (1 - logistic(Xg))^(J - Y) + 
                logistic(Xb) * choose(J, Y - 1) * logistic(Xg)^(Y - 1) * (1 - logistic(Xg))^(J - Y + 1)
          ))))))
        yzetay1 <- yzetay1 + tmp
      }

      yzeta1 <- yzetaJ1 + yzetay1

      yzeta <- ifelse(Tr == 0, yzeta0, yzeta1)

      return(list(xi = xi, eta = eta, yzeta = yzeta))

    }



    uniformM <- function(formula, data, treat, J, xi, eta, yzeta) {

      mf <- model.frame(formula, data)
      n  <- nrow(mf)

      Y  <- model.response(mf)
      X  <- model.matrix(formula, data)
      Tr <- data[row.names(mf), treat]
      k  <- ncol(X)

      # intercept-only models
      if (length(attr(attr(mf, "terms"), "term.labels")) == 0) {

        betaX <- rbind(X[Tr == 1, ], X[Tr == 1, ])
        betaY <- rep(0:1, each = sum(Tr))

        beta.fit <- glm(cbind(betaY, 1 - betaY) ~ 1, weights = c(1 - eta[Tr == 1], eta[Tr == 1]), 
          family = binomial("logit"), control = list(maxit = 100))

        gammaX <- rbind(X, X)
        gammaY <- rep(0:1, each = n)

        gamma.fit <- glm(cbind(gammaY, 1 - gammaY) ~ 1, weights = c(J - yzeta, yzeta), 
          family = binomial("logit"), control = list(maxit = 100))

      # models with covariates
      } else {

        betaX <- rbind(X[Tr == 1, ], X[Tr == 1, ])
        betaY <- rep(0:1, each = sum(Tr))

        beta.fit <- glm(cbind(betaY, 1 - betaY) ~ betaX - 1, weights = c(1 - eta[Tr == 1], eta[Tr == 1]), 
          family = binomial("logit"), control = list(maxit = 100))

        gammaX <- rbind(X, X)
        gammaY <- rep(0:1, each = n)

        gamma.fit <- glm(cbind(gammaY, 1 - gammaY) ~ gammaX - 1, weights = c(J - yzeta, yzeta), 
          family = binomial("logit"), control = list(maxit = 100))
      }

      return(list(beta.fit = beta.fit, gamma.fit = gamma.fit, 
        pars = c(coef(beta.fit), coef(gamma.fit)), p0 = mean(xi[Tr==0]), p1 = mean(xi[Tr==1])))

    }


    uniformEM <- function(formula, data, treat, J, p0 = 0.05, p1 = 0.05, params = NULL, eps = 1e-08, 
      maxIter) {

      mf <- model.frame(formula, data)
      n  <- nrow(mf)

      Y  <- model.response(mf)
      X  <- model.matrix(formula, data)
      Tr <- data[row.names(mf), treat]
      k  <- ncol(X)

      if (is.null(params)) {
        params <- c(coef.treat.start, coef.control.start)
      }

      beta  <- params[1:k]
      gamma <- params[(k+1):(2*k)]

      Xb <- X %*% beta
      Xg <- X %*% gamma

      Estep0 <- uniformE(formula, data, treat, J, p0, p1, params)
      Mstep0 <- uniformM(formula, data, treat, J, 
        eta = Estep0$eta, yzeta = Estep0$yzeta, xi = Estep0$xi)

      Estep <- uniformE(formula, data, treat, J, p0, p1, params)
      Mstep <- uniformM(formula, data, treat, J, 
        eta = Estep$eta, yzeta = Estep$yzeta, xi = Estep$xi)

      par.holder <- c(Mstep$pars, Mstep$p0, Mstep$p1)

      iteration <- 1

      while ((iteration == 1 | max(abs(par.holder - c(Mstep$par, Mstep$p0, Mstep$p1))) > eps) & 
        iteration <= maxIter) {

        par.holder <- c(Mstep$pars, Mstep$p0, Mstep$p1)

        obs.llik0 <- obs.llik.uniform(
          all.pars = c(Mstep$pars, log(Mstep$p0/(1-Mstep$p0)), log(Mstep$p1/(1-Mstep$p1))), 
          formula, data, treat, J)

        Estep <- uniformE(formula, data, treat, J, 
          p0 = Mstep$p0, p1 = Mstep$p1, params = Mstep$pars)
        Mstep <- uniformM(formula, data, treat, J, 
          eta = Estep$eta, yzeta = Estep$yzeta, xi = Estep$xi)

        obs.llik1 <- obs.llik.uniform(
          all.pars = c(Mstep$pars, log(Mstep$p0/(1-Mstep$p0)), log(Mstep$p1/(1-Mstep$p1))), 
          formula, data, treat, J)

        if (signif(obs.llik1) < signif(obs.llik0)) warning("Observed-data likelihood is not monotonically increasing.")
        iteration <- iteration + 1

      }

      if (iteration == maxIter) warning("Maximum number of iterations reached.")

      optim.out <- optim(fn = obs.llik.uniform, hessian = TRUE, 
        par = c(Mstep$pars, log(Mstep$p0/(1-Mstep$p0)), log(Mstep$p1/(1-Mstep$p1))), 
        formula = formula, data = data, treat = treat, J = J, 
        control = list(fnscale = -1))

      vcov <- vcov.flip <- solve(-optim.out$hessian)
      std.errors <- sqrt(diag(vcov))
      vcov.flip[(1:k), (1:k)] <- vcov[(k+1):(2*k), (k+1):(2*k)] 
      vcov.flip[(k+1):(2*k), (k+1):(2*k)] <- vcov[(1:k), (1:k)]

      return(
        list(
          par.treat   = Mstep$par[1:k], 
          par.control = Mstep$par[(k+1):(2*k)],
          se.treat    = std.errors[1:k],
          se.control  = std.errors[(k+1):(2*k)],
          vcov        = vcov.flip, 
          p0.est      = Mstep$p0,
          p0.ci       = c("lwr" = logistic(optim.out$par[2*k+1] - qnorm(0.975)*std.errors[2*k+1]),
            "upr" = logistic(optim.out$par[2*k+1] + qnorm(0.975)*std.errors[2*k+1])), 
          p1.est      = Mstep$p1,
          p1.ci       = c("lwr" = logistic(optim.out$par[2*k+2] - qnorm(0.975)*std.errors[2*k+2]),
            "upr" = logistic(optim.out$par[2*k+2] + qnorm(0.975)*std.errors[2*k+2])), 
          lp0.est     = optim.out$par[2*k+1], 
          lp0.se      = std.errors[2*k+1], 
          lp1.est     = optim.out$par[2*k+2], 
          lp1.se      = std.errors[2*k+2], 
          iterations  = iteration, 
          obs.llik    = optim.out$value
          )
        )

    }

    uniform.out <- uniformEM(formula = formula, data = data, treat = treat, J = J, maxIter = maxIter)
    uniform.par.treat <- uniform.out$par.treat
    uniform.par.control <- uniform.out$par.control
    uniform.vcov <- uniform.out$vcov
    uniform.se.treat <- uniform.out$se.treat
    uniform.se.control <- uniform.out$se.control
    llik.uniform <- uniform.out$obs.llik   
    uniform.p0 <- uniform.out$p0.est
    uniform.p0.ci <- uniform.out$p0.ci
    uniform.lp0 <- uniform.out$lp0.est
    uniform.lp0.se <- uniform.out$lp0.se
    uniform.p1 <- uniform.out$p1.est
    uniform.p1.ci <- uniform.out$p1.ci
    uniform.lp1 <- uniform.out$lp1.est
    uniform.lp1.se <- uniform.out$lp1.se
    uniform.iter <- uniform.out$iterations

  }
  # end of measurement error models 

  # robust ml functionality 

  if (robust) { 

    sweepCross <- function(M, x) t(M) %*% sweep(M, MAR = 1, STATS = x, "*")

    loglik <- function(params, J, Y, T, X) {
      
      n     <- length(Y)
      K     <- length(params)
      gamma <- params[1:(K/2)]
      beta  <- params[(K/2 + 1):K]

      Xb    <- X %*% beta
      Xg    <- X %*% gamma
     
      llik <- c(
         -J * log(1 + exp(Xg)) - log(1 + exp(Xb)) + 
          ifelse(Y == J + 1, Xb + J * Xg, 0) +  
          ifelse(T == 0, log(choose(J, Y)) + Y * Xg + log(1 + exp(Xb)), 0) + 
          ifelse(T == 1 & Y %in% 1:J, (Y - 1) * Xg + log(choose(J, Y - 1) * exp(Xb) + choose(J, Y) * exp(Xg)), 0)
        )

      return(sum(llik))

    }

    MLGMM <- function(params, cW, J, Y, T, X, robust) {
      
      n     <- length(Y)
      K     <- length(params)
      beta  <- params[1:(K/2)]
      gamma <- params[(K/2 + 1):K]

      Xb    <- X %*% beta
      Xg    <- X %*% gamma
      Xbg   <- X %*% (beta + gamma)

      # identification of beta
      beta.coef <- c(
        -logistic(Xb) + 
        ifelse(Y == J + 1, 1, 0) + 
        ifelse(T == 0, logistic(Xb), 0) + 
        ifelse(T == 1 & Y %in% 1:J, choose(J, Y - 1) * exp(Xb)/(choose(J, Y - 1) * exp(Xb) + choose(J, Y) * exp(Xg)), 0)
       )

      beta.foc <- colMeans(X*beta.coef)

      # identification of gamma
      gamma.coef <- c(
        -J*logistic(Xg) + 
        ifelse(Y == J + 1, J, 0) + 
        ifelse(T == 0, Y, 0) + 
        ifelse(T == 1 & Y %in% 1:J, (Y - 1) + choose(J, Y) * exp(Xg)/(choose(J, Y - 1) * exp(Xb) + choose(J, Y) * exp(Xg)), 0)
      )
      
      gamma.foc <- colMeans(X*gamma.coef)

      # dim moment
      if (robust == TRUE) {
        aux.vec <- (logistic(Xb) - mean(Y[T == 1]) + mean(Y[T == 0]))
        aux.mom <- mean(aux.vec)
        cG <- c(beta.foc, gamma.foc, aux.mom) 
      } else {
        cG <- c(beta.foc, gamma.foc)
      }

      # gmm objective
      gmm.objective <- as.numeric(t(cG) %*% ginv(cW) %*% cG)

      return(gmm.objective)

    }

    MLGMM.Grad <- function(params, cW, J, Y, T, X, robust) {
      
      n     <- length(Y)
      K     <- length(params)
      beta  <- params[1:(K/2)]
      gamma <- params[(K/2 + 1):K]

      Xb    <- X %*% beta
      Xg    <- X %*% gamma
      Xbg   <- X %*% (beta + gamma)

      # identification of beta
      beta.coef <- c(
        -logistic(Xb) + 
        ifelse(Y == J + 1, 1, 0) + 
        ifelse(T == 0, logistic(Xb), 0) + 
        ifelse(T == 1 & Y %in% 1:J, choose(J, Y - 1) * exp(Xb)/(choose(J, Y - 1) * exp(Xb) + choose(J, Y) * exp(Xg)), 0)
       )

      beta.foc <- colMeans(X*beta.coef)

      # identification of gamma
      gamma.coef <- c(
        -J*logistic(Xg) + 
        ifelse(Y == J + 1, J, 0) + 
        ifelse(T == 0, Y, 0) + 
        ifelse(T == 1 & Y %in% 1:J, (Y - 1) + choose(J, Y) * exp(Xg)/(choose(J, Y - 1) * exp(Xb) + choose(J, Y) * exp(Xg)), 0)
      )
      
      gamma.foc <- colMeans(X*gamma.coef)

      # dim moment
      if (robust == TRUE) {
        aux.vec <- (logistic(Xb) - mean(Y[T == 1]) + mean(Y[T == 0]))
        aux.mom <- mean(aux.vec)
        cG <- c(beta.foc, gamma.foc, aux.mom) 
      } else {
        cG <- c(beta.foc, gamma.foc)
      }

      # jacobian
      tmp1 <- -logistic(Xb)/(1 + exp(Xb)) + 
        ifelse(T == 0, logistic(Xb)/(1 + exp(Xb)), 0) + 
          ifelse(T == 1 & Y %in% 1:J, choose(J, Y) * choose(J, Y - 1) * exp(Xbg)/((choose(J, Y - 1) * exp(Xb) + choose(J, Y) * exp(Xg))^2), 0)
      tmp2 <- -J*logistic(Xg)/(1+exp(Xg)) + 
        ifelse(T == 1 & Y %in% 1:J, choose(J, Y) * choose(J, Y - 1) * exp(Xbg)/((choose(J, Y - 1) * exp(Xb) + choose(J, Y) * exp(Xg))^2), 0)
      tmp3 <- -ifelse(T == 1 & Y %in% 1:J, choose(J, Y) * choose(J, Y - 1) * exp(Xbg)/((choose(J, Y - 1) * exp(Xb) + choose(J, Y) * exp(Xg))^2), 0)
      tmp4 <- logistic(Xb)/(1+exp(Xb))

      if (robust == TRUE) {
        dcG <- 1/n * rbind(
          cbind(sweepCross(M = X, x = tmp1), sweepCross(M = X, x = tmp3)), 
          cbind(t(sweepCross(M = X, x = tmp3)), sweepCross(M = X, x = tmp2)), 
          cbind(t(colSums(X*c(tmp4))), matrix(0, nc = K/2, nr = 1))
        )
      } else {
        dcG <- 1/n * rbind(
          cbind(sweepCross(M = X, x = tmp1), sweepCross(M = X, x = tmp3)), 
          cbind(t(sweepCross(M = X, x = tmp3)), sweepCross(M = X, x = tmp2))
        )
      }

      return.vec <- c(t(cG) %*% (ginv(cW) + t(ginv(cW))) %*% dcG)
      
      return(return.vec)

    }

    weightMatrix <- function(params, J, Y, T, X, robust) {
      
      n     <- length(Y)
      K     <- length(params)
      beta  <- params[1:(K/2)]
      gamma <- params[(K/2 + 1):K]

      Xb    <- X %*% beta
      Xg    <- X %*% gamma
      Xbg   <- X %*% (beta + gamma)

      # identification of beta
      beta.coef <- c(
        -logistic(Xb) + 
        ifelse(Y == J + 1, 1, 0) + 
        ifelse(T == 0, logistic(Xb), 0) + 
        ifelse(T == 1 & Y %in% 1:J, choose(J, Y - 1) * exp(Xb)/(choose(J, Y - 1) * exp(Xb) + choose(J, Y) * exp(Xg)), 0)
       )

      beta.mat <- X*beta.coef

      # identification of gamma
      gamma.coef <- c(
        -J*logistic(Xg) + 
        ifelse(Y == J + 1, J, 0) + 
        ifelse(T == 0, Y, 0) + 
        ifelse(T == 1 & Y %in% 1:J, (Y - 1) + choose(J, Y) * exp(Xg)/(choose(J, Y - 1) * exp(Xb) + choose(J, Y) * exp(Xg)), 0)
      )
      
      gamma.mat <- X*gamma.coef

      # dim moment
      if (robust == TRUE) {
        aux.vec <- logistic(Xb) - mean(Y[T == 1]) + mean(Y[T == 0])
        Wtmp <- cbind(beta.mat, gamma.mat, aux.vec) # N by (2K + 1)
      } else {
        Wtmp <- cbind(beta.mat, gamma.mat)
      }

      # weight matrix
      cW <- (t(Wtmp) %*% Wtmp)/n # (2K + 1) by (2K + 1)

      return(cW)

    }

    MLGMM.var <- function(params, J, Y, T, X, robust) {
      
      n     <- length(Y)
      K     <- length(params)
      beta  <- params[1:(K/2)]
      gamma <- params[(K/2 + 1):K]

      Xb    <- X %*% beta
      Xg    <- X %*% gamma
      Xbg   <- X %*% (beta + gamma)

      # identification of beta
      beta.coef <- c(
        -logistic(Xb) + 
        ifelse(Y == J + 1, 1, 0) + 
        ifelse(T == 0, logistic(Xb), 0) + 
        ifelse(T == 1 & Y %in% 1:J, choose(J, Y - 1) * exp(Xb)/(choose(J, Y - 1) * exp(Xb) + choose(J, Y) * exp(Xg)), 0)
       )

      beta.foc <- colMeans(X*beta.coef)

      # identification of gamma
      gamma.coef <- c(
        -J*logistic(Xg) + 
        ifelse(Y == J + 1, J, 0) + 
        ifelse(T == 0, Y, 0) + 
        ifelse(T == 1 & Y %in% 1:J, (Y - 1) + choose(J, Y) * exp(Xg)/(choose(J, Y - 1) * exp(Xb) + choose(J, Y) * exp(Xg)), 0)
      )
      
      gamma.foc <- colMeans(X*gamma.coef)

      # dim moment
      if (robust == TRUE) {
        aux.vec <- (logistic(Xb) - mean(Y[T == 1]) + mean(Y[T == 0]))
        aux.mom <- mean(aux.vec)
        cG <- c(beta.foc, gamma.foc, aux.mom) 
      } else {
        cG <- c(beta.foc, gamma.foc)
      }

      # jacobian
      tmp1 <- -logistic(Xb)/(1 + exp(Xb)) + 
        ifelse(T == 0, logistic(Xb)/(1 + exp(Xb)), 0) + 
          ifelse(T == 1 & Y %in% 1:J, choose(J, Y) * choose(J, Y - 1) * exp(Xbg)/((choose(J, Y - 1) * exp(Xb) + choose(J, Y) * exp(Xg))^2), 0)
      tmp2 <- -J*logistic(Xg)/(1+exp(Xg)) + 
        ifelse(T == 1 & Y %in% 1:J, choose(J, Y) * choose(J, Y - 1) * exp(Xbg)/((choose(J, Y - 1) * exp(Xb) + choose(J, Y) * exp(Xg))^2), 0)
      tmp3 <- -ifelse(T == 1 & Y %in% 1:J, choose(J, Y) * choose(J, Y - 1) * exp(Xbg)/((choose(J, Y - 1) * exp(Xb) + choose(J, Y) * exp(Xg))^2), 0)
      tmp4 <- logistic(Xb)/(1+exp(Xb))

      if (robust == TRUE) {
        dcG <- 1/n * rbind(
          cbind(sweepCross(M = X, x = tmp1), sweepCross(M = X, x = tmp3)), 
          cbind(t(sweepCross(M = X, x = tmp3)), sweepCross(M = X, x = tmp2)), 
          cbind(t(colSums(X*c(tmp4))), matrix(0, nc = K/2, nr = 1))
        )
      } else {
        dcG <- 1/n * rbind(
          cbind(sweepCross(M = X, x = tmp1), sweepCross(M = X, x = tmp3)), 
          cbind(t(sweepCross(M = X, x = tmp3)), sweepCross(M = X, x = tmp2))
        )
      }

      cW <- weightMatrix(params, J, Y, T, X, robust)

      return.mat <- ginv(t(dcG) %*% ginv(cW) %*% dcG)

      return(return.mat)

    }

    ictrobust <- function(formula, data, treat, J, robust) {
      
      mf <- model.frame(formula, data)
      n  <- nrow(mf)

      Y  <- model.response(mf)
      X  <- model.matrix(formula, data)
      T  <- data[row.names(mf), treat]
      k  <- ncol(X)

      try(nls.fit <- ictreg(formula, data = data, treat = treat, J = J, method = "nls"))
      if (exists("nls.fit")) {
        par0   <- c(nls.fit$par.treat, nls.fit$par.control)
      } else {
        par0   <- rep(0, 2*k)
      }

      cW0   <- weightMatrix(params = par0, J = J, Y = Y, T = T, X = X, robust = robust)  
      step1 <- optim(par = par0, fn = MLGMM, gr = MLGMM.Grad, robust, 
        J = J, Y = Y, T = T, X = X, cW = cW0, 
        method = "L-BFGS-B", lower = -10, upper = 10, control = list(maxit = 5000))
      par1  <- step1$par
      cW1   <- weightMatrix(params = par1, J = J, Y = Y, T = T, X = X, robust = robust)

      step2 <- optim(par = par1, fn = MLGMM, gr = MLGMM.Grad, robust, 
        J = J, Y = Y, T = T, X = X, cW = cW1, 
        method = "L-BFGS-B", lower = -10, upper = 10, control = list(maxit = 5000))
      par2  <- step2$par
      cW2   <- weightMatrix(params = par2, J = J, Y = Y, T = T, X = X, robust = robust)

      step3 <- optim(par = par2, fn = MLGMM, gr = MLGMM.Grad, robust, 
        J = J, Y = Y, T = T, X = X, cW = cW2, 
        method = "L-BFGS-B", lower = -10, upper = 10, control = list(maxit = 5000))

      vcov  <- MLGMM.var(params = step2$par, J = J, Y = Y, T = T, X = X, robust = FALSE)/n

      return(list(
        par = step3$par, 
        vcov = vcov, 
        se = sqrt(diag(vcov)), 
        converge = 1 - step3$convergence, 
        J.stat = ifelse(robust, n*step3$value, NA),
        p.val = ifelse(robust, 1 - pchisq(q = step3$value, df = 2*k), NA),
        est.mean = mean(logistic(X %*% step3$par[1:k])), 
        diff.mean = mean(Y[T == 1]) - mean(Y[T == 0]), 
        robust = robust, 
        message = step3$message
        )
      )

    }


    ictrobust.out <- ictrobust(formula = formula, data = data, treat = treat, J = J, robust = robust)
    if (ictrobust.out$converge != 1) warning("Optimization routine did not converge.")
    ictrobust.par <- ictrobust.out$par
    ictrobust.vcov <- ictrobust.out$vcov
    ictrobust.se <- sqrt(diag(ictrobust.vcov))

    names(ictrobust.par) <- rep(names(x.all), 2)
    names(ictrobust.se)  <- rep(names(x.all), 2)

    par.control.robust <- ictrobust.par[1:ncol(x.all)]
    par.treat.robust <- ictrobust.par[(ncol(x.all)+1):(ncol(x.all)*2)]
    llik.robust <- loglik(params = ictrobust.par, J = J, Y = y.all, T = treat, X = x.all)   

  }

  # auxiliary data functionality
  if (aux.check) {

    # this function returns the weighting matrix for the GMM estimator
    invF <- function(pars, J, y, treat, x, w = NULL, h, group, matrixMethod = c("efficient", "princomp", "cue")) {

        n <- length(y)

        if (missing(w)) w <- rep(1, n)
        w <- (n/sum(w)) * w
        
        y1 <- y[treat == 1]
        y0 <- y[treat == 0]
        x1 <- x[treat == 1, , drop = FALSE]
        x0 <- x[treat == 0, , drop = FALSE]
        w1 <- w[treat == 1]
        w0 <- w[treat == 0]
        group1 <- group[treat == 1]
        delta <- pars[1:ncol(x)]
        gamma <- pars[(ncol(x) + 1):(ncol(x) * 2)]

        m1 <- c((y1 - J * logistic(x1 %*% gamma) - logistic(x1 %*% delta)) * 
          logistic(x1 %*% delta)/(1 + exp(x1 %*% delta))) * x1
        m0 <- c((y0 - J * logistic(x0 %*% gamma)) * 
          J * logistic(x0 %*% gamma)/(1 + exp(x0 %*% gamma))) * x0
        Em1 <- t(m1) %*% m1/n
        Em0 <- t(m0) %*% m0/n
        F <- adiag(Em1, Em0)

        if (length(h) > 0) {
          group.labels <- names(h)
          for (label in group.labels) {
            aux.mom <- h[label] - logistic(x[group == label, , drop = FALSE] %*% delta)
            F <- adiag(F, t(aux.mom) %*% aux.mom/n)
          }
        }
        
        efficientW <- solve(F, tol = 1e-20)

        if (matrixMethod == "efficient" | matrixMethod == "cue") {
          efficientW
        } else if (matrixMethod == "princomp") {
          k <- ncol(F)        
          decomp <- eigen(F)

          threshold <- .95 * sum(decomp$values)
          curr.sum <- r <- 0
          while (curr.sum <= threshold) {
            r <- r + 1
            curr.sum <- sum(decomp$values[1:r]) 
          }

          princompW <- matrix(0, nrow = k, ncol = k)
          for (value in 1:r) {
            princompW <- princompW + 1/decomp$values[value] * decomp$vectors[, value] %*% t(decomp$vectors[, value])
          }
          princompW
        }

    }

    # GMM objective function
    gmm.nls <- function(pars, J, y, treat, x, w = NULL, h, group, W = NULL) {

        n <- length(y)

        if (missing(w)) w <- rep(1, n)
        if (missing(W)) W <- diag(nrow = length(pars) + length(h)) 
        w <- (n/sum(w)) * w
        
        y1 <- y[treat == 1]
        y0 <- y[treat == 0]
        x1 <- x[treat == 1, , drop = FALSE]
        x0 <- x[treat == 0, , drop = FALSE]
        w1 <- w[treat == 1]
        w0 <- w[treat == 0]
        group1 <- group[treat == 1]
        delta <- pars[1:ncol(x)]
        gamma <- pars[(ncol(x) + 1):(2 * ncol(x))]

        # NLS moments
        m1 <- c((y1 - J * logistic(x1 %*% gamma) - logistic(x1 %*% delta)) * 
          logistic(x1 %*% delta)/(1 + exp(x1 %*% delta))) * x1
        m0 <- c((y0 - J * logistic(x0 %*% gamma)) * 
          J * logistic(x0 %*% gamma)/(1 + exp(x0 %*% gamma))) * x0

        # AUX moments
        g <- c()
        if (length(h) > 0) {
          group.labels <- names(h)
          for (label in group.labels) {
              g.tmp <- h[label] - logistic(x[group == label, , drop = FALSE] %*% delta)
              g <- c(g, sum(g.tmp))
          }
        }

        M <- c(colSums(m1), colSums(m0), g)/n
        as.numeric(t(M) %*% W %*% M)

    }

    # GMM gradient
    gmm.grad <- function(pars, J, y, treat, x, w = NULL, h, group, W = NULL) {

        n <- length(y)

        if (missing(w)) w <- rep(1, n)
        if (missing(W)) W <- diag(nrow = length(pars) + length(h)) 
        w <- (n/sum(w)) * w

        y1 <- y[treat == 1]
        y0 <- y[treat == 0]
        x1 <- x[treat == 1, , drop = FALSE]
        x0 <- x[treat == 0, , drop = FALSE]
        w1 <- w[treat == 1]
        w0 <- w[treat == 0]
        group1 <- group[treat == 1]
        delta <- pars[1:ncol(x)]
        gamma <- pars[(ncol(x) + 1):(2 * ncol(x))]

        # NLS moments
        m1 <- c((y1 - J * logistic(x1 %*% gamma) - logistic(x1 %*% delta)) * 
          logistic(x1 %*% delta)/(1 + exp(x1 %*% delta))) * x1
        m0 <- c((y0 - J * logistic(x0 %*% gamma)) * 
          J * logistic(x0 %*% gamma)/(1 + exp(x0 %*% gamma))) * x0

        # AUX moments
        g <- dh <- c()
        if (length(h) > 0) {
          group.labels <- names(h)
          for (label in group.labels) {
              g.tmp <- h[label] - logistic(x[group == label, , drop = FALSE] %*% delta)
              g <- c(g, sum(g.tmp))
          # AUX Jacobian
              dh.tmp <- -c(logistic(x[group == label, , drop = FALSE] %*% delta)/
                (1 + exp(x[group == label, , drop = FALSE] %*% delta))) * 
                  x[group == label, , drop = FALSE]
              dh <- rbind(dh, t(colSums(dh.tmp))/n)
          }
        } 
        # NLS Jacobian
        Gtmp <- c(logistic(x1 %*% delta)/(1 + exp(x1 %*%
            delta))) * x1
        G1 <- -t(Gtmp * w1) %*% Gtmp/n
        Gtmp <- c(sqrt(J * logistic(x1 %*% delta) * logistic(x1 %*%
            gamma)/((1 + exp(x1 %*% delta)) * (1 + exp(x1 %*%
            gamma))))) * x1
        G2 <- -t(Gtmp * w1) %*% Gtmp/n
        Gtmp <- c(J * logistic(x0 %*% gamma)/(1 + exp(x0 %*%
            gamma))) * x0
        G3 <- -t(Gtmp * w0) %*% Gtmp/n
        G <- rbind(cbind(G1, G2), cbind(matrix(0, ncol = ncol(G1), nrow = nrow(G3)), G3))
        
        # complete moment vector
        M <- c(colSums(m1), colSums(m0), g)/n
        # complete Jacobian
        if (length(h) > 0) G <- rbind(G, cbind(dh, matrix(0, nrow = length(h), ncol = ncol(G3))))

        t(M) %*% (W + t(W)) %*% G

    }

    gmm.cue <- function(pars, J, y, treat, x, w, h, group, matrixMethod = c("efficient", "princomp", "cue")) {

        n <- length(y)
        w <- (n/sum(w)) * w
        y1 <- y[treat == 1]
        y0 <- y[treat == 0]
        x1 <- x[treat == 1, , drop = FALSE]
        x0 <- x[treat == 0, , drop = FALSE]
        w1 <- w[treat == 1]
        w0 <- w[treat == 0]
        group1 <- group[treat == 1]
        delta <- pars[1:ncol(x)]
        gamma <- pars[(ncol(x) + 1):(2 * ncol(x))]

        # NLS moments
        m1 <- c((y1 - J * logistic(x1 %*% gamma) - logistic(x1 %*% delta)) * 
          logistic(x1 %*% delta)/(1 + exp(x1 %*% delta))) * x1
        m0 <- c((y0 - J * logistic(x0 %*% gamma)) * 
          J * logistic(x0 %*% gamma)/(1 + exp(x0 %*% gamma))) * x0

        # AUX moments
        g <- c()
        group.labels <- names(h)
        for (label in group.labels) {
            g.tmp <- c(h[label]) - logistic(x[group == label, , drop = FALSE] %*% delta)
            g <- c(g, sum(g.tmp * w[group == label]))
        }
        M <- c(colSums(m1), colSums(m0), g)/n

        Em1 <- t(m1) %*% m1/n
        Em0 <- t(m0) %*% m0/n
        F <- adiag(Em1, Em0)

        for (label in group.labels) {
          aux.mom <- h[label] - logistic(x[group == label, , drop = FALSE] %*% delta)
          aux.mom.var <- t(aux.mom) %*% aux.mom/n
          F <- adiag(F, aux.mom.var)
        }
        W <- solve(F, tol = 1e-20)

        as.numeric(t(M) %*% W %*% M)

    }

    list.gmm <- function(pars, J, y, treat, x, w, h, group, matrixMethod = c("efficient", "princomp", "cue")) {

        n <- length(y)
        w <- (n/sum(w)) * w

        y1 <- y[treat == 1]
        y0 <- y[treat == 0]
        x1 <- x[treat == 1, , drop = FALSE]
        x0 <- x[treat == 0, , drop = FALSE]
        w1 <- w[treat == 1] 
        w0 <- w[treat == 0]
        group1 <- group[treat == 1]

        if (matrixMethod == "efficient" | matrixMethod == "princomp") {
            # step 1
            W0 <- invF(pars = pars, J = J, y = y, treat = treat, 
                x = x, w = w, h = h, group = group, matrixMethod = matrixMethod)
            fit0 <- optim(par = pars, fn = gmm.nls, gr = gmm.grad,
                J = J, y = y, treat = treat, x = x, w = w, h = h, 
                group = group, W = W0, method = "BFGS", 
                control = list(maxit = maxIter))
            
            # step 2
            W1 <- invF(fit0$par, J = J, y = y, treat = treat, 
                x = x, w = w, h = h, group = group, matrixMethod = matrixMethod)
            fit1 <- optim(par = fit0$par, fn = gmm.nls, gr = gmm.grad,
                J = J, y = y, treat = treat, x = x, w = w, h = h, 
                group = group, W = W1, method = "BFGS", 
                control = list(maxit = maxIter))

        } else if (matrixMethod == "cue") {
            fit1 <- optim(par = pars, fn = gmm.cue, J = J, y = y, 
              treat = treat, x = x, w = w, h = h, group = group, matrixMethod = matrixMethod, 
              method = "BFGS", control = list(maxit = maxIter))

        }
        
        if (fit1$convergence != 0) stop("Optimization routine did not converge.")

        coef <- fit1$par
        delta <- coef[1:ncol(x)]
        gamma <- coef[(ncol(x) + 1):(2 * ncol(x))]

        # AUX Jacobian
        dh <- c()
        if (length(h) > 0) {
          group.labels <- names(h)
          for (label in group.labels) {
              dh.tmp <- -c(logistic(x[group == label, , drop = FALSE] %*% delta)/
                (1 + exp(x[group == label, , drop = FALSE] %*% delta))) * 
                  x[group == label, , drop = FALSE]
              dh <- rbind(dh, t(colSums(dh.tmp))/n)
          }
        }

        # NLS Jacobian
        Gtmp <- c(logistic(x1 %*% delta)/(1 + exp(x1 %*%
            delta))) * x1
        G1 <- -t(Gtmp) %*% Gtmp/n
        Gtmp <- c(sqrt(J * logistic(x1 %*% delta) * logistic(x1 %*%
            gamma)/((1 + exp(x1 %*% delta)) * (1 + exp(x1 %*%
            gamma))))) * x1
        G2 <- -t(Gtmp) %*% Gtmp/n
        Gtmp <- c(J * logistic(x0 %*% gamma)/(1 + exp(x0 %*%
            gamma))) * x0
        G3 <- -t(Gtmp) %*% Gtmp/n
        G <- rbind(cbind(G1, G2), cbind(matrix(0, ncol = ncol(G1), nrow = nrow(G3)), G3))

        # complete Jacobian
        if (length(h) > 0) G <- rbind(G, cbind(dh, matrix(0, nrow = length(h), ncol = ncol(G3))))

        weightMatrix <- invF(pars = coef, J = J, y = y, treat = treat, x = x, w = w, 
          h = h, group = group, matrixMethod = matrixMethod)
        if (matrixMethod == "princomp") 
          vcov.aux <- ginv(t(G) %*% weightMatrix %*% G)/n
        else 
          vcov.aux <- solve(t(G) %*% weightMatrix %*% G, tol = 1e-20)/n
        
        list(coef = coef, vcov = vcov.aux, val = fit1$value)
    
    }
      g.all <- group[na.x == 0 & na.y == 0 & na.w == 0]
      par.nls.std <- c(par.treat.nls.std, par.control.nls.std)
      fit.nls.aux <- list.gmm(pars = par.nls.std, J = J, y = y.all, treat = t, 
        x = x.all, w = w.all, h = h, group = g.all, matrixMethod = matrixMethod)
      par.nls.aux <- fit.nls.aux$coef
      vcov.nls <- fit.nls.aux$vcov
      se.twostep <- sqrt(diag(vcov.nls))
      par.treat.nls.std <- par.nls.aux[1:ncol(x.all)]
      par.control.nls.std <- par.nls.aux[(ncol(x.all)+1):(ncol(x.all)*2)]

      # Sargan-Hansen overidentification test
      J.stat <- fit.nls.aux$val * n
      overid.p <- round(1 - pchisq(J.stat, df = length(h)), 4)
     
  }

  ## 
  ## Set up return object
  
  if (method == "nls") {
    
    if (error == "topcode") {
      return.object <-
        list(
          par.treat = par.treat,
          se.treat = se.treat,
          par.control = par.control,
          se.control = se.control,
          vcov = vcov.nls,
          p.est = p.est,
          p.ci = p.ci,
          coef.names = coef.names,
          J = J,
          design = design,
          method = method,
          fit.start = fit.start,
          overdispersed = overdispersed,
          boundary = boundary,
          multi = multi,
          error = error,
          data = data,
          x = x.all,
          y = y.all,
          treat = t,
          call = match.call()
        )
    } else if (error == "uniform"){
      return.object <-
        list(
          par.treat = par.treat,
          se.treat = se.treat,
          par.control = par.control,
          se.control = se.control,
          vcov = vcov.nls,
          p0.est = p0.est,
          p0.ci = p0.ci,
          p1.est = p1.est,
          p1.ci = p1.ci,
          coef.names = coef.names,
          J = J,
          design = design,
          method = method,
          fit.start = fit.start,
          overdispersed = overdispersed,
          boundary = boundary,
          multi = multi,
          error = error,
          data = data,
          x = x.all,
          y = y.all,
          treat = t,
          call = match.call()
        )
      
    } else {
      
      if (multi == FALSE) {
        
        if(design=="standard")
          par.treat <- par.treat.nls.std
        if(design=="modified")
          par.treat <- par.treat.nls.mod
        se.treat <- se.twostep[1:(length(par.treat))]
        
        if(design=="standard")
          par.control <- par.control.nls.std
        if(design=="modified")
          par.control <- par.control.nls.mod
        se.control <- se.twostep[(length(par.treat)+1):(length(se.twostep))]
        
        names(par.treat) <- names(se.treat) <- coef.names
        
        if (design=="standard")
          names(par.control) <- names(se.control) <- coef.names
        if (design=="modified")
          names(par.control) <- names(se.control) <- rep(coef.names, J)
        
        sum.fit.treat <- summary(fit.treat)
        
        resid.se <- sum.fit.treat$sigma
        resid.df <- sum.fit.treat$df[2]
        
        if(design=="standard") {
          return.object <- list(par.treat=par.treat, se.treat=se.treat, par.control=par.control, se.control=se.control, vcov=vcov.nls, resid.se=resid.se, resid.df=resid.df, coef.names=coef.names,  J=J, design = design, method = method, fit.start = fit.start, overdispersed=overdispersed, boundary = boundary, multi = multi, data = data, x = x.all, y = y.all, treat = t, call = match.call())
          if(weighted == TRUE)
            return.object$weights <- w.all
          
        } else if (design=="modified") {
          return.object <- list(par.treat=par.treat, se.treat=se.treat, par.control=par.control, se.control=se.control, vcov=vcov.twostep, resid.se=resid.se, resid.df=resid.df, coef.names=coef.names, J=J, design = design, method = method, fit.nonsensitive = fit.nonsensitive, data = data, x = x.all, y = y.all, treat = t, boundary = FALSE, multi = FALSE, call = match.call())
          if(weighted == TRUE)
            return.object$weights <- w.all
          
        }
        
      } else if (multi == TRUE) {
        
        par.treat <- par.treat.nls.std
        
        se.treat <- list()
        for (m in 1:length(treatment.values))
          se.treat[[m]] <- se.twostep[[m]][1:(length(par.treat[[m]]))]
        
        par.control <- par.control.nls.std
        se.control <- se.twostep[[1]][(length(par.treat[[1]])+1):(length(se.twostep[[1]]))]
        
        for (m in 1:length(treatment.values)) {
          names(par.treat[[m]]) <- coef.names
          names(se.treat[[m]]) <- coef.names
        }
        
        names(par.control) <- names(se.control) <- coef.names
        
        
        resid.se <- resid.df <- rep(NA, length(treatment.values))
        for (m in 1:length(treatment.values)) {
          sum.fit.treat <- summary(fit.treat[[m]])
          resid.se[m] <- sum.fit.treat$sigma
          resid.df[m] <- sum.fit.treat$df[2]
        }
        
        return.object <- list(par.treat=par.treat, se.treat=se.treat, par.control=par.control, se.control=se.control, vcov=vcov.nls, treat.values = treatment.values, treat.labels = treatment.labels, control.label = control.label, resid.se=resid.se, resid.df=resid.df, J=J,  coef.names=coef.names, design = design, method = method, overdispersed=overdispersed, boundary = boundary, multi = multi, data = data, x = x.all, y = y.all, treat = t, call = match.call())
        if(weighted == TRUE)
          return.object$weights <- w.all
        
      }
      
    }
  }
  
  if (method == "ml" & robust == FALSE) {
    
    if(design == "standard") {
      
      if (multi == FALSE) {
        
        if (constrained == T) {

          par.control <- MLEfit$par[1:(nPar)]

          if (overdispersed == T){
            par.treat <- MLEfit$par[(nPar+2):(nPar*2+1)]
            se.treat <- se.mle[(nPar+2):(nPar*2+1)]
            se.control <- se.mle[1:(nPar)]
            par.overdispersion <- MLEfit$par[nPar+1]
            se.overdispersion <- se.mle[nPar+1]
            names(par.overdispersion) <- names(se.overdispersion) <- "overdispersion"
          } else {
            par.treat <- MLEfit$par[(nPar+1):(nPar*2)]
            se.treat <- se.mle[(nPar+1):(nPar*2)]
            se.control <- se.mle[1:(nPar)]
          }
          
          if(floor==TRUE)
            par.floor <- par[(nPar*2+1):(nPar*2 + nPar.floor)]
          if(floor==TRUE & ceiling==TRUE)
            par.ceiling <- par[(nPar*2 + nPar.floor + 1):(nPar*2 + nPar.floor + nPar.ceiling)]
          if(floor==FALSE & ceiling==TRUE)
            par.ceiling <- par[(nPar*2+1):(nPar*2 + nPar.ceiling)]
          
          if(floor==TRUE)
            se.floor <- se.mle[(nPar*2+1):(nPar*2+ nPar.floor)]
          if(floor==TRUE & ceiling==TRUE)
            se.ceiling <- se.mle[(nPar*2 + nPar.floor + 1):(nPar*2 + nPar.floor + nPar.ceiling)]
          if(floor==FALSE & ceiling==TRUE)
            se.ceiling <- se.mle[(nPar*2+1):(nPar*2 + nPar.ceiling)]
          
          names(par.treat) <- names(se.treat) <- names(par.control) <- names(se.control)  <- coef.names
          
          if(floor==TRUE)
            names(par.floor) <- names(se.floor) <- coef.names.floor
          if(ceiling==TRUE)
            names(par.ceiling) <- names(se.ceiling) <- coef.names.floor
          
          if(boundary == TRUE | multi == TRUE)
            llik.const <- llik
          
          if(boundary==F)
            if (overdispersed == T)
              return.object <- list(par.treat=par.treat, se.treat=se.treat, par.control=par.control, se.control=se.control, par.overdispersion=par.overdispersion, se.overdispersion=se.overdispersion, vcov=vcov.mle, pred.post = w, treat.labels = treatment.labels, control.label = control.label, llik=llik.const, J=J,  coef.names=coef.names, design = design, method = method, overdispersed=overdispersed, constrained=constrained, boundary = boundary, multi = multi, ceiling = ceiling, floor = floor, call = match.call(), data = data, x = x.all, y = y.all, treat = t)
            else
              return.object <- list(par.treat=par.treat, se.treat=se.treat, par.control=par.control, se.control=se.control, vcov=vcov.mle, pred.post = w, treat.labels = treatment.labels, control.label = control.label, llik=llik.const, J=J,  coef.names=coef.names, design = design, method = method, overdispersed=overdispersed, constrained=constrained, boundary = boundary, multi = multi, ceiling = ceiling, floor = floor, call = match.call(), data = data, x = x.all, y = y.all, treat = t)
          
          if(floor==FALSE & ceiling==TRUE)
            return.object <- list(par.treat=par.treat, se.treat=se.treat, par.control=par.control, se.control=se.control, par.ceiling = par.ceiling, se.ceiling = se.ceiling, vcov=vcov.mle, pred.post = w, treat.labels = treatment.labels, control.label = control.label, llik=llik.const, J=J,  coef.names=coef.names, coef.names.ceiling = coef.names.ceiling, design = design, method = method, overdispersed=overdispersed, constrained=constrained, boundary = boundary, multi = multi, ceiling = ceiling, floor = floor, call = match.call(), data = data, x = x.all, y = y.all, treat = t)
          
          if(floor==TRUE & ceiling==FALSE)
            return.object <- list(par.treat=par.treat, se.treat=se.treat, par.control=par.control, se.control=se.control, par.floor = par.floor, se.floor = se.floor, vcov=vcov.mle, pred.post = w, llik=llik.const, treat.labels = treatment.labels, control.label = control.label, J=J, coef.names=coef.names, coef.names.floor = coef.names.floor, design = design, method = method, overdispersed=overdispersed, constrained=constrained, boundary = boundary, multi = multi, ceiling = ceiling, floor = floor, call = match.call(), data = data, x = x.all, y = y.all, treat = t)
          
          if(floor==TRUE & ceiling==TRUE)
            return.object <- list(par.treat=par.treat, se.treat=se.treat, par.control=par.control, se.control=se.control, par.floor = par.floor, se.floor = se.floor, par.ceiling = par.ceiling, se.ceiling = se.ceiling, pred.post = w, vcov=vcov.mle, treat.labels = treatment.labels, control.label = control.label, llik=llik.const, J=J,  coef.names=coef.names, coef.names.floor = coef.names.floor, coef.names.ceiling = coef.names.ceiling, design = design, method = method, overdispersed=overdispersed, constrained=constrained, boundary = boundary, multi = multi, ceiling = ceiling, floor = floor, call = match.call(), data = data, x = x.all, y = y.all, treat = t)
          
        } else if (constrained == FALSE) { 
          
          par.treat <- MLEfit$par[(nPar*2+1):(nPar*3)]
          par.control.psi0 <- MLEfit$par[1:(nPar)]
          par.control.psi1 <- MLEfit$par[(nPar+1):(nPar*2)]
          
          if (overdispersed == T){
            se.treat <- se.mle[(nPar*2+3):(nPar*3 + 2)]
            se.control.psi0 <- se.mle[1:(nPar)]
            se.control.psi1 <- se.mle[(nPar*2+2):(nPar*2+1)]
            par.overdispersion <- MLEfit$par[nPar*2+2]
            se.overdispersion <- se.mle[nPar*2+2]
            names(par.overdispersion) <- names(se.overdispersion) <- "overdispersion"
          } else {
            se.treat <- se.mle[(nPar*2+1):(nPar*3)]
            se.control.psi0 <- se.mle[1:(nPar)]
            se.control.psi1 <- se.mle[(nPar+1):(nPar*2)]
          }
          
          names(par.treat) <- names(se.treat) <- names(par.control.psi0) <- names(se.control.psi0) <- names(par.control.psi1) <- names(se.control.psi1) <- coef.names

          if(overdispersed==T)
            return.object <- list(par.treat=par.treat, se.treat=se.treat, par.control.psi0=par.control.psi0, se.control.psi0=se.control.psi0, par.control.psi1=par.control.psi1, se.control.psi1=se.control.psi1, par.overdispersion=par.overdispersion, se.overdispersion=se.overdispersion, vcov=vcov.mle, pred.post = w, treat.labels = treatment.labels, control.label = control.label, llik=llik, J=J,  coef.names=coef.names, design = design, method = method, overdispersed=overdispersed, constrained=constrained, boundary = boundary, multi = multi, call = match.call(), data=data, x = x.all, y = y.all, treat = t)
          else
            return.object <- list(par.treat=par.treat, se.treat=se.treat, par.control.psi0=par.control.psi0, se.control.psi0=se.control.psi0, par.control.psi1=par.control.psi1, se.control.psi1=se.control.psi1, vcov=vcov.mle, pred.post = w, treat.labels = treatment.labels, control.label = control.label, llik=llik, J=J,  coef.names=coef.names, design = design, method = method, overdispersed=overdispersed, constrained=constrained, boundary = boundary, multi = multi, call = match.call(), data=data, x = x.all, y = y.all, treat = t)
	
        }

        if(weighted == TRUE)
            return.object$weights <- w.all

      } else if (multi == TRUE) {

        par.control <- MLEfit$par[1:(nPar)]
        se.control <- se.mle[1:(nPar)]

        if (multi.condition == "none") {
          se.treat <- list()
          for (m in 1:length(treatment.values))
            se.treat[[m]] <- se.mle[(nPar + (m-1) * nPar + 1) : (nPar + m * nPar)]
          
          for (m in 1:length(treatment.values)) {
            names(par.treat[[m]]) <- coef.names
            names(se.treat[[m]]) <- coef.names
          }
        } else if (multi.condition == "level") {
          se.treat <- list()
          for (m in 1:length(treatment.values))
            se.treat[[m]] <- se.mle[(nPar + (m-1) * (nPar+1) + 1) :
                                    (nPar + m * (nPar+1))]
          
          for (m in 1:length(treatment.values)) {
            names(par.treat[[m]]) <- c(coef.names, "y_i(0)")
            names(se.treat[[m]]) <- c(coef.names, "y_i(0)")
          }
        }
        
        names(par.control) <- names(se.control) <- coef.names
        
        return.object <- list(par.treat=par.treat, se.treat=se.treat, par.control=par.control, se.control=se.control, vcov=vcov.mle, pred.post = w, treat.values = treatment.values, treat.labels = treatment.labels, control.label = control.label, multi.condition = multi.condition, llik=llik, J=J,  coef.names=coef.names, design = design, method = method, overdispersed=overdispersed, constrained=constrained, boundary = boundary, multi = multi, call = match.call(), data = data, x = x.all, y = y.all, treat = t)
        if(weighted == TRUE)
            return.object$weights <- w.all

      }
      
    } else if (design == "modified") {
      
      par.treat <- MLEfit$par[(nPar*J+1):(nPar*(J+1))]
      par.control <- MLEfit$par[1:(nPar*J)]
      
      se.treat <- se.mle[(nPar*J+1):(nPar*(J+1))]
      se.control <- se.mle[1:(nPar*J)]
      
      names(par.treat) <- names(se.treat) <- coef.names
      
      names(par.control) <- names(se.control) <- rep(coef.names, J)
      
      return.object <- list(par.treat=par.treat, se.treat=se.treat, par.control=par.control, se.control=se.control, vcov=vcov.mle, llik=llik, treat.labels = treatment.labels, control.label = control.label, coef.names=coef.names,  J=J, design = design, method = method, boundary = FALSE, multi = FALSE, call = match.call(), data=data, x = x.all, y = y.all, treat = t)
      if(weighted == TRUE)
          return.object$weights <- w.all
    }
  }

  if (robust) {

    par.treat <- ictrobust.par[1:nPar]
    par.control <- ictrobust.par[(nPar+1):(nPar*2)]
    se.treat <- ictrobust.se[1:(nPar)]
    se.control <- ictrobust.se[(nPar+1):(nPar*2)]
  
    # flipping the vcov matrix
    ictrobust.vcov.flip <- ictrobust.vcov
    ictrobust.vcov.flip[1:nPar, 1:nPar] <- ictrobust.vcov[(nPar+1):(nPar*2), (nPar+1):(nPar*2)]
    ictrobust.vcov.flip[(nPar+1):(nPar*2), (nPar+1):(nPar*2)] <- ictrobust.vcov[1:nPar, 1:nPar]

    return.object <- list(par.treat=par.treat, se.treat=se.treat, par.control=par.control, se.control=se.control, 
      vcov=ictrobust.vcov.flip, treat.labels = treatment.labels, control.label = control.label, 
        llik=llik.robust, J=J, coef.names=coef.names, design = design, method = method, robust = robust, 
          overdispersed=overdispersed, constrained=constrained, boundary = boundary, multi = multi, 
            ceiling = ceiling, floor = floor, call = match.call(), data = data, x = x.all, y = y.all, treat = t)
  
  }      

  # measurement error models -- setting up return objects
  if (error == "topcode" & method == "ml") {

    par.treat <- topcode.par.treat
    par.control <- topcode.par.control
    se.treat <- topcode.se.treat
    se.control <- topcode.se.control
    p.est <- topcode.p

    return.object <- list(par.treat=par.treat, se.treat=se.treat, par.control=par.control, se.control=se.control, 
      vcov=topcode.vcov, treat.labels = treatment.labels, control.label = control.label, 
        llik=llik.topcode, J=J, coef.names=coef.names, design = design, method = method, robust = robust, 
          overdispersed=overdispersed, constrained=constrained, boundary = boundary, multi = multi, 
            ceiling = ceiling, floor = floor, call = match.call(), data = data, x = x.all, y = y.all, treat = t, 
              p.est = topcode.p, p.ci = topcode.p.ci, lp.est = topcode.lp, lp.se = topcode.lp.se, iters = topcode.iter)
  
  }      

  if (error == "uniform" & method == "ml") {

    par.treat <- uniform.par.treat
    par.control <- uniform.par.control
    se.treat <- uniform.se.treat
    se.control <- uniform.se.control
    p0.est <- uniform.p0
    p1.est <- uniform.p1
    p0.est <- uniform.p0
    p1.est <- uniform.p1

    return.object <- list(par.treat=par.treat, se.treat=se.treat, par.control=par.control, se.control=se.control, 
      vcov=uniform.vcov, treat.labels = treatment.labels, control.label = control.label, 
        llik=llik.uniform, J=J, coef.names=coef.names, design = design, method = method, robust = robust, 
          overdispersed=overdispersed, constrained=constrained, boundary = boundary, multi = multi, 
            ceiling = ceiling, floor = floor, call = match.call(), data = data, x = x.all, y = y.all, treat = t, 
      p0.est = uniform.p0, p0.ci = uniform.p0.ci, lp0.est = uniform.lp0, lp0.se = uniform.lp0.se, 
        p1.est = uniform.p1, p1.ci = uniform.p1.ci, lp1.est = uniform.lp1, lp1.se = uniform.lp1.se, iters = uniform.iter)

  
  }      

  # auxiliary data functionality -- setting up return object
  return.object$aux <- aux.check 

  if (aux.check) {
    return.object$nh <- length(h)
    return.object$wm <- ifelse(matrixMethod == "cue", "continuously updating", 
      ifelse(matrixMethod == "princomp", "principal components", "efficient"))
    return.object$J.stat <- round(J.stat, 4)
    return.object$overid.p <- overid.p
  }


  return.object$error <- error
  
  if (error == "topcode") {
    return.object$p.est <- p.est
  }

  class(return.object) <- "ictreg"
  
  return.object
  
}


print.ictreg <- function(x, ...){
  
  cat("\nItem Count Technique Regression \n\nCall: ")
  
  dput(x$call)

  cat("\nCoefficient estimates\n")

  tb <- as.matrix(x$par.treat)
  colnames(tb) <- "est."
  
  cat("\n")

  print(coef(x))

  cat("\n")

  treat.print <- c()
  for (i in 1:length(x$treat.labels)) {
    treat.print <- c(treat.print, "'", x$treat.labels[i], "'", sep = "")
    if (i != length(x$treat.labels))
      treat.print <- c(treat.print, " and ")
  }
  
  cat("Number of control items J set to ", x$J, ". Treatment groups were indicated by ", sep = "")
  cat(treat.print, sep ="")
  cat(" and the control group by '", x$control.label, "'.\n\n", sep = "")
     
  # auxiliary data functionality -- print details
  if (x$aux) cat("Incorporating ", x$nh, " auxiliary moment(s). Weighting method: ", x$wm, ".\n", 
    "The overidentification test statistic was: ", x$J.stat, " (p < ", x$overid.p, ")", ".\n", sep = "")

  # measurement error models -- print details
  if (x$error == "topcode") cat("Estimated proportion of top-coded respondents: ", x$p.est, ". 95% CI: (", 
    round(x$p.ci[1], 6), ", ", round(x$p.ci[2], 6), ").\n", sep = "")

  if (x$error == "uniform") cat("Estimated proportion of respondents with uniform error (control): ", x$p0.est, ". 95% CI: (", 
    round(x$p0.ci[1], 6), ", ", round(x$p0.ci[2], 6), ").\n", 
    "Estimated proportion of respondents with uniform error (treated): ", x$p1.est, ". 95% CI: (", 
    round(x$p1.ci[1], 6), ", ", round(x$p1.ci[2], 6), ").\n", sep = "")

  invisible(x)
  
}

print.predict.ictreg <- function(x, ...){
  
  cat("\nList Experiment Prediction\n")
  
  cat("\nProportion of affirmative responses to the sensitive item\n")

  if (class(x$fit) == "data.frame")
    n <- nrow(x$fit)
  else
    n <- length(x$fit)

  for (i in 1:n) {
    cat(paste("\nPrediction:", rownames(as.matrix(x$fit))[i], "\nEst. = ",
              round(as.matrix(x$fit)[i,1], 4)))
    if (is.null(x$se.fit) == FALSE) {
        cat(paste(" (s.e. = ", round(as.matrix(x$se.fit)[i], 4), ")", sep = ""))
    }
    if (class(x$fit) == "data.frame") {
      cat(paste("\n95% confidence interval: (", round(as.matrix(x$fit)[i,2], 4), ", ",
                round(as.matrix(x$fit)[i,3], 4), ")", sep = ""))
    }
    cat("\n")
  }
  
  cat("\n")

  invisible(x)
  
}



#' Predict Method for Item Count Technique
#' 
#' Function to calculate predictions and uncertainties of predictions from
#' estimates from multivariate regression analysis of survey data with the item
#' count technique.
#' 
#' \code{predict.ictreg} produces predicted values, obtained by evaluating the
#' regression function in the frame newdata (which defaults to
#' \code{model.frame(object)}. If the logical \code{se.fit} is \code{TRUE},
#' standard errors of the predictions are calculated. Setting \code{interval}
#' specifies computation of confidence intervals at the specified level or no
#' intervals.
#' 
#' If \code{avg} is set to \code{TRUE}, the mean prediction across all
#' observations in the dataset will be calculated, and if the \code{se.fit}
#' option is set to \code{TRUE} a standard error for this mean estimate will be
#' provided. The \code{interval} option will output confidence intervals
#' instead of only the point estimate if set to \code{TRUE}.
#' 
#' Two additional types of mean prediction are also available. The first, if a
#' \code{newdata.diff} data frame is provided by the user, calculates the mean
#' predicted values across two datasets, as well as the mean difference in
#' predicted value. Standard errors and confidence intervals can also be added.
#' For difference prediction, \code{avg} must be set to \code{TRUE}.
#' 
#' The second type of prediction, triggered if a \code{direct.glm} object is
#' provided by the user, calculates the mean difference in prediction between
#' predictions based on an \code{ictreg} fit and a \code{glm} fit from a direct
#' survey item on the sensitive question. This is defined as the revealed
#' social desirability bias in Blair and Imai (2010).
#' 
#' @param object Object of class inheriting from "ictreg"
#' @param newdata An optional data frame containing data that will be used to
#' make predictions from. If omitted, the data used to fit the regression are
#' used.
#' @param newdata.diff An optional data frame used to compare predictions with
#' predictions from the data in the provided newdata data frame.
#' @param direct.glm A glm object from a logistic binomial regression
#' predicting responses to a direct survey item regarding the sensitive item.
#' The predictions from the ictreg object are compared to the predictions based
#' on this glm object.
#' @param se.fit A switch indicating if standard errors are required.
#' @param interval Type of interval calculation.
#' @param level Significance level for confidence intervals.
#' @param avg A switch indicating if the mean prediction and associated
#' statistics across all obserations in the dataframe will be returned instead
#' of predictions for each observation.
#' @param sensitive.item For multiple sensitive item design list experiments,
#' specify which sensitive item fits to use for predictions. Default is the
#' first sensitive item.
#' @param ... further arguments to be passed to or from other methods.
#' @return \code{predict.ictreg} produces a vector of predictions or a matrix
#' of predictions and bounds with column names fit, lwr, and upr if interval is
#' set. If se.fit is TRUE, a list with the following components is returned:
#' 
#' \item{fit}{vector or matrix as above} \item{se.fit}{standard error of
#' prediction}
#' @author Graeme Blair, UCLA, \email{graeme.blair@ucla.edu}
#' and Kosuke Imai, Princeton University, \email{kimai@princeton.edu}
#' @seealso \code{\link{ictreg}} for model fitting
#' @references Blair, Graeme and Kosuke Imai. (2012) ``Statistical Analysis of
#' List Experiments."  Political Analysis, Vol. 20, No 1 (Winter). available at
#' \url{http://imai.princeton.edu/research/listP.html}
#' 
#' Imai, Kosuke. (2011) ``Multivariate Regression Analysis for the Item Count
#' Technique.'' Journal of the American Statistical Association, Vol. 106, No.
#' 494 (June), pp. 407-416. available at
#' \url{http://imai.princeton.edu/research/list.html}
#' @keywords models regression
#' @examples
#' 
#' data(race)
#' 
#' race.south <- race.nonsouth <- race
#' 
#' race.south[, "south"] <- 1
#' race.nonsouth[, "south"] <- 0
#' 
#' \dontrun{
#' 
#' # Fit EM algorithm ML model with constraint with no covariates
#' 
#' ml.results.south.nocov <- ictreg(y ~ 1, 
#'    data = race[race$south == 1, ], method = "ml", treat = "treat", 
#'    J = 3, overdispersed = FALSE, constrained = TRUE)
#' ml.results.nonsouth.nocov <- ictreg(y ~ 1, 
#'    data = race[race$south == 0, ], method = "ml", treat = "treat", 
#'    J = 3, overdispersed = FALSE, constrained = TRUE)
#' 
#' # Calculate average predictions for respondents in the South 
#' # and the the North of the US for the MLE no covariates 
#' # model, replicating the estimates presented in Figure 1, 
#' # Imai (2010)
#' 
#' avg.pred.south.nocov <- predict(ml.results.south.nocov,
#'    newdata = as.data.frame(matrix(1, 1, 1)), se.fit = TRUE, 
#'    avg = TRUE)
#' avg.pred.nonsouth.nocov <- predict(ml.results.nonsouth.nocov,
#'    newdata = as.data.frame(matrix(1, 1, 1)), se.fit = TRUE, 
#'    avg = TRUE)
#' 
#' # Fit linear regression
#' 
#' lm.results <- ictreg(y ~ south + age + male + college, 
#'    data = race, treat = "treat", J=3, method = "lm")
#' 
#' # Calculate average predictions for respondents in the 
#' # South and the the North of the US for the lm model, 
#' # replicating the estimates presented in Figure 1, Imai (2010)
#' 
#' avg.pred.south.lm <- predict(lm.results, newdata = race.south, 
#'    se.fit = TRUE, avg = TRUE)
#' 
#' avg.pred.nonsouth.lm <- predict(lm.results, newdata = race.nonsouth, 
#'    se.fit = TRUE, avg = TRUE)
#' 
#' # Fit two-step non-linear least squares regression
#' 
#' nls.results <- ictreg(y ~ south + age + male + college, 
#'    data = race, treat = "treat", J=3, method = "nls")
#' 
#' # Calculate average predictions for respondents in the South 
#' # and the the North of the US for the NLS model, replicating
#' # the estimates presented in Figure 1, Imai (2010)
#' 
#' avg.pred.nls <- predict(nls.results, newdata = race.south, 
#'    newdata.diff = race.nonsouth, se.fit = TRUE, avg = TRUE)
#' 
#' # Fit EM algorithm ML model with constraint
#' 
#' ml.constrained.results <- ictreg(y ~ south + age + male + college, 
#'    data = race, treat = "treat", J=3, method = "ml", 
#'    overdispersed = FALSE, constrained = TRUE)
#' 
#' # Calculate average predictions for respondents in the South 
#' # and the the North of the US for the MLE model, replicating the 
#' # estimates presented in Figure 1, Imai (2010)
#' 
#' avg.pred.diff.mle <- predict(ml.constrained.results, 
#'    newdata = race.south, newdata.diff = race.nonsouth,
#'    se.fit = TRUE, avg = TRUE)
#' 
#' # Calculate average predictions from the item count technique
#' # regression and from a direct sensitive item modeled with
#' # a logit.
#' 
#' # Estimate logit for direct sensitive question
#' 
#' data(mis)
#' 
#' mis.list <- subset(mis, list.data == 1)
#' 
#' mis.sens <- subset(mis, sens.data == 1)
#' 
#' # Fit EM algorithm ML model
#' 
#' fit.list <- ictreg(y ~ age + college + male + south,
#'    J = 4, data = mis.list, method = "ml")
#' 
#' # Fit logistic regression with directly-asked sensitive question
#' 
#' fit.sens <- glm(sensitive ~ age + college + male + south, 
#'    data = mis.sens, family = binomial("logit"))
#' 
#' # Predict difference between response to sensitive item
#' # under the direct and indirect questions (the list experiment).
#' # This is an estimate of the revealed social desirability bias
#' # of respondents. See Blair and Imai (2010).
#' 
#' avg.pred.social.desirability <- predict(fit.list, 
#'    direct.glm = fit.sens, se.fit = TRUE)
#' 
#' }
#' 
#' 
predict.ictreg <- function(object, newdata, newdata.diff, direct.glm, se.fit = FALSE,
                           interval = c("none","confidence"), level = .95, avg = FALSE, sensitive.item, ...){

  ##if(object$design != "standard")
  ##  stop("The predict function only currently works for standard design, single sensitive item experiments without ceiling or floor effects.")
  
  if(missing(interval)) interval <- "none"
  
  nPar <- length(object$coef.names)
  
  ## dummy value for constraint
  if(object$method!="ml") object$constrained <- F
  
  logistic <- function(object) exp(object)/(1+exp(object))

  if(missing(sensitive.item)) {
    sensitive.item <- 1
    if(object$multi==TRUE)
      warning("Using the first sensitive item for predictions. Change with the sensitive.item option.")
  }

  if (class(object$par.treat) != "list") {
    ## extract only treatment coef's, var/covar
    beta <- object$par.treat ## was coef(object)[1:nPar]
    var.beta <- vcov(object)[1:nPar, 1:nPar]
  } else {
    ## extract only treatment coef's, var/covar
    beta <- object$par.treat[[sensitive.item]] ## was coef(object)[1:nPar]
    var.beta <- vcov(object)[(sensitive.item*nPar + 1):((sensitive.item+1)*nPar), (sensitive.item*nPar + 1):((sensitive.item+1)*nPar)]
  }
  
  
  if(missing(newdata)) {
    xvar <- object$x
  } else { 
    if(nrow(newdata)==0)
      stop("No data in the provided data frame.")
    ##formula <- do.call(c, as.list(as.character(print(object$call$formula[[3]]))))
    xvar <- model.matrix(as.formula(paste("~", c(object$call$formula[[3]]))), newdata)
  }

  if (missing(direct.glm) == TRUE) {
    
    if (!missing(newdata.diff)) {
      data.list <- list(xvar,
                        model.matrix(as.formula(paste("~", c(object$call$formula[[3]]))),
                                     newdata.diff))
    } else {
      data.list <- list(xvar)
    }
    
    k <- 1
    
    multi.return.object <- list()
    
    for (xvar in data.list) {
      
      n <- nrow(xvar)
      
      if (object$method == "lm")
        pix <- xvar %*% beta
      else
        pix <- logistic(xvar %*% beta)
      
      v <- pix/(1+exp(xvar %*% beta))
      
      if(avg==FALSE){        
          if (object$method == "lm"){
            var.pred <- diag(xvar %*% var.beta %*% t(xvar))
          } else {
            var.pred <- v^2 * diag(xvar %*% var.beta %*% t(xvar))
          }
        
        se <- sqrt(var.pred)
        
      } else {
        
        pix <- mean(pix)
        
        if (object$method == "lm"){
          var.pred <- sum(xvar %*% var.beta %*% t(xvar))
        } else {
          var.pred <- sum(v %*% t(v) * xvar %*% var.beta %*% t(xvar))
        }
        
        se <- c(sqrt(var.pred)/n)
        
      }
      
      return.object <- list(fit=as.vector(pix))
      
      if(interval=="confidence"){
        
        critical.value <- qt(1-(1-level)/2, df = n)
        ci.upper <- pix + critical.value * se
        ci.lower <- pix - critical.value * se
        
        return.object <- list(fit=data.frame(fit = pix, lwr = ci.lower, upr = ci.upper))
        
      }
      
      if(se.fit==T)
        return.object$se.fit <- as.vector(se)
      
      if (length(data.list) == 2)
        multi.return.object[[k]] <- return.object
      
      k <- k + 1
      
    }
    
    if (length(data.list) == 2) {
      
      if (object$method == "lm")
        stop("Currently the difference option does not work with the lm method. Try ml.")
      
      xa <- data.list[[1]]
      xb <- data.list[[2]]
      
      n <- nrow(xa)
      var <- est <- 0
      pixa <- logistic(xa %*% beta)
      pixb <- logistic(xb %*% beta)
      va <- pixa/(1+exp(xa %*% beta))
      vb <- pixb/(1+exp(xb %*% beta))

      var <- sum( va %*% t(va) * xa %*% var.beta %*% t(xa)) + sum(vb %*% t(vb) * xb %*% var.beta %*% t(xb)) - 
        2 * sum(va %*% t(vb) * xa %*% var.beta %*% t(xb))
        
      est <- mean(pixa-pixb)
      
      se <- as.vector(sqrt(var)) / n
      
      return.object <- list(fit = est)
      
      if(interval=="confidence"){
        
        critical.value <- qt(1-(1-level)/2, df = n)
        ci.upper <- est + critical.value * se
        ci.lower <- est - critical.value * se
        
        return.object <- list(fit=data.frame(fit = est, lwr = ci.lower, upr = ci.upper))
        
      }
      
      if(se.fit==T)
        return.object$se.fit <- as.vector(se)
      
      multi.return.object[[3]] <- return.object
      
      ## now put them together
      
      fit <- se <- list()
      for (j in 1:3) {
        fit[[j]] <- multi.return.object[[j]]$fit
        se[[j]] <- multi.return.object[[j]]$se.fit
      }
      
      return.object <- c()
      return.object$fit <- as.data.frame(do.call(rbind, fit))

      rownames(return.object$fit) <- c("newdata", "newdata.diff", "Difference (newdata - newdata.diff)")
      if (se.fit == T)
        return.object$se.fit <- as.vector(do.call(c, se))

      attr(return.object, "concat") <- TRUE
      
    }


  } else {

    ## direct versus list Monte Carlo prediction

    beta.direct <- coef(direct.glm)
    vcov.direct <- vcov(direct.glm)

    n.draws <- 10000
        
    xvar.direct <- model.matrix(as.formula(paste("~", c(object$call$formula[[3]]))),
                                direct.glm$data)
    
    if (ncol(xvar.direct) != nPar)
      stop("Different number of covariates in direct and list regressions.")
     
    draws.list <- mvrnorm(n = n.draws, mu = beta, Sigma = var.beta)
    draws.direct <- mvrnorm(n = n.draws, mu = beta.direct, Sigma = vcov.direct)

    pred.list.mean <- pred.direct.mean <- pred.diff.mean <- rep(NA, n.draws)
    
    for (d in 1:n.draws) {
      
      par.g <- draws.list[d, ]

      if (object$method == "lm")   
        pred.list <- xvar %*% par.g
      else
        pred.list <- logistic(xvar %*% par.g)

      pred.direct <- logistic(xvar.direct %*% draws.direct[d,])
      
      pred.list.mean[d] <- mean(pred.list)
      pred.direct.mean[d] <- mean(pred.direct)
      pred.diff.mean[d] <- pred.list.mean[d] - pred.direct.mean[d]
      
    }

    est.list <- mean(pred.list.mean)
    se.list <- sd(pred.list.mean)
    est.direct <- mean(pred.direct.mean)
    se.direct <- sd(pred.direct.mean)
    est.diff <- mean(pred.diff.mean)
    se.diff <- sd(pred.diff.mean)
    
    critical.value <- qt(1-(1-level)/2, df = nrow(xvar))
    ci.upper.list <- est.list + critical.value * se.list
    ci.lower.list <- est.list - critical.value * se.list
    ci.upper.direct <- est.direct + critical.value * se.direct
    ci.lower.direct <- est.direct - critical.value * se.direct
    ci.upper.diff <- est.diff + critical.value * se.diff
    ci.lower.diff <- est.diff - critical.value * se.diff
       
    fit.matrix <- as.data.frame(rbind(c(est.list, ci.lower.list, ci.upper.list),
                                         c(est.direct, ci.lower.direct, ci.upper.direct),
                                         c(est.diff, ci.lower.diff, ci.upper.diff)))
    names(fit.matrix) <- c("fit","lwr","upr")
    rownames(fit.matrix) <- c("List", "Direct", "Difference (list - direct)")

    return.object <- c()
    return.object$fit <- fit.matrix
    if (se.fit == T)
      return.object$se.fit <- c(se.list, se.direct, se.diff)

    attr(return.object, "concat") <- TRUE
    
  }
  
  class(return.object) <- "predict.ictreg"
  
  return.object
  
}

coef.ictreg <- function(object, ...){
  
  nPar <- length(object$coef.names)
  
  if((object$method=="lm" | object$method=="nls" | object$design=="modified") & object$boundary == FALSE & object$multi == FALSE){
    coef <- c(object$par.treat,object$par.control)
    
    if(object$design=="standard") {
      names(coef) <- c(paste("sensitive.",object$coef.names,sep=""), paste("control.",object$coef.names,sep=""))
    } else if(object$design=="modified") {
      coef.names <- paste("sensitive.",object$coef.names,sep="")
      for(j in 1:object$J)
        coef.names <- c(coef.names, paste("control.", j, ".", object$coef.names, sep=""))
      names(coef) <- coef.names
    }
    
  } else if(object$method=="ml" & object$boundary == FALSE & object$multi == FALSE) {
    
    if(object$design=="standard"){
      if(object$constrained==F) {
        coef <- c(object$par.treat,object$par.control.psi0,object$par.control.psi1)
        names(coef) <- c(paste("sensitive.",object$coef.names,sep=""), paste("control.psi0.",object$coef.names,sep=""),paste("control.psi1.",object$coef.names,sep=""))
      } else if (object$design=="modified") {
        coef <- c(object$par.treat,object$par.control)
        names(coef) <- c(paste("sensitive.",object$coef.names,sep=""), paste("control.",object$coef.names,sep=""))
      } else if (object$constrained == T) {
        coef <- c(object$par.treat,object$par.control)
        names(coef) <- c(paste("sensitive.",object$coef.names,sep=""),
                         paste("control.",object$coef.names,sep=""))
      } 
    } 
  }

  if (object$boundary == TRUE) {
    if (object$floor == TRUE & object$ceiling == TRUE) {
      coef <- c(object$par.control, object$par.treat,
                object$par.floor, object$par.ceiling)
      coef.names <- paste("control.",object$coef.names,sep="")
      coef.names <- c(coef.names, paste("sensitive.",object$coef.names,sep=""))
      coef.names <- c(coef.names, paste("floor.",object$coef.names.floor,sep=""))
      coef.names <- c(coef.names, paste("ceiling.",object$coef.names.ceiling,sep=""))
      names(coef) <- coef.names
    } else if (object$floor == TRUE){
      coef <- c(object$par.control, object$par.treat, object$par.floor)
      coef.names <- paste("control.",object$coef.names,sep="")
      coef.names <- c(coef.names, paste("sensitive.",object$coef.names,sep=""))
      coef.names <- c(coef.names, paste("floor.",object$coef.names.floor,sep=""))
      names(coef) <- coef.names
    } else  if (object$ceiling == TRUE) {
      coef <- c(object$par.control, object$par.treat, object$par.ceiling)
      coef.names <- paste("control.",object$coef.names,sep="")
      coef.names <- c(coef.names, paste("sensitive.",object$coef.names,sep=""))
      coef.names <- c(coef.names, paste("ceiling.",object$coef.names.ceiling,sep=""))
      names(coef) <- coef.names
    }
  }

  if(object$multi == TRUE) {

    if (object$method == "nls" | object$method == "lm")
      object$multi.condition <- "none"
    
    coef <- c(do.call(c, object$par.treat), object$par.control)
    coef.names <- c() 
    if(object$multi.condition == "none")
      for(j in 1:length(object$treat.labels)) coef.names <- c(coef.names, paste("sensitive.", object$treat.labels[j], ".", object$coef.names,sep=""))
    if(object$multi.condition == "level")
      for(j in 1:length(object$treat.labels)) coef.names <- c(coef.names, paste("sensitive.", object$treat.labels[j], ".", c(object$coef.names, "y_i(0)"),sep=""))
    coef.names <- c(coef.names, paste("control.",object$coef.names,sep=""))
    names(coef) <- coef.names

  }
  
  return(coef)
  
}

vcov.ictreg <- function(object, ...){
  
  vcov <- object$vcov
  
  nPar <- length(object$coef.names)
  
  ## dummy value for constraint
  if(object$method=="nls" | object$design=="modified" | object$method == "lm")
    object$constrained <- F
  
  if (object$method == "lm" | (object$method=="ml" & object$constrained==T &
     object$boundary == F & object$multi == F)) {
    vcov <- rbind( cbind( vcov[(nPar+1):(nPar*2), (nPar+1):(nPar*2)],
                         vcov[(nPar+1):(nPar*2), 1:nPar]  ),
                  cbind(vcov[1:nPar, (nPar+1):(nPar*2)] ,vcov[1:nPar, 1:nPar])  )
  } else if (object$method=="ml" & object$constrained==F & object$boundary == F & object$multi == F) {
    vcov <- rbind( cbind(vcov[(nPar*2+1):(nPar*3),(nPar*2+1):(nPar*3)],
                         vcov[(nPar*2+1):(nPar*3),1:nPar],
                         vcov[(nPar*2+1):(nPar*3), (nPar+1):(nPar*2)] ),
                  cbind(vcov[1:nPar, (nPar*2+1):(nPar*3)] , vcov[1:nPar,1:nPar] ,
                        vcov[1:nPar, (nPar+1):(nPar*2)]),
                  cbind(vcov[(nPar+1):(nPar*2), (nPar*2+1):(nPar*3)] ,
                        vcov[(nPar+1):(nPar*2), 1:nPar] ,
                        vcov[(nPar+1):(nPar*2),(nPar+1):(nPar*2)]) )
  }

  if((object$method=="nls" | object$method == "lm") &
     object$boundary == FALSE & object$multi == FALSE){
    if(object$design=="standard") rownames(vcov) <-
      colnames(vcov)<- c(paste("sensitive.",object$coef.names,sep=""),
                         paste("control.",object$coef.names,sep=""))
    else if(object$design=="modified") {
      coef.names<- paste("sensitive.",object$coef.names,sep="")
      for(j in 1:object$J) coef.names <-
        c(coef.names, paste("control.", j, ".", object$coef.names, sep=""))
      rownames(vcov) <- colnames(vcov) <- coef.names
    }
  } else if(object$method=="ml" & object$boundary == FALSE & object$multi == FALSE) {
    if(object$constrained==F) {
      rownames(vcov) <- colnames(vcov) <- c(paste("sensitive.",object$coef.names,sep=""),
                                            paste("control.psi0.",object$coef.names,
                                                  sep=""),
                                            paste("control.psi1.",object$coef.names,
                                                  sep=""))
    } else {
      rownames(vcov) <- colnames(vcov) <-
        c(paste("sensitive.",object$coef.names,sep=""),
          paste("control.",object$coef.names,sep=""))
    }
  }

  if (object$boundary == TRUE) {
    if (object$floor == TRUE & object$ceiling == TRUE) {
      coef.names <- paste("control.",object$coef.names,sep="")
      coef.names <- c(coef.names, paste("sensitive.",object$coef.names,sep=""))
      coef.names <- c(coef.names, paste("floor.",object$coef.names.floor,sep=""))
      coef.names <- c(coef.names, paste("ceiling.",object$coef.names.ceiling,sep=""))
      rownames(vcov) <- colnames(vcov) <- coef.names
    } else if (object$floor == TRUE){
      coef.names <- paste("control.",object$coef.names,sep="")
      coef.names <- c(coef.names, paste("sensitive.",object$coef.names,sep=""))
      coef.names <- c(coef.names, paste("floor.",object$coef.names.floor,sep=""))
      rownames(vcov) <- colnames(vcov) <- coef.names
    } else  if (object$ceiling == TRUE) {
      coef.names <- paste("control.",object$coef.names,sep="")
      coef.names <- c(coef.names, paste("sensitive.",object$coef.names,sep=""))
      coef.names <- c(coef.names, paste("ceiling.",object$coef.names.ceiling,sep=""))
      rownames(vcov) <- colnames(vcov) <- coef.names
    }
  }

  if(object$multi == TRUE) {

    ##if (object$method == "nls")
    ##  stop("This function does not yet work for the NLS estimator for the multiple sensitive item design.")

    coef.names <- paste("control.",object$coef.names,sep="")
    if(object$method != "ml") {
      for(j in 1:length(object$treat.values)) coef.names <-
        c(coef.names, paste("sensitive.", object$treat.labels[j], ".", object$coef.names,sep=""))
    } else {
      if(object$multi.condition == "none")
        for(j in 1:length(object$treat.values)) coef.names <-
          c(coef.names, paste("sensitive.", object$treat.labels[j], ".", object$coef.names,sep=""))
      else if(object$multi.condition == "level")
        for(j in 1:length(object$treat.values)) coef.names <-
          c(coef.names, paste("sensitive.", object$treat.labels[j], ".", c(object$coef.names, "y_i(0)"),sep=""))
      rownames(vcov) <- colnames(vcov) <- coef.names
    }
  }
  
  return(vcov)
  
}

c.predict.ictreg <- function(...){

  x <- list(...)
  
  for (j in 1:length(x)) 
    x[[j]] <- x[[j]]$fit
  
  return.object <- c()
  return.object$fit <- as.matrix(do.call(rbind, x))

  class(return.object) <- "predict.ictreg"

  return.object
  
}



#' Plot Method for the Item Count Technique
#' 
#' Function to plot predictions and confidence intervals of predictions from
#' estimates from multivariate regression analysis of survey data with the item
#' count technique.
#' 
#' \code{plot.predict.ictreg} produces plots with estimated population
#' proportions of respondents answering the sensitive item in a list experiment
#' in the affirmative, with confidence intervals.
#' 
#' The function accepts a set of \code{predict.ictreg} objects calculated in
#' the following manner:
#' 
#' \code{predict(ictreg.object, avg = TRUE, interval = "confidence")}
#' 
#' For each average prediction, a point estimate and its confidence interval is
#' plotted at equally spaced intervals. The x location of the points can be
#' manipulated with the \code{xvec} option.
#' 
#' Either a single predict object can be plotted, or a group of them combined
#' with \code{c(predict.object1, predict.object2)}. Predict objects with the
#' \code{newdata.diff} option, which calculates the mean difference in
#' probability between two datasets, and the \code{direct.glm} option, which
#' calculates the mean difference between the mean predicted support for the
#' sensitive item in the list experiment and in a direct survey item, can also
#' be plotted in the same way as other \code{predict} objects.
#' 
#' @param x object or set of objects of class inheriting from "predict.ictreg".
#' Either a single object from an \code{ictreg()} model fit or multiple
#' \code{predict} objects combined with the c() function.
#' @param labels a vector of labels for each prediction, plotted at the x axis.
#' @param axes.ict a switch indicating if custom plot axes are to be used with
#' the user-provided estimate \code{labels}.
#' @param xlim a title for the y axis.
#' @param ylim a title for the y axis.
#' @param xlab a title for the x axis.
#' @param ylab a title for the y axis.
#' @param axes an indicator for whether default plot axes are included.
#' @param pch either an integer specifying a symbol or a single character to be
#' used as the default in plotting points.
#' @param xvec a vector of x values at which the proportions will be printed.
#' @param ... Other graphical parameters to be passed to the \code{plot()}
#' command are accepted.
#' @author Graeme Blair, UCLA, \email{graeme.blair@ucla.edu}
#' and Kosuke Imai, Princeton University, \email{kimai@princeton.edu}
#' @seealso \code{\link{ictreg}} for model fitting and
#' \code{\link{predict.ictreg}} for predictions based on the model fits.
#' @references Blair, Graeme and Kosuke Imai. (2012) ``Statistical Analysis of
#' List Experiments."  Political Analysis, Vol. 20, No 1 (Winter). available at
#' \url{http://imai.princeton.edu/research/listP.html}
#' 
#' Imai, Kosuke. (2011) ``Multivariate Regression Analysis for the Item Count
#' Technique.'' Journal of the American Statistical Association, Vol. 106, No.
#' 494 (June), pp. 407-416. available at
#' \url{http://imai.princeton.edu/research/list.html}
#' @keywords models regression
#' @examples
#' 
#' data(race)
#' race.south <- race.nonsouth <- race
#' race.south[, "south"] <- 1
#' race.nonsouth[, "south"] <- 0
#' 
#' \dontrun{
#' 
#' # Fit EM algorithm ML model with constraint
#' ml.constrained.results <- ictreg(y ~ south + age + male + college, 
#'    data = race, treat = "treat", J=3, method = "ml", 
#'    overdispersed = FALSE, constrained = TRUE)
#' 
#' # Calculate average predictions for respondents in the South 
#' # and the the North of the US for the MLE model, replicating the 
#' # estimates presented in Figure 1, Imai (2011)
#' avg.pred.south.mle <- predict(ml.constrained.results, 
#'    newdata = race.south, avg = TRUE, interval = "confidence")
#' avg.pred.nonsouth.mle <- predict(ml.constrained.results, 
#'    newdata = race.nonsouth, avg = TRUE, interval = "confidence")
#' 
#' # A plot of a single estimate and its confidence interval
#' plot(avg.pred.south.mle, labels = "South")
#' 
#' # A  plot of the two estimates and their confidence intervals
#' # use c() to combine more than one predict object for plotting
#' plot(c(avg.pred.south.mle, avg.pred.nonsouth.mle), labels = c("South", "Non-South"))
#' 
#' # The difference option can also be used to simultaneously
#' # calculate separate estimates of the two sub-groups
#' # and the estimated difference. This can also be plotted.
#' 
#' avg.pred.diff.mle <- predict(ml.constrained.results, 
#'    newdata = race.south, newdata.diff = race.nonsouth,
#'    se.fit = TRUE, avg = TRUE, interval="confidence")
#' 
#' plot(avg.pred.diff.mle, labels = c("South", "Non-South", "Difference"))
#' 
#' # Social desirability bias plots
#' 
#' # Estimate logit for direct sensitive question
#' 
#' data(mis)
#' 
#' mis.list <- subset(mis, list.data == 1)
#' 
#' mis.sens <- subset(mis, sens.data == 1)
#' 
#' # Fit EM algorithm ML model
#' 
#' fit.list <- ictreg(y ~ age + college + male + south,
#'    J = 4, data = mis.list, method = "ml")
#' 
#' # Fit logistic regression with directly-asked sensitive question
#' 
#' fit.sens <- glm(sensitive ~ age + college + male + south, 
#'    data = mis.sens, family = binomial("logit"))
#' 
#' # Predict difference between response to sensitive item
#' # under the direct and indirect questions (the list experiment).
#' # This is an estimate of the revealed social desirability bias
#' # of respondents. See Blair and Imai (2010).
#' 
#' avg.pred.social.desirability <- predict(fit.list, 
#'    direct.glm = fit.sens, se.fit = TRUE)
#' 
#' plot(avg.pred.social.desirability)
#' 
#' }
#' @method plot predict.ictreg
plot.predict.ictreg <- function(x, labels = NA, axes.ict = TRUE,
                                xlim = NULL, ylim = NULL, xlab = NULL, ylab = "Estimated Proportion",
                                axes = F, pch = 19, xvec = NULL, ...){

  x <- as.matrix(x$fit)

  if (is.null(xlim))
    xlim <- c(0.5,nrow(x)+0.5)
  if (is.null(ylim))
    ylim <- c(min(x[,2]-.05,0), max(x[,3]+.05))  
  if (is.null(xvec))
    xvec <- 1:nrow(x)  
  if (is.null(xlab))
    xlab <- ""

  plot(xvec, x[,1], pch = pch, xlim = xlim, ylim = ylim, 
       xlab = xlab, ylab = ylab, axes = axes, ...)
  
  critical <- abs(qnorm(0.025))
  for (j in 1:nrow(x))
    lines(rep(xvec[j], 2), x[j, 2:3]) 
  
  abline(h = 0, lty = "dashed")
  
  if (axes.ict == TRUE) {
    axis(side = 2, at = seq(from = round(min(x[,2],0),1),
                     to = round(max(x[,3]),1), by = 0.1))
    axis(side = 1, at = xvec, tick = FALSE, labels = labels)
  }
  
}



#' Summary Method for the Item Count Technique
#' 
#' Function to summarize results from list experiment regression based on the
#' ictreg() function, and to produce proportions of liars estimates.
#' 
#' \code{predict.ictreg} produces a summary of the results from an
#' \code{ictreg} object. It displays the coefficients, standard errors, and fit
#' statistics for any model from \code{ictreg}.
#' 
#' \code{predict.ictreg} also produces estimates of the conditional probability
#' of lying and of the population proportion of liars for boundary models from
#' \code{ictreg()} if \code{ceiling = TRUE} or \code{floor = TRUE}.
#' 
#' The conditional probability of lying for the ceiling model is the
#' probability that a respondent with true affirmative views of all the
#' sensitive and non-sensitive items lies and responds negatively to the
#' sensitive item. The conditional probability for the floor model is the
#' probability that a respondent lies to conceal her true affirmative views of
#' the sensitive item when she also holds true negative views of all the
#' non-sensitive items. In both cases, the respondent may believe her privacy
#' is not protected, so may conceal her true affirmative views of the sensitive
#' item.
#' 
#' @param object Object of class inheriting from "ictreg"
#' @param boundary.proportions A switch indicating whether, for models with
#' ceiling effects, floor effects, or both (indicated by the \code{floor =
#' TRUE}, \code{ceiling = TRUE} options in ictreg), the conditional probability
#' of lying and the population proportions of liars are calculated.
#' @param n.draws For quasi-Bayesian approximation based predictions, specify
#' the number of Monte Carlo draws.
#' @param ... further arguments to be passed to or from other methods.
#' @author Graeme Blair, UCLA, \email{graeme.blair@ucla.edu}
#' and Kosuke Imai, Princeton University, \email{kimai@princeton.edu}
#' @seealso \code{\link{ictreg}} for model fitting
#' @references Blair, Graeme and Kosuke Imai. (2012) ``Statistical Analysis of
#' List Experiments."  Political Analysis, Vol. 20, No 1 (Winter). available at
#' \url{http://imai.princeton.edu/research/listP.html}
#' 
#' Imai, Kosuke. (2011) ``Multivariate Regression Analysis for the Item Count
#' Technique.'' Journal of the American Statistical Association, Vol. 106, No.
#' 494 (June), pp. 407-416. available at
#' \url{http://imai.princeton.edu/research/list.html}
#' @keywords models regression
#' @examples
#' 
#' data(race)
#' \dontrun{
#' # Fit standard design ML model with ceiling effects
#' # Replicates Table 7 Columns 3-4 in Blair and Imai (2012)
#' 
#' ceiling.results <- ictreg(y ~ age + college + male + south, treat = "treat", 
#' 		   	  J = 3, data = affirm, method = "ml", fit.start = "nls",
#' 			  ceiling = TRUE, ceiling.fit = "bayesglm",
#' 			  ceiling.formula = ~ age + college + male + south)
#' 
#' # Summarize fit object and generate conditional probability 
#' # of ceiling liars the population proportion of ceiling liars,
#' # both with standard errors.
#' # Replicates Table 7 Columns 3-4 last row in Blair and Imai (2012)
#' 
#' summary(ceiling.results, boundary.proportions = TRUE)
#' }
#' 
#' 
summary.ictreg <- function(object, boundary.proportions = FALSE, n.draws = 10000, ...) {
  object$boundary.proportions <- boundary.proportions
  object$n.draws <- n.draws
  structure(object, class = c("summary.ictreg", class(object)))
}

print.summary.ictreg <- function(x, ...){
  
  cat("\nItem Count Technique Regression \n\nCall: ")
  
  dput(x$call)

  cat("\n")
  
  if(x$method=="nls" | x$design=="modified" | x$method == "lm")
    x$constrained <- T

  if (x$design == "standard") {
  
    if((x$method=="nls" | x$method == "lm") & x$multi == FALSE){
      ##cat(rep(" ", max(nchar(x$coef.names))+8), "Sensitive Item", rep(" ", 5),
      ##    "Control Items \n", sep="")
      tb.treat <- tb.control <- matrix(NA, ncol = 2, nrow = length(x$par.control))
      colnames(tb.treat) <- colnames(tb.control) <- c("Est.", "S.E.")
      if(x$design=="standard") {
        rownames(tb.treat) <- rownames(tb.control) <- x$coef.names
      } else if(x$design=="modified") {
        rownames(tb.treat) <- rownames(tb.control) <- rep(x$coef.names, x$J)
      }
      tb.treat[,1] <- x$par.treat
      tb.treat[,2] <- x$se.treat
      tb.control[,1] <- x$par.control
      tb.control[,2] <- x$se.control
      
      cat("Sensitive item \n")
      print(as.matrix(round(tb.treat,5)))
      cat("\nControl items \n")
      print(as.matrix(round(tb.control,5)))    
      
      summ.stat <- paste("Residual standard error:",signif(x$resid.se, digits=6),
                         "with",x$resid.df,"degrees of freedom")
      
      cat("\n",summ.stat,"\n\n", sep="")
      
    } else if(x$method=="ml" | x$multi == TRUE) {
      
      if (x$multi == FALSE) {
        if(x$constrained==F) {
          ##cat(rep(" ", max(nchar(x$coef.names))+7), "Sensitive Item", rep(" ", 3),
          ##    "Control (psi0)", rep(" ",2), "Control (psi1)\n", sep="")
          tb.treat <- tb.control.psi0 <- tb.control.psi1 <- matrix(NA, ncol = 2, nrow = length(x$par.treat))
          colnames(tb.treat) <- colnames(tb.control.psi0) <- colnames(tb.control.psi1) <- c("Est.", "S.E.")
          rownames(tb.treat) <- rownames(tb.control.psi0) <- rownames(tb.control.psi1) <- x$coef.names
          tb.treat[,1] <- x$par.treat
          tb.treat[,2] <- x$se.treat
          tb.control.psi0[,1] <- x$par.control.psi0
          tb.control.psi0[,2] <- x$se.control.psi0
          tb.control.psi1[,1] <- x$par.control.psi1
          tb.control.psi1[,2] <- x$se.control.psi1

          cat("Sensitive item \n")
          print(as.matrix(round(tb.treat, 5)))
          cat("\nControl items (psi0) \n")
          print(as.matrix(round(tb.control.psi0,5)))
          cat("\nControl items (psi1) \n")
          print(as.matrix(round(tb.control.psi1,5)))
          
        } else {
          #cat(rep(" ", max(nchar(x$coef.names))+8), "Sensitive Item", rep(" ", 5),
          #    "Control Items \n", sep="")
          tb.treat <- tb.control <- matrix(NA, ncol = 2, nrow = length(x$par.treat))
          colnames(tb.treat) <- colnames(tb.control) <- c("Est.", "S.E.")
          rownames(tb.treat) <- rownames(tb.control) <- x$coef.names
          tb.treat[,1] <- x$par.treat
          tb.treat[,2] <- x$se.treat
          tb.control[,1] <- x$par.control
          tb.control[,2] <- x$se.control
          
          cat("Sensitive item \n")
          print(as.matrix(round(tb.treat,5)))
          cat("\nControl items \n")
          print(as.matrix(round(tb.control,5)))
          
        }
        
        if (x$overdispersed == TRUE)
          cat(paste("\nOverdispersion parameter: ",
                    round(x$par.overdispersion, 4), " (s.e. = ", round(x$se.overdispersion, 4), ").\n",
                    sep = ""))
        
      } else {
        ## multi
        
        if (x$method == "nls" | x$method == "lm")
          x$multi.condition <- "none"
        
        ##cat(rep(" ", max(nchar(x$coef.names))+5), x$treat.labels[1], sep="")
        ##for(k in 2:length(x$treat.values))
        ##  cat(rep(" ", 7), x$treat.labels[k], "\n", sep = "")
        
        tb.treat <- tb.control <- matrix(NA, ncol = 2, nrow = length(x$par.treat[[1]]))
        colnames(tb.treat) <- c("Est.", "S.E.")
        if(x$multi.condition == "none")
          rownames(tb.treat) <- x$coef.names
        else if (x$multi.condition == "level")
          rownames(tb.treat) <- c(x$coef.names, "y_i(0)")
        else if (x$multi.condition == "indicators")
          rownames(tb.treat) <- rep("varnames", length(x$par.treat[[1]]))

        for (k in 1:length(x$treat.values)) {

          tb.treat[,1] <- x$par.treat[[k]]
          tb.treat[,2] <- x$se.treat[[k]]

          cat(paste("Sensitive item (", x$treat.labels[k], ")", "\n", sep = ""))
          print(as.matrix(round(tb.treat,5)))
   
          cat("\n")

        }
                
        ##cat("\n",rep(" ", max(nchar(x$coef.names))+16), "Control\n", sep="")
        cat("Control items\n")
        tb.control <- matrix(NA, ncol = 2, nrow = length(x$par.control))
        colnames(tb.control) <- rep(c("Est.", "S.E."), 1)
        rownames(tb.control) <- x$coef.names
        tb.control[,1] <- x$par.control
        tb.control[,2] <- x$se.control
        
        print(as.matrix(round(tb.control, 5)))
        
        if (x$method == "lm")
          summ.stat <- paste("Residual standard error:",signif(x$resid.se, digits=6),
                             "with",x$resid.df,"degrees of freedom")
        
      }
      
      if(x$method == "ml")
        summ.stat <- paste("Log-likelihood:", round(x$llik, 3))
      else if (x$method == "nls") {
        summ.stat <- "Residual standard error of "
        for (m in 1:length(x$treat.values)) {
          if(m==1)
            summ.stat <- paste(summ.stat, signif(x$resid.se[m], digits=6), " with ",x$resid.df[m]," degrees of freedom for sensitive item ",m,"", sep = "")
          else
            summ.stat <- paste(summ.stat, "; residual standard error of ", signif(x$resid.se[m], digits=6), " with ",x$resid.df[m]," degrees of freedom for sensitive item ",m,"", sep ="")
        }
        summ.stat <- paste(summ.stat, ".", sep = "")
      }
      
      if(x$boundary == FALSE)
        cat("\n",summ.stat,"\n\n", sep="")
      
    }
    
    if (x$boundary == TRUE & x$method == "ml" & x$design == "standard") {
      
      tb.ceiling <- tb.floor <- matrix(NA, ncol = 2, nrow = length(x$par.treat))
      colnames(tb.ceiling) <- colnames(tb.floor) <- rep(c("Est.", "S.E."), 1)
      rownames(tb.ceiling) <- rownames(tb.floor) <- x$coef.names
      
      if (x$ceiling == TRUE) {
        cat("\nCeiling\n", sep="")
        tb.ceiling[,1] <- x$par.ceiling
        tb.ceiling[,2] <- x$se.ceiling
        print(as.matrix(round(tb.ceiling, 5)))
      }
      if (x$floor == TRUE) {
        cat("\nFloor\n", sep="")
        tb.floor[,1] <- x$par.floor
        tb.floor[,2] <- x$se.floor
        print(as.matrix(round(tb.floor, 5)))
     }
      
      cat("\n",summ.stat,"\n\n", sep="")
      
      if (x$boundary.proportions == TRUE) {
        
        logistic <- function(x) exp(x)/(1+exp(x))
        logit <- function(x) return(log(x)-log(1-x))
        
        nPar <- ncol(x$x)
        n <- nrow(x$x)
                
        mu.par <- c(x$par.control, x$par.treat)
        if (x$floor == TRUE)
          mu.par <- c(mu.par, x$par.floor)
        if (x$ceiling == TRUE)
          mu.par <- c(mu.par, x$par.ceiling)
        
        try(draws <- mvrnorm(n = x$n.draws, mu = mu.par, Sigma = x$vcov, tol = 1e-1))
        
        if (exists("draws")) {
          
          coef.h.draws <- draws[,1:nPar, drop = FALSE]
          coef.g.draws <- draws[,(nPar+1):(2*nPar), drop = FALSE]
          
          if(x$floor==TRUE) {
            coef.ql.draws <- draws[,(2*nPar+1):(3*nPar), drop = FALSE]
            if(x$ceiling==TRUE){
              coef.qu.draws <- draws[,(3*nPar+1):(4*nPar), drop = FALSE] 
            } 
          } else {
            if(x$ceiling==TRUE){
              coef.qu.draws <- draws[,(2*nPar+1):(3*nPar), drop = FALSE] 
            } 
          }
          
          liar.ceiling.sims <- liar.floor.sims <-
            liar.ceiling.pop.sims <- liar.floor.pop.sims <- rep(NA, x$n.draws)
          
          for (i in 1:x$n.draws) {
            if (x$floor==TRUE)
              qlX <- logistic(x$x %*% as.vector(coef.ql.draws[i, , drop = FALSE]))
            if (x$ceiling==TRUE)
              quX <- logistic(x$x %*% as.vector(coef.qu.draws[i, , drop = FALSE]))
            
            hX <- logistic(x$x %*% as.vector(coef.h.draws[i,, drop = FALSE]))
            gX <- logistic(x$x %*% as.vector(coef.g.draws[i,, drop = FALSE]))
            
            hX0 <- dbinom(x = 0, size = x$J, prob = hX, log = FALSE)
            hXJ <- dbinom(x = x$J, size = x$J, prob = hX, log = FALSE)
            
            if (x$ceiling==TRUE) {
              liar.ceiling.sims[i] <- sum(quX * hXJ * gX) / sum(hXJ * gX)
              liar.ceiling.pop.sims[i] <- (1/n) * sum(quX * hXJ * gX)
            }
            
            if (x$floor==TRUE) {
              liar.floor.sims[i] <- sum(qlX * hX0 * gX) / sum(hX0 * gX)
              liar.floor.pop.sims[i] <- (1/n) * sum(qlX * hX0 * gX)
            }
            
          }
          
          if(x$ceiling==TRUE | x$floor == TRUE)
            cat("Quasi-Bayesian approximation estimates\n")
          
          if(x$ceiling==TRUE) {
            liar.ceiling <- mean(liar.ceiling.sims)
            liar.ceiling.se <- sd(liar.ceiling.sims)/sqrt(nrow(x$x))
            liar.ceiling.pop <- mean(liar.ceiling.pop.sims)
            liar.ceiling.pop.se <- sd(liar.ceiling.pop.sims)/sqrt(nrow(x$x))
            cat("Ceiling liar cond. prob. est. = ", round(liar.ceiling,5), " (s.e. = ",
                round(liar.ceiling.se,5), ")\nCeiling liar pop. prop. est. = ", round(liar.ceiling.pop,5), " (s.e. = ",
                round(liar.ceiling.pop.se,5), ")\n", sep = "")
          } 
          
          if(x$floor==TRUE) {
            liar.floor <- mean(liar.floor.sims)
            liar.floor.se <- sd(liar.floor.sims)/sqrt(nrow(x$x))
            liar.floor.pop <- mean(liar.floor.pop.sims)
            liar.floor.pop.se <- sd(liar.floor.pop.sims)/sqrt(nrow(x$x))
            cat("Floor liar cond. prob. est. = ",round(liar.floor,5), " (s.e. = ",
                round(liar.floor.se,5), ")\nFloor liar pop. prop. est. = ",round(liar.floor.pop,5), " (s.e. = ",
                round(liar.floor.pop.se,5), ")\n", sep = "")
          } 
          
        ##} else {
        }
        
          ## cannot use mvrnorm()
          
          coef.h <- x$par.control
          coef.g <- x$par.treat
          
          if(x$floor == TRUE) {
            coef.ql <- x$par.floor
            qlX <- logistic(x$x %*% coef.ql)
          }
          if(x$ceiling == TRUE) {
            coef.qu <- x$par.ceiling
            quX <- logistic(x$x %*% coef.qu)
          }
          
          hX <- logistic(x$x %*% coef.h)
          gX <- logistic(x$x %*% coef.g)
          
          hX0 <- dbinom(x = 0, size = x$J, prob = hX, log = FALSE)
          hXJ <- dbinom(x = x$J, size = x$J, prob = hX, log = FALSE)
          
           if(x$ceiling == TRUE | x$floor == TRUE)
          cat("\nMaximum likelihood estimates\n")
    
          if (x$ceiling == TRUE) {
            liar.ceiling <- sum(quX * hXJ * gX) / sum(hXJ * gX)
            liar.ceiling.pop <- (1/n) * sum(quX * hXJ * gX)
            cat("Ceiling liar cond. prob. est. =", round(liar.ceiling,5), "\n")
            cat("Ceiling liar pop. prop. est. =", round(liar.ceiling.pop,5),"\n")
          }
          if (x$floor == TRUE) {
            liar.floor <- sum(qlX * hX0 * gX) / sum(hX0 * gX)
            liar.floor.pop <- (1/n) * sum(qlX * hX0 * gX)
            cat("Floor liar cond. prob. est. =", round(liar.floor,5), "\n")
            cat("Floor liar pop. prop. est. =", round(liar.floor.pop,5),"\n")
          }
          if(!exists("draws")) {
          cat("\nNote: covariance matrix was not computationally singular,\n")
          cat("so std. errors for the proportion estimates could not be calculated.\n\n")
          } else {
          cat("\n")
          }
        }
        
      }
    
  } else {

    ## modified design
    
    cat("Sensitive Item \n", sep = "")
    tb.control <- matrix(NA, ncol = 2, nrow = length(x$par.control))
    tb.treat <- matrix(NA, ncol = 2, nrow = length(x$par.treat))
    colnames(tb.control) <- colnames(tb.treat) <- c("Est.", "S.E.")
    rownames(tb.control) <- rep(x$coef.names, x$J)
    rownames(tb.treat) <- x$coef.names
    
    tb.control[,1] <- x$par.control
    tb.control[,2] <- x$se.control

    tb.treat[,1] <- x$par.treat
    tb.treat[,2] <- x$se.treat
    
    print(as.matrix(round(tb.treat,5)))
    
    cat("Control Items \n", sep="")
    print(as.matrix(round(tb.control,5)))

    if (x$method == "nls")
      summ.stat <- paste("Residual standard error:",signif(x$resid.se, digits=6),
                         "with",x$resid.df,"degrees of freedom")
    else if (x$method == "ml")
      summ.stat <- paste("Log-likelihood:", round(x$llik,5))
    
    cat("\n",summ.stat,"\n\n", sep="")
      
  }


  treat.print <- c()
  for (i in 1:length(x$treat.labels)) {
    treat.print <- c(treat.print, "'", x$treat.labels[i], "'", sep = "")
    if (i != length(x$treat.labels))
      treat.print <- c(treat.print, " and ")
  }
  
  cat("Number of control items J set to ", x$J, ". Treatment groups were indicated by ", sep = "")
  cat(treat.print, sep ="")
  cat(" and the control group by '", x$control.label, "'.\n\n", sep = "")
  
  # auxiliary data functionality
  if (x$aux) cat("Incorporating ", x$nh, " auxiliary moment(s). Weighting method: ", x$wm, ".\n", 
    "The overidentification test statistic was: ", x$J.stat, " (p < ", x$overid.p, ")", ".\n", sep = "")

  # measurement error models
  if (x$error == "topcode") cat("Estimated proportion of top-coded respondents: ", x$p.est, ". 95% CI: (", 
    round(x$p.ci[1], 6), ", ", round(x$p.ci[2], 6), ").\n", sep = "")

  if (x$error == "uniform") cat("Estimated proportion of respondents with uniform error (control): ", x$p0.est, ". 95% CI: (", 
    round(x$p0.ci[1], 6), ", ", round(x$p0.ci[2], 6), ").\n", 
    "Estimated proportion of respondents with uniform error (treated): ", x$p1.est, ". 95% CI: (", 
    round(x$p1.ci[1], 6), ", ", round(x$p1.ci[2], 6), ").\n", sep = "")

  invisible(x)
  
}
