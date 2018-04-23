#' Regression Models for Misreporting
#' 
#' Note: This function is in development.  The ratio, NLS, and ML estimators 
#' for list experiments with the standard wording (i.e., the treatment item 
#' directly gauges the sensitive item prevalence) without floor effects are 
#' implemented. 
#' 
#' This function implements the methods introduced in Chou (N.d.) and 
#' Eady (2017) for estimating multiple regression models of misreporting. 
#' The methods can be applied to survey designs that combine direct 
#' questioning with the list experiment, a popular tool for eliciting 
#' sensitive attitudes (see Blair and Imai 2011).  

#' @param data data frame containing all variables
#' @param rhs formula with predictor variables
#' @param y vector of indirect responses
#' @param dir vector of direct responses. should be coded such that 1 implies that the respondent has the sensitive trait.  can contain missing values.
#' @param J number of innocuous control itmes
#' @param treat name of list experiment treatment indicator in data 
#' @param sensitive.treatment indicates whether the treatment item is sensitive (I have the sensitive trait) or not (I do not have the sensitive trait). TRUE/FALSE
#' @param method one of lm, nls, ml, or ratio
#' @param model.extremes TRUE/FALSE
#' @param two.way TRUE/FALSE
#' @param bayes indicates whether to regularize logistic regression coefficients using the method of Gelman et al.  Note that continuous coefficients should be standardized by two standard deviations.  TRUE/FALSE, default is FALSE
#' @param spec.test default is TRUE if method = "ml", otherwise FALSE
#' @param infer.design, default is TRUE
#' @param design an N by 2 matrix that indicates if the respondent is assigned to receive the direct question (column 1) and the list experiment (column 2)
#' @param nls.options a list of options that specifies whether to iterate the fit or to estimate in one-step
#' @param ml.options can specifiy tolerance and maximum number of iterations

#' @author Winston Chou, Princeton University, \email{wchou@princeton.edu}
#' @references Chou, Winston. (2018) ``Lying on Surveys.'' Technical report, Princeton University.

#' @return \code{misreg} returns an object of class "misreg". 

logistic  <- function(x) exp(x)/(1+exp(x))
dlogistic <- function(x) exp(x)/(1+exp(x))^2

misreg <- function(data, rhs = NULL, y, dir, treat, J, sensitive.treatment = TRUE, 
  method = "ratio", model.extremes = FALSE, two.way = FALSE, bayes = FALSE, 
  spec.test = FALSE, infer.design = TRUE, design, nls.options, ml.options) {

  Y  <- data[, y]
  D  <- data[, dir]
  Tr <- data[, treat]

  # if formula specified
  if (!is.null(rhs)) {
    
    list.formula <- as.formula(paste0(y, paste0(as.character(rhs), collapse = "")))
    dir.formula  <- as.formula(paste0(dir, paste0(as.character(rhs), collapse = "")))
    X <- model.matrix(rhs)

  }

  # direct estimate of sensitive trait
  dir.est <- mean(d[t==0])
  dir.var <- var(d[t==0])/sum(t==0)
  dir.se  <- sqrt(dir.var)

  # list estimate of sensitive trait
  tau <- ifelse(sensitive.treatment, 
    mean(Y[Tr==1])-mean(Y[Tr==0]), 
      1-(mean(Y[Tr==1])-mean(Y[Tr==0])))
  tau.var <- var(Y[Tr==1])/sum(Tr==1) + var(Y[Tr==0])/sum(Tr==0)
  tau.se  <- sqrt(tau.var)

  if (method == "ratio") {

    # begin ratio estimator
    ratio.est <- function(y, d, t, eZ) {

      r <- mean(d[t==0])/eZ
      r

    }

    # ratio estimator variance
    ratio.var <- function(y, d, t, eZ) {

      n   <- length(y)
      r   <- mean(d[t==0])/eZ

      Vz  <- var(y[t==1]) + var(y[t==0])
      Vd  <- var(d[t==0])
      rho <- cor(y[t==0], d[t==0])

      1/(n*eZ^2) * (Vz + 2*rho*r + Vd * r^2)

    }

    # compute ratio
    r.est <- 1-ratio.est(y=Y, d=D, t=Tr, eZ = tau)
    r.var  <- ratio.var(y=Y, d=D, t=Tr, eZ = tau)

    # return object
    return.obj <- list("s.est" = r.est, "s.var" = r.var, 
      "z.est" = tau, "z.var" = tau.var, 
      "d.est" = dir.est, "d.var" = dir.var)
    class(return.obj) <- c("misreg", "ratio")

    # diagnostics
    return.obj$s.flag <- 0
    return.obj$z.flag <- 0

    if (r.est < 0 | r.est > 1) {
      return.obj$s.flag <- 1
    }

    if ((tau-qnorm(0.975)*tau.se) < 0) {
      return.obj$z.flag <- 1
    }

  } # end ratio estimator

  if (method == "nls") {

    # begin nls estimator
    # nls functions
    {

      ss0 <- function(gamma, y, J, X) {
        resid <- y - J * logistic(X %*% gamma)
        sum(resid^2)
      }

      ss1 <- function(delta, gamma, z, J, X) {
        resid <- z - logistic(X %*% delta)
        sum(resid^2)
      }

      ss2 <- function(beta, delta, d, J, X) {
        resid <- (1-d/logistic(X %*% delta)) - logistic(X %*% beta)
        sum(resid^2)
      }

      m1 <- function(y, J, delta, gamma, beta, X, d, t) {
        n <- nrow(X)
        Xd <- X %*% delta
        Xg <- X %*% gamma
        Xb <- X %*% beta
        Xbd <- X %*% (delta + beta)

        resid1 <- c(y - J*logistic(Xg) - logistic(Xb))[t==1]

        m1 <- c(resid1 * dlogistic(Xd[t==1, ])) * X[t==1, ]
        m1
      }

      m0 <- function(y, J, delta, gamma, beta, X, d, t) {
        n <- nrow(X)
        Xd <- X %*% delta
        Xg <- X %*% gamma
        Xb <- X %*% beta
        Xbd <- X %*% (delta + beta)

        resid0 <- c(y - J*logistic(Xg))[t==0]
        m0 <- (resid0 * J * dlogistic(Xg[t==0, ])) * X[t==0, ]
        m0
      }

      m2 <- function(y, J, delta, gamma, beta, X, d, t) {
        n <- nrow(X)
        Xd <- X %*% delta
        Xg <- X %*% gamma
        Xb <- X %*% beta
        Xbd <- X %*% (delta + beta)

        resid2 <- c(d - logistic(Xd)*logistic(Xb))[t==0]
        m2 <- (resid2 * logistic(Xd[t==0, ]) * dlogistic(Xb[t==0, ])) * X[t==0, ]
        m2
      }

      G1 <- function(delta, gamma, beta, X, d, t) {
        n <- nrow(X)
        Xd <- X %*% delta
        Xg <- X %*% gamma
        Xb <- X %*% beta
        Xbd <- X %*% (delta + beta)
        -t(t*c(exp(2*Xd)/((1+exp(Xd))^4))*X) %*% X/n
      }

      G2 <- function(delta, gamma, beta, X, d, t) {
        n <- nrow(X)
        Xd <- X %*% delta
        Xg <- X %*% gamma
        Xb <- X %*% beta
        Xdg <- X %*% (delta + gamma)
        Xbd <- X %*% (delta + beta)
        -t(t*c(J*exp(Xdg)/((1+exp(Xd))^2*(1+exp(Xg))^2)) * X) %*% X/n
      }

      G3 <- function(delta, gamma, beta, X, d, t) {
        n <- nrow(X)
        Xd <- X %*% delta
        Xg <- X %*% gamma
        Xb <- X %*% beta
        Xdg <- X %*% (delta + gamma)
        Xbd <- X %*% (delta + beta)
        -t((1-t)*c(J^2*exp(2*Xg)/((1+exp(Xg))^4)) * X) %*% X/n
      }

      G4 <- function(delta, gamma, beta, X, d, t) {
        n <- nrow(X)
        Xd <- X %*% delta
        Xg <- X %*% gamma
        Xb <- X %*% beta
        Xdg <- X %*% (delta + gamma)
        Xbd <- X %*% (delta + beta)

        -t(c((d - (2*exp(Xbd)/((1+exp(Xb))*(1+exp(Xd))))) * 
          exp(Xbd)/((1+exp(Xb))^2*(1+exp(Xd))^2))[t==0] * X[t==0, ]) %*% 
            X[t==0, ]/n
      }

      G5 <- function(delta, gamma, beta, X, d, t) {
        n <- nrow(X)
        Xd <- X %*% delta
        Xg <- X %*% gamma
        Xb <- X %*% beta
        Xdg <- X %*% (delta + gamma)
        Xbd <- X %*% (delta + beta)
        -t(c((d*(exp(Xb)-1) - 
            (exp(X %*% (delta+2*beta)) - 2*exp(Xbd))/
              ((1+exp(Xb))*(1+exp(Xd)))) * 
                exp(Xbd)/
                  ((1+exp(Xb))^3*(1+exp(Xd))))[t==0] * X[t==0, ]) %*% 
                    X[t==0, ]/n
      }

      bread <- function(delta, gamma, beta, X, d, t) {

        n <- nrow(X)
        k <- ncol(X)
        z.mat <- matrix(0, nr = k, nc = k)

        G1 <- G1(delta, gamma, beta, X, d, t)
        G2 <- G2(delta, gamma, beta, X, d, t)
        G3 <- G3(delta, gamma, beta, X, d, t)
        G4 <- G4(delta, gamma, beta, X, d, t)
        G5 <- G5(delta, gamma, beta, X, d, t)

        rbind(cbind(G1, G2, z.mat), 
          cbind(z.mat, G3, z.mat), 
            cbind(G4, z.mat, G5))
      }

      br <- bread(delta, gamma, beta, X, d, t)

      meat <- function(y, J, delta, gamma, beta, X, d, t) {

        n <- nrow(X)
        k <- ncol(X)
        z.mat <- matrix(0, nr = k, nc = k)

        m0 <- m0(y, J, delta, gamma, beta, X, d, t)
        m1 <- m1(y, J, delta, gamma, beta, X, d, t)
        m2 <- m2(y, J, delta, gamma, beta, X, d, t)

        rbind(
          cbind(crossprod(m1), z.mat, z.mat), 
          cbind(z.mat, crossprod(m0), t(m0) %*% m2),  
            cbind(z.mat, t(m0) %*% m2, crossprod(m2)))/n

      }

      nls.var <- function(y, J, delta, gamma, beta, X, d, t) {

        n <- nrow(X)
        br <- bread(delta, gamma, beta, X, d, t)
        mt <- meat(y, J, delta, gamma, beta, X, d, t)
        var <- 1/n * solve(br) %*% mt %*% t(solve(br))
        var

      }
    }

    # initial values 
    init.fit    <-  ictreg(formula = list.formula, data, treat, J, method = "lm")
    gamma.start <-  1/5 * init.fit$par.control
    delta.start <-  1/5 * init.fit$par.treat
    beta.start  <- -1/5 * coef(lm(formula = dir.formula, data = data))

    fit.control <- optim(fn = ss0, par = gamma.start, y = Y[Tr == 0], J, X = X[Tr == 0, ])
    y.hat <- J*logistic(X %*% fit.control$par)

    fit.treat <- optim(fn = ss1, par = delta.start, gamma = fit.control$par, 
      z = Y[Tr == 1] - y.hat[Tr == 1], J , X = X[Tr == 1, ])
    z.hat <- logistic(X %*% fit.treat$par)

    fit.lying <- optim(fn = ss2, par = beta.start, delta = fit.treat$par, 
      d = D[Tr == 0], J, X = X[Tr == 0, ])

    return.obj <- list()

    return.obj$par.treat      <- fit.treat$par
    return.obj$par.control    <- fit.control$par
    return.obj$par.misreport  <- fit.lying$par
    
    return.obj$vcov <- nls.var(delta = fit.treat$par, gamma = fit.control$par, beta = fit.lying$par, 
      y = Y, J, X, d = D, t = Tr)
    return.obj$se   <- sqrt(diag(return.obj$vcov))  

  } # end nls estimator

  if (method == "ml") {

    # begin ml estimator
    # ml functions 
    obs.llik <- function(par, X, treat, dir, y, J, bayes = FALSE) {

      n     <- nrow(X)
      k     <- ncol(X)

      delta <- par[1:k]
      gamma <- par[(k+1):(2*k)]
      beta  <- par[(2*k+1):(3*k)]

      f <- logistic(X %*% gamma)
      g <- logistic(X %*% delta)
      h <- logistic(X %*% beta)

      lik <- rep(NA, n)

      lik[treat ==  0 & dir == 0] <- ((1-g*(1-h))*choose(J,y)*f^y*(1-f)^(J-y))[treat ==  0 & dir == 0]
      lik[treat ==  0 & dir == 1] <- (g*(1-h)*choose(J,y)*f^y*(1-f)^(J-y))[treat ==  0 & dir == 1]
      lik[treat ==  1 & y == 0]   <- ((1-g)*(1-f)^J)[treat ==  1 & y == 0]
      lik[treat ==  1 & y > 0 & y <= J] <- (g*choose(J, y-1)*f^(y-1)*(1-f)^(J-y+1) + 
        (1-g)*choose(J, y)*f^y*(1-f)^(J-y))[treat ==  1 & y > 0 & y <= J]
      lik[treat ==  1 & y == J+1] <- (g*f^J)[treat ==  1 & y == J+1]

      intercept.index <- c(1, 1 + k, 1 + 2*k)
      bayes.adjust <- ifelse(bayes, 
        sum(dcauchy(x = par[intercept.index], scale = rep(10, 3), log = TRUE)) +  
        sum(dcauchy(x = par[-intercept.index], scale = rep(2.5, length(par)-3), log = TRUE)), 0)

      return(-sum(log(lik)) - bayes.adjust)

    }

    eZ.step <- function(par, X, treat, dir, y, J) {

      n     <- nrow(X)
      k     <- ncol(X)

      delta <- par[1:k]
      gamma <- par[(k+1):(2*k)]
      beta  <- par[(2*k+1):(3*k)]

      f <- logistic(X %*% gamma)
      g <- logistic(X %*% delta)
      h <- logistic(X %*% beta)

      z.update <- rep(NA, n)

      z.update[treat ==  0 & dir == 0] <- ((g*h)/(1-g*(1-h)))[treat ==  0 & dir == 0]
      z.update[treat ==  0 & dir == 1] <- 1
      z.update[treat ==  1 & y == 0]   <- 0
      z.update[treat ==  1 & y > 0 & y <= J] <- 
        ((g * choose(J, y-1) * f^(y-1) * (1-f)^(J-y+1))/(g * choose(J, y-1) * f^(y-1) * (1-f)^(J-y+1) + 
          (1-g) * choose(J, y) * f^y * (1-f)^(J-y)))[treat ==  1 & y > 0 & y <= J]
      z.update[treat ==  1 & y == J+1] <- 1

      return(z.update)

    }

    eS.step <- function(par, X, treat, dir, y, J) {

      n     <- nrow(X)
      k     <- ncol(X)

      delta <- par[1:k]
      gamma <- par[(k+1):(2*k)]
      beta  <- par[(2*k+1):(3*k)]

      f <- logistic(X %*% gamma)
      g <- logistic(X %*% delta)
      h <- logistic(X %*% beta)
      
      s.update <- rep(NA, n)

      s.update[treat ==  0 & dir == 0] <- 1
      s.update[treat ==  0 & dir == 1] <- 0
      s.update[treat ==  1 & y == 0]   <- 0 # can be anything
      s.update[treat ==  1 & y > 0] <- h[treat ==  1 & y > 0]

      return(s.update)

    }


    # start of algorithm
    fit0.list <- ictreg(list.formula, data, J, treat = treat, method = "lm")
    fit0.dir  <- glm(dir.formula, data, family = binomial("logit"))

    delta.update <- 1/5 * fit0.list$par.treat
    gamma.update <- 1/5 * fit0.list$par.control
    beta.update  <- -coef(fit0.dir)

    init.llik <- obs.llik(par = c(delta.update, gamma.update, beta.update), 
      X = X, treat = Tr, dir = D, y = Y, J, bayes = bayes)

    lliks <- rep(init.llik, 2)

    iter <- 1

    # initial variables

    dat          <- rbind(data, data)
    dat$y        <- rep(1:0, each = nrow(data))
    dat$y.hat    <- c(ifelse(Tr == 1, Y-1, Y), Y)

    while((iter < 10) | abs(lliks[iter] - lliks[iter+1]) > 1e-06) {

      eZ <- eZ.step(par = c(delta.update, gamma.update, beta.update), 
        X = X, treat = Tr, dir = D, y = Y, J = J)
      eS <- eS.step(par = c(delta.update, gamma.update, beta.update), 
        X = X, treat = Tr, dir = D, y = Y, J = J)

      dat$z.weight <- c(eZ, 1-eZ)
      dat$s.weight <- c(eZ, eZ) * c(eS, 1-eS)

      delta.update <- coef(glm(as.formula(paste0("cbind(y, 1-y)", 
          paste0(as.character(rhs), collapse = ""))), 
        family = binomial("logit"), data = dat, weights = z.weight))
      gamma.update <- coef(glm(as.formula(paste0("cbind(y.hat, J-y.hat)", 
          paste0(as.character(rhs), collapse = ""))), 
        family = binomial("logit"), data = subset(dat, z.weight > 0), 
        weights = z.weight))
      beta.update  <- coef(glm(as.formula(paste0("cbind(y, 1-y)", 
          paste0(as.character(rhs), collapse = ""))), 
        family = binomial("logit"), data = dat, weights = s.weight))

      post.llik <- obs.llik(par = c(delta.update, gamma.update, beta.update), 
        X = X, treat = Tr, dir = D, y = Y, J = J, bayes = bayes)

      if (!(signif(post.llik, 6) <= signif(init.llik, 6))) {
        warning("Observed data log-likelihood is not monotonically increasing.")
      }

      iter <- iter + 1

      lliks <- c(lliks, post.llik)

    }

    optim.out <- optim(fn = obs.llik, 
      par = c(delta.update, gamma.update, beta.update), 
      X = X, treat = Tr, dir = D, y = Y, J = J, bayes = bayes, 
      method = "BFGS", hessian = TRUE, control = list(maxit = 8))

    return.obj <- list()

    return.obj$par.treat      <- delta.update
    return.obj$par.control    <- gamma.update
    return.obj$par.misreport  <- beta.update
    
    return.obj$vcov <- solve(optim.out$hessian)
    return.obj$se   <- sqrt(diag(return.obj$vcov))  

  }

  return(return.obj)

}

print.misreg <- function(obj) {
  
  if (class(obj)[1] != "misreg") {
    stop("Not a misreg object.")
  } 

  if (class(obj)[2] == "ratio") {
    s.se <- sqrt(obj$s.var)
    z.se <- sqrt(obj$z.var)
    d.se <- sqrt(obj$d.var)

    cat(paste0(
      "\n", 
      "    Ratio Estimator of the Misreporting Rate", "\n", 
      "    =============================================", "\n", 
      "                                Est.    ", "Std. Err.", "\n", 
      "    Sensitive Trait (direct)    ", 
      sprintf("%0.4f", round(obj$d.est, 4)), 
        "  (", sprintf("%0.4f", round(d.se, 4)), ")", "\n", 
      "    Sensitive Trait (indirect)  ", 
      sprintf("%0.4f", round(obj$z.est, 4)), 
        "  (", sprintf("%0.4f", round(z.se, 4)), ")", "\n", 
      "    Misreporting Rate           ", 
      sprintf("%0.4f", round(obj$s.est, 4)), 
        "  (", sprintf("%0.4f", round(s.se, 4)), ")", 
      "\n"), 
        paste0(ifelse(obj$s.flag == 1, 
          "\n    WARNING: Ratio estimate does not lie between 0 and 1.", "")), 
        paste0(ifelse(obj$z.flag == 1, 
          "\n    WARNING: List estimate does not significantly differ from 0.", "")), 
      "\n\n"
      )
  }

  # if (class(obj)[2] == "nls") {

  #   row.labels <- paste0("\t", names(obj$par.treat))

  # }

}

# cat(paste(paste(c("Sensitive Item", 
#   paste0("  ", names(m.out$par.treat))), collapse = "\n"), "\n"))