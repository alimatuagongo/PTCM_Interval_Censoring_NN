ptm <- proc.time()
library(e1071)
library(xtable)
library(survival)
library(dplyr)
library(writexl)
library(readxl)
library(PRROC)
library(MASS)
library(caTools)
library(splines)


# Creating interval censoring
LR_int=function(y1,len1,l1){ # function to find the left and right end points that contain the event time, provided the event takes place
  if (len1 < 0) stop("len1 must be positive.")
  if (y1 < 0) stop("y1 must be positive.")
  if(y1>0 & y1<=l1){
    a=c(.Machine$double.eps,l1)
  }else{
    k=as.integer((y1-l1)/len1)+1 # i changed from as.integer to round
    a=c(l1+((k-1)*len1),l1+(k*len1))
  }
  return(a)
}

#Generating Data
data.pois.ic <- function(n, alpha, beta, delta, setting) {
  
  n_betas_needed <- ifelse(setting %in% 1:2, 2, 5)
  beta <- c(beta, rep(0, max(0, n_betas_needed - length(beta))))[1:n_betas_needed]
  
  # Generate covariates
  if (setting == 3) {
    p <- 10
    Sigma <- matrix(0, nrow = p, ncol = p)
    for (i in 1:p) {
      for (j in 1:p) {
        Sigma[i, j] <- 0.8^abs(i - j)
      }
    }
    
    Z <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
    Z <- as.data.frame(Z)
    names(Z) <- paste0("z", 1:10)
    
  } else {
    Z <- data.frame(
      z1 = rnorm(n),
      z2 = rnorm(n)
    )
  }
  
  z1 <- Z$z1
  z2 <- Z$z2
  
  # True uncured probability
  if (setting == 1) {
    eta_inc <- 0.3 - 5*z1 - 3*z2
    piz <- exp(eta_inc) / (1 + exp(eta_inc))
  }
  
  if (setting == 2) {
    eta_inc <- 0.5 +4*z1^2 - 2*z2^2
    piz <- exp(eta_inc) / (1 + exp(eta_inc))
  }
  
  if (setting == 3) {
    eta_inc <- (-0.8 * (Z$z1 * Z$z2)) +
      (1.1 * (Z$z2 * Z$z4)) +
      (0.5 * Z$z3) +
      (0.2 * Z$z7^2) -
      (1.3 * sin(Z$z5 * Z$z6)) +
      (1.9 * cos(Z$z7 * Z$z8)) -
      (1.5 * exp( Z$z5 * Z$z6 * Z$z7)) -
      (1.6 * Z$z7 * Z$z8 * Z$z9 * Z$z10) +
      (0.8 * Z$z6 * Z$z7 * Z$z8^2 * Z$z9^2) +
      ( 1.8 * cos(Z$z5 * Z$z6 * Z$z7 * Z$z8 * Z$z9)) +
      (1.2 * sqrt(abs(Z$z6 * Z$z7 * Z$z8 * Z$z9 * Z$z10))) -
      1.2
    
    
    piz <- exp(-exp(eta_inc))
  }
  
  # Keep probabilities away from 0 and 1
  piz <- pmin(pmax(piz, .Machine$double.eps), 1 - .Machine$double.eps)
  
  # Latency linear predictor
  if (setting == 3) {
    X.lat <- as.matrix(Z[, 1:5])
  } else {
    X.lat <- as.matrix(Z[, c("z1", "z2")])
  }
  
  eta_lat <- drop(X.lat %*% beta)
  
  # Generate censoring and cure indicators
  C = rexp(n,rate=delta) # censoring time 
  U = runif(n,0,1)
  L = rep(NA,n) # left end point 
  R = rep(NA,n) # right end point
  D = rep(NA,n) # censoring indicator
  J = rep(NA,n) # cured indicator (J=0 implies cured)
  Sp = rep(NA,n) # overall survival function
  S1 = rep(NA,n) # susceptible survival function
  S = rep(NA,n) # survival function of promotion times
  Y = rep(NA,n)
  
  for (i in 1:n) {
    
    if (U[i] <= 1 - piz[i]) {
      # cured subject
      L[i] <- C[i]
      R[i] <- Inf
      D[i] <- 0
      J[i] <- 0
      Y[i] <- L[i]
      
    } else {
      # uncured subject
      U1 <- runif(1, 1 - piz[i], 1)
      
      T1 <- qweibull(
        log(U1) / log(1 - piz[i]),
        shape = alpha,
        scale = exp(-eta_lat[i] / alpha)
      )
      
      J[i] <- 1
      
      if (C[i] < T1) {
        # censored before event
        D[i] <- 0
        L[i] <- C[i]
        R[i] <- Inf
        Y[i] <- L[i]
        
      } else {
        # interval-observed event
        len <- runif(1, 0.2, 0.7)
        l <- runif(1, 0, 1)
        ans <- LR_int(T1, len, l)
        
        L[i] <- ans[1]
        R[i] <- ans[2]
        D[i] <- 1
        Y[i] <- (L[i] + R[i]) / 2
      }
    }
  }
  
  # True survival quantities evaluated at Y
  S <- exp(-(Y^alpha) * exp(eta_lat))
  Sp <- (1 - piz)^(1 - S)
  S1 <- (Sp - (1 - piz)) / piz
  
  return(data.frame(
    L = L,
    R = R,
    D = D,
    Z,
    J = J,
    uncure = piz,
    Sp = Sp,
    S1 = S1,
    S = S
  ))
}

#LOGIT
em.logit.ic <- function(TimeL.train, TimeR.train, Status.train, X.train, Z.train,
                        TimeL.test, TimeR.test, Status.test, X.test, Z.test,
                        offsetvar, b, alpha, beta, emmax, eps, data) {
  
  n <- length(Status.train)
  m <- length(Status.test)
  
  TimeY.train <- ifelse(Status.train == 0, TimeL.train, TimeR.train)
  
  convergence <- 1000
  i <- 1
  
  while (convergence > eps && i < emmax) {
    
    ## Current incidence estimate on training data
    eta_inc.train <- drop(Z.train %*% b)
    uncureprob.train <- plogis(eta_inc.train)
    uncureprob.train <- pmin(pmax(uncureprob.train, .Machine$double.eps),
                             1 - .Machine$double.eps)
    
    ## Current latency survival on training data
    eta_lat.train <- drop(X.train[, -1, drop = FALSE] %*% beta)
    H.train <- (TimeY.train^alpha) * exp(eta_lat.train)
    survival.train <- exp(-H.train)
    
    ## E-step
    M.train <- Status.train - survival.train * log(1 - uncureprob.train)
    
    ## M-step: incidence/logit part
    Q1 <- function(par) {
      eta1 <- drop(Z.train %*% par)
      uncure1 <- plogis(eta1)
      uncure1 <- pmin(pmax(uncure1, 1e-6), 1 - 1e-6)
      
      theta1 <- -log(1 - uncure1)
      theta1 <- pmax(theta1, 1e-10)
      
      out1 <- sum(M.train * log(theta1)) +
        sum(log(1 - uncure1))
      
      if (!is.finite(out1)) return(1e20)
      
      return(-out1)
    }
    
    update_b <- optim(par = b, fn = Q1, method = "Nelder-Mead")$par
    
    ## M-step: latency/Weibull part
    Q2 <- function(par) {
      alpha_new <- par[1]
      beta_new <- par[-1]
      
      if (alpha_new <= 0) return(1e20)
      
      eta2 <- drop(X.train[, -1, drop = FALSE] %*% beta_new)
      
      HY <- (TimeY.train^alpha_new) * exp(eta2)
      S.Y <- exp(-HY)
      S.Y <- pmin(pmax(S.Y, 1e-12), 1)
      
      out1 <- sum((M.train - Status.train) * log(S.Y))
      
      event_idx <- which(Status.train == 1)
      
      if (length(event_idx) > 0) {
        HL <- (TimeL.train[event_idx]^alpha_new) * exp(eta2[event_idx])
        HR <- (TimeR.train[event_idx]^alpha_new) * exp(eta2[event_idx])
        
        S.L <- exp(-HL)
        S.R <- exp(-HR)
        
        diff_SR <- S.L - S.R
        diff_SR <- pmax(diff_SR, 1e-12)
        
        out2 <- sum(log(diff_SR))
      } else {
        out2 <- 0
      }
      
      out <- out1 + out2
      
      if (!is.finite(out)) return(1e20)
      
      return(-out)
    }
    
    update_latency <- optim(
      par = c(alpha, beta),
      fn = Q2,
      method = "Nelder-Mead"
    )$par
    
    convergence <- sum(c(update_b - b, update_latency - c(alpha, beta))^2)
    
    b <- update_b
    alpha <- update_latency[1]
    beta <- update_latency[-1]
    
    i <- i + 1
  }
  
  ## =========================
  ## Training predictions
  ## =========================
  
  TimeK.train <- ifelse(is.infinite(TimeR.train), TimeL.train, TimeR.train)
  Time.mid.train <- (TimeL.train + TimeK.train) / 2
  
  eta_inc.train <- drop(Z.train %*% b)
  uncureprob.train <- plogis(eta_inc.train)
  uncureprob.train <- pmin(pmax(uncureprob.train, .Machine$double.eps),
                           1 - .Machine$double.eps)
  
  eta_lat.train <- drop(X.train[, -1, drop = FALSE] %*% beta)
  H.mid.train <- (Time.mid.train^alpha) * exp(eta_lat.train)# cumulative hazard of the progression times for train
  S.train <- exp(-H.mid.train)# survival function of the progression times for train
  
  Sp.train <- (1 - uncureprob.train)^(1 - S.train)# survival prob. of the population
  S1.train <- (Sp.train - (1 - uncureprob.train)) / uncureprob.train# survival prob. of the susceptible group
  
  ## =========================
  ## Test predictions
  ## =========================
  
  TimeK.test <- ifelse(is.infinite(TimeR.test), TimeL.test, TimeR.test)
  Time.mid.test <- (TimeL.test + TimeK.test) / 2
  
  eta_inc.test <- drop(Z.test %*% b)
  uncureprob.test <- plogis(eta_inc.test)
  uncureprob.test <- pmin(pmax(uncureprob.test, .Machine$double.eps),
                          1 - .Machine$double.eps)
  
  eta_lat.test <- drop(X.test[, -1, drop = FALSE] %*% beta)
  H.mid.test <- (Time.mid.test^alpha) * exp(eta_lat.test)
  S.test <- exp(-H.mid.test)
  
  Sp.test <- (1 - uncureprob.test)^(1 - S.test)
  S1.test <- (Sp.test - (1 - uncureprob.test)) / uncureprob.test
  
  
  
  return(list(
    b = b,
    latencyfit = c(alpha, beta),
    
    train = list(
      Time.mid = Time.mid.train,
      Uncureprob = uncureprob.train,
      Sp = Sp.train,
      S1 = S1.train,
      S = S.train
    ),
    
    test = list(
      Time.mid = Time.mid.test,
      Uncureprob = uncureprob.test,
      Sp = Sp.test,
      S1 = S1.test,
      S = S.test
    ),
    
    tau = convergence,
    iter = i
  ))
}

smcure.logit.ic <- function(formula, cureform,
                            alpha.init, beta.init,
                            train_data, test_data,
                            offset = NULL,
                            na.action = na.omit,
                            Var = F,
                            emmax = 1000,
                            eps = 1e-3,
                            nboot = 100) {
  
  call <- match.call()
  cat("Program is running..be patient...")
  
  ## =========================
  ## Prepare training data
  ## =========================
  
  train_data <- na.action(train_data)
  test_data  <- na.action(test_data)
  
  n <- nrow(train_data)
  m <- nrow(test_data)
  
  mf.train <- model.frame(formula, train_data)
  mf.test  <- model.frame(formula, test_data)
  
  cvars <- all.vars(cureform)
  
  Z.train <- as.matrix(cbind(1, train_data[, cvars, drop = FALSE]))
  Z.test  <- as.matrix(cbind(1, test_data[, cvars, drop = FALSE]))
  
  colnames(Z.train) <- c("(Intercept)", cvars)
  colnames(Z.test)  <- c("(Intercept)", cvars)
  
  if (!is.null(offset)) {
    offsetvar <- all.vars(offset)
    offsetvar <- train_data[, offsetvar]
  } else {
    offsetvar <- NULL
  }
  
  Y.train <- model.extract(mf.train, "response")
  Y.test  <- model.extract(mf.test, "response")
  
  X.train <- model.matrix(attr(mf.train, "terms"), mf.train)
  X.test  <- model.matrix(attr(mf.test, "terms"), mf.test)
  
  if (!inherits(Y.train, "Surv")) stop("Training response must be a survival object")
  if (!inherits(Y.test, "Surv")) stop("Test response must be a survival object")
  
  TimeL.train <- Y.train[, 1]
  TimeR.train <- Y.train[, 2]
  Status.train <- train_data$D
  Status.test  <- test_data$D
  
  TimeL.test <- Y.test[, 1]
  TimeR.test <- Y.test[, 2]
  
  
  bnm <- colnames(Z.train)
  nb <- ncol(Z.train)
  
  betanm <- colnames(X.train)[-1]
  nbeta <- ncol(X.train) - 1
  
  ## Initial incidence parameters
  b <- rep(0.5, nb)
  
  ## =========================
  ## Fit EM on training data
  ## =========================
  
  emfit <- em.logit.ic(
    TimeL.train = TimeL.train,
    TimeR.train = TimeR.train,
    Status.train = Status.train,
    X.train = X.train,
    Z.train = Z.train,
    TimeL.test = TimeL.test,
    TimeR.test = TimeR.test,
    Status.test = Status.test,
    X.test = X.test,
    Z.test = Z.test,
    offsetvar = offsetvar,
    b = b,
    alpha = alpha.init,
    beta = beta.init,
    emmax = emmax,
    eps = eps,
    data = train_data
  )
  
  latency <- emfit$latencyfit
  b.est <- emfit$b
  
  alpha.est <- latency[1]
  beta.est <- latency[-1]
  
  conv <- emfit$tau
  
  ## =========================
  ## Bootstrap variance
  ## =========================
  
  if (Var) {
    
    latency_boot <- matrix(NA, nrow = nboot, ncol = nbeta + 1)
    
    tempdata <- data.frame(
      TimeL = TimeL.train,
      TimeR = TimeR.train,
      Status = Status.train,
      X.train,
      Z.train
    )
    
    data1 <- subset(tempdata, Status == 1)
    data0 <- subset(tempdata, Status == 0)
    
    n1 <- nrow(data1)
    n0 <- nrow(data0)
    
    i <- 1
    
    while (i <= nboot) {
      
      id1 <- sample(1:n1, n1, replace = TRUE)
      id0 <- sample(1:n0, n0, replace = TRUE)
      
      bootdata <- rbind(data1[id1, ], data0[id0, ])
      
      bootZ <- as.matrix(bootdata[, bnm, drop = FALSE])
      bootX <- as.matrix(cbind(1, bootdata[, betanm, drop = FALSE]))
      
      bootfit <- em.logit.ic(
        TimeL.train = bootdata[, "TimeL"],
        TimeR.train = bootdata[, "TimeR"],
        Status.train = bootdata[, "Status"],
        X.train = bootX,
        Z.train = bootZ,
        TimeL.test = TimeL.test,
        TimeR.test = TimeR.test,
        Status.test = Status.test,
        X.test = X.test,
        Z.test = Z.test,
        offsetvar = offsetvar,
        b = b.est,
        alpha = alpha.est,
        beta = beta.est,
        emmax = emmax,
        eps = eps,
        data = train_data
      )
      
      if (bootfit$tau < eps) {
        latency_boot[i, ] <- bootfit$latencyfit
        i <- i + 1
      }
    }
    
    latency_var <- apply(latency_boot, 2, var, na.rm = TRUE)
    latency_sd <- sqrt(latency_var)
  }
  
  ## =========================
  ## Return fitted object
  ## =========================
  
  fit <- list()
  class(fit) <- "smcure.logit.ic"
  
  fit$call <- call
  fit$b <- b.est
  fit$latency <- latency
  fit$tau <- conv
  fit$iter <- emfit$iter
  
  fit$train <- emfit$train
  fit$test <- emfit$test
  
  if (Var) {
    fit$latency_var <- latency_var
    fit$latency_sd <- latency_sd
    fit$latency_zvalue <- fit$latency / latency_sd
    fit$latency_pvalue <- 2 * (1 - pnorm(abs(fit$latency_zvalue)))
  }
  
  cat(" done.\n")
  
  return(fit)
}



#SPLINE
em.spline.ic <- function(TimeL.train, TimeR.train, Status.train, X.train, Z.train,
                         TimeL.test, TimeR.test, Status.test, X.test, Z.test,
                         offsetvar, b, alpha, beta, emmax, eps, data) {
  
  n <- length(Status.train)
  m <- length(Status.test)
  
  TimeY.train <- ifelse(Status.train == 0, TimeL.train, TimeR.train)
  
  convergence <- 1000
  iter <- 1
  
  while (convergence > eps && iter < emmax) {
    
    ## =========================
    ## Current incidence estimate
    ## =========================
    
    eta_inc.train <- drop(Z.train %*% b)
    uncureprob.train <- plogis(eta_inc.train)
    uncureprob.train <- pmin(pmax(uncureprob.train, .Machine$double.eps),
                             1 - .Machine$double.eps)
    
    ## =========================
    ## Current latency survival
    ## =========================
    
    eta_lat.train <- drop(X.train[, -1, drop = FALSE] %*% beta)
    H.train <- (TimeY.train^alpha) * exp(eta_lat.train)
    survival.train <- exp(-H.train)
    
    ## =========================
    ## E-step
    ## =========================
    
    M.train <- Status.train - survival.train * log(1 - uncureprob.train)
    
    ## =========================
    ## M-step: spline incidence part
    ## =========================
    
    Q1 <- function(par) {
      eta1 <- drop(Z.train %*% par)
      uncure1 <- plogis(eta1)
      uncure1 <- pmin(pmax(uncure1, 1e-6), 1 - 1e-6)
      
      theta1 <- -log(1 - uncure1)
      theta1 <- pmax(theta1, 1e-10)
      
      out1 <- sum(M.train * log(theta1)) +
        sum(log(1 - uncure1))
      
      if (!is.finite(out1)) return(1e20)
      
      return(-out1)
    }
    
    update_b <- optim(
      par = b,
      fn = Q1,
      method = "Nelder-Mead"
    )$par
    
    ## =========================
    ## M-step: Weibull latency part
    ## =========================
    
    Q2 <- function(par) {
      alpha_new <- par[1]
      beta_new <- par[-1]
      
      if (alpha_new <= 0) return(1e20)
      
      eta2 <- drop(X.train[, -1, drop = FALSE] %*% beta_new)
      
      HY <- (TimeY.train^alpha_new) * exp(eta2)
      S.Y <- exp(-HY)
      S.Y <- pmin(pmax(S.Y, 1e-12), 1)
      
      out1 <- sum((M.train - Status.train) * log(S.Y))
      
      event_idx <- which(Status.train == 1)
      
      if (length(event_idx) > 0) {
        HL <- (TimeL.train[event_idx]^alpha_new) * exp(eta2[event_idx])
        HR <- (TimeR.train[event_idx]^alpha_new) * exp(eta2[event_idx])
        
        S.L <- exp(-HL)
        S.R <- exp(-HR)
        
        diff_SR <- S.L - S.R
        diff_SR <- pmax(diff_SR, 1e-12)
        
        out2 <- sum(log(diff_SR))
      } else {
        out2 <- 0
      }
      
      out <- out1 + out2
      
      if (!is.finite(out)) return(1e20)
      
      return(-out)
    }
    
    update_latency <- optim(
      par = c(alpha, beta),
      fn = Q2,
      method = "Nelder-Mead"
    )$par
    
    convergence <- sum(c(update_b - b, update_latency - c(alpha, beta))^2)
    
    b <- update_b
    alpha <- update_latency[1]
    beta <- update_latency[-1]
    
    iter <- iter + 1
  }
  
  ## =========================
  ## Training predictions
  ## =========================
  
  TimeK.train <- ifelse(is.infinite(TimeR.train), TimeL.train, TimeR.train)
  Time.mid.train <- (TimeL.train + TimeK.train) / 2
  
  eta_inc.train <- drop(Z.train %*% b)
  uncureprob.train <- plogis(eta_inc.train)
  uncureprob.train <- pmin(pmax(uncureprob.train, .Machine$double.eps),
                           1 - .Machine$double.eps)
  
  eta_lat.train <- drop(X.train[, -1, drop = FALSE] %*% beta)
  H.mid.train <- (Time.mid.train^alpha) * exp(eta_lat.train)
  S.train <- exp(-H.mid.train)
  
  Sp.train <- (1 - uncureprob.train)^(1 - S.train)
  S1.train <- (Sp.train - (1 - uncureprob.train)) / uncureprob.train
  
  ## =========================
  ## Test predictions
  ## =========================
  
  TimeK.test <- ifelse(is.infinite(TimeR.test), TimeL.test, TimeR.test)
  Time.mid.test <- (TimeL.test + TimeK.test) / 2
  
  eta_inc.test <- drop(Z.test %*% b)
  uncureprob.test <- plogis(eta_inc.test)
  uncureprob.test <- pmin(pmax(uncureprob.test, .Machine$double.eps),
                          1 - .Machine$double.eps)
  
  eta_lat.test <- drop(X.test[, -1, drop = FALSE] %*% beta)
  H.mid.test <- (Time.mid.test^alpha) * exp(eta_lat.test)
  S.test <- exp(-H.mid.test)
  
  Sp.test <- (1 - uncureprob.test)^(1 - S.test)
  S1.test <- (Sp.test - (1 - uncureprob.test)) / uncureprob.test
  
  return(list(
    b = b,
    latencyfit = c(alpha, beta),
    
    train = list(
      Time.mid = Time.mid.train,
      Uncureprob = uncureprob.train,
      Sp = Sp.train,
      S1 = S1.train,
      S = S.train
    ),
    
    test = list(
      Time.mid = Time.mid.test,
      Uncureprob = uncureprob.test,
      Sp = Sp.test,
      S1 = S1.test,
      S = S.test
    ),
    
    tau = convergence,
    iter = iter
  ))
}
smcure.spline.ic <- function(formula, cureform,
                             alpha.init, beta.init,
                             train_data, test_data,
                             df_spline = 1,
                             offset = NULL,
                             na.action = na.omit,
                             Var = F,
                             emmax = 1000,
                             eps = 1e-3,
                             nboot = 100) {
  
  call <- match.call()
  cat("Program is running..be patient...")
  
  train_data <- na.action(train_data)
  test_data  <- na.action(test_data)
  
  n <- nrow(train_data)
  m <- nrow(test_data)
  
  mf.train <- model.frame(formula, train_data)
  mf.test  <- model.frame(formula, test_data)
  
  cvars <- all.vars(cureform)
  
  ## =========================
  ## Build spline incidence design matrices
  ## =========================
  
  spline_train_list <- list()
  spline_test_list  <- list()
  
  for (v in cvars) {
    basis_train <- splines::ns(train_data[[v]], df = df_spline)
    basis_test  <- predict(basis_train, test_data[[v]])
    
    colnames(basis_train) <- paste0(v, "_s", seq_len(ncol(basis_train)))
    colnames(basis_test)  <- paste0(v, "_s", seq_len(ncol(basis_test)))
    
    spline_train_list[[v]] <- basis_train
    spline_test_list[[v]]  <- basis_test
  }
  
  Z.train <- as.matrix(cbind(1, do.call(cbind, spline_train_list)))
  Z.test  <- as.matrix(cbind(1, do.call(cbind, spline_test_list)))
  
  colnames(Z.train)[1] <- "(Intercept)"
  colnames(Z.test)[1]  <- "(Intercept)"
  
  if (!is.null(offset)) {
    offsetvar <- all.vars(offset)
    offsetvar <- train_data[, offsetvar]
  } else {
    offsetvar <- NULL
  }
  
  Y.train <- model.extract(mf.train, "response")
  Y.test  <- model.extract(mf.test, "response")
  
  X.train <- model.matrix(attr(mf.train, "terms"), mf.train)
  X.test  <- model.matrix(attr(mf.test, "terms"), mf.test)
  
  if (!inherits(Y.train, "Surv")) stop("Training response must be a survival object")
  if (!inherits(Y.test, "Surv")) stop("Test response must be a survival object")
  
  TimeL.train <- Y.train[, 1]
  TimeR.train <- Y.train[, 2]
  Status.train <- train_data$D
  Status.test  <- test_data$D
  TimeL.test <- Y.test[, 1]
  TimeR.test <- Y.test[, 2]
  
  
  nb <- ncol(Z.train)
  b <- rep(0.5, nb)
  
  betanm <- colnames(X.train)[-1]
  nbeta <- ncol(X.train) - 1
  
  ## =========================
  ## Fit EM
  ## =========================
  
  emfit <- em.spline.ic(
    TimeL.train = TimeL.train,
    TimeR.train = TimeR.train,
    Status.train = Status.train,
    X.train = X.train,
    Z.train = Z.train,
    TimeL.test = TimeL.test,
    TimeR.test = TimeR.test,
    Status.test = Status.test,
    X.test = X.test,
    Z.test = Z.test,
    offsetvar = offsetvar,
    b = b,
    alpha = alpha.init,
    beta = beta.init,
    emmax = emmax,
    eps = eps,
    data = train_data
  )
  
  latency <- emfit$latencyfit
  b.est <- emfit$b
  alpha.est <- latency[1]
  beta.est <- latency[-1]
  conv <- emfit$tau
  
  ## =========================
  ## Bootstrap variance
  ## =========================
  
  if (Var) {
    
    latency_boot <- matrix(NA, nrow = nboot, ncol = nbeta + 1)
    
    tempdata <- data.frame(
      TimeL = TimeL.train,
      TimeR = TimeR.train,
      Status = Status.train,
      X.train,
      train_data[, cvars, drop = FALSE]
    )
    
    data1 <- subset(tempdata, Status == 1)
    data0 <- subset(tempdata, Status == 0)
    
    n1 <- nrow(data1)
    n0 <- nrow(data0)
    
    boot_count <- 1
    
    while (boot_count <= nboot) {
      
      id1 <- sample(1:n1, n1, replace = TRUE)
      id0 <- sample(1:n0, n0, replace = TRUE)
      
      bootdata <- rbind(data1[id1, ], data0[id0, ])
      
      ## rebuild spline basis inside bootstrap
      boot_spline_train_list <- list()
      boot_spline_test_list  <- list()
      
      for (v in cvars) {
        basis_boot <- splines::ns(bootdata[[v]], df = df_spline)
        basis_test <- predict(basis_boot, test_data[[v]])
        
        colnames(basis_boot) <- paste0(v, "_s", seq_len(ncol(basis_boot)))
        colnames(basis_test) <- paste0(v, "_s", seq_len(ncol(basis_test)))
        
        boot_spline_train_list[[v]] <- basis_boot
        boot_spline_test_list[[v]]  <- basis_test
      }
      
      bootZ <- as.matrix(cbind(1, do.call(cbind, boot_spline_train_list)))
      bootZ.test <- as.matrix(cbind(1, do.call(cbind, boot_spline_test_list)))
      
      colnames(bootZ)[1] <- "(Intercept)"
      colnames(bootZ.test)[1] <- "(Intercept)"
      
      bootX <- as.matrix(cbind(1, bootdata[, betanm, drop = FALSE]))
      
      bootfit <- em.spline.ic(
        TimeL.train = bootdata[, "TimeL"],
        TimeR.train = bootdata[, "TimeR"],
        Status.train = bootdata[, "Status"],
        X.train = bootX,
        Z.train = bootZ,
        TimeL.test = TimeL.test,
        TimeR.test = TimeR.test,
        Status.test = Status.test,
        X.test = X.test,
        Z.test = bootZ.test,
        offsetvar = offsetvar,
        b = b.est,
        alpha = alpha.est,
        beta = beta.est,
        emmax = emmax,
        eps = eps,
        data = train_data
      )
      
      if (bootfit$tau < eps) {
        latency_boot[boot_count, ] <- bootfit$latencyfit
        boot_count <- boot_count + 1
      }
    }
    
    latency_var <- apply(latency_boot, 2, var, na.rm = TRUE)
    latency_sd <- sqrt(latency_var)
  }
  
  fit <- list()
  class(fit) <- "smcure.spline.ic"
  
  fit$call <- call
  fit$b <- b.est
  fit$latency <- latency
  fit$tau <- conv
  fit$iter <- emfit$iter
  
  fit$train <- emfit$train
  fit$test <- emfit$test
  
  fit$df_spline <- df_spline
  fit$cure_covariates <- cvars
  
  if (Var) {
    fit$latency_var <- latency_var
    fit$latency_sd <- latency_sd
    fit$latency_zvalue <- fit$latency / latency_sd
    fit$latency_pvalue <- 2 * (1 - pnorm(abs(fit$latency_zvalue)))
  }
  
  cat(" done.\n")
  
  return(fit)
}


#NEURAL NETWORK
em.nn.ic <- function(TimeL.train, TimeR.train, Status.train, X.train, Z.train,
                     TimeL.test, TimeR.test, Status.test, X.test, Z.test,
                     offsetvar, uncureprob, alpha, beta,
                     emmax, eps, data, best_params) {
  
  n <- length(Status.train)
  m <- length(Status.test)
  M.imp <- 5
  
  
  TimeY.train <- ifelse(
    Status.train == 0,
    TimeL.train,
    TimeR.train
  )
  
  ## NN incidence covariates
  NN.train <- as.data.frame(Z.train[, -1, drop = FALSE])
  NN.test  <- as.data.frame(Z.test[, -1, drop = FALSE])
  
  colnames(NN.train) <- paste0("Input_", seq_len(ncol(NN.train)))
  colnames(NN.test)  <- colnames(NN.train)
  
  nn_formula <- as.formula(
    paste("y ~", paste(colnames(NN.train), collapse = " + "))
  )
  
  convergence <- 1000
  iter <- 1
  
  final_models <- vector("list", M.imp)
  
  while (convergence > eps && iter < emmax) {
    
    ## ==================================================
    ## Current Weibull latency survival
    ## ==================================================
    
    eta_lat.train <- drop(
      X.train[, -1, drop = FALSE] %*% beta
    )
    
    H.train <- (TimeY.train^alpha) * exp(eta_lat.train)
    survival.train <- exp(-H.train)
    
    uncureprob <- pmin(
      pmax(uncureprob, 1e-6),
      1 - 1e-6
    )
    
    ## ==================================================
    ## E-step
    ## ==================================================
    
    w <- Status.train +
      (1 - Status.train) *
      (1 - (1 - uncureprob)^survival.train)
    
    w <- pmin(
      pmax(w, 1e-6),
      1 - 1e-6
    )
    
    M.train <- Status.train -
      survival.train * log1p(-uncureprob)
    
    ## ==================================================
    ## M-step: NN incidence component
    ## ==================================================
    
    imputed_labels <- matrix(
      NA_integer_,
      nrow = n,
      ncol = M.imp
    )
    
    for (j in seq_len(M.imp)) {
      imputed_labels[, j] <- rbinom(
        n,
        size = 1,
        prob = w
      )
    }
    
    train_probs <- matrix(
      NA_real_,
      nrow = n,
      ncol = M.imp
    )
    
    current_models <- vector("list", M.imp)
    
    for (j in seq_len(M.imp)) {
      
      yj <- imputed_labels[, j]
      
      ## A binary NN cannot be fitted with one response class
      if (length(unique(yj)) < 2) {
        train_probs[, j] <- uncureprob
        next
      }
      
      datj <- data.frame(
        y = yj,
        NN.train,
        check.names = FALSE
      )
      
      modj <- tryCatch(
        neuralnet::neuralnet(
          formula = nn_formula,
          data = datj,
          hidden = best_params,
          linear.output = FALSE,
          err.fct = "ce",
          lifesign = "none",
          stepmax = 5e4,
          threshold = 0.05
        ),
        error = function(e) NULL
      )
      
      if (is.null(modj)) {
        train_probs[, j] <- uncureprob
        next
      }
      
      predj <- tryCatch(
        as.numeric(
          neuralnet::compute(
            modj,
            covariate = NN.train
          )$net.result
        ),
        error = function(e) rep(NA_real_, n)
      )
      
      bad_pred <- !is.finite(predj)
      
      if (any(bad_pred)) {
        predj[bad_pred] <- uncureprob[bad_pred]
      }
      
      predj <- pmin(
        pmax(predj, 1e-6),
        1 - 1e-6
      )
      
      train_probs[, j] <- predj
      current_models[[j]] <- modj
    }
    
    valid_models <- sum(
      !vapply(current_models, is.null, logical(1))
    )
    
    if (valid_models == 0) {
      stop("All five neural-network fits failed during the EM iteration.")
    }
    
    update_cureprob <- rowMeans(
      train_probs,
      na.rm = TRUE
    )
    
    bad_train <- !is.finite(update_cureprob)
    
    if (any(bad_train)) {
      update_cureprob[bad_train] <- uncureprob[bad_train]
    }
    
    update_cureprob <- pmin(
      pmax(update_cureprob, 1e-6),
      1 - 1e-6
    )
    
    ## ==================================================
    ## M-step: Weibull latency component
    ## ==================================================
    
    Q2 <- function(par) {
      
      alpha_new <- par[1]
      beta_new <- par[-1]
      
      if (!is.finite(alpha_new) || alpha_new <= 0) {
        return(1e20)
      }
      
      eta2 <- drop(
        X.train[, -1, drop = FALSE] %*% beta_new
      )
      
      HY <- (TimeY.train^alpha_new) * exp(eta2)
      
      S.Y <- exp(-HY)
      S.Y <- pmin(pmax(S.Y, 1e-12), 1)
      
      out1 <- sum(
        (M.train - Status.train) * log(S.Y)
      )
      
      event_idx <- which(Status.train == 1)
      
      if (length(event_idx) > 0) {
        
        HL <- (TimeL.train[event_idx]^alpha_new) *
          exp(eta2[event_idx])
        
        HR <- (TimeR.train[event_idx]^alpha_new) *
          exp(eta2[event_idx])
        
        S.L <- exp(-HL)
        S.R <- exp(-HR)
        
        diff_SR <- pmax(
          S.L - S.R,
          1e-12
        )
        
        out2 <- sum(log(diff_SR))
        
      } else {
        out2 <- 0
      }
      
      out <- out1 + out2
      
      if (!is.finite(out)) {
        return(1e20)
      }
      
      return(-out)
    }
    
    optim_latency <- optim(
      par = c(alpha, beta),
      fn = Q2,
      method = "Nelder-Mead"
    )
    
    update_latency <- optim_latency$par
    
    if (any(!is.finite(update_latency))) {
      stop("The Weibull latency update produced non-finite parameters.")
    }
    
    convergence <- mean(
      (update_cureprob - uncureprob)^2
    ) +
      sum(
        (update_latency - c(alpha, beta))^2
      )
    
    uncureprob <- update_cureprob
    alpha <- update_latency[1]
    beta <- update_latency[-1]
    
    final_models <- current_models
    
    iter <- iter + 1
  }
  
  ## ====================================================
  ## Training predictions
  ## ====================================================
  
  Time.mid.train <- ifelse(
    Status.train == 1,
    (TimeL.train + TimeR.train) / 2,
    TimeL.train
  )
  
  eta_lat.train <- drop(
    X.train[, -1, drop = FALSE] %*% beta
  )
  
  H.mid.train <- (Time.mid.train^alpha) *
    exp(eta_lat.train)
  
  S.train <- exp(-H.mid.train)
  
  Sp.train <- (1 - uncureprob)^(1 - S.train)
  
  S1.train <- (
    Sp.train - (1 - uncureprob)
  ) / uncureprob
  
  ## ====================================================
  ## Test incidence predictions
  ## ====================================================
  
  test_probs <- matrix(
    NA_real_,
    nrow = m,
    ncol = M.imp
  )
  
  for (j in seq_len(M.imp)) {
    
    if (is.null(final_models[[j]])) {
      next
    }
    
    pred_test_j <- tryCatch(
      as.numeric(
        neuralnet::compute(
          final_models[[j]],
          covariate = NN.test
        )$net.result
      ),
      error = function(e) rep(NA_real_, m)
    )
    
    pred_test_j[!is.finite(pred_test_j)] <- NA_real_
    
    test_probs[, j] <- pmin(
      pmax(pred_test_j, 1e-6),
      1 - 1e-6
    )
  }
  
  if (all(is.na(test_probs))) {
    stop("No valid final neural-network models are available for test prediction.")
  }
  
  uncureprob.test <- rowMeans(
    test_probs,
    na.rm = TRUE
  )
  
  bad_test <- !is.finite(uncureprob.test)
  
  if (any(bad_test)) {
    uncureprob.test[bad_test] <- mean(
      uncureprob,
      na.rm = TRUE
    )
  }
  
  uncureprob.test <- pmin(
    pmax(uncureprob.test, 1e-6),
    1 - 1e-6
  )
  
  ## ====================================================
  ## Test survival predictions
  ## ====================================================
  
  Time.mid.test <- ifelse(
    Status.test == 1,
    (TimeL.test + TimeR.test) / 2,
    TimeL.test
  )
  
  eta_lat.test <- drop(
    X.test[, -1, drop = FALSE] %*% beta
  )
  
  H.mid.test <- (Time.mid.test^alpha) *
    exp(eta_lat.test)
  
  S.test <- exp(-H.mid.test)
  
  Sp.test <- (1 - uncureprob.test)^(1 - S.test)
  
  S1.test <- (
    Sp.test - (1 - uncureprob.test)
  ) / uncureprob.test
  
  return(list(
    latencyfit = c(alpha, beta),
    
    train = list(
      Time.mid = Time.mid.train,
      Uncureprob = uncureprob,
      Sp = Sp.train,
      S1 = S1.train,
      S = S.train
    ),
    
    test = list(
      Time.mid = Time.mid.test,
      Uncureprob = uncureprob.test,
      Sp = Sp.test,
      S1 = S1.test,
      S = S.test
    ),
    
    tau = convergence,
    iter = iter,
    converged = is.finite(convergence) && convergence <= eps,
    models = final_models,
    best_params = best_params
  ))
}
smcure.nn.ic <- function(formula, cureform,
                         alpha.init, beta.init,
                         train_data, test_data,
                         offset = NULL,
                         na.action = na.omit,
                         Var = FALSE,
                         emmax = 100,
                         eps = 1e-3,
                         nboot = 100) {
  
  call <- match.call()
  cat("Program is running..be patient...")
  
  ## =====================================================
  ## Prepare training and test data
  ## =====================================================
  
  train_data <- na.action(train_data)
  test_data  <- na.action(test_data)
  
  n <- nrow(train_data)
  
  if (!("D" %in% names(train_data))) {
    stop("train_data must contain the status variable D.")
  }
  
  if (!("D" %in% names(test_data))) {
    stop("test_data must contain the status variable D.")
  }
  
  mf.train <- model.frame(formula, data = train_data)
  mf.test  <- model.frame(formula, data = test_data)
  
  Y.train <- model.extract(mf.train, "response")
  Y.test  <- model.extract(mf.test, "response")
  
  if (!inherits(Y.train, "Surv")) {
    stop("Training response must be a Surv object.")
  }
  
  if (!inherits(Y.test, "Surv")) {
    stop("Test response must be a Surv object.")
  }
  
  X.train <- model.matrix(attr(mf.train, "terms"), mf.train)
  X.test  <- model.matrix(attr(mf.test, "terms"), mf.test)
  
  ## Ensure train and test latency matrices have matching columns
  if (!identical(colnames(X.train), colnames(X.test))) {
    stop("Training and test latency design matrices have different columns.")
  }
  
  cvars <- all.vars(cureform)
  
  Z.train <- as.matrix(
    cbind(
      "(Intercept)" = 1,
      train_data[, cvars, drop = FALSE]
    )
  )
  
  Z.test <- as.matrix(
    cbind(
      "(Intercept)" = 1,
      test_data[, cvars, drop = FALSE]
    )
  )
  
  storage.mode(Z.train) <- "double"
  storage.mode(Z.test)  <- "double"
  
  if (!identical(colnames(Z.train), colnames(Z.test))) {
    stop("Training and test incidence matrices have different columns.")
  }
  
  TimeL.train <- Y.train[, 1]
  TimeR.train <- Y.train[, 2]
  Status.train <- as.numeric(train_data$D)
  
  TimeL.test <- Y.test[, 1]
  TimeR.test <- Y.test[, 2]
  Status.test <- as.numeric(test_data$D)
  
  if (any(!Status.train %in% c(0, 1))) {
    stop("Training status D must contain only 0 and 1.")
  }
  
  if (any(!Status.test %in% c(0, 1))) {
    stop("Test status D must contain only 0 and 1.")
  }
  
  if (length(beta.init) != ncol(X.train) - 1) {
    stop(
      "Length of beta.init must equal ",
      ncol(X.train) - 1,
      ", the number of latency covariates."
    )
  }
  
  if (!is.finite(alpha.init) || alpha.init <= 0) {
    stop("alpha.init must be positive and finite.")
  }
  
  offsetvar <- if (!is.null(offset)) {
    train_data[, all.vars(offset), drop = TRUE]
  } else {
    NULL
  }
  
  nbeta <- ncol(X.train) - 1
  
  ## =====================================================
  ## Construct NN incidence data
  ## =====================================================
  
  NN.train <- as.data.frame(
    Z.train[, -1, drop = FALSE]
  )
  
  NN.test <- as.data.frame(
    Z.test[, -1, drop = FALSE]
  )
  
  nn_names <- paste0(
    "Input_",
    seq_len(ncol(NN.train))
  )
  
  colnames(NN.train) <- nn_names
  colnames(NN.test)  <- nn_names
  
  nn_formula <- as.formula(
    paste(
      "y ~",
      paste(nn_names, collapse = " + ")
    )
  )
  
  ## =====================================================
  ## Initial NN response
  ## =====================================================
  
  init_label <- Status.train
  
  if (length(unique(init_label)) < 2) {
    stop(
      "Initial NN labels contain only one class. ",
      "The NN incidence model cannot be fitted."
    )
  }
  
  init_data <- data.frame(
    y = init_label,
    NN.train,
    check.names = FALSE
  )
  
  ## =====================================================
  ## NN architecture and training settings
  ## =====================================================
  
  #hidden_grid <- list(c(2),c(3),c(4),c(5),c(6),c(3, 2),c(4, 2),c(5, 2),c(6, 3))
  
  hidden_grid <-list(c(2, 2), c(3, 2), c(3, 3), c(4, 2), c(4, 3), c(5, 3), c(6, 3))
  nn_stepmax <- 5e4
  
  ## =====================================================
  ## Tune NN architecture once before EM
  ## =====================================================
  
  best_loss.init <- Inf
  best_params.init <- NULL
  best_model.init <- NULL
  
  for (hidden_val in hidden_grid) {
    
    fit_tmp <- tryCatch(
      neuralnet::neuralnet(
        formula = nn_formula,
        data = init_data,
        hidden = hidden_val,
        linear.output = FALSE,
        err.fct = "ce",
        lifesign = "none",
        stepmax = nn_stepmax
      ),
      error = function(e) NULL
    )
    
    if (is.null(fit_tmp)) {
      next
    }
    
    pred_tmp <- tryCatch(
      as.numeric(
        neuralnet::compute(
          fit_tmp,
          covariate = NN.train
        )$net.result
      ),
      error = function(e) rep(NA_real_, n)
    )
    
    if (length(pred_tmp) != n ||
        any(!is.finite(pred_tmp))) {
      next
    }
    
    pred_tmp <- pmin(
      pmax(pred_tmp, 1e-6),
      1 - 1e-6
    )
    
    loss_tmp <- -mean(
      init_label * log(pred_tmp) +
        (1 - init_label) * log1p(-pred_tmp)
    )
    
    if (is.finite(loss_tmp) &&
        loss_tmp < best_loss.init) {
      
      best_loss.init <- loss_tmp
      best_params.init <- hidden_val
      best_model.init <- fit_tmp
    }
  }
  
  
  
  
  
  ## =====================================================
  ## Initial uncured probabilities
  ## =====================================================
  
  if (is.null(best_model.init) || is.null(best_params.init)) {
    stop(
      "All initial NN architecture fits failed. ",
      "Check the incidence predictors or increase nn_stepmax."
    )
  }
  
  best_params <- best_params.init
  mod.init <- best_model.init
  
  
  uncureprob.init <- tryCatch(
    as.numeric(
      neuralnet::compute(
        mod.init,
        covariate = NN.train
      )$net.result
    ),
    error = function(e) rep(NA_real_, n)
  )
  
  if (length(uncureprob.init) != n ||
      any(!is.finite(uncureprob.init))) {
    stop("The initial NN produced invalid uncured probabilities.")
  }
  
  uncureprob.init <- pmin(
    pmax(uncureprob.init, 1e-6),
    1 - 1e-6
  )
  
  ## =====================================================
  ## Run the NN–Weibull interval-censored EM
  ## =====================================================
  
  emfit <- em.nn.ic(
    TimeL.train = TimeL.train,
    TimeR.train = TimeR.train,
    Status.train = Status.train,
    X.train = X.train,
    Z.train = Z.train,
    
    TimeL.test = TimeL.test,
    TimeR.test = TimeR.test,
    Status.test = Status.test,
    X.test = X.test,
    Z.test = Z.test,
    
    offsetvar = offsetvar,
    uncureprob = uncureprob.init,
    alpha = alpha.init,
    beta = beta.init,
    emmax = emmax,
    eps = eps,
    data = train_data,
    best_params = best_params
  )
  
  latency <- emfit$latencyfit
  alpha.est <- latency[1]
  beta.est <- latency[-1]
  
  ## =====================================================
  ## Bootstrap latency variance
  ## =====================================================
  
  if (Var) {
    
    latency_boot <- matrix(
      NA_real_,
      nrow = nboot,
      ncol = nbeta + 1
    )
    
    event_index <- which(Status.train == 1)
    cens_index  <- which(Status.train == 0)
    
    if (length(event_index) == 0 ||
        length(cens_index) == 0) {
      stop(
        "Bootstrap requires both event and censored ",
        "observations in the training data."
      )
    }
    
    boot_count <- 1
    attempt <- 0
    max_attempts <- 10 * nboot
    
    while (boot_count <= nboot &&
           attempt < max_attempts) {
      
      attempt <- attempt + 1
      
      boot_event <- sample(
        event_index,
        size = length(event_index),
        replace = TRUE
      )
      
      boot_cens <- sample(
        cens_index,
        size = length(cens_index),
        replace = TRUE
      )
      
      boot_index <- c(
        boot_event,
        boot_cens
      )
      
      boot_TimeL <- TimeL.train[boot_index]
      boot_TimeR <- TimeR.train[boot_index]
      boot_Status <- Status.train[boot_index]
      
      bootX <- X.train[
        boot_index,
        ,
        drop = FALSE
      ]
      
      bootZ <- Z.train[
        boot_index,
        ,
        drop = FALSE
      ]
      
      bootNN <- as.data.frame(
        bootZ[, -1, drop = FALSE]
      )
      
      colnames(bootNN) <- nn_names
      
      boot_init_data <- data.frame(
        y = boot_Status,
        bootNN,
        check.names = FALSE
      )
      
      boot_initial_model <- tryCatch(
        neuralnet::neuralnet(
          formula = nn_formula,
          data = boot_init_data,
          hidden = best_params,
          linear.output = FALSE,
          err.fct = "ce",
          lifesign = "none",
          stepmax = 5e4
        ),
        error = function(e) NULL
      )
      
      if (is.null(boot_initial_model)) {
        next
      }
      
      boot_uncure.init <- tryCatch(
        as.numeric(
          neuralnet::compute(
            boot_initial_model,
            covariate = bootNN
          )$net.result
        ),
        error = function(e) {
          rep(NA_real_, length(boot_index))
        }
      )
      
      if (length(boot_uncure.init) != length(boot_index) ||
          any(!is.finite(boot_uncure.init))) {
        next
      }
      
      boot_uncure.init <- pmin(
        pmax(boot_uncure.init, 1e-6),
        1 - 1e-6
      )
      
      bootfit <- tryCatch(
        em.nn.ic(
          TimeL.train = boot_TimeL,
          TimeR.train = boot_TimeR,
          Status.train = boot_Status,
          X.train = bootX,
          Z.train = bootZ,
          
          TimeL.test = TimeL.test,
          TimeR.test = TimeR.test,
          Status.test = Status.test,
          X.test = X.test,
          Z.test = Z.test,
          
          offsetvar = offsetvar,
          uncureprob = boot_uncure.init,
          alpha = alpha.est,
          beta = beta.est,
          emmax = emmax,
          eps = eps,
          data = train_data,
          best_params = best_params
        ),
        error = function(e) NULL
      )
      
      if (is.null(bootfit)) {
        next
      }
      
      if (!is.finite(bootfit$tau) ||
          any(!is.finite(bootfit$latencyfit))) {
        next
      }
      
      if (bootfit$tau <= eps) {
        
        latency_boot[boot_count, ] <-
          bootfit$latencyfit
        
        boot_count <- boot_count + 1
      }
    }
    
    valid_boot <- complete.cases(latency_boot)
    nboot_valid <- sum(valid_boot)
    
    if (nboot_valid < 2) {
      warning(
        "Fewer than two valid NN bootstrap fits were obtained. ",
        "Bootstrap standard deviations cannot be computed."
      )
      
      latency_var <- rep(
        NA_real_,
        nbeta + 1
      )
      
      latency_sd <- rep(
        NA_real_,
        nbeta + 1
      )
      
    } else {
      
      latency_var <- apply(
        latency_boot[
          valid_boot,
          ,
          drop = FALSE
        ],
        2,
        var
      )
      
      latency_sd <- sqrt(latency_var)
    }
    
  } else {
    
    latency_boot <- NULL
    latency_var <- NULL
    latency_sd <- NULL
    nboot_valid <- NULL
  }
  
  ## =====================================================
  ## Return fitted object
  ## =====================================================
  
  fit <- list()
  class(fit) <- "smcure.nn.ic"
  
  fit$call <- call
  
  fit$latency <- latency
  fit$alpha <- alpha.est
  fit$beta <- beta.est
  
  fit$tau <- emfit$tau
  fit$iter <- emfit$iter
  fit$converged <- emfit$converged
  
  fit$train <- emfit$train
  fit$test <- emfit$test
  
  fit$models <- emfit$models
  
  fit$best_params <- best_params
  fit$best_hidden <- paste(best_params, collapse = "-")
  fit$best_loss <- best_loss.init
  
  fit$initial_model <- best_model.init
  fit$initial_uncureprob <- uncureprob.init
  
  if (Var) {
    
    fit$latency_boot <- latency_boot
    fit$nboot_valid <- nboot_valid
    
    fit$latency_var <- latency_var
    fit$latency_sd <- latency_sd
    
    fit$latency_zvalue <- fit$latency /
      fit$latency_sd
    
    fit$latency_pvalue <- 2 * (
      1 - pnorm(
        abs(fit$latency_zvalue)
      )
    )
  }
  
  cat(" done.\n")
  
  return(fit)
}







# =============================================================================
# SIMULATION SCRIPT: INTERVAL-CENSORED PCM WITH TRAIN/TEST
# =============================================================================


MC <- 1
g  <- 400   # train size
m  <- 200   # test size


methods <- c("logit")
#methods <- c("logit", "spline", "nn")
settingop <- c(3)

all_settings_results <- list()

true_props <- vector("list", 3)
for (s in settingop) {
  true_props[[s]] <- data.frame(
    run = integer(),
    Cured_Prop = numeric(),
    Censor_Overall = numeric(),
    Censor_Among_Susceptible = numeric(),
    Event_Overall = numeric()
  )
}

get_simulation_params <- function(setting) {
  switch(
    as.character(setting),
    "1" = list(beta = c(1, 0.5), alpha = 0.5, delta = 0.2),
    "2" = list(beta = c(1, 0.5), alpha = 0.5, delta = 0.2),
    "3" = list(beta = c(-1.5, -0.5, -1.3, 0.6, 1.4), alpha = 1.5, delta = 0.2),
    stop("Unsupported setting")
  )
}

safe_auc <- function(response, predictor) {
  if (length(unique(response)) < 2) return(NA_real_)
  as.numeric(pROC::roc(response = factor(response),
                       predictor = as.numeric(predictor),
                       auc = TRUE,
                       quiet = TRUE)$auc)
}

extract_auc <- function(results, model, type) {
  sapply(results, function(x) {
    value <- x$auc[[model]][[type]]
    if (is.null(value)) return(NA_real_)
    value
  })
}

extract_metric <- function(results, model, metric_type, train_test, metric_name) {
  sapply(results, function(x) {
    value <- x[[metric_type]][[model]][[train_test]][[metric_name]]
    if (is.null(value)) return(NA_real_)
    value
  })
}

for (current_setting in settingop) {
  
  cat("\n====================================================\n")
  cat("STARTING SIMULATION FOR SETTING:", current_setting, "\n")
  cat("====================================================\n")
  
  results_per_setting <- list()
  
  for (r in 1:MC) {
    
    cat("\nSimulation run:", r, "for setting:", current_setting, "\n")
    set.seed(250 + r)
    
    params <- get_simulation_params(current_setting)
    
    data_all <- data.pois.ic(
      n = g + m,
      alpha = params$alpha,
      beta = params$beta,
      delta = params$delta,
      setting = current_setting
    )
    

    
    cured_overall <- mean(data_all$J == 0)
    censor_overall <- mean(data_all$D == 0)
    susc_prop <- mean(data_all$J == 1)
    censor_given_J1 <- if (susc_prop > 0) {
      mean(data_all$D == 0 & data_all$J == 1) / susc_prop
    } else {
      NA_real_
    }
    event_overall <- mean(data_all$D == 1)
    
    true_props[[current_setting]] <- rbind(
      true_props[[current_setting]],
      data.frame(
        run = r,
        Cured_Prop = cured_overall,
        Censor_Overall = censor_overall,
        Censor_Among_Susceptible = censor_given_J1,
        Event_Overall = event_overall
      )
    )
    
    split <- caTools::sample.split(data_all$J, SplitRatio = g / (g + m))
    trainingdata <- subset(data_all, split == TRUE)
    test <- subset(data_all, split == FALSE)
    
    
    if (current_setting == 1) {
      surv_formula <- as.formula("Surv(L, R,  type = 'interval2') ~ z1 + z2")
      cure_formula <- as.formula("~ z1 + z2")
    }
    
    if (current_setting == 2) {
      surv_formula <- as.formula("Surv(L, R, type = 'interval2') ~ z1 + z2")
      cure_formula <- as.formula("~ z1 + z2")
    }
    
    if (current_setting == 3) {
      surv_formula <- as.formula("Surv(L, R, type = 'interval2') ~ z1 + z2 + z3 + z4 + z5")
      cure_formula <- as.formula("~ z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8 + z9 + z10")
    }
    
    auc_results <- list()
    bias_results <- list()
    absolute_bias_results <- list()
    mse_results <- list()
    
    for (method_name in methods) {
      
      cat("Running method:", method_name, "\n")
      
      func <- get(paste0("smcure.", method_name, ".ic"))
      
      out <- tryCatch({
        func(
          formula = surv_formula,
          cureform = cure_formula,
          alpha.init = params$alpha,
          beta.init = params$beta,
          train_data = trainingdata,
          test_data = test,
          Var = FALSE,
          emmax = 1000,
          eps = 1e-3
        )
      }, error = function(e) {
        warning(sprintf(
          "Method %s failed in run %d for setting %d: %s",
          method_name, r, current_setting, e$message
        ))
        NULL
      })
      
      
      
      
      
      if (!is.null(out)) {
       
        
        auc_results[[method_name]] <- list(
          trainingdata = safe_auc(trainingdata$J, out$train$Uncureprob),
          test = safe_auc(test$J, out$test$Uncureprob)
        )
        
        bias_test <- list(
          uncure = mean(out$train$Uncureprob - trainingdata$uncure, na.rm = TRUE),
          Sp     = mean(out$train$Sp - trainingdata$Sp, na.rm = TRUE),
          S1     = mean(out$train$S1 - trainingdata$S1, na.rm = TRUE),
          S      = mean(out$train$S - trainingdata$S, na.rm = TRUE)
        )
        
        abs_bias_test <- list(
          uncure = mean(abs(out$train$Uncureprob - trainingdata$uncure), na.rm = TRUE),
          Sp     = mean(abs(out$train$Sp - trainingdata$Sp), na.rm = TRUE),
          S1     = mean(abs(out$train$S1 - trainingdata$S1), na.rm = TRUE),
          S      = mean(abs(out$train$S - trainingdata$S), na.rm = TRUE)
        )
        
        mse_test <- list(
          uncure = mean((out$train$Uncureprob - trainingdata$uncure)^2, na.rm = TRUE),
          Sp     = mean((out$train$Sp - trainingdata$Sp)^2, na.rm = TRUE),
          S1     = mean((out$train$S1 - trainingdata$S1)^2, na.rm = TRUE),
          S      = mean((out$train$S - trainingdata$S)^2, na.rm = TRUE)
        )
        
        bias_results[[method_name]] <- list(testresult = bias_test)
        absolute_bias_results[[method_name]] <- list(testresult = abs_bias_test)
        mse_results[[method_name]] <- list(testresult = mse_test)
      }
    }
    
    results_per_setting[[r]] <- list(
      auc = auc_results,
      bias = bias_results,
      absoluteBias = absolute_bias_results,
      mse = mse_results
    )
  }
  
  all_settings_results[[current_setting]] <- results_per_setting
}

# =============================================================================
# SUMMARIZE RESULTS
# =============================================================================

for (current_setting in settingop) {
  
  cat("\n\n--- RESULTS FOR SETTING", current_setting, "---\n")
  
  results <- all_settings_results[[current_setting]]
  
  results_table <- data.frame(
    Metric = c(
      "Bias : Uncure", "Abs. Bias : Uncure", "MSE : Uncure",
      "Bias : Sp(t)", "Abs. Bias : Sp(t)", "MSE : Sp(t)",
      "Bias : S1(t)", "Abs. Bias : S1(t)", "MSE : S1(t)",
      "Bias : S(t)", "Abs. Bias : S(t)", "MSE : S(t)",
      "AUC (Train)", "AUC (Test)"
    )
  )
  
  for (model in methods) {
    results_table[[toupper(model)]] <- c(
      mean(extract_metric(results, model, "bias", "testresult", "uncure"), na.rm = TRUE),
      mean(extract_metric(results, model, "absoluteBias", "testresult", "uncure"), na.rm = TRUE),
      mean(extract_metric(results, model, "mse", "testresult", "uncure"), na.rm = TRUE),
      
      mean(extract_metric(results, model, "bias", "testresult", "Sp"), na.rm = TRUE),
      mean(extract_metric(results, model, "absoluteBias", "testresult", "Sp"), na.rm = TRUE),
      mean(extract_metric(results, model, "mse", "testresult", "Sp"), na.rm = TRUE),
      
      mean(extract_metric(results, model, "bias", "testresult", "S1"), na.rm = TRUE),
      mean(extract_metric(results, model, "absoluteBias", "testresult", "S1"), na.rm = TRUE),
      mean(extract_metric(results, model, "mse", "testresult", "S1"), na.rm = TRUE),
      
      mean(extract_metric(results, model, "bias", "testresult", "S"), na.rm = TRUE),
      mean(extract_metric(results, model, "absoluteBias", "testresult", "S"), na.rm = TRUE),
      mean(extract_metric(results, model, "mse", "testresult", "S"), na.rm = TRUE),
      
      mean(extract_auc(results, model, "trainingdata"), na.rm = TRUE),
      mean(extract_auc(results, model, "test"), na.rm = TRUE)
    )
  }
  
  for (col_name in toupper(methods)) {
    results_table[[col_name]] <- round(results_table[[col_name]], 4)
  }
  
  print(results_table)
  
  tp <- true_props[[current_setting]]
  
  se <- function(x) {
    sd(x, na.rm = TRUE) / sqrt(length(na.omit(x)))
  }
  
  summary_props <- data.frame(
    Setting = current_setting,
    Cured_Prop_mean = mean(tp$Cured_Prop, na.rm = TRUE),
    Cured_Prop_se = se(tp$Cured_Prop),
    Censor_Overall_mean = mean(tp$Censor_Overall, na.rm = TRUE),
    Censor_Overall_se = se(tp$Censor_Overall),
    Censor_Among_Susceptible_mean = mean(tp$Censor_Among_Susceptible, na.rm = TRUE),
    Censor_Among_Susceptible_se = se(tp$Censor_Among_Susceptible),
    Event_Overall_mean = mean(tp$Event_Overall, na.rm = TRUE),
    Event_Overall_se = se(tp$Event_Overall)
  )
  
  summary_props[, -1] <- lapply(summary_props[, -1], function(x) round(x, 4))
  
  cat("\nTrue generation-time proportions for setting", current_setting, ":\n")
  print(summary_props)
}

time_taken <- proc.time() - ptm
time_taken