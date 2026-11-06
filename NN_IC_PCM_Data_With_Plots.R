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
library(rpart)
library(randomForest)
library(xgboost)
library(neuralnet)


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


#LOGIT
em.LOGIT.ic <- function(TimeL.train, TimeR.train, Status.train, X.train, Z.train,TimeL.test, TimeR.test,X.test, Z.test,b, alpha, beta,emmax = 1000,eps = 1e-3) {
  
  ## =====================================================
  ## Basic quantities
  ## =====================================================
  
  TimeY.train <- ifelse(Status.train == 0, TimeL.train, TimeR.train)
  convergence <- Inf
  iter <- 1
  
  ## =====================================================
  ## EM algorithm
  ## =====================================================
  
  while (convergence > eps && iter < emmax) {
    
    eta_inc.train <- drop(Z.train %*% b)
    uncureprob.train <- plogis(eta_inc.train)
    uncureprob.train <- pmin(pmax(uncureprob.train, 1e-6), 1 - 1e-6)
      eta_lat.train <- drop(X.train[, -1, drop = FALSE] %*% beta)
    H.train <- (TimeY.train^alpha) * exp(eta_lat.train)
    survival.train <- exp(-H.train)
    ## ===================================================
    ## E-step
    ## ===================================================
    
    M.train <- Status.train - survival.train * log1p(-uncureprob.train)
    
    ## ===================================================
    ## M-step: incidence/LOGIT component
    ## ===================================================
    Q1 <- function(par) {
      eta1 <- drop(Z.train %*% par)
      uncure1 <- plogis(eta1)
      uncure1 <- pmin(pmax(uncure1, 1e-6), 1 - 1e-6)
      theta1 <- -log1p(-uncure1)
      theta1 <- pmax(theta1, 1e-10)
      out <- sum(M.train * log(theta1)) + sum(log1p(-uncure1))
      if (!is.finite(out)) {
        return(1e20)
      }
      return(-out)
    }
    update_b <- optim(par = b, fn = Q1, method = "Nelder-Mead")$par
    ## ===================================================
    ## M-step: Weibull latency component
    ## ===================================================
    Q2 <- function(par) {
      alpha_new <- par[1]
      beta_new <- par[-1]
      if (!is.finite(alpha_new) || alpha_new <= 0) {
        return(1e20)
      }
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
        diff_SR <- pmax(S.L - S.R, 1e-12)
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
    update_latency <- optim(par = c(alpha, beta), fn = Q2, method = "Nelder-Mead")$par
    ## ===================================================
    ## Convergence
    ## ===================================================
    convergence <- sum(c(update_b - b, update_latency - c(alpha, beta))^2)
    ## Update parameters
    b <- update_b
    alpha <- update_latency[1]
    beta <- update_latency[-1]
    iter <- iter + 1
  }
  ## =====================================================
  ## Training predictions
  ## =====================================================
  TimeK.train <- ifelse(is.infinite(TimeR.train), TimeL.train, TimeR.train)
  Time.mid.train <- (TimeL.train + TimeK.train) / 2
  eta_inc.train <- drop(Z.train %*% b)
  uncureprob.train <- plogis(eta_inc.train)
  uncureprob.train <- pmin(pmax(uncureprob.train, 1e-6), 1 - 1e-6)
  eta_lat.train <- drop(X.train[, -1, drop = FALSE] %*% beta)
  H.mid.train <- (Time.mid.train^alpha) * exp(eta_lat.train)
  S.train <- exp(-H.mid.train)
  Sp.train <- (1 - uncureprob.train)^(1 - S.train)
  S1.train <- (Sp.train - (1 - uncureprob.train)) / uncureprob.train
  ## =====================================================
  ## Test predictions
  ## =====================================================
  TimeK.test <- ifelse(is.infinite(TimeR.test), TimeL.test, TimeR.test)
  Time.mid.test <- (TimeL.test + TimeK.test) / 2
  eta_inc.test <- drop(Z.test %*% b)
  uncureprob.test <- plogis(eta_inc.test)
  uncureprob.test <- pmin(pmax(uncureprob.test, 1e-6), 1 - 1e-6)
  eta_lat.test <- drop(X.test[, -1, drop = FALSE] %*% beta)
  H.mid.test <- (Time.mid.test^alpha) * exp(eta_lat.test)
  S.test <- exp(-H.mid.test)
  Sp.test <- (1 - uncureprob.test)^(1 - S.test)
  S1.test <- (Sp.test - (1 - uncureprob.test)) / uncureprob.test
  ## =====================================================
  ## Return
  ## =====================================================
  return(list(
    b = b,
    latencyfit = c(alpha, beta),
    train = list(Time.mid = Time.mid.train, Uncureprob = uncureprob.train, Sp = Sp.train, S1 = S1.train, S = S.train),
    test = list(Time.mid = Time.mid.test, Uncureprob = uncureprob.test, Sp = Sp.test, S1 = S1.test, S = S.test),
    tau = convergence,
    iter = iter
  ))
  
}
smcure.LOGIT.ic <- function(formula, cureform, alpha.init, beta.init, train_data, test_data, na.action = na.omit, Var = FALSE, emmax = 1000, eps = 1e-3, nboot = 100) {
  call <- match.call()
  ## =====================================================
  ## Prepare data
  ## =====================================================
  train_data <- na.action(train_data)
  test_data <- na.action(test_data)
  if (!("D" %in% names(train_data))) stop("train_data must contain status column D.")
  if (!("D" %in% names(test_data))) stop("test_data must contain status column D.")
  ## =====================================================
  ## Build latency model matrices: X
  ## =====================================================
  mf.train <- model.frame(formula, data = train_data)
  mf.test <- model.frame(formula, data = test_data)
  Y.train <- model.extract(mf.train, "response")
  Y.test <- model.extract(mf.test, "response")
  if (!inherits(Y.train, "Surv")) stop("Training response must be a survival object.")
  if (!inherits(Y.test, "Surv")) stop("Test response must be a survival object.")
  X.train <- model.matrix(attr(mf.train, "terms"), mf.train)
  X.test <- model.matrix(attr(mf.test, "terms"), mf.test)
  ## =====================================================
  ## Build incidence model matrices: Z
  ## =====================================================
  Z.train <- model.matrix(cureform, data = train_data)
  Z.test <- model.matrix(cureform, data = test_data)
  ## =====================================================
  ## Survival quantities
  ## =====================================================
  TimeL.train <- Y.train[, 1]
  TimeR.train <- Y.train[, 2]
  Status.train <- train_data$D
  TimeL.test <- Y.test[, 1]
  TimeR.test <- Y.test[, 2]
  ## =====================================================
  ## Parameter dimensions
  ## =====================================================
  nb <- ncol(Z.train)
  nbeta <- ncol(X.train) - 1
  if (length(beta.init) != nbeta) stop("Length of beta.init must equal the number of latency covariates: ", nbeta)
  b <- rep(0.5, nb)
  ## =====================================================
  ## Fit EM
  ## =====================================================
  emfit <- em.LOGIT.ic(TimeL.train = TimeL.train, TimeR.train = TimeR.train, Status.train = Status.train, X.train = X.train, Z.train = Z.train, TimeL.test = TimeL.test, TimeR.test = TimeR.test, X.test = X.test, Z.test = Z.test, b = b, alpha = alpha.init, beta = beta.init, emmax = emmax, eps = eps)
  latency <- emfit$latencyfit
  b.est <- emfit$b
  alpha.est <- latency[1]
  beta.est <- latency[-1]
  conv <- emfit$tau
  ## =====================================================
  ## Bootstrap inference
  ## =====================================================
  if (Var) {
    latency_boot <- matrix(NA_real_, nrow = nboot, ncol = nbeta + 1)
    incidence_boot <- matrix(NA_real_, nrow = nboot, ncol = nb)
    colnames(latency_boot) <- c("alpha", colnames(X.train)[-1])
    colnames(incidence_boot) <- colnames(Z.train)
    event_idx <- which(Status.train == 1)
    cens_idx <- which(Status.train == 0)
    n1 <- length(event_idx)
    n0 <- length(cens_idx)
    if (n1 == 0 || n0 == 0) stop("Bootstrap requires both event and censored observations.")
    boot_count <- 1
    attempt <- 0
    max_attempts <- 10 * nboot
    while (boot_count <= nboot && attempt < max_attempts) {
      attempt <- attempt + 1
      id1 <- sample(event_idx, size = n1, replace = TRUE)
      id0 <- sample(cens_idx, size = n0, replace = TRUE)
      boot_id <- c(id1, id0)
      boot_TimeL <- TimeL.train[boot_id]
      boot_TimeR <- TimeR.train[boot_id]
      boot_Status <- Status.train[boot_id]
      bootX <- X.train[boot_id, , drop = FALSE]
      bootZ <- Z.train[boot_id, , drop = FALSE]
      bootfit <- tryCatch(em.LOGIT.ic(TimeL.train = boot_TimeL, TimeR.train = boot_TimeR, Status.train = boot_Status, X.train = bootX, Z.train = bootZ, TimeL.test = TimeL.test, TimeR.test = TimeR.test, X.test = X.test, Z.test = Z.test, b = b.est, alpha = alpha.est, beta = beta.est, emmax = emmax, eps = eps), error = function(e) NULL)
      if (!is.null(bootfit) && is.finite(bootfit$tau) && bootfit$tau < eps && all(is.finite(bootfit$latencyfit)) && all(is.finite(bootfit$b))) {
        latency_boot[boot_count, ] <- bootfit$latencyfit
        incidence_boot[boot_count, ] <- bootfit$b
        boot_count <- boot_count + 1
      }
    }
    nboot_valid <- boot_count - 1
    if (nboot_valid < 2) {
      warning("Fewer than two valid bootstrap fits were obtained.")
      latency_boot <- latency_boot[seq_len(max(nboot_valid, 1)), , drop = FALSE]
      incidence_boot <- incidence_boot[seq_len(max(nboot_valid, 1)), , drop = FALSE]
      latency_var <- rep(NA_real_, nbeta + 1)
      latency_sd <- rep(NA_real_, nbeta + 1)
    } else {
      latency_boot <- latency_boot[seq_len(nboot_valid), , drop = FALSE]
      incidence_boot <- incidence_boot[seq_len(nboot_valid), , drop = FALSE]
      latency_var <- apply(latency_boot, 2, var)
      latency_sd <- sqrt(latency_var)
    }
  }
  ## =====================================================
  ## Return fitted object
  ## =====================================================
  fit <- list()
  class(fit) <- "smcure.LOGIT.ic"
  fit$call <- call
  fit$b <- b.est
  fit$latency <- latency
  fit$tau <- conv
  fit$iter <- emfit$iter
  fit$train <- emfit$train
  fit$test <- emfit$test
  fit$incidence_names <- colnames(Z.train)
  fit$latency_names <- c("alpha", colnames(X.train)[-1])
  if (Var) {
    fit$latency_boot <- latency_boot
    fit$incidence_boot <- incidence_boot
    fit$nboot_valid <- nboot_valid
    fit$latency_var <- latency_var
    fit$latency_sd <- latency_sd
    fit$latency_zvalue <- fit$latency / latency_sd
    fit$latency_pvalue <- 2 * (1 - pnorm(abs(fit$latency_zvalue)))
  }
  return(fit)
}
#SPLINE
em.SPLINE.ic <- function(TimeL.train, TimeR.train, Status.train,X.train, Z.train,TimeL.test, TimeR.test,X.test, Z.test,b, alpha, beta,emmax = 1000, eps = 1e-3) {
  
  ## =====================================================
  ## Basic quantities
  ## =====================================================
  TimeY.train <- ifelse(Status.train == 0, TimeL.train, TimeR.train)
  convergence <- Inf
  iter <- 1
  ## =====================================================
  ## EM algorithm
  ## =====================================================
  while (convergence > eps && iter < emmax) {
    # ===================================================
    eta_inc.train <- drop(Z.train %*% b)
    uncureprob.train <- plogis(eta_inc.train)
    uncureprob.train <- pmin(pmax(uncureprob.train, 1e-6), 1 - 1e-6)
    eta_lat.train <- drop(X.train[, -1, drop = FALSE] %*% beta)
    H.train <- (TimeY.train^alpha) * exp(eta_lat.train)
    survival.train <- exp(-H.train)
    ## ===================================================
    ## E-step
    ## ===================================================
    M.train <- Status.train - survival.train * log1p(-uncureprob.train)
    ## ===================================================
    ## M-step: SPLINE incidence component
    ## ===================================================
    Q1 <- function(par) {
      eta1 <- drop(Z.train %*% par)
      uncure1 <- plogis(eta1)
      uncure1 <- pmin(pmax(uncure1, 1e-6), 1 - 1e-6)
      theta1 <- -log1p(-uncure1)
      theta1 <- pmax(theta1, 1e-10)
      out <- sum(M.train * log(theta1)) + sum(log1p(-uncure1))
      if (!is.finite(out)) {
        return(1e20)
      }
      return(-out)
    }
    update_b <- optim(par = b, fn = Q1, method = "Nelder-Mead")$par
    ## ===================================================
    ## M-step: Weibull latency component
    ## ===================================================
    Q2 <- function(par) {
      alpha_new <- par[1]
      beta_new <- par[-1]
      if (!is.finite(alpha_new) || alpha_new <= 0) {
        return(1e20)
      }
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
        diff_SR <- pmax(S.L - S.R, 1e-12)
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
    update_latency <- optim(par = c(alpha, beta), fn = Q2, method = "Nelder-Mead")$par
    ## ===================================================
    ## Convergence
    ## ===================================================
    convergence <- sum(c(update_b - b, update_latency - c(alpha, beta))^2)
    ## ===================================================
    ## Update parameters
    ## ===================================================
    b <- update_b
    alpha <- update_latency[1]
    beta <- update_latency[-1]
    iter <- iter + 1
  }
  ## =====================================================
  ## Training predictions
  ## =====================================================
  TimeK.train <- ifelse(is.infinite(TimeR.train), TimeL.train, TimeR.train)
  Time.mid.train <- (TimeL.train + TimeK.train) / 2
  eta_inc.train <- drop(Z.train %*% b)
  uncureprob.train <- plogis(eta_inc.train)
  uncureprob.train <- pmin(pmax(uncureprob.train, 1e-6), 1 - 1e-6)
  eta_lat.train <- drop(X.train[, -1, drop = FALSE] %*% beta)
  H.mid.train <- (Time.mid.train^alpha) * exp(eta_lat.train)
  S.train <- exp(-H.mid.train)
  Sp.train <- (1 - uncureprob.train)^(1 - S.train)
  S1.train <- (Sp.train - (1 - uncureprob.train)) / uncureprob.train
  ## =====================================================
  ## Test predictions
  ## =====================================================
  TimeK.test <- ifelse(is.infinite(TimeR.test), TimeL.test, TimeR.test)
  Time.mid.test <- (TimeL.test + TimeK.test) / 2
  eta_inc.test <- drop(Z.test %*% b)
  uncureprob.test <- plogis(eta_inc.test)
  uncureprob.test <- pmin(pmax(uncureprob.test, 1e-6), 1 - 1e-6)
  eta_lat.test <- drop(X.test[, -1, drop = FALSE] %*% beta)
  H.mid.test <- (Time.mid.test^alpha) * exp(eta_lat.test)
  S.test <- exp(-H.mid.test)
  Sp.test <- (1 - uncureprob.test)^(1 - S.test)
  S1.test <- (Sp.test - (1 - uncureprob.test)) / uncureprob.test
  ## =====================================================
  ## Return
  ## =====================================================
  return(list(
    b = b,
    latencyfit = c(alpha, beta),
    train = list(Time.mid = Time.mid.train, Uncureprob = uncureprob.train, Sp = Sp.train, S1 = S1.train, S = S.train),
    test = list(Time.mid = Time.mid.test, Uncureprob = uncureprob.test, Sp = Sp.test, S1 = S1.test, S = S.test),
    tau = convergence,
    iter = iter
  ))
}

smcure.SPLINE.ic <- function(formula, cureform, alpha.init, beta.init, train_data, test_data, df_SPLINE = 2, na.action = na.omit, Var = TRUE, emmax = 1000, eps = 1e-3, nboot = 100) {
  call <- match.call()
  ## =====================================================
  ## Prepare data
  ## =====================================================
  train_data <- na.action(train_data)
  test_data <- na.action(test_data)
  if (!("D" %in% names(train_data))) stop("train_data must contain status column D.")
  if (!("D" %in% names(test_data))) stop("test_data must contain status column D.")
  ## =====================================================
  ## Latency model: X
  ## =====================================================
  mf.train <- model.frame(formula, data = train_data)
  mf.test <- model.frame(formula, data = test_data)
  Y.train <- model.extract(mf.train, "response")
  Y.test <- model.extract(mf.test, "response")
  if (!inherits(Y.train, "Surv")) stop("Training response must be a survival object.")
  if (!inherits(Y.test, "Surv")) stop("Test response must be a survival object.")
  X.train <- model.matrix(attr(mf.train, "terms"), mf.train)
  X.test <- model.matrix(attr(mf.test, "terms"), mf.test)
  ## =====================================================
  ## Incidence variables
  ## =====================================================
  cvars <- all.vars(cureform)
  ## =====================================================
  ## Build SPLINE incidence matrices: Z
  ## =====================================================
  Z.train.parts <- list()
  Z.test.parts <- list()
  SPLINE_bases <- list()
  for (v in cvars) {
    x.train <- train_data[[v]]
    x.test <- test_data[[v]]
    if (is.numeric(x.train) && length(unique(x.train)) > 2) {
      basis.train <- splines::ns(x.train, df = df_SPLINE)
      basis.test <- predict(basis.train, x.test)
      colnames(basis.train) <- paste0(v, "_s", seq_len(ncol(basis.train)))
      colnames(basis.test) <- colnames(basis.train)
      Z.train.parts[[v]] <- basis.train
      Z.test.parts[[v]] <- basis.test
      SPLINE_bases[[v]] <- basis.train
    } else {
      temp_formula <- as.formula(paste0("~ ", v))
      temp.train <- model.matrix(temp_formula, data = train_data)[, -1, drop = FALSE]
      temp.test <- model.matrix(temp_formula, data = test_data)[, -1, drop = FALSE]
      Z.train.parts[[v]] <- temp.train
      Z.test.parts[[v]] <- temp.test
    }
  }
  Z.train <- as.matrix(cbind("(Intercept)" = 1, do.call(cbind, Z.train.parts)))
  Z.test <- as.matrix(cbind("(Intercept)" = 1, do.call(cbind, Z.test.parts)))
  ## =====================================================
  ## Survival quantities
  ## =====================================================
  TimeL.train <- Y.train[, 1]
  TimeR.train <- Y.train[, 2]
  Status.train <- train_data$D
  TimeL.test <- Y.test[, 1]
  TimeR.test <- Y.test[, 2]
  ## =====================================================
  ## Parameter dimensions
  ## =====================================================
  nb <- ncol(Z.train)
  nbeta <- ncol(X.train) - 1
  if (length(beta.init) != nbeta) stop("Length of beta.init must equal the number of latency covariates: ", nbeta)
  b <- rep(0.5, nb)
  ## =====================================================
  ## Fit EM
  ## =====================================================
  emfit <- em.SPLINE.ic(TimeL.train = TimeL.train, TimeR.train = TimeR.train, Status.train = Status.train, X.train = X.train, Z.train = Z.train, TimeL.test = TimeL.test, TimeR.test = TimeR.test, X.test = X.test, Z.test = Z.test, b = b, alpha = alpha.init, beta = beta.init, emmax = emmax, eps = eps)
  latency <- emfit$latencyfit
  b.est <- emfit$b
  alpha.est <- latency[1]
  beta.est <- latency[-1]
  conv <- emfit$tau
  ## =====================================================
  ## Bootstrap inference
  ## =====================================================
  if (Var) {
    latency_boot <- matrix(NA_real_, nrow = nboot, ncol = nbeta + 1)
    incidence_boot <- matrix(NA_real_, nrow = nboot, ncol = nb)
    SPLINE_boot <- vector("list", nboot)
    colnames(latency_boot) <- c("alpha", colnames(X.train)[-1])
    colnames(incidence_boot) <- colnames(Z.train)
    event_idx <- which(Status.train == 1)
    cens_idx <- which(Status.train == 0)
    n1 <- length(event_idx)
    n0 <- length(cens_idx)
    if (n1 == 0 || n0 == 0) stop("Bootstrap requires both event and censored observations.")
    boot_count <- 1
    attempt <- 0
    max_attempts <- 10 * nboot
    while (boot_count <= nboot && attempt < max_attempts) {
      attempt <- attempt + 1
      ## ===================================================
      ## Stratified bootstrap
      ## ===================================================
      id1 <- sample(event_idx, size = n1, replace = TRUE)
      id0 <- sample(cens_idx, size = n0, replace = TRUE)
      boot_id <- c(id1, id0)
      bootdata <- train_data[boot_id, , drop = FALSE]
      boot_TimeL <- TimeL.train[boot_id]
      boot_TimeR <- TimeR.train[boot_id]
      boot_Status <- Status.train[boot_id]
      ## ===================================================
      ## Bootstrap latency matrix
      ## ===================================================
      boot_mf <- model.frame(formula, data = bootdata)
      bootX <- model.matrix(attr(boot_mf, "terms"), boot_mf)
      ## ===================================================
      ## Bootstrap SPLINE incidence matrices
      ## ===================================================
      boot.Z.parts <- list()
      test.Z.parts <- list()
      boot.SPLINE.bases <- list()
      for (v in cvars) {
        x.boot <- bootdata[[v]]
        x.test <- test_data[[v]]
        if (is.numeric(x.boot) && length(unique(x.boot)) > 2) {
          basis.boot <- splines::ns(x.boot, df = df_SPLINE)
          basis.test <- predict(basis.boot, x.test)
          colnames(basis.boot) <- paste0(v, "_s", seq_len(ncol(basis.boot)))
          colnames(basis.test) <- colnames(basis.boot)
          boot.Z.parts[[v]] <- basis.boot
          test.Z.parts[[v]] <- basis.test
          boot.SPLINE.bases[[v]] <- basis.boot
        } else {
          temp_formula <- as.formula(paste0("~ ", v))
          temp.boot <- model.matrix(temp_formula, data = bootdata)[, -1, drop = FALSE]
          temp.test <- model.matrix(temp_formula, data = test_data)[, -1, drop = FALSE]
          boot.Z.parts[[v]] <- temp.boot
          test.Z.parts[[v]] <- temp.test
        }
      }
      bootZ <- as.matrix(cbind("(Intercept)" = 1, do.call(cbind, boot.Z.parts)))
      bootZ.test <- as.matrix(cbind("(Intercept)" = 1, do.call(cbind, test.Z.parts)))
      ## ===================================================
      ## Bootstrap EM
      ## ===================================================
      bootfit <- tryCatch(em.SPLINE.ic(TimeL.train = boot_TimeL, TimeR.train = boot_TimeR, Status.train = boot_Status, X.train = bootX, Z.train = bootZ, TimeL.test = TimeL.test, TimeR.test = TimeR.test, X.test = X.test, Z.test = bootZ.test, b = b.est, alpha = alpha.est, beta = beta.est, emmax = emmax, eps = eps), error = function(e) NULL)
      ## ===================================================
      ## Keep successful bootstrap fits
      ## ===================================================
      if (!is.null(bootfit) && is.finite(bootfit$tau) && bootfit$tau < eps && all(is.finite(bootfit$latencyfit)) && all(is.finite(bootfit$b))) {
        latency_boot[boot_count, ] <- bootfit$latencyfit
        incidence_boot[boot_count, ] <- bootfit$b
        SPLINE_boot[[boot_count]] <- boot.SPLINE.bases
        boot_count <- boot_count + 1
      }
    }
    ## ===================================================
    ## Keep successful bootstrap fits only
    ## ===================================================
    nboot_valid <- boot_count - 1
    if (nboot_valid < 2) {
      warning("Fewer than two valid bootstrap fits were obtained.")
      latency_boot <- latency_boot[seq_len(max(nboot_valid, 1)), , drop = FALSE]
      incidence_boot <- incidence_boot[seq_len(max(nboot_valid, 1)), , drop = FALSE]
      SPLINE_boot <- SPLINE_boot[seq_len(max(nboot_valid, 1))]
      latency_var <- rep(NA_real_, nbeta + 1)
      latency_sd <- rep(NA_real_, nbeta + 1)
    } else {
      latency_boot <- latency_boot[seq_len(nboot_valid), , drop = FALSE]
      incidence_boot <- incidence_boot[seq_len(nboot_valid), , drop = FALSE]
      SPLINE_boot <- SPLINE_boot[seq_len(nboot_valid)]
      latency_var <- apply(latency_boot, 2, var)
      latency_sd <- sqrt(latency_var)
    }
  }
  ## =====================================================
  ## Return fitted object
  ## =====================================================
  fit <- list()
  class(fit) <- "smcure.SPLINE.ic"
  fit$call <- call
  fit$b <- b.est
  fit$latency <- latency
  fit$tau <- conv
  fit$iter <- emfit$iter
  fit$train <- emfit$train
  fit$test <- emfit$test
  fit$df_SPLINE <- df_SPLINE
  fit$cure_covariates <- cvars
  fit$SPLINE_bases <- SPLINE_bases
  fit$incidence_names <- colnames(Z.train)
  fit$latency_names <- c("alpha", colnames(X.train)[-1])
  if (Var) {
    fit$latency_boot <- latency_boot
    fit$incidence_boot <- incidence_boot
    fit$SPLINE_boot <- SPLINE_boot
    fit$nboot_valid <- nboot_valid
    fit$latency_var <- latency_var
    fit$latency_sd <- latency_sd
    fit$latency_zvalue <- fit$latency / latency_sd
    fit$latency_pvalue <- 2 * (1 - pnorm(abs(fit$latency_zvalue)))
  }
  return(fit)
}







#NEURAL NETWORK
em.NN.ic <- function(TimeL.train, TimeR.train, Status.train,
                     X.train, Z.train,
                     TimeL.test, TimeR.test,
                     X.test, Z.test,
                     uncureprob, alpha, beta,
                     best_params,
                     emmax = 1000, eps = 1e-3, M.imp = 5) {
  ## =====================================================
  ## Basic quantities
  ## =====================================================
  n <- length(Status.train)
  m <- nrow(Z.test)
  TimeY.train <- ifelse(Status.train == 0, TimeL.train, TimeR.train)
  NN.train <- as.data.frame(Z.train[, -1, drop = FALSE])
  NN.test <- as.data.frame(Z.test[, -1, drop = FALSE])
  colnames(NN.train) <- paste0("Input_", seq_len(ncol(NN.train)))
  colnames(NN.test) <- colnames(NN.train)
  NN_formula <- as.formula(paste("y ~", paste(colnames(NN.train), collapse = " + ")))
  convergence <- Inf
  iter <- 1
  final_models <- vector("list", M.imp)
  ## =====================================================
  ## EM algorithm
  ## =====================================================
  while (convergence > eps && iter < emmax) {
    ## ===================================================
    ## Current Weibull latency survival
    ## ===================================================
    eta_lat.train <- drop(X.train[, -1, drop = FALSE] %*% beta)
    H.train <- (TimeY.train^alpha) * exp(eta_lat.train)
    survival.train <- exp(-H.train)
    uncureprob <- pmin(pmax(uncureprob, 1e-6), 1 - 1e-6)
    ## ===================================================
    ## E-step
    ## ===================================================
    w <- Status.train + (1 - Status.train) * (1 - (1 - uncureprob)^survival.train)
    w <- pmin(pmax(w, 1e-6), 1 - 1e-6)
    M.train <- Status.train - survival.train * log1p(-uncureprob)
    ## ===================================================
    ## M-step: NN incidence component
    ## ===================================================
    train_probs <- matrix(NA_real_, nrow = n, ncol = M.imp)
    current_models <- vector("list", M.imp)
    for (j in seq_len(M.imp)) {
      imputed_label <- rbinom(n = n, size = 1, prob = w)
      if (length(unique(imputed_label)) < 2) {
        train_probs[, j] <- uncureprob
        next
      }
      datj <- data.frame(y = imputed_label, NN.train, check.names = FALSE)
      modj <- tryCatch(
        neuralnet::neuralnet(
          formula = NN_formula,
          data = datj,
          hidden = best_params,
          linear.output = FALSE,
          err.fct = "ce",
          lifesign = "none",
          stepmax = 5e4,
          threshold = 0.02
        ),
        error = function(e) {
          NULL
        }
      )
      if (is.null(modj)) {
        train_probs[, j] <- uncureprob
        next
      }
      predj <- tryCatch(
        as.numeric(neuralnet::compute(modj, covariate = NN.train)$net.result),
        error = function(e) {
          rep(NA_real_, n)
        }
      )
      if (length(predj) != n) {
        train_probs[, j] <- uncureprob
        next
      }
      bad_pred <- !is.finite(predj)
      if (any(bad_pred)) {
        predj[bad_pred] <- uncureprob[bad_pred]
      }
      predj <- pmin(pmax(predj, 1e-6), 1 - 1e-6)
      train_probs[, j] <- predj
      current_models[[j]] <- modj
    }
    valid_models <- sum(!vapply(current_models, is.null, logical(1)))
    if (valid_models == 0) {
      stop("All neural-network fits failed during the EM iteration.")
    }
    update_uncureprob <- rowMeans(train_probs, na.rm = TRUE)
    bad_train <- !is.finite(update_uncureprob)
    if (any(bad_train)) {
      update_uncureprob[bad_train] <- uncureprob[bad_train]
    }
    update_uncureprob <- pmin(pmax(update_uncureprob, 1e-6), 1 - 1e-6)
    ## ===================================================
    ## M-step: Weibull latency component
    ## ===================================================
    Q2 <- function(par) {
      alpha_new <- par[1]
      beta_new <- par[-1]
      if (!is.finite(alpha_new) || alpha_new <= 0) {
        return(1e20)
      }
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
        diff_SR <- pmax(S.L - S.R, 1e-12)
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
    update_latency <- optim(par = c(alpha, beta), fn = Q2, method = "Nelder-Mead")$par
    if (any(!is.finite(update_latency))) {
      stop("The Weibull latency update produced non-finite parameters.")
    }
    ## ===================================================
    ## Convergence
    ## ===================================================
    convergence <- sum(c(mean(update_uncureprob) - mean(uncureprob), update_latency - c(alpha, beta))^2)
    uncureprob <- update_uncureprob
    alpha <- update_latency[1]
    beta <- update_latency[-1]
    final_models <- current_models
    iter <- iter + 1
  }
  ## =====================================================
  ## Training predictions
  ## =====================================================
  TimeK.train <- ifelse(is.infinite(TimeR.train), TimeL.train, TimeR.train)
  Time.mid.train <- (TimeL.train + TimeK.train) / 2
  eta_lat.train <- drop(X.train[, -1, drop = FALSE] %*% beta)
  H.mid.train <- (Time.mid.train^alpha) * exp(eta_lat.train)
  S.train <- exp(-H.mid.train)
  Sp.train <- (1 - uncureprob)^(1 - S.train)
  S1.train <- (Sp.train - (1 - uncureprob)) / uncureprob
  ## =====================================================
  ## Test incidence predictions
  ## =====================================================
  test_probs <- matrix(NA_real_, nrow = m, ncol = M.imp)
  for (j in seq_len(M.imp)) {
    if (is.null(final_models[[j]])) {
      next
    }
    pred_test_j <- tryCatch(
      as.numeric(neuralnet::compute(final_models[[j]], covariate = NN.test)$net.result),
      error = function(e) {
        rep(NA_real_, m)
      }
    )
    if (length(pred_test_j) != m) {
      next
    }
    pred_test_j[!is.finite(pred_test_j)] <- NA_real_
    test_probs[, j] <- pmin(pmax(pred_test_j, 1e-6), 1 - 1e-6)
  }
  if (all(is.na(test_probs))) {
    stop("No valid final neural-network models are available for test prediction.")
  }
  uncureprob.test <- rowMeans(test_probs, na.rm = TRUE)
  bad_test <- !is.finite(uncureprob.test)
  if (any(bad_test)) {
    uncureprob.test[bad_test] <- mean(uncureprob, na.rm = TRUE)
  }
  uncureprob.test <- pmin(pmax(uncureprob.test, 1e-6), 1 - 1e-6)
  ## =====================================================
  ## Test survival predictions
  ## =====================================================
  TimeK.test <- ifelse(is.infinite(TimeR.test), TimeL.test, TimeR.test)
  Time.mid.test <- (TimeL.test + TimeK.test) / 2
  eta_lat.test <- drop(X.test[, -1, drop = FALSE] %*% beta)
  H.mid.test <- (Time.mid.test^alpha) * exp(eta_lat.test)
  S.test <- exp(-H.mid.test)
  Sp.test <- (1 - uncureprob.test)^(1 - S.test)
  S1.test <- (Sp.test - (1 - uncureprob.test)) / uncureprob.test
  ## =====================================================
  ## Return
  ## =====================================================
  return(list(
    latencyfit = c(alpha, beta),
    train = list(Time.mid = Time.mid.train, Uncureprob = uncureprob, Sp = Sp.train, S1 = S1.train, S = S.train),
    test = list(Time.mid = Time.mid.test, Uncureprob = uncureprob.test, Sp = Sp.test, S1 = S1.test, S = S.test),
    tau = convergence,
    iter = iter,
    converged = is.finite(convergence) && convergence <= eps,
    models = final_models,
    best_params = best_params
  ))
}
smcure.NN.ic <- function(formula, cureform,
                         alpha.init, beta.init,
                         train_data, test_data,
                         na.action = na.omit,
                         Var = FALSE,
                         emmax = 1000,
                         eps = 1e-3,
                         nboot = 100,
                         hidden_grid = list(c(2), c(3), c(4), c(5), c(6), c(8), c(10), c(3, 2), c(4, 2), c(5, 2), c(5, 3), c(6, 3), c(8, 4)),
                         M.imp = 5,
                         NN_stepmax = 5e4,
                         NN_threshold = 0.02) {
  call <- match.call()
  ## =====================================================
  ## Prepare training and test data
  ## =====================================================
  train_data <- na.action(train_data)
  test_data <- na.action(test_data)
  if (!("D" %in% names(train_data))) {
    stop("train_data must contain status column D.")
  }
  if (!("D" %in% names(test_data))) {
    stop("test_data must contain status column D.")
  }
  ## =====================================================
  ## Build latency model matrices: X
  ## =====================================================
  mf.train <- model.frame(formula, data = train_data)
  mf.test <- model.frame(formula, data = test_data)
  Y.train <- model.extract(mf.train, "response")
  Y.test <- model.extract(mf.test, "response")
  if (!inherits(Y.train, "Surv")) {
    stop("Training response must be a survival object.")
  }
  if (!inherits(Y.test, "Surv")) {
    stop("Test response must be a survival object.")
  }
  X.train <- model.matrix(attr(mf.train, "terms"), mf.train)
  X.test <- model.matrix(attr(mf.test, "terms"), mf.test)
  if (!identical(colnames(X.train), colnames(X.test))) {
    stop("Training and test latency design matrices have different columns.")
  }
  ## =====================================================
  ## Build incidence model matrices: Z
  ## =====================================================
  Z.train <- model.matrix(cureform, data = train_data)
  Z.test <- model.matrix(cureform, data = test_data)
  if (!identical(colnames(Z.train), colnames(Z.test))) {
    stop("Training and test incidence design matrices have different columns.")
  }
  ## =====================================================
  ## Survival quantities
  ## =====================================================
  TimeL.train <- Y.train[, 1]
  TimeR.train <- Y.train[, 2]
  Status.train <- as.numeric(train_data$D)
  TimeL.test <- Y.test[, 1]
  TimeR.test <- Y.test[, 2]
  if (any(!Status.train %in% c(0, 1))) {
    stop("Training status D must contain only 0 and 1.")
  }
  ## =====================================================
  ## Parameter dimensions
  ## =====================================================
  nbeta <- ncol(X.train) - 1
  if (length(beta.init) != nbeta) {
    stop("Length of beta.init must equal the number of latency covariates: ", nbeta)
  }
  if (!is.finite(alpha.init) || alpha.init <= 0) {
    stop("alpha.init must be positive and finite.")
  }
  ## =====================================================
  ## Construct NN incidence data
  ## =====================================================
  NN.train <- as.data.frame(Z.train[, -1, drop = FALSE])
  NN.test <- as.data.frame(Z.test[, -1, drop = FALSE])
  NN_names <- paste0("Input_", seq_len(ncol(NN.train)))
  colnames(NN.train) <- NN_names
  colnames(NN.test) <- NN_names
  NN_formula <- as.formula(paste("y ~", paste(NN_names, collapse = " + ")))
  ## =====================================================
  ## Initial NN response
  ## =====================================================
  init_label <- Status.train
  if (length(unique(init_label)) < 2) {
    stop("Initial NN labels contain only one class. The NN incidence model caNNot be fitted.")
  }
  init_data <- data.frame(y = init_label, NN.train, check.names = FALSE)
  ## =====================================================
  ## Tune NN architecture once before EM
  ## =====================================================
  best_loss.init <- Inf
  best_params.init <- NULL
  best_model.init <- NULL
  for (hidden_val in hidden_grid) {
    fit_tmp <- tryCatch(
      neuralnet::neuralnet(
        formula = NN_formula,
        data = init_data,
        hidden = hidden_val,
        linear.output = FALSE,
        err.fct = "ce",
        lifesign = "none",
        stepmax = NN_stepmax,
        threshold = NN_threshold
      ),
      error = function(e) {
        NULL
      }
    )
    if (is.null(fit_tmp)) {
      next
    }
    pred_tmp <- tryCatch(
      as.numeric(neuralnet::compute(fit_tmp, covariate = NN.train)$net.result),
      error = function(e) {
        rep(NA_real_, nrow(train_data))
      }
    )
    if (length(pred_tmp) != nrow(train_data) || any(!is.finite(pred_tmp))) {
      next
    }
    pred_tmp <- pmin(pmax(pred_tmp, 1e-6), 1 - 1e-6)
    loss_tmp <- -mean(init_label * log(pred_tmp) + (1 - init_label) * log1p(-pred_tmp))
    if (is.finite(loss_tmp) && loss_tmp < best_loss.init) {
      best_loss.init <- loss_tmp
      best_params.init <- hidden_val
      best_model.init <- fit_tmp
    }
  }
  if (is.null(best_model.init) || is.null(best_params.init)) {
    stop("All initial NN architecture fits failed. Check the incidence predictors or increase NN_stepmax.")
  }
  ## =====================================================
  ## Initial uncured probabilities
  ## =====================================================
  best_params <- best_params.init
  mod.init <- best_model.init
  uncureprob.init <- tryCatch(
    as.numeric(neuralnet::compute(mod.init, covariate = NN.train)$net.result),
    error = function(e) {
      rep(NA_real_, nrow(train_data))
    }
  )
  if (length(uncureprob.init) != nrow(train_data) || any(!is.finite(uncureprob.init))) {
    stop("The initial NN produced invalid uncured probabilities.")
  }
  uncureprob.init <- pmin(pmax(uncureprob.init, 1e-6), 1 - 1e-6)
  ## =====================================================
  ## Run EM algorithm
  ## =====================================================
  emfit <- em.NN.ic(
    TimeL.train = TimeL.train,
    TimeR.train = TimeR.train,
    Status.train = Status.train,
    X.train = X.train,
    Z.train = Z.train,
    TimeL.test = TimeL.test,
    TimeR.test = TimeR.test,
    X.test = X.test,
    Z.test = Z.test,
    uncureprob = uncureprob.init,
    alpha = alpha.init,
    beta = beta.init,
    best_params = best_params,
    emmax = emmax,
    eps = eps,
    M.imp = M.imp
  )
  latency <- emfit$latencyfit
  alpha.est <- latency[1]
  beta.est <- latency[-1]
  conv <- emfit$tau
  ## =====================================================
  ## Bootstrap latency variance
  ## =====================================================
  if (Var) {
    latency_boot <- matrix(NA_real_, nrow = nboot, ncol = nbeta + 1)
    NN_models_boot <- vector("list", nboot)
    event_index <- which(Status.train == 1)
    cens_index <- which(Status.train == 0)
    if (length(event_index) == 0 || length(cens_index) == 0) {
      stop("Bootstrap requires both event and censored observations in the training data.")
    }
    boot_count <- 1
    attempt <- 0
    max_attempts <- 10 * nboot
    while (boot_count <= nboot && attempt < max_attempts) {
      attempt <- attempt + 1
      boot_event <- sample(event_index, size = length(event_index), replace = TRUE)
      boot_cens <- sample(cens_index, size = length(cens_index), replace = TRUE)
      boot_index <- c(boot_event, boot_cens)
      boot_TimeL <- TimeL.train[boot_index]
      boot_TimeR <- TimeR.train[boot_index]
      boot_Status <- Status.train[boot_index]
      bootX <- X.train[boot_index, , drop = FALSE]
      bootZ <- Z.train[boot_index, , drop = FALSE]
      ## =================================================
      ## Initial NN for bootstrap sample
      ## =================================================
      if (length(unique(boot_Status)) < 2) {
        next
      }
      bootNN <- as.data.frame(bootZ[, -1, drop = FALSE])
      colnames(bootNN) <- NN_names
      boot_init_data <- data.frame(y = boot_Status, bootNN, check.names = FALSE)
      boot_initial_model <- tryCatch(
        neuralnet::neuralnet(
          formula = NN_formula,
          data = boot_init_data,
          hidden = best_params,
          linear.output = FALSE,
          err.fct = "ce",
          lifesign = "none",
          stepmax = NN_stepmax,
          threshold = NN_threshold
        ),
        error = function(e) {
          NULL
        }
      )
      if (is.null(boot_initial_model)) {
        next
      }
      boot_uncure.init <- tryCatch(
        as.numeric(neuralnet::compute(boot_initial_model, covariate = bootNN)$net.result),
        error = function(e) {
          rep(NA_real_, length(boot_index))
        }
      )
      if (length(boot_uncure.init) != length(boot_index) || any(!is.finite(boot_uncure.init))) {
        next
      }
      boot_uncure.init <- pmin(pmax(boot_uncure.init, 1e-6), 1 - 1e-6)
      ## =================================================
      ## Bootstrap EM
      ## =================================================
      bootfit <- tryCatch(
        em.NN.ic(
          TimeL.train = boot_TimeL,
          TimeR.train = boot_TimeR,
          Status.train = boot_Status,
          X.train = bootX,
          Z.train = bootZ,
          TimeL.test = TimeL.test,
          TimeR.test = TimeR.test,
          X.test = X.test,
          Z.test = Z.test,
          uncureprob = boot_uncure.init,
          alpha = alpha.est,
          beta = beta.est,
          best_params = best_params,
          emmax = emmax,
          eps = eps,
          M.imp = M.imp
        ),
        error = function(e) {
          NULL
        }
      )
      if (is.null(bootfit)) {
        next
      }
      if (!is.finite(bootfit$tau) || any(!is.finite(bootfit$latencyfit))) {
        next
      }
      if (bootfit$tau <= eps && !is.null(bootfit$models) && any(!vapply(bootfit$models, is.null, logical(1)))) {
        latency_boot[boot_count, ] <- bootfit$latencyfit
        NN_models_boot[[boot_count]] <- bootfit$models
        boot_count <- boot_count + 1
      }
    }
    nboot_valid <- boot_count - 1
    if (nboot_valid < 2) {
      warning("Fewer than two valid NN bootstrap fits were obtained. Bootstrap standard deviations caNNot be computed.")
      latency_boot <- latency_boot[seq_len(max(nboot_valid, 1)), , drop = FALSE]
      NN_models_boot <- NN_models_boot[seq_len(max(nboot_valid, 1))]
      latency_var <- rep(NA_real_, nbeta + 1)
      latency_sd <- rep(NA_real_, nbeta + 1)
    } else {
      latency_boot <- latency_boot[seq_len(nboot_valid), , drop = FALSE]
      NN_models_boot <- NN_models_boot[seq_len(nboot_valid)]
      latency_var <- apply(latency_boot, 2, var)
      latency_sd <- sqrt(latency_var)
    }
  } else {
    latency_boot <- NULL
    NN_models_boot <- NULL
    latency_var <- NULL
    latency_sd <- NULL
    nboot_valid <- NULL
  }
  ## =====================================================
  ## Return fitted object
  ## =====================================================
  fit <- list()
  class(fit) <- "smcure.NN.ic"
  fit$call <- call
  fit$latency <- latency
  fit$alpha <- alpha.est
  fit$beta <- beta.est
  fit$tau <- conv
  fit$iter <- emfit$iter
  fit$converged <- emfit$converged
  fit$train <- emfit$train
  fit$test <- emfit$test
  fit$models <- emfit$models
  fit$best_params <- best_params
  fit$best_hidden <- paste(best_params, collapse = "-")
  fit$best_loss <- best_loss.init
  fit$initial_model <- mod.init
  fit$initial_uncureprob <- uncureprob.init
  ## =====================================================
  ## Store generic covariate names
  ## =====================================================
  fit$incidence_names <- colnames(Z.train)
  fit$latency_names <- c("alpha", colnames(X.train)[-1])
  ## =====================================================
  ## Bootstrap inference
  ## =====================================================
  if (Var) {
    fit$latency_boot <- latency_boot
    fit$NN_models_boot <- NN_models_boot
    fit$nboot_valid <- nboot_valid
    fit$latency_var <- latency_var
    fit$latency_sd <- latency_sd
    fit$latency_zvalue <- fit$latency / fit$latency_sd
    fit$latency_pvalue <- 2 * (1 - pnorm(abs(fit$latency_zvalue)))
  }
  return(fit)
}




## =====================================================
## 1. Load and prepare smoking data
## =====================================================
smoke <- read.table("~/Desktop/Desktop/SmokingDataBiocs.txt", header = TRUE)
#smoke <- read.table("/cloud/project/SmokingDataBiocs.txt", header = TRUE)
smoke$x1 <- smoke$SexF
smoke$Duration_original <- smoke$Duration
smoke$F10Cigs_original <- smoke$F10Cigs
smoke$x2 <- smoke$Duration
smoke$x3 <- smoke$F10Cigs
smoke$D <- smoke$Relapse
smoke$L <- smoke$Timept1
smoke$R <- ifelse(is.na(smoke$Timept2), Inf, smoke$Timept2)
#smoke$x2 = (smoke$x2- mean(smoke$x2))/(sd(smoke$x2))
#smoke$x3 = (smoke$x3- mean(smoke$x3))/(sd(smoke$x3))
duration_mean <- mean(smoke$x2)
duration_sd <- sd(smoke$x2)
f10cigs_mean <- mean(smoke$x3)
f10cigs_sd <- sd(smoke$x3)
smoke$x2 <- (smoke$x2 - duration_mean) / duration_sd
smoke$x3 <- (smoke$x3 - f10cigs_mean) / f10cigs_sd
smoke.data <- smoke[, c("L", "R", "D", "x1", "x2", "x3")]

## =====================================================
## INCIDENCE PREDICTION GRIDS
## Male = 0, Female = 1
## =====================================================
duration_seq <- seq(
  min(smoke$Duration_original),
  max(smoke$Duration_original),
  length.out = 40
)

f10cigs_seq <- seq(
  min(smoke$F10Cigs_original),
  max(smoke$F10Cigs_original),
  length.out = 40
)

pred_grid_male <- expand.grid(
  Duration = duration_seq,
  F10Cigs = f10cigs_seq
)

pred_grid_male$x1 <- 0
pred_grid_male$x2 <- (pred_grid_male$Duration - duration_mean) / duration_sd
pred_grid_male$x3 <- (pred_grid_male$F10Cigs - f10cigs_mean) / f10cigs_sd

pred_grid_female <- pred_grid_male
pred_grid_female$x1 <- 1
## =====================================================
## 2. Settings
## ====================================================
set.seed(2025)
#methods <- c("LOGIT")
methods <- c("LOGIT","SPLINE", "NN")#
K <- 10
Var <- FALSE
nboot <- 100
surv_formula <- Surv(L, R, type = "interval2") ~ x1 + x2 + x3
cure_formula <- ~ x1 + x2 + x3
alpha.init <- 1.5
beta.init <- rep(0, 3)
n <- nrow(smoke.data)
## =====================================================
## 3. Create stratified 10-fold CV
## =====================================================
fold_id <- rep(NA_integer_, n)
for (dval in unique(smoke.data$D)) {
  idx <- which(smoke.data$D == dval)
  idx <- sample(idx)
  fold_id[idx] <- rep(1:K, length.out = length(idx))
}
## =====================================================
## 4. 10-fold cross-validation
## =====================================================
cv_results <- list()
cv_fits <- list()
for (method_name in methods) {
  cat("\nRunning", method_name, "\n")
  pred_uncure <- rep(NA_real_, n)
  pred_S <- rep(NA_real_, n)
  pred_Sp <- rep(NA_real_, n)
  pred_S1 <- rep(NA_real_, n)
  method_fits <- vector("list", K)
  func <- get(paste0("smcure.", method_name, ".ic"))
  for (k in 1:K) {
    train_data <- smoke.data[fold_id != k, , drop = FALSE]
    test_data <- smoke.data[fold_id == k, , drop = FALSE]
    test_index <- which(fold_id == k)
    fit <- tryCatch({
      func(
        formula = surv_formula,
        cureform = cure_formula,
        alpha.init = alpha.init,
        beta.init = beta.init,
        train_data = train_data,
        test_data = test_data,
        Var = Var,
        emmax = 1000,
        eps = 1e-3,
        nboot = nboot
      )
    }, error = function(e) {
      cat("FAILED:", method_name, "fold", k, "\n")
      print(e$message)
      NULL
    })
    method_fits[[k]] <- fit
    if (!is.null(fit)) {
      pred_uncure[test_index] <- fit$test$Uncureprob
      pred_S[test_index] <- fit$test$S
      pred_Sp[test_index] <- fit$test$Sp
      pred_S1[test_index] <- fit$test$S1
    }
  }
  cv_results[[method_name]] <- list(
    Uncureprob = pred_uncure,
    S = pred_S,
    Sp = pred_Sp,
    S1 = pred_S1
  )
  cv_fits[[method_name]] <- method_fits
}

## =====================================================
## 5. ROC AND AUC USING 500 CURE-STATUS IMPUTATIONS
## =====================================================
set.seed(2025)
B <- 500
roc_results <- list()
auc_results <- data.frame(Method = methods, Mean_AUC = NA_real_, SD_AUC = NA_real_)

for (method_name in methods) {
  pred_pi <- cv_results[[method_name]]$Uncureprob
  pred_S <- cv_results[[method_name]]$S
  valid_idx <- complete.cases(pred_pi, pred_S)
  pred_pi_valid <- pred_pi[valid_idx]
  pred_S_valid <- pred_S[valid_idx]
  D_valid <- smoke.data$D[valid_idx]
  w_post <- D_valid + (1 - D_valid) * (1 - (1 - pred_pi_valid)^pred_S_valid)
  w_post <- pmin(pmax(w_post, 1e-6), 1 - 1e-6)
  auc_vec <- rep(NA_real_, B)
  roc_list <- vector("list", B)
  cens_idx <- which(D_valid == 0)
  
  for (b in 1:B) {
    V_imp <- D_valid
    V_imp[cens_idx] <- rbinom(length(cens_idx), size = 1, prob = w_post[cens_idx])
    if (length(unique(V_imp)) < 2) next
    roc_b <- pROC::roc(response = V_imp, predictor = pred_pi_valid, quiet = TRUE, direction = "<")
    auc_vec[b] <- as.numeric(pROC::auc(roc_b))
    roc_list[[b]] <- roc_b
  }
  
  roc_results[[method_name]] <- list(
    auc = auc_vec,
    roc = roc_list,
    mean_auc = mean(auc_vec, na.rm = TRUE),
    sd_auc = sd(auc_vec, na.rm = TRUE),
    w_post = w_post
  )
  
  auc_results[auc_results$Method == method_name, "Mean_AUC"] <- mean(auc_vec, na.rm = TRUE)
  auc_results[auc_results$Method == method_name, "SD_AUC"] <- sd(auc_vec, na.rm = TRUE)
}

print(auc_results)

## =====================================================
## 6. FINAL FULL-DATA FITS
## =====================================================
set.seed(2026)
final_nboot <- 100

final_fits <- list()

for (method_name in methods) {
  cat("\nFinal full-data fit:", method_name, "\n")
  func <- get(paste0("smcure.", method_name, ".ic"))
  
  final_fits[[method_name]] <- func(
    formula = surv_formula,
    cureform = cure_formula,
    alpha.init = alpha.init,
    beta.init = beta.init,
    train_data = smoke.data,
    test_data = smoke.data,
    Var = T,
    emmax = 1000,
    eps = 1e-3,
    nboot = final_nboot
  )
}

## =====================================================
## 7. LATENCY ESTIMATES, SEs AND P-VALUES
## =====================================================
latency_results <- do.call(rbind, lapply(names(final_fits), function(method_name) {
  fit <- final_fits[[method_name]]
  
  data.frame(
    Method = method_name,
    Parameter = fit$latency_names,
    Estimate = fit$latency,
    SE = fit$latency_sd,
    P_value = fit$latency_pvalue
  )
}))

latency_results$Parameter <- c("alpha", "Sex", "Duration", "F10Cigs")[
  match(latency_results$Parameter, c("alpha", "x1", "x2", "x3"))
]

rownames(latency_results) <- NULL
print(latency_results)

## =====================================================
## 8. LOGIT UNCURED-PROBABILITY SURFACE
## =====================================================
get_LOGIT_surface <- function(fit, pred_grid) {
  
  Z.grid <- model.matrix(cure_formula, data = pred_grid)
  
  estimate <- plogis(drop(Z.grid %*% fit$b))
  
  boot_pred <- t(apply(
    fit$incidence_boot,
    1,
    function(b) plogis(drop(Z.grid %*% b))
  ))
  
  lower <- apply(boot_pred, 2, quantile, probs = 0.025, na.rm = TRUE)
  upper <- apply(boot_pred, 2, quantile, probs = 0.975, na.rm = TRUE)
  
  data.frame(
    pred_grid,
    Estimate = estimate,
    Lower = lower,
    Upper = upper
  )
}

fit <- final_fits$LOGIT

surface_LOGIT_male <- get_LOGIT_surface(
  fit,
  pred_grid_male
)

surface_LOGIT_female <- get_LOGIT_surface(
  fit,
  pred_grid_female
)

## =====================================================
## 9. SPLINE PREDICTION FUNCTION
## =====================================================
make_SPLINE_grid <- function(newdata, fit, SPLINE_bases) {
  Z.parts <- list()
  
  for (v in fit$cure_covariates) {
    x.new <- newdata[[v]]
    
    if (!is.null(SPLINE_bases[[v]])) {
      basis.new <- predict(SPLINE_bases[[v]], x.new)
      colnames(basis.new) <- paste0(v, "_s", seq_len(ncol(basis.new)))
      Z.parts[[v]] <- basis.new
    } else {
      temp_formula <- as.formula(paste0("~ ", v))
      Z.parts[[v]] <- model.matrix(
        temp_formula,
        data = newdata
      )[, -1, drop = FALSE]
    }
  }
  
  as.matrix(cbind(
    "(Intercept)" = 1,
    do.call(cbind, Z.parts)
  ))
}


## =====================================================
## 10. SPLINE UNCURED-PROBABILITY SURFACE
## =====================================================

get_SPLINE_surface <- function(fit, pred_grid) {
  
  Z.grid <- make_SPLINE_grid(
    pred_grid,
    fit,
    fit$SPLINE_bases
  )
  
  estimate <- plogis(drop(Z.grid %*% fit$b))
  
  boot_pred <- matrix(
    NA_real_,
    nrow = fit$nboot_valid,
    ncol = nrow(pred_grid)
  )
  
  for (b in seq_len(fit$nboot_valid)) {
    
    Z.boot <- make_SPLINE_grid(
      pred_grid,
      fit,
      fit$SPLINE_boot[[b]]
    )
    
    boot_pred[b, ] <- plogis(
      drop(Z.boot %*% fit$incidence_boot[b, ])
    )
  }
  
  lower <- apply(boot_pred, 2, quantile, probs = 0.025, na.rm = TRUE)
  upper <- apply(boot_pred, 2, quantile, probs = 0.975, na.rm = TRUE)
  
  data.frame(
    pred_grid,
    Estimate = estimate,
    Lower = lower,
    Upper = upper
  )
}

fit <- final_fits$SPLINE

surface_SPLINE_male <- get_SPLINE_surface(
  fit,
  pred_grid_male
)

surface_SPLINE_female <- get_SPLINE_surface(
  fit,
  pred_grid_female
)
## =====================================================
## 11. NN PREDICTION FUNCTION
## =====================================================
predict_NN_grid <- function(models, newdata) {
  NN.grid <- as.data.frame(
    newdata[, c("x1", "x2", "x3"), drop = FALSE]
  )
  
  colnames(NN.grid) <- paste0(
    "Input_",
    seq_len(ncol(NN.grid))
  )
  
  pred.mat <- sapply(models, function(mod) {
    if (is.null(mod)) {
      return(rep(NA_real_, nrow(NN.grid)))
    }
    
    tryCatch(
      as.numeric(
        neuralnet::compute(
          mod,
          covariate = NN.grid
        )$net.result
      ),
      error = function(e) {
        rep(NA_real_, nrow(NN.grid))
      }
    )
  })
  
  if (is.null(dim(pred.mat))) {
    pred.mat <- matrix(pred.mat, ncol = 1)
  }
  
  rowMeans(pred.mat, na.rm = TRUE)
}

## =====================================================
## 12. NN UNCURED-PROBABILITY SURFACE
## =====================================================
get_NN_surface <- function(fit, pred_grid) {
  
  estimate <- predict_NN_grid(
    fit$models,
    pred_grid
  )
  
  boot_pred <- matrix(
    NA_real_,
    nrow = fit$nboot_valid,
    ncol = nrow(pred_grid)
  )
  
  for (b in seq_len(fit$nboot_valid)) {
    
    boot_pred[b, ] <- predict_NN_grid(
      fit$NN_models_boot[[b]],
      pred_grid
    )
  }
  
  lower <- apply(boot_pred, 2, quantile, probs = 0.025, na.rm = TRUE)
  upper <- apply(boot_pred, 2, quantile, probs = 0.975, na.rm = TRUE)
  
  data.frame(
    pred_grid,
    Estimate = estimate,
    Lower = lower,
    Upper = upper
  )
}

fit <- final_fits$NN

surface_NN_male <- get_NN_surface(
  fit,
  pred_grid_male
)

surface_NN_female <- get_NN_surface(
  fit,
  pred_grid_female
)

## =====================================================
## 13. CHECK INCIDENCE SURFACE RESULTS
## =====================================================
cat("\nLogit surface range:\n")
print(range(surface_logit$Estimate, na.rm = TRUE))

cat("\nSPLINE surface range:\n")
print(range(surface_SPLINE$Estimate, na.rm = TRUE))

cat("\nNN surface range:\n")
print(range(surface_NN$Estimate, na.rm = TRUE))

## =====================================================
## 14. CONVERT SURFACE DATA TO MATRICES FOR 3D PLOTS
## =====================================================
surface_to_matrix <- function(surface_data, variable) {
  matrix(
    surface_data[[variable]],
    nrow = length(duration_seq),
    ncol = length(f10cigs_seq)
  )
}

LOGIT_male_est <- surface_to_matrix(surface_LOGIT_male, "Estimate")
LOGIT_male_low <- surface_to_matrix(surface_LOGIT_male, "Lower")
LOGIT_male_up <- surface_to_matrix(surface_LOGIT_male, "Upper")

LOGIT_female_est <- surface_to_matrix(surface_LOGIT_female, "Estimate")
LOGIT_female_low <- surface_to_matrix(surface_LOGIT_female, "Lower")
LOGIT_female_up <- surface_to_matrix(surface_LOGIT_female, "Upper")

SPLINE_male_est <- surface_to_matrix(surface_SPLINE_male, "Estimate")
SPLINE_male_low <- surface_to_matrix(surface_SPLINE_male, "Lower")
SPLINE_male_up <- surface_to_matrix(surface_SPLINE_male, "Upper")

SPLINE_female_est <- surface_to_matrix(surface_SPLINE_female, "Estimate")
SPLINE_female_low <- surface_to_matrix(surface_SPLINE_female, "Lower")
SPLINE_female_up <- surface_to_matrix(surface_SPLINE_female, "Upper")

NN_male_est <- surface_to_matrix(surface_NN_male, "Estimate")
NN_male_low <- surface_to_matrix(surface_NN_male, "Lower")
NN_male_up <- surface_to_matrix(surface_NN_male, "Upper")

NN_female_est <- surface_to_matrix(surface_NN_female, "Estimate")
NN_female_low <- surface_to_matrix(surface_NN_female, "Lower")
NN_female_up <- surface_to_matrix(surface_NN_female, "Upper")

## =====================================================
## CHECK 95% CONFIDENCE BOUNDS
## =====================================================
check_bounds <- function(surface_data, method_name) {
  below <- sum(surface_data$Estimate < surface_data$Lower, na.rm = TRUE)
  above <- sum(surface_data$Estimate > surface_data$Upper, na.rm = TRUE)
  
  cat(
    method_name,
    ": estimate below lower =", below,
    "| estimate above upper =", above,
    "\n"
  )
}

check_bounds(surface_logit, "Logit")
check_bounds(surface_SPLINE, "SPLINE")
check_bounds(surface_NN, "NN")
## =====================================================
## 15. INTERACTIVE 3D PLOT FUNCTION
## Estimate + 95% pointwise bootstrap confidence bounds
## =====================================================
library(plotly)

plot_incidence_surface <- function(est, lower, upper, title) {
  
  plot_ly() %>%
    
    add_surface(
      x = duration_seq,
      y = f10cigs_seq,
      z = t(est),
      name = "Estimate",
      showscale = TRUE,
      opacity = 1
    ) %>%
    
    add_surface(
      x = duration_seq,
      y = f10cigs_seq,
      z = t(lower),
      name = "Lower 95% CI",
      showscale = FALSE,
      opacity = 0.25
    ) %>%
    
    add_surface(
      x = duration_seq,
      y = f10cigs_seq,
      z = t(upper),
      name = "Upper 95% CI",
      showscale = FALSE,
      opacity = 0.25
    ) %>%
    
    layout(
      title = title,
      scene = list(
        xaxis = list(title = "Duration"),
        yaxis = list(title = "Average Cigarettes Smoked"),
        zaxis = list(
          title = "Uncured Probability",
          range = c(0, 1)
        )
      )
    )
}


## =====================================================
## 16. DRAW 3D INCIDENCE PLOTS
## =====================================================
plot_LOGIT_male <- plot_incidence_surface(
  LOGIT_male_est,
  LOGIT_male_low,
  LOGIT_male_up,
  "Model: LOGIT   Sex: male"
)

plot_LOGIT_female <- plot_incidence_surface(
  LOGIT_female_est,
  LOGIT_female_low,
  LOGIT_female_up,
  "Model: LOGIT   Sex: female"
)

plot_SPLINE_male <- plot_incidence_surface(
  SPLINE_male_est,
  SPLINE_male_low,
  SPLINE_male_up,
  "Model: SPLINE   Sex: male"
)

plot_SPLINE_female <- plot_incidence_surface(
  SPLINE_female_est,
  SPLINE_female_low,
  SPLINE_female_up,
  "Model: SPLINE   Sex: female"
)

plot_NN_male <- plot_incidence_surface(
  NN_male_est,
  NN_male_low,
  NN_male_up,
  "Model: NN   Sex: male"
)

plot_NN_female <- plot_incidence_surface(
  NN_female_est,
  NN_female_low,
  NN_female_up,
  "Model: NN   Sex: female"
)

plot_LOGIT_male
plot_LOGIT_female

plot_SPLINE_male
plot_SPLINE_female

plot_NN_male
plot_NN_female
# =====================================================
## 17. AVERAGE ROC CURVES
## =====================================================
grid_fpr <- seq(0, 1, length.out = 101)

plot(
  NULL,
  xlim = c(0, 1),
  ylim = c(0, 1),
  xlab = "False Positive Rate",
  ylab = "True Positive Rate",
  main = ""
)

abline(0, 1, lty = 2)

cols <- c(
 LOGIT = "black",
  SPLINE = "blue",
  NN = "brown"
)

for (method_name in methods) {
  tpr_mat <- matrix(
    NA_real_,
    nrow = B,
    ncol = length(grid_fpr)
  )
  
  for (b in 1:B) {
    roc_b <- roc_results[[method_name]]$roc[[b]]
    
    if (is.null(roc_b)) next
    
    fpr_b <- 1 - roc_b$specificities
    tpr_b <- roc_b$sensitivities
    ord <- order(fpr_b)
    
    tpr_mat[b, ] <- approx(
      x = fpr_b[ord],
      y = tpr_b[ord],
      xout = grid_fpr,
      ties = mean,
      rule = 2
    )$y
  }
  
  mean_tpr <- colMeans(
    tpr_mat,
    na.rm = TRUE
  )
  
  lines(
    grid_fpr,
    mean_tpr,
    col = cols[method_name],
    lwd = 2
  )
}

legend(
  "bottomright",
  legend = paste0(
    methods,
    " AUC = ",
    sprintf("%.4f", auc_results$Mean_AUC)
  ),
  col = cols[methods],
  lwd = 2
)










time_taken <- proc.time() - ptm
time_taken