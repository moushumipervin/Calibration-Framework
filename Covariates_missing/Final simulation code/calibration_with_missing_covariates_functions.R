###############################################################################
# kfold_functions.R
#
# Purpose:
# This script contains helper functions for simulation and estimation in a
# missing-covariate setting. The workflow is:
#
# 1. Generate simulated data with outcome model (OR) and missingness model (PS)
# 2. Construct K-fold splits that preserve labeled and unlabeled observations
# 3. Use cross-fitting to predict the missing covariate z via z.hat
# 4. Estimate regression parameters using calibration-based weighting methods
#    under empirical tilting (ET) and Hellinger distance (HD)
#
# Main functions:
#   - generate_data()
#   - variance_theta_diag()
#   - variance_theta_diag1()
#   - build_SU_folds()
#   - k_fold_function()
#   - estimate_theta_EM_kfold_CVXR_ET()
#   - estimate_theta_EM_kfold_CVXR_HD()
#   - estimate_theta_nested_kfold_CVXR_ET()
#   - estimate_theta_nested_kfold_CVXR_HD()
#   - estimate_theta_EM_NOfold_CVXR()
###############################################################################

library(MASS)
library(mgcv)
library(CVXR)
library(caret)
library(dplyr)

###############################################################################
# Function: generate_data
#
# Description:
# This function generates simulated data for the regression setting with a
# partially observed covariate z. The argument OR controls the outcome model,
# and PS controls the missingness/selection model for D.
#
# Output:
# A data frame containing x, y, z, D, D.hat, and ID.
###############################################################################
generate_data <- function(n, OR, PS) {
  
  x <- rnorm(n, 0, 1)
  z <- rbinom(n, 1, 0.5)
  
  if (OR == 1) {
    beta1 <- -1
    beta2 <-  1
    beta3 <-  2
    y <- beta1 + beta2 * x + beta3 * z + rnorm(n, 0, 1)
  } else {
    b1 <-  0.5
    b2 <-  2
    b3 <- -1.5
    b4 <-  0.25
    b5 <-  1
    
    y <- b1 +
      b2 * sin(pi * x) +
      b3 * cos(2 * pi * x) +
      b4 * (x^3) +
      b5 * z +
      rnorm(n, 0, 2)
  }
  
  if (PS == 1) {
    eta1 <- -1
    eta2 <-  0.5
    eta3 <-  0.5
    px <- exp(eta1 + eta2 * y + eta3 * x) /
      (1 + exp(eta1 + eta2 * y + eta3 * x))
    D <- rbinom(n, 1, px)
  } else {
    eta1 <- -1
    eta2 <-  0
    eta3 <-  0
    px <- exp(eta1 + eta2 * y + eta3 * x) /
      (1 + exp(eta1 + eta2 * y + eta3 * x))
    D <- rbinom(n, 1, px)
  }
  
  D.hat <- fitted(glm(D ~ y + x, family = binomial))
  
  dat <- data.frame(
    x = x,
    y = y,
    z = z,
    D = D,
    D.hat = D.hat
  )
  
  dat$ID <- seq_len(nrow(dat))
  return(dat)
}

###############################################################################
# Function: variance_theta_diag
#
# Description:
# This function computes a sandwich-style variance estimate for theta_hat based
# on weighted estimating equations. It returns standard errors and confidence
# intervals for the regression coefficients.
#
# Note:
# This version rescales weights to sum approximately to the number of observed
# units used in the estimating equations.
###############################################################################
variance_theta_diag <- function(y, X, w_hat, theta_hat, ci_level = 0.95) {
  
  N <- length(y)
  p <- length(theta_hat)
  
  W <- w_hat
  n_eff <- length(W)
  W_rescaled <- W * n_eff
  
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("X", 1:p)
  }
  
  resid <- as.numeric(y - X %*% theta_hat)
  
  D_mat <- matrix(NA_real_, nrow = N, ncol = p)
  for (i in 1:N) {
    x_i <- as.numeric(X[i, ])
    D_mat[i, ] <- W_rescaled[i] * resid[i] * x_i
  }
  
  d_bar <- colMeans(D_mat)
  
  tau_hat <- matrix(0, p, p)
  for (i in 1:N) {
    x_i <- as.numeric(X[i, ])
    tau_hat <- tau_hat + W_rescaled[i] * (x_i %o% x_i)
  }
  tau_hat <- tau_hat / N
  
  middle <- matrix(0, p, p)
  for (i in 1:N) {
    diff_i <- D_mat[i, ] - d_bar
    middle <- middle + diff_i %o% diff_i
  }
  middle <- middle / (N * (N - 1))
  
  tau_inv <- solve(tau_hat)
  var_theta <- tau_inv %*% middle %*% t(tau_inv)
  
  variances <- diag(var_theta)
  names(variances) <- colnames(X)
  
  se <- sqrt(variances)
  alpha <- 1 - ci_level
  z_val <- qnorm(1 - alpha / 2)
  
  lower <- theta_hat - z_val * se
  upper <- theta_hat + z_val * se
  width <- upper - lower
  
  result <- data.frame(
    Variable = colnames(X),
    Estimate = theta_hat,
    Variance = variances,
    SE = se,
    CI_Lower = lower,
    CI_Upper = upper,
    CI_Width = width,
    row.names = NULL
  )
  
  return(result)
}

###############################################################################
# Function: variance_theta_diag1
#
# Description:
# This is a second variance estimator where the weights are explicitly rescaled
# to sum to the full sample size N. This is useful when inference should be
# reported relative to the full stacked data size rather than the observed
# labeled sample size alone.
###############################################################################
variance_theta_diag1 <- function(y, X, w_hat, theta_hat, N, ci_level = 0.95) {
  
  n1 <- length(y)
  p  <- length(theta_hat)
  
  w_star <- N * w_hat / sum(w_hat)
  resid <- as.numeric(y - X %*% theta_hat)
  
  D_mat <- matrix(NA_real_, nrow = n1, ncol = p)
  for (i in 1:n1) {
    x_i <- as.numeric(X[i, ])
    D_mat[i, ] <- w_star[i] * resid[i] * x_i
  }
  
  d_bar <- colSums(D_mat) / N
  
  tau_hat <- matrix(0, p, p)
  for (i in 1:n1) {
    x_i <- as.numeric(X[i, ])
    tau_hat <- tau_hat + w_star[i] * (x_i %o% x_i)
  }
  tau_hat <- tau_hat / N
  
  middle <- matrix(0, p, p)
  for (i in 1:n1) {
    diff_i <- D_mat[i, ] - d_bar
    middle <- middle + diff_i %o% diff_i
  }
  middle <- middle / (N * (N - 1))
  
  tau_inv <- solve(tau_hat)
  var_theta <- tau_inv %*% middle %*% t(tau_inv)
  
  variances <- diag(var_theta)
  names(variances) <- colnames(X)
  
  se <- sqrt(variances)
  alpha <- 1 - ci_level
  z_val <- qnorm(1 - alpha / 2)
  
  lower <- theta_hat - z_val * se
  upper <- theta_hat + z_val * se
  width <- upper - lower
  
  result <- data.frame(
    Variable = colnames(X),
    Estimate = theta_hat,
    Variance = variances,
    SE = se,
    CI_Lower = lower,
    CI_Upper = upper,
    CI_Width = width,
    row.names = NULL
  )
  
  return(result)
}

###############################################################################
# Function: build_SU_folds
#
# Description:
# This function splits the data into K folds while keeping labeled observations
# (D = 1) and unlabeled observations (D = 0) separate. For each fold:
#   - S_k contains labeled validation observations
#   - U_k contains all validation observations in that fold
#
# Output:
# A list containing fold-specific indices and fold membership vectors.
###############################################################################
build_SU_folds <- function(df, K, id_col = "ID", delta_col = "D", seed = 1234) {
  
  stopifnot(all(c(id_col, delta_col) %in% names(df)))
  
  set.seed(seed)
  df[[id_col]] <- seq_len(nrow(df))
  
  id <- df[[id_col]]
  delta <- df[[delta_col]]
  
  idx_S  <- which(delta == 1)
  idx_U0 <- which(delta == 0)
  
  nS  <- length(idx_S)
  nU0 <- length(idx_U0)
  
  if (nS < K || nU0 < K) {
    stop("Need at least K labeled and K unlabeled observations.")
  }
  
  split_K <- function(idx, K) {
    if (length(idx) == 0) {
      return(rep(list(integer(0)), K))
    }
    idx <- sample(idx)
    split(idx, rep(1:K, length.out = length(idx)))
  }
  
  S_parts  <- split_K(idx_S, K)
  U0_parts <- split_K(idx_U0, K)
  
  folds <- vector("list", K)
  fold_id_S <- integer(nrow(df))
  fold_id_U <- integer(nrow(df))
  
  for (k in 1:K) {
    S_k_idx  <- S_parts[[k]]
    U0_k_idx <- U0_parts[[k]]
    U_k_idx  <- c(S_k_idx, U0_k_idx)
    
    folds[[k]] <- list(
      S_k_ids = id[S_k_idx],
      U_k_ids = id[U_k_idx],
      S_k_idx = S_k_idx,
      U_k_idx = U_k_idx
    )
    
    fold_id_S[S_k_idx] <- k
    fold_id_U[U_k_idx] <- k
  }
  
  list(
    folds = folds,
    fold_id_S = fold_id_S,
    fold_id_U = fold_id_U
  )
}

###############################################################################
# Function: k_fold_function
#
# Description:
# This function performs cross-fitting for the missing covariate z. In each fold,
# it trains a logistic regression model for z using only labeled training data
# and then predicts z.hat for the validation fold.
#
# Output:
# A list of fold-specific validation data frames, each containing z.hat.
###############################################################################
k_fold_function <- function(df, K, seed) {
  
  res <- build_SU_folds(df, K, id_col = "ID", delta_col = "D", seed = seed)
  lab_idx_all <- which(df$D == 1)
  
  data_unlabeled <- list()
  
  for (k in 1:K) {
    train_idx_lab <- setdiff(lab_idx_all, res$folds[[k]]$S_k_idx)
    
    validation_data_unlabeled <- df[res$folds[[k]]$U_k_idx, ]
    train_data_labeled <- df[train_idx_lab, ]
    
    predictor_names <- colnames(
      subset(train_data_labeled, select = -c(D, z, ID, D.hat))
    )
    
    formula.glm <- as.formula(
      paste0("z ~ ", paste0(predictor_names, collapse = " + "))
    )
    
    glm_model <- glm(
      formula.glm,
      data = as.data.frame(subset(train_data_labeled, select = -c(D, ID, D.hat))),
      family = binomial(link = "logit")
    )
    
    z_hat <- predict(
      glm_model,
      newdata = as.data.frame(validation_data_unlabeled),
      type = "response"
    )
    
    validation_data_unlabeled$z.hat <- z_hat
    data_unlabeled[[k]] <- validation_data_unlabeled
  }
  
  return(data_unlabeled)
}

###############################################################################
# Function: estimate_theta_EM_kfold_CVXR_ET
#
# Description:
# This function estimates theta using an EM-style iterative procedure with
# cross-fitted z.hat and entropy-type calibration weights based on empirical
# tilting (ET). The calibration step is solved using CVXR.
#
# Output:
# A list containing theta, weights, and an inference table.
###############################################################################
estimate_theta_EM_kfold_CVXR_ET <- function(th, data_full, K, seed,
                                            max.iter = 50, eps = 1e-6) {
  
  theta <- as.matrix(th)
  fold_hat <- k_fold_function(df = data_full, K = K, seed = seed)
  data_all <- do.call(rbind, fold_hat)
  n_k <- sapply(fold_hat, nrow)
  ofs <- c(0, cumsum(n_k))
  
  new.data <- as.matrix(data_all[, c("x", "z.hat"), drop = FALSE])
  new.data1 <- cbind(1, new.data)
  
  D  <- data_all$D
  I1 <- which(D == 1)
  y  <- data_all$y
  g  <- log(data_all$D.hat)
  
  iter <- 0
  
  repeat {
    iter <- iter + 1
    
    error <- y - as.numeric(new.data1 %*% theta)
    H <- new.data1 * error
    A <- H
    
    w <- CVXR::Variable(length(I1), pos = TRUE)
    
    constr <- list(
      sum(w) == 1,
      sum(w * g[I1]) == mean(g)
    )
    
    for (k in seq_len(K)) {
      idx_k_all <- (ofs[k] + 1):ofs[k + 1]
      idx_k <- intersect(idx_k_all, I1)
      if (length(idx_k) == 0) next
      
      A_k <- A[idx_k, , drop = FALSE]
      map <- match(idx_k, I1)
      w_k <- w[map]
      A_k_all <- A[idx_k_all, , drop = FALSE]
      
      constr <- c(constr, list(
        t(A_k) %*% w_k == matrix(colSums(A_k_all) / nrow(A), ncol = 1)
      ))
    }
    
    ones <- rep(1, length(I1))
    objective <- CVXR::Minimize(sum(CVXR::kl_div(w, ones)))
    
    prob <- CVXR::Problem(objective, constr)
    sol <- CVXR::solve(prob)
    
    if (sol$status != "optimal") {
      return(list(theta = matrix(NA_real_, nrow = 1, ncol = ncol(new.data1))))
    }
    
    W <- as.numeric(sol$getValue(w))
    
    new.dat <- cbind(
      data.frame(y = data_all$y, x = data_all$x, z = data_all$z)[I1, ],
      w2 = W
    )
    
    th.new <- as.matrix(as.numeric(coefficients(lm(y ~ x + z, data = new.dat, weights = w2))))
    X1 <- cbind(1, as.matrix(data_all[, c("x", "z"), drop = FALSE]))
    
    if (max(abs(th.new - as.numeric(theta))) < eps) {
      return(list(
        theta = matrix(th.new, ncol = 1),
        w = W,
        estimate.table = variance_theta_diag1(
          y = y[I1],
          X = X1[I1, ],
          w_hat = W,
          theta_hat = th.new,
          N = nrow(data_all)
        )
      ))
    }
    
    if (iter >= max.iter) {
      message("Maximum iterations reached.")
      return(list(
        theta = matrix(th.new, ncol = 1),
        w = W,
        estimate.table = variance_theta_diag(
          y = y[I1],
          X = X1[I1, ],
          w_hat = W,
          theta_hat = th.new
        )
      ))
    }
    
    theta <- as.matrix(th.new)
  }
}

###############################################################################
# Function: estimate_theta_EM_kfold_CVXR_HD
#
# Description:
# This function is the HD counterpart of the ET estimator above. It uses the
# same EM-style structure and cross-fitted data, but solves the calibration
# problem under a Hellinger-distance-type objective.
###############################################################################
estimate_theta_EM_kfold_CVXR_HD <- function(th, data_full, K, seed,
                                            max.iter = 50, eps = 1e-6) {
  
  fold_hat <- k_fold_function(df = data_full, K = K, seed = seed)
  data_all <- do.call(rbind, fold_hat)
  n_k <- sapply(fold_hat, nrow)
  ofs <- c(0, cumsum(n_k))
  
  new.data <- as.matrix(data_all[, c("x", "z.hat"), drop = FALSE])
  new.data1 <- cbind(1, new.data)
  
  D  <- data_all$D
  I1 <- which(D == 1)
  y  <- data_all$y
  g  <- -sqrt(data_all$D.hat)/2
  
  theta <- as.matrix(th)
  iter <- 0
  
  repeat {
    iter <- iter + 1
    
    error <- y - as.numeric(new.data1 %*% theta)
    H <- new.data1 * error
    A <- H
    
    w <- CVXR::Variable(length(I1), pos = TRUE)
    
    constr <- list(
      sum(w) == 1,
      sum(w * g[I1]) == mean(g)
    )
    
    for (k in seq_len(K)) {
      idx_k_all <- (ofs[k] + 1):ofs[k + 1]
      idx_k <- intersect(idx_k_all, I1)
      if (length(idx_k) == 0) next
      
      A_k <- A[idx_k, , drop = FALSE]
      map <- match(idx_k, I1)
      w_k <- w[map]
      A_k_all <- A[idx_k_all, , drop = FALSE]
      
      constr <- c(constr, list(
        t(A_k) %*% w_k == matrix(colSums(A_k_all) / nrow(A), ncol = 1)
      ))
    }
    
    objective <- CVXR::Minimize(-sum(sqrt(w)))
    prob <- CVXR::Problem(objective, constr)
    sol <- CVXR::solve(prob, solver = "ECOS")
    
    if (sol$status != "optimal") {
      return(list(theta = matrix(NA_real_, nrow = 1, ncol = ncol(new.data1))))
    }
    
    W <- as.numeric(sol$getValue(w))
    
    new.dat <- cbind(
      data.frame(y = data_all$y, x = data_all$x, z = data_all$z)[I1, ],
      w2 = W
    )
    
    th.new <- as.matrix(as.numeric(coefficients(lm(y ~ x + z, data = new.dat, weights = w2))))
    X1 <- cbind(1, as.matrix(data_all[, c("x", "z"), drop = FALSE]))
    
    if (max(abs(th.new - as.numeric(theta))) < eps) {
      message("convergence achieved.")
      return(list(
        theta = matrix(th.new, ncol = 1),
        w = W,
        estimate.table = variance_theta_diag1(
          y = y[I1],
          X = X1[I1, ],
          w_hat = W,
          theta_hat = th.new,
          N = nrow(data_all)
        )
      ))
    }
    
    if (iter >= max.iter) {
      message("Maximum iterations reached.")
      return(list(
        theta = matrix(th.new, ncol = 1),
        w = W,
        estimate.table = variance_theta_diag(
          y = y[I1],
          X = X1[I1, ],
          w_hat = W,
          theta_hat = th.new
        )
      ))
    }
    
    theta <- as.matrix(th.new)
  }
}

###############################################################################
# Function: estimate_theta_nested_kfold_CVXR_ET
#
# Description:
# This function implements a nested optimization version of the ET estimator.
# For a given theta, it solves the inner calibration problem and then optimizes
# the resulting objective over theta using BFGS.
###############################################################################
estimate_theta_nested_kfold_CVXR_ET <- function(theta_init, data_full, K, seed) {
  
  fold_hat <- k_fold_function(df = data_full, K = K, seed = seed)
  data_all <- do.call(rbind, fold_hat)
  n_k <- sapply(fold_hat, nrow)
  ofs <- c(0, cumsum(n_k))
  
  new.data <- as.matrix(data_all[, c("x", "z.hat"), drop = FALSE])
  new.data1 <- cbind(1, new.data)
  
  D  <- data_all$D
  I1 <- which(D == 1)
  y  <- data_all$y
  g  <- log(data_all$D.hat)
  
  L_fn <- function(theta) {
    error <- y - as.numeric(new.data1 %*% theta)
    H <- new.data1 * error
    A <- H
    
    w <- CVXR::Variable(length(I1), pos = TRUE)
    
    new.data2 <- cbind(1, as.matrix(data_all[, c("x", "z"), drop = FALSE]))
    resid <- (y - as.matrix(new.data2) %*% theta)[I1]
    U_mat1 <- as.matrix(new.data2[I1, ]) * as.vector(resid)
    
    constr <- list(
      sum(w) == 1,
      sum(w * g[I1]) == mean(g),
      t(U_mat1) %*% w == rep(0, ncol(U_mat1))
    )
    
    mu <- c()
    for (k in seq_len(K)) {
      idx_k_all <- (ofs[k] + 1):ofs[k + 1]
      idx_k <- intersect(idx_k_all, I1)
      if (length(idx_k) == 0) next
      
      A_k <- A[idx_k, , drop = FALSE]
      map <- match(idx_k, I1)
      w_k <- w[map]
      A_k_all <- A[idx_k_all, , drop = FALSE]
      
      constr <- c(constr, list(
        t(A_k) %*% w_k == matrix(colSums(A_k_all) / nrow(A), ncol = 1)
      ))
      
      mu <- c(mu, colSums(A_k_all) / nrow(A))
    }
    
    U1 <- c(1, mean(g), mu)
    ones <- rep(1, length(I1))
    
    objective <- CVXR::Minimize(sum(CVXR::kl_div(w, ones)))
    prob <- CVXR::Problem(objective, constr)
    sol <- CVXR::solve(prob)
    
    if (sol$status != "optimal") {
      return(NA_real_)
    }
    
    lamda <- lapply(constr, sol$getDualValue)
    lamda <- lamda[-3]
    lambda_hat <- -unlist(lamda)
    
    W <- as.numeric(sol$getValue(w))
    val <- -(sum(W)) + as.vector(U1 %*% lambda_hat)
    return(val)
  }
  
  optim_result <- optim(
    par = theta_init,
    fn = L_fn,
    method = "BFGS",
    control = list(reltol = 1e-8, maxit = 300)
  )
  
  list(
    theta_hat = optim_result$par,
    value = optim_result$value,
    convergence = optim_result$convergence,
    message = optim_result$message
  )
}

###############################################################################
# Function: estimate_theta_nested_kfold_CVXR_HD
#
# Description:
# This function is the nested optimization version of the HD estimator. It
# follows the same idea as the ET version, but uses the HD objective in the
# inner calibration problem.
###############################################################################
estimate_theta_nested_kfold_CVXR_HD <- function(theta_init, data_full, K, seed) {
  
  L_fn <- function(theta) {
    fold_hat <- k_fold_function(df = data_full, K = K, seed = seed)
    data_all <- do.call(rbind, fold_hat)
    n_k <- sapply(fold_hat, nrow)
    ofs <- c(0, cumsum(n_k))
    
    new.data <- as.matrix(data_all[, c("x", "z.hat"), drop = FALSE])
    new.data1 <- cbind(1, new.data)
    
    D  <- data_all$D
    I1 <- which(D == 1)
    y  <- data_all$y
    g  <- -sqrt(data_all$D.hat)/2
    
    error <- y - as.numeric(new.data1 %*% theta)
    H <- new.data1 * error
    A <- H
    
    w <- CVXR::Variable(length(I1), pos = TRUE)
    
    new.data2 <- cbind(1, as.matrix(data_all[, c("x", "z"), drop = FALSE]))
    resid <- (y - as.matrix(new.data2) %*% theta)[I1]
    U_mat1 <- as.matrix(new.data2[I1, ]) * as.vector(resid)
    
    constr <- list(
      sum(w) == 1,
      sum(w * g[I1]) == mean(g),
      t(U_mat1) %*% w == rep(0, ncol(U_mat1))
    )
    
    mu <- c()
    for (k in seq_len(K)) {
      idx_k_all <- (ofs[k] + 1):ofs[k + 1]
      idx_k <- intersect(idx_k_all, I1)
      if (length(idx_k) == 0) next
      
      A_k <- A[idx_k, , drop = FALSE]
      map <- match(idx_k, I1)
      w_k <- w[map]
      A_k_all <- A[idx_k_all, , drop = FALSE]
      
      constr <- c(constr, list(
        t(A_k) %*% w_k == matrix(colSums(A_k_all) / nrow(A), ncol = 1)
      ))
      
      mu <- c(mu, colSums(A_k_all) / nrow(A))
    }
    
    U1 <- c(1, mean(g), mu)
    
    objective <- CVXR::Minimize(-sum(sqrt(w)))
    prob <- CVXR::Problem(objective, constr)
    sol <- CVXR::solve(prob)
    
    if (sol$status != "optimal") {
      return(NA_real_)
    }
    
    lamda <- lapply(constr, sol$getDualValue)
    lamda <- lamda[-3]
    lambda_hat <- -unlist(lamda)
    
    W <- as.numeric(sol$getValue(w))
    val <- -(sum(W)) + as.vector(U1 %*% lambda_hat)
    return(val)
  }
  
  optim_result <- optim(
    par = theta_init,
    fn = L_fn,
    method = "BFGS",
    control = list(reltol = 1e-8, maxit = 300)
  )
  
  list(
    theta_hat = optim_result$par,
    value = optim_result$value,
    convergence = optim_result$convergence,
    message = optim_result$message
  )
}

###############################################################################
# Function: estimate_theta_EM_NOfold_CVXR
#
# Description:
# This function estimates theta without K-fold cross-fitting. It uses labeled
# and unlabeled samples directly in a single calibration problem and updates
# theta iteratively until convergence.
###############################################################################
estimate_theta_EM_NOfold_CVXR <- function(th, y.hat, data_labelled,
                                          data_unlabelled,
                                          max.iter = 50, eps = 1e-6) {
  
  iter <- 0
  
  while (TRUE) {
    iter <- iter + 1
    theta <- as.matrix(th)
    
    new.data <- rbind(data_labelled[, -1], data_unlabelled)
    new.data1 <- cbind(1, new.data)
    
    error <- y.hat - as.numeric(as.matrix(new.data1) %*% theta)
    H1 <- apply(new.data1, 2, function(a) error * a)
    
    dat.new <- as.data.frame(cbind(
      D = c(rep(1, nrow(data_labelled)), rep(0, nrow(data_unlabelled))),
      new.data
    ))
    
    pi.hat <- fitted(glm(D ~ ., data = as.data.frame(dat.new), family = binomial()))
    H <- cbind(H1, log(pi.hat))
    
    p1 <- nrow(data_labelled)
    w <- Variable(p1, pos = TRUE)
    
    n <- nrow(data_labelled)
    A <- cbind(1, H[1:n, ])
    b <- c(1, colMeans(H))
    
    constraints <- list((t(A) %*% w) == b)
    objective <- Minimize(sum(kl_div(w, rep(1, p1))))
    
    problem <- Problem(objective, constraints = constraints)
    result3 <- solve(problem)
    
    if (result3$status == "optimal") {
      w2 <- as.numeric(result3$getValue(w))
      Q <- as.matrix(cbind(1, data_labelled[, -1]))
      th.new <- solve(t(Q) %*% diag(w2) %*% Q) %*%
        (t(Q) %*% diag(w2) %*% data_labelled[, 1])
      
      if (max(abs(theta - th.new)) < eps) {
        return(th.new)
      }
      
      if (iter >= max.iter) {
        message("Maximum iterations reached. Returning current estimate.")
        return(as.matrix(rep(NA, ncol(data_unlabelled) + 1)))
      }
      
      th <- th.new
    } else {
      return(as.matrix(rep(NA, ncol(data_unlabelled) + 1)))
    }
  }
}
