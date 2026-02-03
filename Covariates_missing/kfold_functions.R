
library(MASS)
library(mgcv)
library(CVXR)
library(caret)
library(dplyr)
library(mgcv)
library(caret)

generate_data<-function(n,PS,MAR){
  
  x<-rnorm(n,0,1)
  #z<-rexp(n,1)
  z<-rbinom(n,1,0.5)
  #z<-rnorm(n,0,1)
  if(MAR==1){
    
    #beta1<--2;beta2<-1;beta3<-2
    beta1<--1;beta2<-1;beta3<-2
    y<-beta1+beta2*x+beta3*z+rnorm(n,0,1)
  }else{
    b1 <- .5      # intercept
    b2 <- 2      # sin(pi x)
    b3 <- -1.5   # cos(2 pi x)
    b4 <- 0.25    # x^3
    b5 <- 1   # z * exp(0.5 x)
    b6 <- 1   # z * x^2
    
    y <- b1 +
      b2 * sin(pi * x) +
      b3 * cos(2 * pi * x) +
      b4 * (x^3) +b5*z+
      rnorm(n, 0, 2)
  }
  
  
 
  if(PS==1){
    eta1<--1;eta2<-.5;eta3<-0.5##with this choice 70% data is missing. ##PS1 with this choice 65% data is missing.final choice with exponential
    #eta1<--1;eta2<-0.5;eta3<-0.5 ###45% missing
    #px<-exp(eta1+eta2*y+eta3*x)/(1+exp(eta1+eta2*y+eta3*x))
    px<-exp(eta1+eta2*y+eta3*x)/(1+exp(eta1+eta2*y+eta3*x))
    D<-rbinom(n,1,px)
    mean(D)
  } else {
    #eta1<-0.5;eta2<-1;eta3<--.25 ###about 60% is missing
    #eta1<--.5;eta2<-.5;eta3<-1###20% IS MISSING FInal choice
    # eta1<--1.5;eta2<-.5;eta3<-.5
    #px<-exp(eta1+eta2*x+eta3*(y-1)^2)/(1+exp(eta1+eta2*x+eta3*(y-1)^2))
    #D<-rbinom(n,1,px)
    #mean(D)
    eta1<--1;eta2<-0;eta3<-0##with this choice 70% data is missing. ##PS1 with this choice 65% data is missing.final choice with exponential
    #eta1<--1;eta2<-0.5;eta3<-0.5 ###45% missing
    #px<-exp(eta1+eta2*y+eta3*x)/(1+exp(eta1+eta2*y+eta3*x))
    px<-exp(eta1+eta2*y+eta3*x)/(1+exp(eta1+eta2*y+eta3*x))
    D<-rbinom(n,1,px)
    mean(D)
  }
  
  
  
  D.hat<-fitted(glm(D~y+x,family=binomial))
  
  
  dat<-data.frame(x=x,y=y,z=z,D=D,D.hat=D.hat)
  dat$ID <- seq_len(nrow(dat))
  #dat$z[dat$D==0]<-NA
  return(dat)
}



variance_theta_diag <- function(y, X, w_hat, theta_hat, ci_level = 0.95) {
  N <- length(y)
  p <- length(theta_hat)
  # current:
  W <- w_hat              # sum(W) = 1
  
  # better: rescale to sum ≈ n_eff
  n_eff <- length(W)
  W_rescaled <- W * n
  # Ensure column names exist
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("X", 1:p)
  }
  
  # Residuals
  resid <- as.numeric(y - X %*% theta_hat)
  
  # Compute d_i
  D <- matrix(NA, nrow = N, ncol = p)
  for (i in 1:N) {
    x_i <- as.numeric(X[i, ])
    D[i, ] <- W_rescaled[i] * resid[i] * x_i
  }
  d_bar <- colMeans(D)
  
  # τ_hat
  tau_hat <- matrix(0, p, p)
  for (i in 1:N) {
    x_i <- as.numeric(X[i, ])
    tau_hat <- tau_hat + W_rescaled[i] * (x_i %o% x_i)
  }
  tau_hat <- tau_hat / N
  
  # Middle term
  middle <- matrix(0, p, p)
  for (i in 1:N) {
    diff_i <- D[i, ] - d_bar
    middle <- middle + diff_i %o% diff_i
  }
  middle <- middle / (N * (N - 1))
  
  # Sandwich variance
  tau_inv <- solve(tau_hat)
  var_theta <- tau_inv %*% middle %*% t(tau_inv)
  
  # Extract diagonals (variances)
  variances <- diag(var_theta)
  names(variances) <- colnames(X)
  
  # Standard errors and confidence intervals
  se <- sqrt(variances)
  alpha <- 1 - ci_level
  z_val <- qnorm(1 - alpha / 2)
  lower <- theta_hat - z_val * se
  upper <- theta_hat + z_val * se
  width=upper-lower
  # Output table
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


variance_theta_diag1 <- function(y, X, w_hat, theta_hat, N,ci_level = 0.95) {
  
  n1 <- length(y)
  p  <- length(theta_hat)
  
  # rescale weights so they sum to N
  w_star <- N * w_hat / sum(w_hat)   # if sum(w_hat)=1, this is just N * w_hat
  
  resid <- as.numeric(y - X %*% theta_hat)
  
  D <- matrix(NA_real_, nrow = n1, ncol = p)
  for (i in 1:n1) {
    x_i   <- as.numeric(X[i, ])
    D[i,] <- w_star[i] * resid[i] * x_i
  }
  d_bar <- colSums(D)/N
  
  tau_hat <- matrix(0, p, p)
  for (i in 1:n1) {
    x_i    <- as.numeric(X[i, ])
    tau_hat <- tau_hat + w_star[i] * (x_i %o% x_i)
  }
  tau_hat <- tau_hat / N            # now this is O(1)
  
  middle <- matrix(0, p, p)
  for (i in 1:n1) {
    diff_i <- D[i, ] - d_bar
    middle <- middle + diff_i %o% diff_i
  }
  middle <- middle / (N * (N - 1))  # empirical cov of D_i is O(1/N)
  
  tau_inv   <- solve(tau_hat)
  var_theta <- tau_inv %*% middle %*% t(tau_inv)
  
  # Extract diagonals (variances)
  variances <- diag(var_theta)
  names(variances) <- colnames(X)
  
  # Standard errors and confidence intervals
  se <- sqrt(variances)
  alpha <- 1 - ci_level
  z_val <- qnorm(1 - alpha / 2)
  lower <- theta_hat - z_val * se
  upper <- theta_hat + z_val * se
  width=upper-lower
  # Output table
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



build_SU_folds <- function(df, K , id_col = "ID", delta_col = "D", seed = seed) {
  stopifnot(all(c(id_col, delta_col) %in% names(df)))
  set.seed(seed)
  df$ID <- seq_len(nrow(df))
  
  id <- df[[id_col]]
  delta <- df[[delta_col]]
  
  idx_S  <- which(delta == 1)  # labeled indices
  idx_U0 <- which(delta == 0)  # unlabeled-only indices
  
  nS  <- length(idx_S)
  nU0 <- length(idx_U0)
  if (nS < K || nU0 < K) stop("Need at least K labeled and K unlabeled observations.")
  
  # helper: random, roughly-equal K-way split
  split_K <- function(idx, K) {
    if (length(idx) == 0) return(rep(list(integer(0)), K))
    idx <- sample(idx)  # shuffle
    split(idx, rep(1:K, length.out = length(idx)))
  }
  
  S_parts  <- split_K(idx_S,  K)  # disjoint labeled parts S_k
  U0_parts <- split_K(idx_U0, K)  # disjoint unlabeled-only parts U0_k
  
  folds <- vector("list", K)
  # fold assignment vectors (same length as df rows), 0 = not in that split
  fold_id_S  <- integer(nrow(df))
  fold_id_U  <- integer(nrow(df))
  
  for (k in 1:K) {
    S_k_idx  <- S_parts[[k]]
    U0_k_idx <- U0_parts[[k]]
    U_k_idx  <- c(S_k_idx, U0_k_idx)  # ensure S_k ⊆ U_k
    
    folds[[k]] <- list(
      S_k_ids  = id[S_k_idx],
      U_k_ids  = id[U_k_idx],
      S_k_idx  = S_k_idx,
      U_k_idx  = U_k_idx
    )
    
    fold_id_S[S_k_idx] <- k
    fold_id_U[U_k_idx] <- k
  }
  
  list(
    folds = folds,
    fold_id_S = fold_id_S,  # for labeled rows (0 if unlabeled)
    fold_id_U = fold_id_U   # for all rows in U_k (0 if outside that fold)
  )
}

# ---------------- Example ----------------
# df must have ID and delta (1=labeled, 0=unlabeled)

k_fold_function<-function(df,K,seed){
  res <- build_SU_folds(df, K , id_col = "ID", delta_col = "D", seed = seed)
  
  lab_idx_all <- which(df$D == 1)
  
  data_unlabeled<-list()
  for(k in 1:K){
    train_idx_lab <- setdiff(lab_idx_all, res$folds[[k]]$S_k_idx)
    validation_data_labeled <- df[res$folds[[k]]$S_k_idx,]
    validation_data_unlabeled<-df[res$folds[[k]]$U_k_idx,]
    
    train_data_labeled <-df[train_idx_lab,] ###get the data frame for trained data
    formula.glm<-as.formula(paste0("z~",paste0(colnames( subset(train_data_labeled, select = -c(D, z, ID,D.hat))),collapse = "+")))
    #formula<-as.formula(paste0("z~",paste0("s(",colnames( subset(train_data_labeled, select = -c(D, z, ID,D.hat))),")",collapse = "+")))
   # gam_model<-gam(formula,data=as.data.frame( subset(train_data_labeled, select = -c(D,ID,D.hat))))
    #z_hat<-predict( gam_model, newdata=as.data.frame(validation_data_unlabeled))
    glm_model<-glm(formula.glm,data=as.data.frame( subset(train_data_labeled, select = -c(D,ID,D.hat))),family=binomial(link="logit"))
    z_hat<-predict( glm_model, newdata=as.data.frame(validation_data_unlabeled), type = "response")
    
    
    #lm_model<-lm(formula.glm,data=as.data.frame( subset(train_data_labeled, select = -c(D,ID,D.hat))))
    
    #gam_model<-gam(formula,data=as.data.frame( subset(train_data_labeled, select = -c(D,ID,D.hat))),method="REML")
   # z_hat<-predict( gam_model, newdata=as.data.frame(validation_data_unlabeled))
    
    validation_data_unlabeled$z.hat<-z_hat
    
    data_unlabeled[[k]]<-validation_data_unlabeled
  }
  return(data_unlabeled)
  
}


#############without pi in the cv





estimate_theta_EM_kfold_CVXR_ET <- function(th, data_full, K, seed, max.iter=50, eps=1e-6) {
  theta <- as.matrix(th)
  fold_hat <- k_fold_function(df = data_full, K = K, seed = seed)
  data_all <- do.call(rbind, fold_hat)
  n_k      <- sapply(fold_hat, nrow)
  ofs      <- c(0, cumsum(n_k))
  
  # features
  # define p properly:
  new.data  <- as.matrix(data_all[, c("x","z.hat"), drop=FALSE])
  p      <- ncol(new.data)
  
  new.data1 <- cbind(1, new.data)
  
  z.hat <- data_all$z.hat
  D     <- data_all$D
  I1    <- which(D == 1)
  y<-data_all$y
  g<-log(data_all$D.hat)
  
  iter  <- 0
  
  repeat {
    iter <- iter + 1
    # --- prepare cross-fitted data, stacked in fold order ---
    
    # residual-based moments
    error <- y - as.numeric(new.data1 %*% theta)
    # H: N x (p+1)
    H <- new.data1 * error  # vectorized (faster & clearer than apply)
    
    # Moment matrix for calibration (choose what you want to balance)
    A <- H # N x r, r = p+1
    
   
    
    # --- CVXR problem on treated only ---
    r <- ncol(A)
    w <- CVXR::Variable(length(I1), pos = TRUE)
    a <- rep(1/length(I1), length(I1))  # base (uniform)
    
    
    constr <- list(sum(w) == 1,
                   sum(w*g[I1]) == mean(g))
    
    
    # per-fold constraints on treated rows (means match)
    for (k in seq_len(K)) {
      idx_k_all <- (ofs[k]+1):ofs[k+1]
      idx_k     <- intersect(idx_k_all, I1)
      if (length(idx_k) == 0) next
      
      A_k <- A[idx_k, , drop=FALSE]
      map <- match(idx_k, I1)
      w_k <- w[map]
      A_k_all<-A[ idx_k_all, , drop=FALSE]
      constr <- c(constr, list(
        t(A_k) %*% w_k ==  matrix(colSums(A_k_all)/nrow(A), ncol=1)
      ))
    }
    
    ones <- rep(1,length(I1))
    #objective <- CVXR::Minimize( sum(CVXR::kl_div(w, a)) )
    objective <- CVXR::Minimize( sum(CVXR::kl_div(w, ones)) )
    
    
    prob <- CVXR::Problem(objective, constr)
    sol  <- CVXR::solve(prob)
    
    
    if (sol$status != "optimal") {
      return(list(theta=matrix(NA_real_, nrow=1, ncol=ncol(new.data1))))
    }
    
    W  <- as.numeric(sol$getValue(w))
    
    new.dat<-cbind(data.frame(y=data_all$y,x=data_all$x,z=data_all$z)[I1,],w2=W)
    th.new<-as.matrix(as.numeric(coefficients(lm(y~x+z,data=new.dat,weights = w2))))
    
    X1<-cbind(1,as.matrix(data_all[, c("x","z"), drop=FALSE]))
    #print(th.new)
    if (max(abs(th.new - as.numeric(theta))) < eps) {
      #message("convergence achieved.")
      return(list(theta=matrix(th.new, ncol=1),w=W,estimate.table=variance_theta_diag1(
        y = y,
        X = X1[I1, ],
        w_hat = W,
        theta_hat = th.new,N=nrow(data_all)
      )))
    }
    if (iter >= max.iter) {
      message("Maximum iterations reached.")
      return(list(theta=matrix(th.new, ncol=1),w=W,estimate.table=variance_theta_diag(
        y = y,
        X = X1[I1, ],
        w_hat = W,
        theta_hat = th.new
      )))
    }
    theta <- as.matrix(th.new)
  }
}





estimate_theta_EM_kfold_CVXR_HD <- function(th, data_full, K, seed, max.iter=50, eps=1e-6) {
  # --- prepare cross-fitted data, stacked in fold order ---
  fold_hat <- k_fold_function(df = data_full, K = K, seed = seed)
  data_all <- do.call(rbind, fold_hat)
  n_k      <- sapply(fold_hat, nrow)
  ofs      <- c(0, cumsum(n_k))
  
  # features
  # define p properly:
  new.data  <- as.matrix(data_all[, c("x","z.hat"), drop=FALSE])
  p      <- ncol(new.data)
  
  new.data1 <- cbind(1, new.data)
  
  z.hat <- data_all$z.hat
  D     <- data_all$D
  I1    <- which(D == 1)
  y<-data_all$y
  g<--sqrt(data_all$D.hat)
  
  theta <- as.matrix(th)
  iter  <- 0
  
  repeat {
    iter <- iter + 1
    
    
    # residual-based moments
    error <- y - as.numeric(new.data1 %*% theta)
    # H: N x (p+1)
    H <- new.data1 * error  # vectorized (faster & clearer than apply)
    
    # Moment matrix for calibration (choose what you want to balance)
    A <- H # N x r, r = p+1
    
    
    
    # --- CVXR problem on treated only ---
    r <- ncol(A)
    w <- CVXR::Variable(length(I1), pos = TRUE)
    
   
    
    constr <- list(sum(w) == 1,
                   sum(w*g[I1]) == mean(g))
    
    
    # per-fold constraints on treated rows (means match)
    for (k in seq_len(K)) {
      idx_k_all <- (ofs[k]+1):ofs[k+1]
      idx_k     <- intersect(idx_k_all, I1)
      if (length(idx_k) == 0) next
      
      A_k <- A[idx_k, , drop=FALSE]
      map <- match(idx_k, I1)
      w_k <- w[map]
      A_k_all<-A[ idx_k_all, , drop=FALSE]
      constr <- c(constr, list(
        t(A_k) %*% w_k ==  matrix(colSums(A_k_all)/nrow(A), ncol=1)
      ))
    }
    
    objective <- Minimize(-sum(sqrt(w)))
    
    prob <- CVXR::Problem(objective, constr)
    sol  <- CVXR::solve(prob,solver = "ECOS")
    
    
    if (sol$status != "optimal") {
      return(list(theta=matrix(NA_real_, nrow=1, ncol=ncol(new.data1))))
      
    }
    
    W  <- as.numeric(sol$getValue(w))
    
    new.dat<-cbind(data.frame(y=data_all$y,x=data_all$x,z=data_all$z)[I1,],w2=W)
    th.new<-as.matrix(as.numeric(coefficients(lm(y~x+z,data=new.dat,weights = w2))))
    X1<-cbind(1,as.matrix(data_all[, c("x","z"), drop=FALSE]))
   # print(max(abs(th.new - as.numeric(theta))))
    if (max(abs(th.new - as.numeric(theta))) < eps) {
      message("convergence achieved.")
      return(list(theta=matrix(th.new, ncol=1),w=W,estimate.table=variance_theta_diag1(
        y = y[I1],
        X = X1[I1,],
        w_hat = W,
        theta_hat = th.new,N=nrow(data_all)
      )))
    }
    if (iter >= max.iter) {
      message("Maximum iterations reached.")
      return(list(theta=matrix(th.new, ncol=1),w=W,estimate.table=variance_theta_diag(
        y = y[I1],
        X = X1[I1,],
        w_hat = W,
        theta_hat = th.new
      )))
    }
    theta <- as.matrix(th.new)
  }
}





#####This function applies the two-loop function
estimate_theta_nested_kfold_CVXR_ET <- function(theta_init, data_full, K , seed ) {
  # ------------------------------
  # Define the objective function L(theta)
  # ------------------------------
  fold_hat <- k_fold_function(df = data_full, K = K, seed = seed)
  data_all <- do.call(rbind, fold_hat)
  n_k      <- sapply(fold_hat, nrow)
  ofs      <- c(0, cumsum(n_k))
  
  # features
  # define p properly:
  new.data  <- as.matrix(data_all[, c("x","z.hat"), drop=FALSE])
  p      <- ncol(new.data)
  
  new.data1 <- cbind(1, new.data)
  
  z.hat <- data_all$z.hat
  D     <- data_all$D
  I1    <- which(D == 1)
  y<-data_all$y
  g <- log(data_all$D.hat)
  L_fn <- function(theta) {
    # --- prepare cross-fitted data, stacked in fold order ---
    
    # residual-based moments
    error <- y - as.numeric(new.data1 %*% theta)
    # H: N x (p+1)
    H <- new.data1 * error  # vectorized (faster & clearer than apply)
    
    # Moment matrix for calibration (choose what you want to balance)
    A <- H # N x r, r = p+1
    
    
    
    # --- CVXR problem on treated only ---
    r <- ncol(A)
    
    # --- CVXR setup ---
    w <- CVXR::Variable(length(I1), pos = TRUE)
    
    
    new.data2 <- cbind(1, as.matrix(data_all[, c("x","z"), drop=FALSE]))
    resid <- (y - as.matrix(new.data2) %*% theta)[I1]
    U_mat1 <- as.matrix(new.data2[I1, ]) * as.vector(resid)
    
    # Constraints: global + per-fold
    constr <- list(
      sum(w) == 1,
      sum(w * g[I1]) == mean(g),
      t(U_mat1) %*% w == rep(0, ncol(U_mat1))
    )
    
    mu <- c()
    for (k in seq_len(K)) {
      idx_k_all <- (ofs[k] + 1):ofs[k + 1]
      idx_k     <- intersect(idx_k_all, I1)
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
    
    ones <- rep(1, length(I1))
    U1 <- c(1, mean(g), mu)
    
    objective <- CVXR::Minimize(sum(CVXR::kl_div(w, ones)))
    prob <- CVXR::Problem(objective, constr)
    sol  <- CVXR::solve(prob)
    
    if (sol$status != "optimal") {
      return(NA_real_)
    }
    
    # Extract duals (lambda)
    lamda <- lapply(constr, sol$getDualValue)
    lamda <- lamda[-3]                    # drop constraint #3 if needed
    lambda_hat <- -unlist(lamda)
    
    W <- as.numeric(sol$getValue(w))
    val <- -(sum(W)) + as.vector(U1 %*% lambda_hat)
    return(val)
  }
  
  # ------------------------------
  # Run optimization over theta
  # ------------------------------
  optim_result <- optim(
    par = theta_init,
    fn = L_fn,
    method = "BFGS",
    control = list(reltol = 1e-8, maxit = 300)
  )
  
  # Return estimated theta and value
  list(
    theta_hat = optim_result$par,
    value = optim_result$value,
    convergence = optim_result$convergence,
    message = optim_result$message
  )
}





estimate_theta_nested_kfold_CVXR_HD <- function(theta_init, data_full, K , seed ) {
  # ------------------------------
  # Define the objective function L(theta)
  # ------------------------------
  L_fn <- function(theta) {
    # --- prepare cross-fitted data, stacked in fold order ---
    fold_hat <- k_fold_function(df = data_full, K = K, seed = seed)
    data_all <- do.call(rbind, fold_hat)
    n_k      <- sapply(fold_hat, nrow)
    ofs      <- c(0, cumsum(n_k))
    
    # features
    # define p properly:
    new.data  <- as.matrix(data_all[, c("x","z.hat"), drop=FALSE])
    p      <- ncol(new.data)
    
    new.data1 <- cbind(1, new.data)
    
    z.hat <- data_all$z.hat
    D     <- data_all$D
    I1    <- which(D == 1)
    y<-data_all$y
    # residual-based moments
    error <- y - as.numeric(new.data1 %*% theta)
    # H: N x (p+1)
    H <- new.data1 * error  # vectorized (faster & clearer than apply)
    
    # Moment matrix for calibration (choose what you want to balance)
    A <- H # N x r, r = p+1
    
    
    
    # --- CVXR problem on treated only ---
    r <- ncol(A)
    
    # --- CVXR setup ---
    w <- CVXR::Variable(length(I1), pos = TRUE)
    g<--sqrt(data_all$D.hat)
    
    new.data2 <- cbind(1, as.matrix(data_all[, c("x","z"), drop=FALSE]))
    resid <- (y - as.matrix(new.data2) %*% theta)[I1]
    U_mat1 <- as.matrix(new.data2[I1, ]) * as.vector(resid)
    
    # Constraints: global + per-fold
    constr <- list(
      sum(w) == 1,
      sum(w * g[I1]) == mean(g),
      t(U_mat1) %*% w == rep(0, ncol(U_mat1))
    )
    
    mu <- c()
    for (k in seq_len(K)) {
      idx_k_all <- (ofs[k] + 1):ofs[k + 1]
      idx_k     <- intersect(idx_k_all, I1)
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
    
    objective <- Minimize(-sum(sqrt(w)))
    prob <- CVXR::Problem(objective, constr)
    sol  <- CVXR::solve(prob)
    
    if (sol$status != "optimal") {
      return(NA_real_)
    }
    
    # Extract duals (lambda)
    lamda <- lapply(constr, sol$getDualValue)
    lamda <- lamda[-3]                    # drop constraint #3 if needed
    lambda_hat <- -unlist(lamda)
    
    W <- as.numeric(sol$getValue(w))
    val <- -(sum(W)) + as.vector(U1 %*% lambda_hat)
    return(val)
  }
  
  # ------------------------------
  # Run optimization over theta
  # ------------------------------
  optim_result <- optim(
    par = theta_init,
    fn = L_fn,
    method = "BFGS",
    control = list(reltol = 1e-8, maxit = 300)
  )
  
  # Return estimated theta and value
  list(
    theta_hat = optim_result$par,
    value = optim_result$value,
    convergence = optim_result$convergence,
    message = optim_result$message
  )
}













estimate_theta_EM_NOfold_CVXR <- function(th, y.hat,data_labelled,data_unlabelled, max.iter = 50, eps = 1e-6)
{
  
  iter <- 0
  
  while (TRUE) {
    iter <- iter + 1
    
    theta <- as.matrix(th)
    
    ###Primal
    new.data<-(rbind(data_labelled[,-1],data_unlabelled))
    new.data1<-cbind(1,new.data)
    
    error<-y.hat-as.numeric(as.matrix(new.data1)%*%theta)
    #max(error^2)
    H1<-apply(new.data1,2,function(a)((error)*a))
    
    
    #pi.hat1<-fitted(glm(D~.,data=as.data.frame(data_full),family=binomial()))
    #pi.hat<-pi.hat1[order(-D)]
    
    dat.new<-as.data.frame(cbind(D=c(rep(1,nrow(data_labelled)),rep(0,nrow(data_unlabelled))),new.data))
    pi.hat<-fitted(glm(D~.,data=as.data.frame(dat.new),family=binomial()))
    
    H<-cbind(H1,log(pi.hat))
    
    p1<-sum(nrow(data_labelled))#+nrow(data_unlabelled)
    
    w <- Variable(p1,pos=TRUE)
    
    ones <- rep(1, p1)
    n<-nrow(data_labelled);N<-nrow(data_unlabelled)
    #A <- cbind(rep(1, n), R)
    A <- cbind(1,H[1:n,])
    
    
    
    
    #y<-data_labelled[,1]
    #resid <- y -  as.matrix(new.data1[1:length(y),] )%*% theta
    #U_mat1 <-  as.matrix(new.data1[1:length(y),] ) * as.vector(resid)  # N x p matrix
    #A <- cbind(1,H[1:n,],  U_mat1)
    
    
    #b <- c(1, colMeans(H),rep(0,ncol(H)-1))
    b <- c(1, colMeans(H))
    
    D=c(rep(1,nrow(data_labelled)),rep(0,nrow(data_unlabelled)))
    #constraints <- list((t(A) %*% (D*w)) == b)
    constraints <- list((t(A) %*% (w)) == b)
    
    
    # Define the objective function using kl_div
    
    
    #objective <- Minimize(sum(kl_div(w, ones)-ones)) 
    objective <- Minimize(sum(kl_div(w, ones))) 
    
    
    problem <- Problem(objective, constraints = constraints)
    result3 <- solve(problem)
    ############get dual value#######
    lamda<--result3$getDualValue(constraints[[1]])
    ####Dr. kim's dual formulation is not correct, in CVXR it is plus sign in lagrange multiplier
    w3=as.numeric(exp(A%*%lamda))
    
    if(result3$status=="optimal"){
      w2=as.numeric(result3$getValue(w))
      #w2<-w4[D==1]
      Q<-as.matrix(cbind(1,data_labelled[,-1]))
      th.new<-solve(t(Q)%*%diag(w2)%*%Q)%*%(t(Q)%*%diag(w2)%*%data_labelled[,1])
      
      #print(max(abs(theta - th.new)))
      if (max(abs(theta - th.new)) < eps) {
        #message("Convergence achieved.")
        return(th.new)
      }
      
      if (iter >= max.iter) {
        message("Maximum iterations reached. Returning current estimate.")
        return(as.matrix(rep(NA,ncol(data_unlabelled)+1)))
      }
      
      th <- th.new
    }else{
      return(as.matrix(rep(NA,ncol(data_unlabelled)+1)))
    }
    
  }
}

