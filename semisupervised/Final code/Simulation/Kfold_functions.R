
library(MASS)
library(mgcv)
library(CVXR)
library(caret)
library(dplyr)
build_SU_folds <- function(df, K , id_col = "ID", delta_col = "D", seed = seed) {
  #stopifnot(all(c(id_col, delta_col) %in% names(df)))
  set.seed(seed)
  #df$ID <- seq_len(nrow(df))
  
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
    
    formula<-as.formula(paste0("Y~",paste0("s(",colnames( subset(train_data_labeled, select = -c(D, Y, ID,pi.hat))),")",collapse = "+")))
    gam_model<-gam(formula,data=as.data.frame( subset(train_data_labeled, select = -c(D, ID,pi.hat))))
    y_hat<-predict( gam_model, newdata=as.data.frame(validation_data_unlabeled))
    validation_data_unlabeled$y.hat<-y_hat
    
    
    ##############################partition for pi.hat
    #train_data_labeled_pi<-df[-res$folds[[k]]$U_k_idx,]
    #train_data_pi<-as.data.frame(subset(train_data_labeled_pi, select = -c(Y,ID,pi.hat)))
    #validation_data_pi<-as.data.frame(subset(df[res$folds[[k]]$U_k_idx,], select = -c(Y,ID,pi.hat)))
    #formula.glm<-as.formula(paste0("D~",paste0(colnames( train_data_pi[,-1]),collapse = "+")))
    #glm_model<-glm(formula.glm,data=train_data_pi,family=binomial())
    #pi_hat_cv<-predict( glm_model, newdata=as.data.frame(validation_data_pi[,-1]))
    #validation_data_unlabeled$pi.hat.cv<-pi_hat_cv
    data_unlabeled[[k]]<-validation_data_unlabeled
  }
  return(data_unlabeled)
  
}


#############without pi in the cv





estimate_theta_EM_kfold_CVXR_ET <- function(th, data_full, K, seed, max.iter=50, eps=1e-6) {
  theta <- as.matrix(th)
  iter  <- 0
  
  repeat {
    iter <- iter + 1
    
    # --- prepare cross-fitted data, stacked in fold order ---
    fold_hat <- k_fold_function(df = data_full, K = K, seed = seed)
    data_all <- do.call(rbind, fold_hat)
    n_k      <- sapply(fold_hat, nrow)
    ofs      <- c(0, cumsum(n_k))
    
    # features
    # define p properly:
    x_cols <- grep("^X\\d+$", names(data_all), value = TRUE)
    p      <- length(x_cols)
    new.data  <- as.matrix(data_all[, x_cols, drop=FALSE])
    new.data1 <- cbind(1, new.data)
    
    y.hat <- data_all$y.hat
    D     <- data_all$D
    I1    <- which(D == 1)
    
    # residual-based moments
    error <- y.hat - as.numeric(new.data1 %*% theta)
    # H: N x (p+1)
    H <- new.data1 * error  # vectorized (faster & clearer than apply)
    
    # Moment matrix for calibration (choose what you want to balance)
    A <- H # N x r, r = p+1
    
    y<-data_all$Y[I1]
    resid <- y -  as.numeric(new.data1[I1,] %*% theta)
    U_mat1 <-  new.data1[I1,] * as.vector(resid)  # N x p matrix
    
    
    # --- CVXR problem on treated only ---
    r <- ncol(A)
    w <- CVXR::Variable(length(I1), pos = TRUE)
    a <- rep(1/length(I1), length(I1))  # base (uniform)
    g<-log(data_all$pi.hat)
    
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
      return(matrix(NA_real_, nrow=1, ncol=ncol(new.data1)))
    }
    
    W  <- as.numeric(sol$getValue(w))
    Q  <- new.data1[I1, , drop=FALSE]
    Y1 <- data_all$Y[I1]
    
    th.new <- solve(t(Q) %*% diag(W) %*% Q, t(Q) %*% diag(W) %*% Y1)
    #print(th.new)
    if (max(abs(th.new - as.numeric(theta))) < eps) {
      return(matrix(th.new, ncol=1))
    }
    if (iter >= max.iter) {
      message("Maximum iterations reached.")
      return(matrix(th.new, ncol=1))
    }
    theta <- as.matrix(th.new)
  }
}



estimate_theta_EM_kfold_CVXR_HD <- function(th, data_full, K, seed, max.iter=50, eps=1e-6) {
  theta <- as.matrix(th)
  iter  <- 0
  
  repeat {
    iter <- iter + 1
    
    # --- prepare cross-fitted data, stacked in fold order ---
    fold_hat <- k_fold_function(df = data_full, K = K, seed = seed)
    data_all <- do.call(rbind, fold_hat)
    n_k      <- sapply(fold_hat, nrow)
    ofs      <- c(0, cumsum(n_k))
    
    # features
    # define p properly:
    x_cols <- grep("^X\\d+$", names(data_all), value = TRUE)
    p      <- length(x_cols)
    new.data  <- as.matrix(data_all[, x_cols, drop=FALSE])
    new.data1 <- cbind(1, new.data)
    
    y.hat <- data_all$y.hat
    D     <- data_all$D
    I1    <- which(D == 1)
    
    # residual-based moments
    error <- y.hat - as.numeric(new.data1 %*% theta)
    # H: N x (p+1)
    H <- new.data1 * error  # vectorized (faster & clearer than apply)
    
    # Moment matrix for calibration (choose what you want to balance)
    A <- H # N x r, r = p+1
    
    y<-data_all$Y[I1]
    resid <- y -  as.numeric(new.data1[I1,] %*% theta)
    U_mat1 <-  new.data1[I1,] * as.vector(resid)  # N x p matrix
    
    
    # --- CVXR problem on treated only ---
    r <- ncol(A)
    w <- CVXR::Variable(length(I1), pos = TRUE)
    a <- rep(1/length(I1), length(I1))  # base (uniform)
    g<--sqrt(data_all$pi.hat)
    
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
    sol  <- CVXR::solve(prob)
    
    
    if (sol$status != "optimal") {
      return(matrix(NA_real_, nrow=1, ncol=ncol(new.data1)))
    }
    
    W  <- as.numeric(sol$getValue(w))
    Q  <- new.data1[I1, , drop=FALSE]
    Y1 <- data_all$Y[I1]
    
    th.new <- solve(t(Q) %*% diag(W) %*% Q, t(Q) %*% diag(W) %*% Y1)
    #print(th.new)
    if (max(abs(th.new - as.numeric(theta))) < eps) {
      return(matrix(th.new, ncol=1))
    }
    if (iter >= max.iter) {
      message("Maximum iterations reached.")
      return(matrix(th.new, ncol=1))
    }
    theta <- as.matrix(th.new)
  }
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




estimate_theta_nested_kfold_CVXR <- function(theta_init, data_full, K = 5, seed ) {
  # ------------------------------
  # Define the objective function L(theta)
  # ------------------------------
  L_fn <- function(theta) {
    # Create folds
    fold_hat <- k_fold_function(df = data_full, K = K, seed = seed)
    data_all <- do.call(rbind, fold_hat)
    n_k      <- sapply(fold_hat, nrow)
    ofs      <- c(0, cumsum(n_k))
    
    # Features
    x_cols <- grep("^X\\d+$", names(data_all), value = TRUE)
    p      <- length(x_cols)
    new.data  <- as.matrix(data_all[, x_cols, drop=FALSE])
    new.data1 <- cbind(1, new.data)
    
    y.hat <- data_all$y.hat
    D     <- data_all$D
    I1    <- which(D == 1)
    
    # residual-based moments
    error <- y.hat - as.numeric(new.data1 %*% theta)
    H <- new.data1 * error
    A <- H
    r <- ncol(A)
    
    # --- CVXR setup ---
    w <- CVXR::Variable(length(I1), pos = TRUE)
    g <- log(data_all$pi.hat)
    y <- data_all[I1, "Y"]
    resid <- y - as.matrix(new.data1[I1, ]) %*% theta
    U_mat1 <- as.matrix(new.data1[I1, ]) * as.vector(resid)
    
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







solve_lambda_dual <- function(theta,
                              data_full,          # data.frame with Y, D, y.hat, pi.hat, X1..Xp
                              K, seed,            # fold controls
                              lambda0 = NULL,
                              tol = 1e-8, maxit = 200,
                              armijo_c = 1e-4, shrink = 0.5,
                              ridge = 1e-4,            # robust default
                              include_g = TRUE,        # include global g constraint & term
                              include_score = TRUE) {  # include t(U) w = 0 & U-term in eta
  
  # ---- deps & helpers ----
  if (!requireNamespace("numDeriv", quietly = TRUE))
    stop("Please install.packages('numDeriv').")
  
  softmax_vec <- function(z) {
    m <- max(z); v <- exp(z - m); v / sum(v)
  }
  
  # ---- folds ----
  if (!is.function(k_fold_function))
    stop("k_fold_function(df, K, seed) must exist.")
  fold_hat <- k_fold_function(df = data_full, K = K, seed = seed)
  K        <- length(fold_hat)
  data_all <- do.call(rbind, fold_hat)
  n_k      <- sapply(fold_hat, nrow)
  ofs      <- c(0L, cumsum(n_k))
  
  # ---- features ----
  x_cols <- grep("^X\\d+$", names(data_all), value = TRUE)
  if (length(x_cols) == 0) stop("No X1..Xp columns in data_full.")
  X   <- as.matrix(data_all[, x_cols, drop = FALSE])
  X1  <- cbind(1, X)                      # (p+1)
  rA  <- ncol(X1)
  
  yhat_all <- data_all$y.hat; if (is.null(yhat_all)) stop("Missing y.hat.")
  D        <- data_all$D;      if (is.null(D))       stop("Missing D.")
  I1       <- which(D == 1)
  if (!length(I1)) stop("No treated rows (D==1).")
  
  # residual-based moments A (all rows)
  err_all <- yhat_all - as.numeric(X1 %*% theta)
  A_mat   <- X1 * err_all                # N x rA
  N       <- nrow(A_mat)
  
  # fold targets: match mean of A within each fold (all rows in fold)
  fold_targets <- lapply(seq_len(K), function(k) {
    idx_k_all <- (ofs[k] + 1L):ofs[k + 1L]
    colSums(A_mat[idx_k_all, , drop = FALSE])/N       # rA
  })
  
  # g term
  g_vec  <- log(data_all$pi.hat)
  if (any(!is.finite(g_vec))) stop("pi.hat must be > 0.")
  g_mean <- mean(g_vec)
  
  # treated score moments U on treated only (for constraint and eta if included)
  U_mat1 <- NULL
  if (include_score) {
    y_tr     <- data_all$Y[I1]; if (is.null(y_tr)) stop("Missing Y.")
    resid_tr <- y_tr - as.numeric(X1[I1, ] %*% theta)
    U_mat1   <- X1[I1, , drop = FALSE] * as.vector(resid_tr)   # |I1| x rA
  }
  
  # ---- lambda layout: per-fold A blocks + optional global g + optional global U ----
  r_fold <- rA                              # per-fold A multipliers
  L_tot  <- K * r_fold + as.integer(include_g) + if (include_score) rA else 0L
  
  # index helpers
  idx_lam_A_k <- function(k) {
    start <- (k - 1L) * r_fold + 1L
    start:(start + r_fold - 1L)
  }
  idx_lam_g <- if (include_g) (K * r_fold + 1L) else integer(0)
  idx_lam_U <- if (include_score) {
    (max(c(0, idx_lam_g)) + 1L):(max(c(0, idx_lam_g)) + rA)
  } else integer(0)
  
  lambda <- if (is.null(lambda0)) rep(0, L_tot) else {
    stopifnot(length(lambda0) == L_tot); as.numeric(lambda0)
  }
  
  # ---- g(lambda): stacked moment equations ----
  g_of <- function(lam) {
    # Build global eta for treated only
    eta_all <- numeric(length(I1))
    
    # per-fold A contribution
    for (k in seq_len(K)) {
      idx_k_all <- (ofs[k] + 1L):ofs[k + 1L]
      idx_k_tr  <- intersect(idx_k_all, I1)
      if (length(idx_k_tr) == 0) next
      
      lam_Ak <- lam[idx_lam_A_k(k)]
      A_k    <- A_mat[idx_k_tr, , drop = FALSE]          # |Ik| x rA
      map    <- match(idx_k_tr, I1)
      
      eta_all[map] <- eta_all[map] + as.vector(A_k %*% lam_Ak)
    }
    
    # global g term
    if (include_g) {
      lam_g <- lam[idx_lam_g]
      eta_all <- eta_all + g_vec[I1] * lam_g
    }
    
    # global U term
    if (include_score) {
      lam_U <- lam[idx_lam_U]
      eta_all <- eta_all + as.vector(U_mat1 %*% lam_U)
    }
    
    # global softmax → weights on treated (sum = 1 automatically)
    w_all <- softmax_vec(eta_all)
    
    # per-fold A moment equations
    U_stack <- c()
    for (k in seq_len(K)) {
      idx_k_all <- (ofs[k] + 1L):ofs[k + 1L]
      idx_k_tr  <- intersect(idx_k_all, I1)
      if (length(idx_k_tr) == 0) next
      
      targ_k <- fold_targets[[k]]                         # rA
      map    <- match(idx_k_tr, I1)
      A_k    <- A_mat[idx_k_tr, , drop = FALSE]
      w_k    <- w_all[map]
      
      U_k <- as.vector(t(A_k) %*% w_k) - targ_k           # rA
      U_stack <- c(U_stack, U_k)
    }
    
    # global constraints (no sum(w)=1; softmax enforces it)
    cons <- c()
    if (include_g)     cons <- c(cons, sum(w_all * g_vec[I1]) - g_mean)  # scalar
    if (include_score) cons <- c(cons, as.vector(t(U_mat1) %*% w_all))   # rA
    
    c(U_stack, cons)
  }
  
  merit <- function(lam) { gv <- g_of(lam); 0.5 * sum(gv * gv) }
  
  # ---- Gauss–Newton + Armijo, adaptive ridge ----
  it <- 0L
  repeat {
    gval <- g_of(lambda)
    if (sqrt(sum(gval * gval)) <= tol || it >= maxit) break
    
    J <- numDeriv::jacobian(g_of, lambda)       # m x L_tot
    JTJ <- crossprod(J)
    Jtg <- crossprod(J, gval)
    
    # robust direction
    ridge_try <- ridge
    d <- NA
    for (attempt in 1:6) {
      JTJ_r <- JTJ + diag(ridge_try, ncol(JTJ))
      d <- tryCatch(as.vector(-solve(JTJ_r, Jtg)),
                    error = function(e) NA)
      if (all(is.finite(d))) break
      ridge_try <- ridge_try * 10
    }
    if (!all(is.finite(d))) {
      # steepest descent fallback
      d <- -as.vector(Jtg)
    }
    
    # Armijo backtracking
    m0 <- merit(lambda)
    grad_phi <- as.vector(Jtg)                  # ∇(½||g||²) = Jᵗ g
    alpha <- 1
    while (merit(lambda + alpha * d) > m0 + armijo_c * alpha * sum(grad_phi * d)) {
      alpha <- alpha * shrink
      if (alpha < 1e-12) break
    }
    
    lambda <- lambda + alpha * d
    it <- it + 1L
  }
  
  # ---- final weights (treated) ----
  # rebuild eta_all at solution
  rebuild_eta_all <- function(lam) {
    eta_all <- numeric(length(I1))
    for (k in seq_len(K)) {
      idx_k_all <- (ofs[k] + 1L):ofs[k + 1L]
      idx_k_tr  <- intersect(idx_k_all, I1)
      if (length(idx_k_tr) == 0) next
      lam_Ak <- lam[idx_lam_A_k(k)]
      A_k    <- A_mat[idx_k_tr, , drop = FALSE]
      map    <- match(idx_k_tr, I1)
      eta_all[map] <- eta_all[map] + as.vector(A_k %*% lam_Ak)
    }
    if (include_g) {
      lam_g <- lam[idx_lam_g]
      eta_all <- eta_all + g_vec[I1] * lam_g
    }
    if (include_score) {
      lam_U <- lam[idx_lam_U]
      eta_all <- eta_all + as.vector(U_mat1 %*% lam_U)
    }
    eta_all
  }
  
  eta_all <- rebuild_eta_all(lambda)
  w_all   <- softmax_vec(eta_all)
  
  # per-fold treated weight slices (for inspection)
  w_by_fold <- vector("list", K)
  for (k in seq_len(K)) {
    idx_k_all <- (ofs[k] + 1L):ofs[k + 1L]
    idx_k_tr  <- intersect(idx_k_all, I1)
    if (length(idx_k_tr) == 0) { w_by_fold[[k]] <- numeric(0); next }
    map <- match(idx_k_tr, I1)
    w_by_fold[[k]] <- w_all[map]
  }
  lamda0<-log(1/sum(exp(eta_all)))
  mu<-c(1,fold_targets,g_mean)
  
  list(
    lambda_hat  = as.numeric(c(lamda0,lambda)),
    w_treated   = w_all,          # |I1| vector, sums to 1
    mu_target  = mu      # list by fold (treated only)
    # idx_treated = I1,
    # converged   = (sqrt(sum(g_of(lambda)^2)) <= tol),
    # iters       = it,
  )
}

estimate_theta_nested_kfold_NO_CVXR <- function(theta_init, data_full, K = 5, seed ) {
  
  L_fn <- function(theta) {
    
    sol1=solve_lambda_dual(theta,
                           data_full,          # data.frame with Y, D, y.hat, pi.hat, X1..Xp
                           K, seed,            # fold controls
                           lambda0 = NULL,
                           tol = 1e-8, maxit = 200,
                           armijo_c = 1e-4, shrink = 0.5,
                           ridge = 1e-4,            # robust default
                           include_g = TRUE,        # include global g constraint & term
                           include_score = TRUE) 
    
    r <- p + 1
    lambda_hat<-sol1$lambda_hat[1:(K*r+2)]
    
    U1<-unlist(sol1$mu_target)
    W<-sol1$w_treated
    val <- -(sum(W)) + as.vector(U1 %*% lambda_hat)
    return(val)
  }
  
  
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


#estimate_theta_nested_kfold_NO_CVXR(theta_init, data_full, K = 5, seed )
#estimate_theta_nested_kfold_CVXR(theta_init, data_full, K = 5, seed )
