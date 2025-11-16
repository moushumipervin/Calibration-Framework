





library(MASS)
library(mgcv)
library(CVXR)
library(caret)
library(dplyr)

build_SU_folds <- function(df, K , id_col = "ID",  idx_S,idx_U0,seed = seed) {
  #stopifnot(all(c(id_col, delta_col) %in% names(df)))
  set.seed(seed)
  df$ID <- seq_len(nrow(df))
  
  id <- df[[id_col]]
  
  
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
###For the treatment group idx_U0=which(D==0) and idx_S=which(D==1)
##For the control group idx_U0=which(D==1) and idx_S=which(D==0)
k_fold_function1<-function(df,K,idx_S,idx_U0,seed){
  res <- build_SU_folds(df, K , id_col = "ID",idx_S,idx_U0, seed = seed)
  df$ID <- seq_len(nrow(df))
  lab_idx_all <- which(df$treat == 1)
  
  data_unlabeled<-list()
  for(k in 1:K){
    train_idx_lab <- setdiff(lab_idx_all, res$folds[[k]]$S_k_idx)
    validation_data_labeled <- df[res$folds[[k]]$S_k_idx,]
    validation_data_unlabeled<-df[res$folds[[k]]$U_k_idx,]
    
    train_data_labeled <-df[train_idx_lab,] ###get the data frame for trained data
    
    formula<-as.formula(paste0("re78~",paste0(colnames( subset(train_data_labeled, select = -c(data_id,treat, re78, ID,pi.hat))),collapse = "+")))
    
    lm_model<-lm(formula,data=as.data.frame( subset(train_data_labeled, select = -c(data_id,treat, ID,pi.hat))))
    
    
    y_hat<-predict( lm_model, newdata=as.data.frame(validation_data_unlabeled))
    validation_data_unlabeled$y.hat<-y_hat
    
    
    
    data_unlabeled[[k]]<-validation_data_unlabeled
  }
  return(data_unlabeled)
  
}



k_fold_function<-function(df,K,idx_S,idx_U0,seed){
  res <- build_SU_folds(df, K , id_col = "ID",idx_S,idx_U0, seed = seed)
  df$ID <- seq_len(nrow(df))
  lab_idx_all <- which(df$treat == 1)
  
  data_unlabeled<-list()
  for(k in 1:K){
    train_idx_lab <- setdiff(lab_idx_all, res$folds[[k]]$S_k_idx)
    validation_data_labeled <- df[res$folds[[k]]$S_k_idx,]
    validation_data_unlabeled<-df[res$folds[[k]]$U_k_idx,]
    
    train_data_labeled <-df[train_idx_lab,] ###get the data frame for trained data
    
    #formula<-as.formula(paste0("re78~",paste0("s(",colnames( subset(train_data_labeled, select = -c(data_id,treat, re78, ID,pi.hat))),")",collapse = "+")))
    #gam_model<-gam(re78 ~ s(age) + s(education) + black + hispanic + married + 
    #      nodegree + s(re74) + s(re75),data=as.data.frame( subset(train_data_labeled, select = -c(data_id,treat, ID,pi.hat))))
    
    
    # --- Prepare the training data safely ---
    train_sub <- subset(train_data_labeled, select = -c(data_id, treat, ID, pi.hat))
    train_sub <- as.data.frame(train_sub)
    
    # Drop variables with <= 2 unique values
    train_sub <- train_sub[sapply(train_sub, function(x) length(unique(x)) > 1)]
    
    
    # Convert binary predictors to factors (so mgcv doesn’t smooth them)
    for (v in c("black", "hispanic", "married", "nodegree")) {
      if (v %in% names(train_sub)) train_sub[[v]] <- factor(train_sub[[v]])
    }
    
    # --- Fit a stable GAM model (use cubic regression splines + REML) ---
    gam_model <- gam(
      re78 ~ s(age, bs = "cr", k = 5) +
        s(education, bs = "cr", k = 5) +
        black + hispanic + married + nodegree +
        s(re75, bs = "cr", k = 5),
      data   = train_sub,
      method = "REML"
    )
    
    # --- Predict on validation data safely ---
    val_sub <- as.data.frame(validation_data_unlabeled)
    
    # Ensure same factor levels in validation data
    for (v in intersect(c("black", "hispanic", "married", "nodegree"), names(val_sub))) {
      val_sub[[v]] <- factor(val_sub[[v]], levels = levels(train_sub[[v]]))
    }
    
    # Make predictions
    val_sub$y.hat <- predict(gam_model, newdata = val_sub)
    
    # Update the validation dataset
    validation_data_unlabeled$y.hat <- val_sub$y.hat
    
    
    data_unlabeled[[k]]<-validation_data_unlabeled
  }
  return(data_unlabeled)
  
}

aipw_var <- function(Y, D, ehat, m1x, m0x) {
  n <- length(Y)
  tau_i <- m1x - m0x + D*(Y - m1x)/ehat - (1 - D)*(Y - m0x)/(1 - ehat)
  tau_hat <- mean(tau_i)
  psi <- tau_i - tau_hat
  se <- sqrt(sum(psi^2) / (n*(n-1)))       # same as sd(psi)/sqrt(n)
  ci <- tau_hat + c(-1,1)*1.96*se
  ci_width <- diff(ci)    # upper - lower
  list(tau = tau_hat, se = se, ci = ci,ci_width = ci_width   )
}

variance_proposed<-function(tau_hat,U1,U0,W1,W0){
  ps1<-(W1*U1)-sum(W1*U1)
  ps0<-(W0*U0)-sum(W0*U0)
  
  se<-sqrt((sum(ps1^2)+sum(ps0^2)))
  ci<-tau_hat+c(-1,1)*se*1.96
  ci_width<-diff(ci)
  return(list(se=se,ci=ci,ci_width=ci_width))
  
}

aipw_fit <- function(dat, y="re78", d="treat",
                     covars=c("age","education","black","hispanic","married","nodegree","re74","re75")) {
  X <- intersect(covars, names(dat))
  
  # propensity
  ps <- glm(reformulate(X, d), data=dat, family=binomial())
  ehat <- pmin(pmax(predict(ps, type="response"), 1e-6), 1-1e-6)
  
  # outcome models by arm
  m1 <- lm(reformulate(X, y), data=dat[dat[[d]]==1, , drop=FALSE])
  m0 <- lm(reformulate(X, y), data=dat[dat[[d]]==0, , drop=FALSE])
  
  m1x <- as.numeric(predict(m1, newdata=dat))
  m0x <- as.numeric(predict(m0, newdata=dat))
  
  aipw_var(dat[[y]], dat[[d]], ehat, m1x, m0x)
}



estimate_theta_EM_kfold_CVXR_ET<- function(theta1,theta0,data_full, K, D,seed, max.iter=50, eps=1e-4) {
  #theta <- as.matrix(th)
  iter  <- 0
  
  repeat {
    iter <- iter + 1
    theta<-theta1-theta0
    # --- prepare cross-fitted data for treatment group T = 1, stacked in fold order ---
    fold_T1 <- k_fold_function1(df = data_full, K = K, idx_S = which(D==1), idx_U0 = which(D==0), seed = seed)
    data_T1 <- do.call(rbind, fold_T1)
    n_fold_T1 <- sapply(fold_T1, nrow)
    offsets_T1 <- c(0, cumsum(n_fold_T1))
    
    # --- features ---
    x_cols_T1 <- grep("^x\\d+$", names(data_T1), value = TRUE)
    p_T1 <- length(x_cols_T1)
    X_T1 <- as.matrix(data_T1[, x_cols_T1, drop = FALSE])
    X1_T1 <- cbind(1, X_T1)  # add intercept
    
    # --- fitted values and treatment indicators ---
    yhat_T1 <- data_T1$y.hat-theta1
    T1 <- data_T1$treat
    idx_T1 <- which(T1 == 1)  # indices for treatment group (T = 1)
    
    
    # --- moment matrix for calibration (select what to balance) ---
    A_T1 <- cbind(yhat_T1)
    
    # --- observed outcomes for treatment group ---
    y_T1 <- data_T1$re78[idx_T1]
    
    
    
    # --- CVXR problem on treatment group only ---
    r <- ncol(A_T1)
    w1 <- CVXR::Variable(length(idx_T1 ), pos = TRUE)
    
    g1<-log(data_T1$pi.hat)
    
    constr_T1 <- list(sum(w1) == 1,
                      sum(w1*g1[idx_T1]) == mean(g1))
    
    
    # per-fold constraints on treated rows (means match)
    for (k in seq_len(K)) {
      
      idx_k_all <- ( offsets_T1[k]+1): offsets_T1[k+1]
      idx_k     <- intersect(idx_k_all, idx_T1)
      if (length(idx_k) == 0) next
      
      A_k <- A_T1[idx_k, , drop=FALSE]
      map <- match(idx_k, idx_T1)
      w_k <- w1[map]
      A_k_all<-A_T1[ idx_k_all, , drop=FALSE]
      constr_T1 <- c(constr_T1, list(
        t(A_k) %*% w_k ==  matrix(colSums(A_k_all)/nrow(A_T1), ncol=1)
      ))
    }
    
    ones <- rep(1,length(idx_T1))
    #objective <- CVXR::Minimize( sum(CVXR::kl_div(w, a)) )
    objective_T1 <- CVXR::Minimize( sum(CVXR::kl_div(w1, ones)) )
    
    
    
    prob_T1 <- CVXR::Problem(objective_T1, constr_T1)
    sol_T1  <- CVXR::solve(prob_T1)
    #############################################################################################
    #############################################################################################
    ############################################################################################
    
    # --- prepare cross-fitted data for CONTROL group T = 0, stacked in fold order ---
    fold_T0 <- k_fold_function1(df = data_full, K = K, idx_S = which(D==0), idx_U0 = which(D==1), seed = seed)
    data_T0 <- do.call(rbind, fold_T0)
    n_fold_T0 <- sapply(fold_T0, nrow)
    offsets_T0 <- c(0, cumsum(n_fold_T0))
    
    # --- features ---
    x_cols_T0 <- grep("^x\\d+$", names(data_T0), value = TRUE)
    p_T0 <- length(x_cols_T0)
    X_T0 <- as.matrix(data_T0[, x_cols_T0, drop = FALSE])
    X1_T0 <- cbind(1, X_T0)  # add intercept
    
    # --- fitted values and treatment indicators ---
    yhat_T0 <- data_T0$y.hat-theta0
    T0 <- data_T0$treat
    idx_T0 <- which(T0 == 0)  # indices for control group (T = 0)
    
    
    # --- moment matrix for calibration (select what to balance) ---
    A_T0 <- cbind(yhat_T0)
    
    # --- observed outcomes for control group ---
    y_T0 <- data_T0$re78[idx_T0]
    
    
    
    # --- CVXR problem on control group only ---
    
    w0 <- CVXR::Variable(length(idx_T0 ), pos = TRUE)
    
    g0<-log(1-data_T0$pi.hat)
    
    constr_T0 <- list(sum(w0) == 1,
                      sum(w0*g0[idx_T0]) == mean(g0))
    
    
    # per-fold constraints on control rows (means match)
    for (k in seq_len(K)) {
      
      idx_k_all <- ( offsets_T0[k]+1): offsets_T0[k+1]
      idx_k     <- intersect(idx_k_all, idx_T0)
      if (length(idx_k) == 0) next
      
      A_k <- A_T0[idx_k, , drop=FALSE]
      map <- match(idx_k, idx_T0)
      w_k <- w0[map]
      A_k_all<-A_T0[ idx_k_all, , drop=FALSE]
      constr_T0 <- c(constr_T0, list(
        t(A_k) %*% w_k ==  matrix(colSums(A_k_all)/nrow(A_T0), ncol=1)
      ))
    }
    
    
    ones0 <- rep(1,length(idx_T0))
    #objective <- CVXR::Minimize( sum(CVXR::kl_div(w, a)) )
    objective_T0 <- CVXR::Minimize( sum(CVXR::kl_div(w0, ones0)) )
    
    
    
    prob_T0 <- CVXR::Problem(objective_T0, constr_T0)
    sol_T0  <- CVXR::solve(prob_T0)
    
    
    
    if(sol_T1$status=="optimal" && sol_T0$status=="optimal"){
      w1<-as.numeric(sol_T1$getValue(w1))
      w0<-as.numeric(sol_T0$getValue(w0))
      theta1.new<-sum(w1*y_T1);theta0.new<-sum(w0*y_T0)
      #th.new<-sum(w1*y_T1)-sum(w0*y_T0)
      th.new<-theta1.new-theta0.new
      print(th.new)
      if (max(abs(th.new - as.numeric(theta))) < eps) {
        return(list(theta=matrix(th.new, ncol=1),w1=w1,w0=w0,y1=data_T1$re78,y0=data_T0$re78,yhat_T1= yhat_T1, yhat_T0= yhat_T0,theta1.new=theta1.new,theta0.new=theta0.new,g1=g1,g0=g0,ps1=data_T1$pi.hat,ps0=data_T0$pi.hat,T1=T1,T0=T0))
      }
      if (iter >= max.iter) {
        message("Maximum iterations reached.")
        return(list(theta=matrix(th.new, ncol=1),w1=w1,w0=w0,y1=data_T1$re78,y0=data_T0$re78,yhat_T1= yhat_T1, yhat_T0= yhat_T0,theta1.new=theta1.new,theta0.new=theta0.new,g1=g1,g0=g0,ps1=data_T1$pi.hat,ps0=data_T0$pi.hat,T1=T1,T0=T0))
      }
      
      theta1<-theta1.new;theta0<-theta0.new
      
    }
  }
}



estimate_theta_EM_kfold_CVXR_HD <- function(theta1,theta0,data_full, K, D,seed, max.iter=50, eps=1e-6) {
  #theta <- as.matrix(th)
  iter  <- 0
  
  repeat {
    iter <- iter + 1
    theta<-theta1-theta0
    # --- prepare cross-fitted data for treatment group T = 1, stacked in fold order ---
    fold_T1 <- k_fold_function1(df = data_full, K = K, idx_S = which(D==1), idx_U0 = which(D==0), seed = seed)
    data_T1 <- do.call(rbind, fold_T1)
    n_fold_T1 <- sapply(fold_T1, nrow)
    offsets_T1 <- c(0, cumsum(n_fold_T1))
    
    # --- features ---
    x_cols_T1 <-subset(data_T1,select=-c(data_id,treat,pi.hat,ID,y.hat,re78))
    p_T1 <- length(x_cols_T1)
    X_T1 <- as.matrix(data_T1[, colnames(x_cols_T1), drop = FALSE])
    X1_T1 <- cbind(1, X_T1)  # add intercept
    
    # --- fitted values and treatment indicators ---
    yhat_T1 <- data_T1$y.hat-theta1
    T1 <- data_T1$treat
    idx_T1 <- which(T1 == 1)  # indices for treatment group (T = 1)
    
    
    # --- moment matrix for calibration (select what to balance) ---
    #A_T1 <- cbind(yhat_T1, -sqrt((data_T1$pi.hat)))
    A_T1 <- cbind(yhat_T1)
    
    # --- observed outcomes for treatment group ---
    y_T1 <- data_T1$re78[idx_T1]
    
    
    
    # --- CVXR problem on treatment group only ---
    r <- ncol(A_T1)
    w1 <- CVXR::Variable(length(idx_T1 ), pos = TRUE)
    
    g1<--sqrt((data_T1$pi.hat))
    
    constr_T1 <- list(sum(w1) == 1,
                      sum(w1*g1[idx_T1]) == mean(g1))
    
    
    # per-fold constraints on treated rows (means match)
    for (k in seq_len(K)) {
      
      idx_k_all <- ( offsets_T1[k]+1): offsets_T1[k+1]
      idx_k     <- intersect(idx_k_all, idx_T1)
      if (length(idx_k) == 0) next
      
      A_k <- A_T1[idx_k, , drop=FALSE]
      map <- match(idx_k, idx_T1)
      w_k <- w1[map]
      A_k_all<-A_T1[ idx_k_all, , drop=FALSE]
      constr_T1 <- c(constr_T1, list(
        t(A_k) %*% w_k ==  matrix(colSums(A_k_all)/nrow(A_T1), ncol=1)
      ))
    }
    
    
    objective_T1 <-  Minimize(-sum(2*sqrt(w1)))
    
    
    prob_T1 <- CVXR::Problem(objective_T1, constr_T1)
    sol_T1  <- CVXR::solve(prob_T1)
    lamda1<-sapply(constr_T1,sol_T1$getDualValue)[2]
    #############################################################################################
    #############################################################################################
    ############################################################################################
    
    # --- prepare cross-fitted data for CONTROL group T = 0, stacked in fold order ---
    fold_T0 <- k_fold_function1(df = data_full, K = K, idx_S = which(D==0), idx_U0 = which(D==1), seed = seed)
    data_T0 <- as.data.frame(do.call(rbind, fold_T0))
    n_fold_T0 <- sapply(fold_T0, nrow)
    offsets_T0 <- c(0, cumsum(n_fold_T0))
    
    # --- features ---
    x_cols_T0 <- subset(data_T0,select=-c(data_id,treat,pi.hat,ID,y.hat,re78))
    p_T0 <- length(x_cols_T0)
    X_T0 <- as.matrix(data_T0[,colnames(x_cols_T0), drop = FALSE])
    X1_T0 <- cbind(1, X_T0)  # add intercept
    
    # --- fitted values and treatment indicators ---
    yhat_T0 <- data_T0$y.hat-theta0
    T0 <- data_T0$treat
    idx_T0 <- which(T0 == 0)  # indices for control group (T = 0)
    
    
    # --- moment matrix for calibration (select what to balance) ---
    A_T0 <- cbind(yhat_T0)
    
    # --- observed outcomes for control group ---
    y_T0 <- data_T0$re78[idx_T0]
    
    
    
    # --- CVXR problem on control group only ---
    
    w0 <- CVXR::Variable(length(idx_T0 ), pos = TRUE)
    
    g0<--sqrt((1-data_T0$pi.hat))
    
    constr_T0 <- list(sum(w0) == 1,
                      sum(w0*g0[idx_T0]) == mean(g0))
    
    
    # per-fold constraints on control rows (means match)
    for (k in seq_len(K)) {
      
      idx_k_all <- ( offsets_T0[k]+1): offsets_T0[k+1]
      idx_k     <- intersect(idx_k_all, idx_T0)
      if (length(idx_k) == 0) next
      
      A_k <- A_T0[idx_k, , drop=FALSE]
      map <- match(idx_k, idx_T0)
      w_k <- w0[map]
      A_k_all<-A_T0[ idx_k_all, , drop=FALSE]
      constr_T0 <- c(constr_T0, list(
        t(A_k) %*% w_k ==  matrix(colSums(A_k_all)/nrow(A_T0), ncol=1)
      ))
    }
    
    
    objective_T0 <- Minimize(-sum(2*sqrt(w0)))
    
    
    prob_T0 <- CVXR::Problem(objective_T0, constr_T0)
    sol_T0  <- CVXR::solve(prob_T0)
    
    lamda0<-sapply(constr_T0,sol_T0$getDualValue)[2]
    
    if(sol_T1$status=="optimal" && sol_T0$status=="optimal"){
      w1<-as.numeric(sol_T1$getValue(w1))
      w0<-as.numeric(sol_T0$getValue(w0))
      theta1.new<-sum(w1*y_T1);theta0.new<-sum(w0*y_T0)
      #th.new<-sum(w1*y_T1)-sum(w0*y_T0)
      th.new<-theta1.new-theta0.new
      if (max(abs(th.new - as.numeric(theta))) < eps) {
        return(list(theta=matrix(th.new, ncol=1),w1=w1,w0=w0,y1=data_T1$re78,y0=data_T0$re78,yhat_T1= yhat_T1, yhat_T0= yhat_T0,theta1.new=theta1.new,theta0.new=theta0.new,g1=g1,g0=g0,ps1=data_T1$pi.hat,ps0=data_T0$pi.hat,lamda1=lamda1,lamda0=lamda0,T1=T1,T0=T0))
      }
      if (iter >= max.iter) {
        message("Maximum iterations reached.")
        # return(matrix(NA, ncol=1))
        return(list(theta=matrix(th.new, ncol=1),w1=w1,w0=w0,y1=data_T1$re78,y0=data_T0$re78,yhat_T1= yhat_T1, yhat_T0= yhat_T0,theta1.new=theta1.new,theta0.new=theta0.new,g1=g1,g0=g0,ps1=data_T1$pi.hat,ps0=data_T0$pi.hat,lamda1=lamda1,lamda0=lamda0,T1=T1,T0=T0))
      }
      
      theta1<-theta1.new;theta0<-theta0.new
      
    }
  }
}
