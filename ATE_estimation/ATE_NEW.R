
library(dplyr)

generate_data<-function(n,OR,PS,seed){
  set.seed(seed)
  
  
Z <- mvrnorm(n,rep(0,p),diag(1,p))

if(PS==1){
  
  # px<-1 / (1 + exp(-0.10+Z[,1] - .5* Z[,2] )) ###Final 500
  
  
  # px<-1 / (1 + exp(-1+Z[,1] -0.5*  Z[,2]+.25*Z[,3]+.1*Z[,4] ))###testing
  
  
  px<-1 / (1 + exp(-0.1+Z[,1] - .5* Z[,2]+.5*Z[,3]+.25*Z[,4] )) ###This is the final one
  #px<-1/(1+exp(-(Z%*%c(1,1,1,1))))
  D<-rbinom(n,1,px)
  mean(D)
} else {
  #X <- cbind((Z[,2]/(1+exp(Z[,1])))+10,(Z[,1]*Z[,2]/25+0.6)^3)
  #X2<-scale(X1,center=TRUE,scale=TRUE)
  #px<-1 / (1 + exp(-0.25*(Z[,1]+Z[,2])+Z[,1]^2*Z[,2])) ###Final for 2 covariates
  
  #px<-1 / (1 + exp(-1+0.5*(Z[,1]-1)*(Z[,2]-2)+0.5*Z[,3]^2-0.5*Z[,4]^2+0.5*cos(Z[,4]))) ###Final but producing less bias
  px<-1 / (1 + exp(-1+0.5*(Z[,1]-1)*(Z[,2]-2)+0.5*Z[,3]^2*Z[,4]+0.25*Z[,1]*Z[,2]-0.5*Z[,4]^2+0.5*cos(Z[,4]))) 
  D<-rbinom(n,1,px)
  mean(D)
}

alpha0=1
alpha1=rep(1,p)
alpha2=rep(1,p)

if(OR==1){
  sigma= 1  # the sd of the error term eta       
  alpha0=1
  alpha1=rep(.5,p)
  alpha2=rep(1,p)
  eta=rnorm(n,0,sigma)
  y0<-1+0.5*Z[,1]+Z[,2]+Z[,3]+Z[,4]+rnorm(n,0,1)
  
}else{
  #X <- cbind((Z[,2]/(1+exp(Z[,1])))+10,(Z[,1]*Z[,2]/25+0.6)^3)###Final
  
  #y0<-210 + (27.4*X[,1] + 13.7*X[,2]) + rnorm(n,0,1)###Final
  sigma= 1  # the sd of the error term eta       
  alpha0=1
  alpha1=rep(1,p)
  alpha2=rep(1,p)
  eta=rnorm(n,0,sigma)
  y0<-alpha0+Z%*%alpha1+(Z^3-Z^2+exp(Z)+cos(Z))%*%alpha2+eta ###OM2
  
  #X <- cbind(exp(Z[,1])/2,(Z[,2]/(1+exp(Z[,1])))+10,(Z[,1]*Z[,3]/25+0.6)^3,(Z[,2]+Z[,4]+20)^2)
  #y0<-1 + (0.5*X[,1] +0.5*X[,2]+0.5*X[,3]+0.5*X[,4]) + rnorm(n,0,1)
}



y1<-y0
y<-D*y1+(1-D)*y0

pi.hat<-fitted(glm(D~Z,family=binomial(link="logit")))
dat<-data.frame(y=y,Z,D=D,pi.hat=pi.hat)
#dat<-data.frame(Z,y=y,D=D,y.hat0=y.hat0,y.hat1=y.hat1,pi.hat=pi.hat)
colnames(dat)[2:(p+1)]<-c(paste0("x",1:p))

return(dat)
}






generate_data_kang_schafer<-function(n,OR,PS,seed){
  set.seed(seed)
  
  
  Z <- matrix(rnorm(4*n),ncol=4,nrow=n) 
  
  X <- cbind(exp(Z[,1])/2,(Z[,2]/(1+exp(Z[,1])))+10,(Z[,1]*Z[,3]/25+0.6)^3,(Z[,2]+Z[,4]+20)^2)
  
  
  #X<-scale(X,center=TRUE,scale=TRUE)
  
  if(PS==1){
    
    px<-1 / (1 + exp(Z[,1] - 0.5 * Z[,2] + 0.25*Z[,3] + 0.1 * Z[,4]))
    
    D<-rbinom(n,1,px)
  } else {
    px<-1 / (1 + exp(X[,1] - 0.5 * X[,2]*X[,1] - 1*X[,3]^2 + 0.1 * X[,4])) ###Highly non response 0.2% responses
    
    D<-rbinom(n,1,px)
    }
  
  if(OR==1){
    b<-(27.4*Z[,1] + 13.7*Z[,2] +13.7*Z[,3] + 13.7*Z[,4])
    y1<-210 + b + rnorm(n)
    y0<-200-0.5*b+ rnorm(n)
    y<-D*y1+(1-D)*y0
    dat1<-data.frame(y,Z,D)
    colnames(dat1)<-c("y",paste0("x",1:4),"D")
  }else{
    b1<- (27.4*X[,1] + 13.7*X[,2] +13.7*X[,3] + 13.7*X[,4])
    y1<-210 +b1 + rnorm(n)
    y0<-200-0.5*b1+ rnorm(n)
    y<-D*y1+(1-D)*y0
    #dat1<-data.frame(y,Z,D)
    #colnames(dat1)<-c("y",paste0("x",1:4),"D")
  }
  
  pi.hat<-fitted(glm(D~Z,family=binomial(link="logit")))
  w.hat<-pi.hat/(1-pi.hat)
  
  
 # new.dat0<-dat1%>%filter(D==0,)
  
  #model1<-gam(y~s(x1)+s(x2)+s(x3)+s(x4),data=new.dat0,family = gaussian())
  
  #new.data0<-dat1[,2:5]
  #y.hat0<-predict(model1,newdata=new.data0)
  
  #new.dat1<-dat1%>%filter(D==1,)
  
  
  
  
  #model2<- model1<-gam(y~s(x1)+s(x2)+s(x3)+s(x4),data=new.dat1,family = gaussian())
 
  #y.hat1<-predict(model2,newdata=new.data0)
  
  dat<-data.frame(y=y,Z,D=D,pi.hat=pi.hat)
  #dat<-data.frame(Z,y=y,D=D,y.hat0=y.hat0,y.hat1=y.hat1,pi.hat=pi.hat)
  colnames(dat)[2:(p+1)]<-c(paste0("x",1:p))
  
  return(dat)
}




generate_data_hainmuller<-function(n,OR,PS,seed){
  set.seed(seed)
  
  s<-matrix(c(2,1,-1,1,1,-0.5,-1,-0.5,1),nrow=3,byrow=TRUE)
  Z1 <- mvrnorm(n,rep(0,3),s)
  Z2<-runif(n,-3,3)
  Z3<-rchisq(n,df=1)
  Z4<-rbinom(n,1,0.5)
  Z<-cbind(Z1,Z2,Z3,Z4)
  X <- cbind(exp(Z[,1])/2,(Z[,2]/(1+exp(Z[,1])))+10,(Z[,1]*Z[,3]/25+0.6)^3,(Z[,2]+Z[,4]+20)^2)
  
  
  #X<-scale(X,center=TRUE,scale=TRUE)
  
  if(PS==1){
    
    D=ifelse(Z[,1]+2*Z[,2]-2*Z[,3]-Z[,4]-0.5*Z[,5]+Z[,6]+rnorm(n,0,sqrt(30))>0,1,0)
    
   
  } else {
    px<-1 / (1 + exp(X[,1] - 0.5 * X[,2]*X[,1] - 1*X[,3]^2 + 0.1 * X[,4])) ###Highly non response 0.2% responses
    
    D<-rbinom(n,1,px)
  }
  
  if(OR==1){
    
    y1<-210 + b + rnorm(n)
    y0<-200-0.5*b+ rnorm(n)
    y<-D*y1+(1-D)*y0
    dat1<-data.frame(y,Z,D)
    colnames(dat1)<-c("y",paste0("x",1:4),"D")
  }else{
    y1<-(Z[,1]+Z[,2]+Z[,5])^2+rnorm(n)
    y0<-y1
    y<-D*y1+(1-D)*y0
    #dat1<-data.frame(y,Z,D)
    #colnames(dat1)<-c("y",paste0("x",1:4),"D")
  }
  
  pi.hat<-fitted(glm(D~Z,family=binomial(link="logit")))
  w.hat<-pi.hat/(1-pi.hat)
  
  
  # new.dat0<-dat1%>%filter(D==0,)
  
  #model1<-gam(y~s(x1)+s(x2)+s(x3)+s(x4),data=new.dat0,family = gaussian())
  
  #new.data0<-dat1[,2:5]
  #y.hat0<-predict(model1,newdata=new.data0)
  
  #new.dat1<-dat1%>%filter(D==1,)
  
  
  
  
  #model2<- model1<-gam(y~s(x1)+s(x2)+s(x3)+s(x4),data=new.dat1,family = gaussian())
  
  #y.hat1<-predict(model2,newdata=new.data0)
  
  dat<-data.frame(y=y,Z,D=D,pi.hat=pi.hat)
  #dat<-data.frame(Z,y=y,D=D,y.hat0=y.hat0,y.hat1=y.hat1,pi.hat=pi.hat)
  colnames(dat)[2:(p+1)]<-c(paste0("x",1:p))
  
  return(dat)
}



y.hat<-function(fold0,fold1,q){
  
  y0.fitted<-matrix(nrow = nrow(new.dat0), ncol =1)
  y1.fitted<-matrix(nrow = nrow(new.dat1), ncol =1)
  
  
  #################
  y0.predicted<-matrix(nrow = nrow(new.dat1), ncol =q)
  y1.predicted<-matrix(nrow = nrow(new.dat0), ncol =q)
  
  for(k in 1:q){
    
    ############For the control group
    test_data0 <- new.dat0[fold0[[k]],]
    train_data0 <- new.dat0[-fold0[[k]],]
    
    formula<-as.formula(paste0("y~",paste0("s(",paste0("x",1:p),collapse = "+",")")))
    
   
    model0<-gam(formula,data=train_data0)
    
    
    pred0 <- predict(model0, newdata=as.data.frame(test_data0[,2:(p+1)]))
    
    y0.fitted[fold0[[k]], ] <- pred0
    
    
    
    #############For the treatment group
    test_data1 <- new.dat1[fold1[[k]],]
    train_data1 <- new.dat1[-fold1[[k]],]
    
    

    model1<-gam(formula,data=train_data1)
    
    
    pred1 <- predict(model1, newdata=as.data.frame(test_data1[,-c(1)]))
    
    y1.fitted[fold1[[k]], ] <- pred1 #####predictions for fitted model
    
    
    predU0<-as.numeric(predict(model1, newdata=as.data.frame(new.dat0)))
    y1.predicted[,k]<-predU0
    
    predU1<-as.numeric(predict(model0, newdata=as.data.frame(new.dat1)))
    y0.predicted[,k]<-predU1
  }
  
  y.hat11<-as.numeric(apply(y1.predicted,1,function(a)mean(a)))
  y1.hat<-c(y1.fitted,y.hat11)
  
  y.hat00<-as.numeric(apply(y0.predicted,1,function(a)mean(a)))
  y0.hat<-c(y0.fitted,y.hat00)
  
  
  return(cbind(y1.hat,y0.hat))
}



newton_EL_primal <- function(th,y,y.hat0,y.hat1,pi.hat,D,max.iter=200,eps=1e-6)
{
  iter <- 0
  
  while (TRUE) {
    iter <- iter + 1
    
    theta<-th
    b_theta1<-y.hat1
    g.hat<-1/pi.hat
    H<-cbind(b_theta1,g.hat)
    
    ###Primal
    R=H[D==1,]
    
    p<-sum(D)
    
    w <- Variable(p)
    
    
    objective <- Minimize(-sum(log(w)))
    
    
    A<-cbind(rep(1,nrow(R)),R)
    b<-c(1,colMeans(H))
    
    constraints<-list(t(A)%*%w==b)
    
    
    problem <- Problem(objective, constraints = constraints)
    
    result3 <- solve(problem)
    
    ####control group
    
    b_theta0<-y.hat0
    g.hat0<-1/(1-pi.hat)
    H0<-cbind(b_theta0,g.hat0)
    
    ###Primal
    R0=H0[D==0,]
    
    p0<-sum(1-D)
    
    w0 <- Variable(p0)
    
    
    objective0 <- Minimize(-sum(log(w0)))
    
    
    A0<-cbind(rep(1,nrow(R0)),R0)
    b0<-c(1,colMeans(H0))
    
    constraints_0<-list(t(A0)%*%w0==b0)
    
    
    problem0 <- Problem(objective0, constraints = constraints_0)
    
    result0 <- solve(problem0)
    
    if(result3$status=="optimal" && result0$status=="optimal"){
      w1<-as.numeric(result3$getValue(w))
      w2<-as.numeric(result0$getValue(w0))
      
      th.new<-sum(w1*y[D==1])-sum(w2*y[D==0])
      
      print((max(abs(theta - th.new))))
      if ((max(abs(theta - th.new))< eps)){
        message("Convergence achieved.")
        return(th.new)
      }
      
      if (iter >= max.iter) {
        
        message("Maximum iterations reached. Returning current estimate.")
        return(as.matrix(rep(NA,length(th)))) # Return x if maximum iterations are reached
      } 
      
      
      th<- th.new
    }else{
      message("Status is not optimal. Therefore returining NA values.")
      return(as.matrix(rep(NA,length(th))))
    }
    
    print(iter)
    # print(th.new)
    
  }
  
}


###HD means Hellinger distance
newton_HD_primal <-function(y,pi.hat1,pi.hat0,D,y.hat1,y.hat0,max.iter=200,eps=1e-6)
{
  
  g.hat<--sqrt(pi.hat1)
  H<-cbind(y.hat1,g.hat)
  
  ###Primal
  R=H[1:nrow(new.dat1),]
  
  p<-sum(D)
  
  w <- Variable(p)
  
  
  
  # Define the objective function using kl_div
  objective <- Minimize(-sum(sqrt(w)))
  #objective <- Minimize(-sum((sqrt(w)-1)^2))
  
  
  A<-cbind(rep(1,nrow(R)),R)
  b<-c(1,colMeans(H))
  
  constraints<-list(t(A)%*%w==b)
  
  
  problem <- Problem(objective, constraints = constraints)
  result3 <- solve(problem)
  
  
  
  
  ####control group
  
  
  g.hat0<--sqrt((1-pi.hat0))
  H0<-cbind(y.hat0,g.hat0)
  
  ###Primal
  R0=H0[1:nrow(new.dat0),]
  
  p0<-sum(1-D)
  
  w0 <- Variable(p0)
  
  
  objective0 <- Minimize(-sum(sqrt(w0)))
  #objective0 <- Minimize(-sum((sqrt(w0)-1)^2))
  
  
  A0<-cbind(rep(1,nrow(R0)),R0)
  b0<-c(1,colMeans(H0))
  
  constraints_0<-list(t(A0)%*%w0==b0)
  
  
  problem0 <- Problem(objective0, constraints = constraints_0)
  
  result0 <- solve(problem0)
  
  if(result3$status=="optimal" && result0$status=="optimal"){
    w1<-as.numeric(result3$getValue(w))
    w2<-as.numeric(result0$getValue(w0))
    
    th.new<-sum(w1*y[D==1])-sum(w2*y[D==0])
    
    
    return(th.new)
  }else{
    message("Status is not optimal. Therefore returining NA values.")
    return(NA)
  }
}



newton_ET_primal <-function(y,pi.hat1,pi.hat0,D,y.hat1,y.hat0,max.iter=200,eps=1e-6)
{
  
  g.hat<-log(pi.hat1)
  H<-cbind(y.hat1,g.hat)
  
  ###Primal
  R=H[D==1,]
  
  R=H[1:nrow(new.dat1),]
  
  p<-sum(D)
  
  w <- Variable(p)
  
  
  ones <- rep(1, p)
  
  # Define the objective function using kl_div
  objective <- Minimize(sum(kl_div(w, ones)-ones)) 
  
  A<-cbind(rep(1,nrow(R)),R)
  b<-c(1,colMeans(H))
  
  constraints<-list(t(A)%*%w==b)
  
  
  problem <- Problem(objective, constraints = constraints)
  result3 <- solve(problem)
  
  
  ####control group
  g.hat0<-log((1-pi.hat0))
  H0<-cbind(y.hat0,g.hat0)
  
  ###Primal
  #R0=H0[D==0,]
  
  R0=H0[1:nrow(new.dat0),]
  
  p0<-sum(1-D)
  
  w0 <- Variable(p0)
  ones0 <- rep(1, p0)
  
  objective0 <- Minimize(sum(kl_div(w0, ones0)-ones0)) 
  
  
  A0<-cbind(rep(1,nrow(R0)),R0)
  b0<-c(1,colMeans(H0))
  
  constraints_0<-list(t(A0)%*%w0==b0)
  
  
  problem0 <- Problem(objective0, constraints = constraints_0)
  
  result0 <- solve(problem0)
  
  if(result3$status=="optimal" && result0$status=="optimal"){
    w1<-as.numeric(result3$getValue(w))
    w2<-as.numeric(result0$getValue(w0))
    
    th.new<-sum(w1*y[D==1])-sum(w2*y[D==0])
    return(th.new)
  }else{
    message("Status is not optimal. Therefore returining NA values.")
    return(NA)
  }
  
}


##########################################################################################################
#############################################################################################################
############################################################################################################
######################################################Fold functions#####################################
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

k_fold_function<-function(df,K,idx_S,idx_U0,seed){
  res <- build_SU_folds(df, K , id_col = "ID",idx_S,idx_U0, seed = seed)
  df$ID <- seq_len(nrow(df))
  lab_idx_all <- which(df$D == 1)
  
  data_unlabeled<-list()
  for(k in 1:K){
    train_idx_lab <- setdiff(lab_idx_all, res$folds[[k]]$S_k_idx)
    validation_data_labeled <- df[res$folds[[k]]$S_k_idx,]
    validation_data_unlabeled<-df[res$folds[[k]]$U_k_idx,]
    
    train_data_labeled <-df[train_idx_lab,] ###get the data frame for trained data
    
    formula<-as.formula(paste0("y~",paste0("s(",colnames( subset(train_data_labeled, select = -c(D, y, ID,pi.hat))),")",collapse = "+")))
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





estimate_theta_EM_kfold_CVXR_ET <- function(th, data_full, K, D,seed, max.iter=50, eps=1e-6) {
  theta <- as.matrix(th)
  iter  <- 0
  
  repeat {
    iter <- iter + 1
    
    # --- prepare cross-fitted data for treatment group T = 1, stacked in fold order ---
    fold_T1 <- k_fold_function(df = data_full, K = K, idx_S = which(D==1), idx_U0 = which(D==0), seed = seed)
    data_T1 <- do.call(rbind, fold_T1)
    n_fold_T1 <- sapply(fold_T1, nrow)
    offsets_T1 <- c(0, cumsum(n_fold_T1))
    
    # --- features ---
    x_cols_T1 <- grep("^x\\d+$", names(data_T1), value = TRUE)
    p_T1 <- length(x_cols_T1)
    X_T1 <- as.matrix(data_T1[, x_cols_T1, drop = FALSE])
    X1_T1 <- cbind(1, X_T1)  # add intercept
    
    # --- fitted values and treatment indicators ---
    yhat_T1 <- data_T1$y.hat
    T1 <- data_T1$D
    idx_T1 <- which(T1 == 1)  # indices for treatment group (T = 1)
    
    
    # --- moment matrix for calibration (select what to balance) ---
    A_T1 <- cbind(yhat_T1)
    
    # --- observed outcomes for treatment group ---
    y_T1 <- data_T1$y[idx_T1]
    
    
    
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
    fold_T0 <- k_fold_function(df = data_full, K = K, idx_S = which(D==0), idx_U0 = which(D==1), seed = seed)
    data_T0 <- do.call(rbind, fold_T0)
    n_fold_T0 <- sapply(fold_T0, nrow)
    offsets_T0 <- c(0, cumsum(n_fold_T0))
    
    # --- features ---
    x_cols_T0 <- grep("^x\\d+$", names(data_T0), value = TRUE)
    p_T0 <- length(x_cols_T0)
    X_T0 <- as.matrix(data_T0[, x_cols_T0, drop = FALSE])
    X1_T0 <- cbind(1, X_T0)  # add intercept
    
    # --- fitted values and treatment indicators ---
    yhat_T0 <- data_T0$y.hat
    T0 <- data_T0$D
    idx_T0 <- which(T0 == 0)  # indices for control group (T = 0)
    
    
    # --- moment matrix for calibration (select what to balance) ---
    A_T0 <- cbind(yhat_T0)
    
    # --- observed outcomes for control group ---
    y_T0 <- data_T0$y[idx_T0]
    
    
    
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
      
      th.new<-sum(w1*y_T1)-sum(w0*y_T0)
      
      
      return(th.new)
    }else{
      message("Status is not optimal. Therefore returining NA values.")
      return(NA)
    }
  }
}



estimate_theta_EM_kfold_CVXR_HD <- function(th, data_full, K, D,seed, max.iter=50, eps=1e-6) {
  theta <- as.matrix(th)
  iter  <- 0
  
  repeat {
    iter <- iter + 1
    
    # --- prepare cross-fitted data for treatment group T = 1, stacked in fold order ---
    fold_T1 <- k_fold_function(df = data_full, K = K, idx_S = which(D==1), idx_U0 = which(D==0), seed = seed)
    data_T1 <- do.call(rbind, fold_T1)
    n_fold_T1 <- sapply(fold_T1, nrow)
    offsets_T1 <- c(0, cumsum(n_fold_T1))
    
    # --- features ---
    x_cols_T1 <- grep("^x\\d+$", names(data_T1), value = TRUE)
    p_T1 <- length(x_cols_T1)
    X_T1 <- as.matrix(data_T1[, x_cols_T1, drop = FALSE])
    X1_T1 <- cbind(1, X_T1)  # add intercept
    
    # --- fitted values and treatment indicators ---
    yhat_T1 <- data_T1$y.hat
    T1 <- data_T1$D
    idx_T1 <- which(T1 == 1)  # indices for treatment group (T = 1)
    
    
    # --- moment matrix for calibration (select what to balance) ---
    A_T1 <- cbind(yhat_T1)
    
    # --- observed outcomes for treatment group ---
    y_T1 <- data_T1$y[idx_T1]
    
    
    
    # --- CVXR problem on treatment group only ---
    r <- ncol(A_T1)
    w1 <- CVXR::Variable(length(idx_T1 ), pos = TRUE)
    a <- rep(1/length(I1), length(I1))  # base (uniform)
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
    
    
    objective_T1 <-  Minimize(-sum(sqrt(w1)))
    
    
    prob_T1 <- CVXR::Problem(objective_T1, constr_T1)
    sol_T1  <- CVXR::solve(prob_T1)
    #############################################################################################
    #############################################################################################
    ############################################################################################
    
    # --- prepare cross-fitted data for CONTROL group T = 0, stacked in fold order ---
    fold_T0 <- k_fold_function(df = data_full, K = K, idx_S = which(D==0), idx_U0 = which(D==1), seed = seed)
    data_T0 <- do.call(rbind, fold_T0)
    n_fold_T0 <- sapply(fold_T0, nrow)
    offsets_T0 <- c(0, cumsum(n_fold_T0))
    
    # --- features ---
    x_cols_T0 <- grep("^x\\d+$", names(data_T0), value = TRUE)
    p_T0 <- length(x_cols_T0)
    X_T0 <- as.matrix(data_T0[, x_cols_T0, drop = FALSE])
    X1_T0 <- cbind(1, X_T0)  # add intercept
    
    # --- fitted values and treatment indicators ---
    yhat_T0 <- data_T0$y.hat
    T0 <- data_T0$D
    idx_T0 <- which(T0 == 0)  # indices for control group (T = 0)
    
    
    # --- moment matrix for calibration (select what to balance) ---
    A_T0 <- cbind(yhat_T0)
    
    # --- observed outcomes for control group ---
    y_T0 <- data_T0$y[idx_T0]
    
    
    
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
    
    
    objective_T0 <- Minimize(-sum(sqrt(w0)))
    
    
    prob_T0 <- CVXR::Problem(objective_T0, constr_T0)
    sol_T0  <- CVXR::solve(prob_T0)
    
    
    
    if(sol_T1$status=="optimal" && sol_T0$status=="optimal"){
      w1<-as.numeric(sol_T1$getValue(w1))
      w0<-as.numeric(sol_T0$getValue(w0))
      
      th.new<-sum(w1*y_T1)-sum(w0*y_T0)
      
      
      return(th.new)
    }else{
      message("Status is not optimal. Therefore returining NA values.")
      return(NA)
    }
  }
}



estimate_theta_EM_kfold_CVXR_HD_1 <- function(theta1,theta0,data_full, K, D,seed, max.iter=50, eps=1e-4) {
  #theta <- as.matrix(th)
  iter  <- 0
  
  repeat {
    iter <- iter + 1
    theta<-theta1-theta0
    # --- prepare cross-fitted data for treatment group T = 1, stacked in fold order ---
    fold_T1 <- k_fold_function(df = data_full, K = K, idx_S = which(D==1), idx_U0 = which(D==0), seed = seed)
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
    T1 <- data_T1$D
    idx_T1 <- which(T1 == 1)  # indices for treatment group (T = 1)
    
    
    # --- moment matrix for calibration (select what to balance) ---
    #A_T1 <- cbind(yhat_T1, -sqrt((data_T1$pi.hat)))
    A_T1 <- cbind(yhat_T1)
    
    # --- observed outcomes for treatment group ---
    y_T1 <- data_T1$y[idx_T1]
    
    
    
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
    
    
    objective_T1 <-  Minimize(-sum(sqrt(w1)))
    
    
    prob_T1 <- CVXR::Problem(objective_T1, constr_T1)
    sol_T1  <- CVXR::solve(prob_T1)
    #############################################################################################
    #############################################################################################
    ############################################################################################
    
    # --- prepare cross-fitted data for CONTROL group T = 0, stacked in fold order ---
    fold_T0 <- k_fold_function(df = data_full, K = K, idx_S = which(D==0), idx_U0 = which(D==1), seed = seed)
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
    T0 <- data_T0$D
    idx_T0 <- which(T0 == 0)  # indices for control group (T = 0)
    
    
    # --- moment matrix for calibration (select what to balance) ---
    A_T0 <- cbind(yhat_T0)
    
    # --- observed outcomes for control group ---
    y_T0 <- data_T0$y[idx_T0]
    
    
    
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
    
    
    objective_T0 <- Minimize(-sum(sqrt(w0)))
    
    
    prob_T0 <- CVXR::Problem(objective_T0, constr_T0)
    sol_T0  <- CVXR::solve(prob_T0)
    
    
    
    if(sol_T1$status=="optimal" && sol_T0$status=="optimal"){
      w1<-as.numeric(sol_T1$getValue(w1))
      w0<-as.numeric(sol_T0$getValue(w0))
      theta1.new<-sum(w1*y_T1);theta0.new<-sum(w0*y_T0)
      #th.new<-sum(w1*y_T1)-sum(w0*y_T0)
      th.new<-theta1.new-theta0.new
      if (max(abs(th.new - as.numeric(theta))) < eps) {
        return(matrix(th.new, ncol=1))
      }
      if (iter >= max.iter) {
        message("Maximum iterations reached.")
        return(matrix(NA, ncol=1))
      }
      
      theta1<-theta1.new;theta0<-theta0.new
      
    }
  }
}






estimate_theta_EM_kfold_CVXR_ET_1 <- function(theta1,theta0,data_full, K, D,seed, max.iter=50, eps=1e-4) {
  #theta <- as.matrix(th)
  iter  <- 0
  
  repeat {
    iter <- iter + 1
    theta<-theta1-theta0
    # --- prepare cross-fitted data for treatment group T = 1, stacked in fold order ---
    fold_T1 <- k_fold_function(df = data_full, K = K, idx_S = which(D==1), idx_U0 = which(D==0), seed = seed)
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
    T1 <- data_T1$D
    idx_T1 <- which(T1 == 1)  # indices for treatment group (T = 1)
    
    
    # --- moment matrix for calibration (select what to balance) ---
    A_T1 <- cbind(yhat_T1)
    
    # --- observed outcomes for treatment group ---
    y_T1 <- data_T1$y[idx_T1]
    
    
    
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
    fold_T0 <- k_fold_function(df = data_full, K = K, idx_S = which(D==0), idx_U0 = which(D==1), seed = seed)
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
    T0 <- data_T0$D
    idx_T0 <- which(T0 == 0)  # indices for control group (T = 0)
    
    
    # --- moment matrix for calibration (select what to balance) ---
    A_T0 <- cbind(yhat_T0)
    
    # --- observed outcomes for control group ---
    y_T0 <- data_T0$y[idx_T0]
    
    
    
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
      if (max(abs(th.new - as.numeric(theta))) < eps) {
        return(matrix(th.new, ncol=1))
      }
      if (iter >= max.iter) {
        message("Maximum iterations reached.")
        return(matrix(NA, ncol=1))
      }
      
      theta1<-theta1.new;theta0<-theta0.new
      
    }
  }
}




estimate_theta_EM_kfold_CVXR_ET_2 <- function(theta1,theta0,fold_T1,fold_T0,K, D,max.iter=50, eps=1e-4) {
  #theta <- as.matrix(th)
  iter  <- 0
  
  repeat {
    iter <- iter + 1
    theta<-theta1-theta0
    # --- prepare cross-fitted data for treatment group T = 1, stacked in fold order ---
    
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
    T1 <- data_T1$D
    idx_T1 <- which(T1 == 1)  # indices for treatment group (T = 1)
    
    
    # --- moment matrix for calibration (select what to balance) ---
    A_T1 <- cbind(yhat_T1)
    
    # --- observed outcomes for treatment group ---
    y_T1 <- data_T1$y[idx_T1]
    
    
    
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
    T0 <- data_T0$D
    idx_T0 <- which(T0 == 0)  # indices for control group (T = 0)
    
    
    # --- moment matrix for calibration (select what to balance) ---
    A_T0 <- cbind(yhat_T0)
    
    # --- observed outcomes for control group ---
    y_T0 <- data_T0$y[idx_T0]
    
    
    
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
      if (max(abs(th.new - as.numeric(theta))) < eps) {
        return(matrix(th.new, ncol=1))
      }
      if (iter >= max.iter) {
        message("Maximum iterations reached.")
        return(matrix(NA, ncol=1))
      }
      
      theta1<-theta1.new;theta0<-theta0.new
      
    }
  }
}


































