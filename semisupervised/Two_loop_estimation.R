# Softmax weights for given lambda and theta


# Define moment function: U(theta; z_i)
rm(list=ls())
source("two_loop_function.R")

# Simulation parameters
library(mgcv)
rep=20
p <- 3
theta_hat<-matrix(NA,nrow=rep,ncol=p)
theta_hat_lm<-matrix(NA,nrow=rep,ncol=p)
theta_ET<-matrix(NA,nrow=rep,ncol=p)
system.time(
for(k in 1:rep){
  set.seed(k+123)
  
  
  ## ----- toy data -----
  n <- 1000
  p <- 3
  
  X <- cbind(1, matrix(rnorm(n * (p-1)), n, p-1))  # include intercept
  beta_true <- c(1, 2, -1)
  y <- as.vector(X %*% beta_true+exp(X)%*%beta_true + rnorm(n, sd = 1.0))
  
  
  n1=500
 # D=rbinom(n,1,n1/(n))
  beta1<-rep(.5,p)
  px=1/(1+exp(-(-1.5+X%*%beta1)))
  #px<-1/(1+exp(X[,1]-X[,2]+.25*X[,3]+0.1*X[,4]))
  D=rbinom(n,1,px)
  
  data_labelled=cbind(y[which(D==1)],X[which(D==1),])
  data_unlabelled=cbind(X[which(D==0),])
  colnames(data_labelled)<-c("y",paste0("X",1:(p)))
  colnames(data_unlabelled)<-c(paste0("X",1:(p)))
  formula<-as.formula(paste0("y~",paste0("s(",colnames(data_labelled[,-c(1,2)]),")",collapse = "+")))
  
  gam_model<-gam(formula,data=as.data.frame(data_labelled[,-2]))
  gam.predict<-predict(gam_model,newdata = as.data.frame(data_unlabelled))
  y.hat<-c(fitted(gam_model),gam.predict)
  
  formula1<-as.formula(paste0("y~",paste0(colnames(data_labelled[,-c(1,2)]),collapse = "+")))
  
  lm_model<-lm(formula1,data=as.data.frame(data_labelled[,-2]))
  theta_hat_lm[k,] <-as.numeric(coefficients(lm_model))
  
  # Initial guess for theta
  theta_init <-   theta_hat_lm[k,]
  
  # Run two-loop ET estimator
  theta_hat[k, ] <- get_theta(theta_init, data_labelled, data_unlabelled, y.hat)
  theta_ET[k, ]<-newton_ET_primal_test(theta_init, y.hat,data_labelled,data_unlabelled)
  
})

