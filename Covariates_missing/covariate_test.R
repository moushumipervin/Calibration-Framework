  
    rm(list = ls())
  
    library(dplyr)
    library(foreach)
    library(doParallel)
    library(CVXR)
    library(MASS)
    library(dplyr)
  
   setwd("F:/Research/Functional calibration/Covariate missing/Final simulation code")
    
   source("kfold_functions.R")
      #source("Newton_test.R")
      
   seed=1234; m=500;k=3;p=2
    
   beta_HT<-matrix(NA,nrow=m,ncol=p+1)
    
   beta_AIPW<-matrix(NA,nrow=m,ncol=p+1)
   beta_AIPW1<-matrix(NA,nrow=m,ncol=p+1)
    
   beta_ET_EM<-matrix(NA,nrow=m,ncol=p+1)
   beta_HD_EM<-matrix(NA,nrow=m,ncol=p+1)
    
   beta_HD<-matrix(NA,nrow=m,ncol=p+1)
   beta_cc<-matrix(NA,nrow=m,ncol=p+1)
   beta_full<-matrix(NA,nrow=m,ncol=p+1)
   
   
   HD_low  <- matrix(NA,nrow=m,ncol=p+1)
   HD_high<-matrix(NA,nrow=m,ncol=p+1)
   
   full_low  <- matrix(NA,nrow=m,ncol=p+1)
   full_high<-matrix(NA,nrow=m,ncol=p+1)
   
   cc_low  <- matrix(NA,nrow=m,ncol=p+1)
   cc_high<-matrix(NA,nrow=m,ncol=p+1)
    for(i in 1:m){
      
      
      
      #######################################################  
      ################Data generation########################
      
      ####Complete case analysis
      set.seed(i+seed)
      dat<-generate_data(n=1000,PS=2,MAR=1)
      
      new.dat1<-dat[dat$D==1,]
      new.dat0<-dat[dat$D==0,]
      row.names(new.dat1) <- NULL
      row.names(new.dat0) <- NULL
      ###########################################################
      ##################Complete data model when delta =1 ######################
      model_cc<-lm(y~x+z,data=new.dat1)
     
      S1<-summary( model_cc)
      
      beta_cc[i,]<-as.numeric(S1$coefficients[,"Estimate"])
      cc_low[i,]<-beta_cc[i,]-1.96*as.numeric(S1$coefficients[,"Std. Error"])
      cc_high[i,]<-beta_cc[i,]+1.96*as.numeric(S1$coefficients[,"Std. Error"])
      #######################################################  
      ################Full data model##################
      
      MODEL<-lm(y~x+z,data=dat)
      S<-summary(MODEL)
      
      beta_full[i,]<-as.numeric(S$coefficients[,"Estimate"])
      full_low[i,]<-beta_full[i,]-1.96*as.numeric(S$coefficients[,"Std. Error"])
      full_high[i,]<-beta_full[i,]+1.96*as.numeric(S$coefficients[,"Std. Error"])
      #######################################################  
      ################HT estimator###########################
      
      model_HT<-lm(y~x+z,data=new.dat1,weights = 1/(new.dat1$D.hat))
      beta_HT[i,]<-as.numeric(coefficients(model_HT))
      
    
      
      ######################################################
      #######################Compute z.hat##################
      fold_hat <- k_fold_function(df = dat, K = k, seed =  (i+seed))
      data_all <- do.call(rbind, fold_hat)
     
      z.hat <- data_all$z.hat
      D     <- data_all$D
      I1    <- which(D == 1)
      y<-data_all$y
      n<-nrow(data_all)
      
      dat.new<-data_all[I1,]
      
      Q0<-cbind(1,dat.new$x,dat.new$z)
      Q1<-cbind(1,dat.new$x,dat.new$z.hat)
      
      Q2<-cbind(1,data_all$x,data_all$z.hat)
      #variance_theta_diag(y=y, X=, w_hat=, theta_hat, ci_level = 0.95)
      
      #######################################################################
      #################################AIPW Without N scale##################
      N<-nrow(data_all)
      N_hat<-sum(data_all$D/data_all$D.hat)
      A1<--(t(Q1)%*%diag((1/(N_hat*dat.new$D.hat)))%*%Q1)+(t(Q0)%*%diag((1/(N_hat*dat.new$D.hat)))%*%Q0)+((t(Q2)%*%Q2)/N)
      
      B1<-(t(Q0)%*%diag((1/(N_hat*dat.new$D.hat)))%*%dat.new$y)-(t(Q1)%*%diag(((1/(N_hat*dat.new$D.hat))))%*%dat.new$y)+((t(Q2)%*%data_all$y)/N)
      
      beta_AIPW1[i,]<-c(solve(A1)%*%B1)
      
    
      
    
   #rr=cbind(beta_cc,beta_HT,beta_AIPW1)
    #apply(rr,2,sd)
      
      ####################################################
      ########################Proposed method#############
      ###Empirical Tilting
      
      #beta_ET_EM[i,]<-estimate_theta_EM_kfold_CVXR_ET(th=beta_HT[i,], data_full=dat, K=k, seed=(i+seed), max.iter=30, eps=1e-4)$theta  
      #variance_theta_diag(y=data_all$y[I1], X= new.data1[I1,], w_hat=res_HD$w, theta_hat=res_HD$theta, ci_level = 0.95) 
      #res_HD$estimate.table
      ###HD
      res_HD<-estimate_theta_EM_kfold_CVXR_HD(th=beta_AIPW1[i,], data_full=dat, K=k, seed=(i+seed), max.iter=30, eps=1e-4)
      
      beta_HD_EM[i,]<-res_HD$theta
      
      if (is.null(res_HD$estimate.table)) {
        # no result → fill with NA
        HD_low[i, ]  <- rep(NA_real_, p + 1)
        HD_high[i, ] <- rep(NA_real_, p + 1)
      } else {
        # have result → use CI columns
        HD_low[i, ]  <- res_HD$estimate.table[, "CI_Lower"]
        HD_high[i, ] <- res_HD$estimate.table[, "CI_Upper"]
      }
      #beta_ET_nested[i,]<-estimate_theta_nested_kfold_CVXR_ET (theta_init=beta_HT[i,], data_full=dat, K=k , seed =(i+seed  ))$theta_hat
    }
   
      
  
   
   
   
   
   
   
   
   
   
   
  
  
