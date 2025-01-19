  
    rm(list = ls())
  
    library(dplyr)
    library(foreach)
    library(doParallel)
    library(CVXR)
    library(MASS)
    library(dplyr)
  
    n_cores=10
    core_seeds<-sample(10^5,n_cores)
    covariate_missing_function<-function(seed,n_cores,n_rep,n,OR,PS,MAR){
    
    
      set.seed(seed)
      m<-n_rep/n_cores
    
    
    
    
    source("Newton_test.R")
    
    result_th1<-data.frame()
    
    result_th2<-data.frame()
    
    result_th3<-data.frame()
    
    
    for(i in 1:m){
      
      
      
      #######################################################  
      ################Data generation########################
      
      ####Complete case analysis
      set.seed(i+seed)
      dat<-generate_data(n=500,PS=1,MAR=1)
      
      new.dat1<-dat[dat$D==1,]
      new.dat0<-dat[dat$D==0,]
      row.names(new.dat1) <- NULL
      row.names(new.dat0) <- NULL
      ###########################################################
      ##################Complete data model when delta =1 ######################
      model_cc<-lm(y~x+z,data=new.dat1)
      beta_cc<-as.numeric(coefficients(model_cc))
      th1_cc<-beta_cc[1]
      th2_cc<-beta_cc[2]
      th3_cc<-beta_cc[3]
      
      
      #######################################################  
      ################Full data model##################
      
      MODEL<-lm(y~x+z,data=dat)
      
      beta_full<-as.numeric(coefficients(MODEL))
      b1<-beta_full[1]
      b2<-beta_full[2]
      b3<-beta_full[3]
      
      
      #######################################################  
      ################HT estimator###########################
      
      model_HT<-lm(y~x+z,data=new.dat1,weights = 1/(new.dat1$D.hat))
      beta_HT<-as.numeric(coefficients(model_HT))
      
      
      th1_HT<-beta_HT[1]
      th2_HT<-beta_HT[2]
      th3_HT<-beta_HT[3]
      
      ######################################################
      #######################Compute z.hat##################
      q=5
      fold<-createFolds(1:nrow(new.dat1), k = q, list = TRUE, returnTrain = FALSE)
      
      z.hat<-z.hat.function(fold=fold,q=q,new.dat0=new.dat0,new.dat1=new.dat1)
      
      dat.new=cbind(rbind(new.dat1,new.dat0),z.hat)
      
      
      
      ####################################################
      #####################AIPW with N scaled###########################
      
      Q0<-cbind(1,dat.new$x,dat.new$z)
      Q1<-cbind(1,dat.new$x,dat.new$z.hat)
      
      n=nrow(dat)
      N.scale<-n/sum(dat.new$D/dat.new$D.hat)
      A<-(t(Q1)%*%diag(1-(N.scale*(N.scale*(dat.new$D/dat.new$D.hat))))%*%Q1)+(t(Q0)%*%diag((N.scale*(dat.new$D/dat.new$D.hat)))%*%Q0)
      
      B<-(t(Q0)%*%diag(N.scale*(dat.new$D/dat.new$D.hat))%*%dat.new$y)+(t(Q1)%*%diag(1-(N.scale*(dat.new$D/dat.new$D.hat)))%*%dat.new$y)
      
      beta_AIPW<-solve(A)%*%B
      
      th1_AIPW<-beta_AIPW[1]
      th2_AIPW<-beta_AIPW[2]
      th3_AIPW<-beta_AIPW[3]
      
      
      
      #######################################################################
      #################################AIPW Without N scale##################
      A1<-(t(Q1)%*%diag(1-(dat.new$D/dat.new$D.hat))%*%Q1)+(t(Q0)%*%diag((dat.new$D/dat.new$D.hat))%*%Q0)
      
      B1<-(t(Q0)%*%diag((dat.new$D/dat.new$D.hat))%*%dat.new$y)+(t(Q1)%*%diag(1-((dat.new$D/dat.new$D.hat)))%*%dat.new$y)
      
      beta_AIPW1<-solve(A1)%*%B1
      
      
      th1_AIPW1<-beta_AIPW1[1]
      th2_AIPW1<-beta_AIPW1[2]
      th3_AIPW1<-beta_AIPW1[3]
      
      
      
      
      ####################################################
      ########################Proposed method#############
      ###Empirical Tilting
      
      beta_ET<-newton_ET_primal(th=beta_HT,x=dat.new$x,y=dat.new$y,z=dat.new$z,D.hat=dat.new$D.hat,z.hat=dat.new$z.hat,D=dat.new$D)
      
      th1_ET<-beta_ET[1]
      th2_ET<-beta_ET[2]
      th3_ET<-beta_ET[3]
      
      
      ###Empirical likelihood
      
      #beta_EL<-newton_EL_primal(th=beta_HT,x=dat.new$x,y=dat.new$y,z=dat.new$z,D.hat=dat.new$D.hat,z.hat=dat.new$z.hat,D=dat.new$D)
      
      
      #th1_EL<-beta_EL[1]
      #th2_EL<-beta_EL[2]
      #th3_EL<-beta_EL[3]
      
      
      ###Hellinger distance
      
      beta_HD<-newton_HD_primal(th=beta_HT,x=dat.new$x,y=dat.new$y,z=dat.new$z,D.hat=dat.new$D.hat,z.hat=dat.new$z.hat,D=dat.new$D)
      
      
      th1_HD<-beta_HD[1]
      th2_HD<-beta_HD[2]
      th3_HD<-beta_HD[3]
      
      ######################################################################################################
      #############################################Record the results#######################################
      
      result1<-data.frame(th1_full=b1,th1_cc=th1_cc,th1_HT=th1_HT,th1_AIPWN=th1_AIPW,th1_AIPW=th1_AIPW1,th1_HD=th1_HD,th1_ET=th1_ET)
      
      result2<-data.frame(th2_full=b2,th2_cc=th2_cc,th2_HT=th2_HT,th2_AIPWN=th2_AIPW,th2_AIPW=th2_AIPW1,th2_HD=th2_HD,th2_ET=th2_ET)
      
      result3<-data.frame(th3_full=b3,th3_cc=th3_cc,th3_HT=th3_HT,th3_AIPWN=th3_AIPW,th3_AIPW=th3_AIPW1,th3_HD=th3_HD,th3_ET=th3_ET)
      
      
      
      result_th1<-rbind(result_th1,result1)
      result_th2<-rbind(result_th2,result2)
      result_th3<-rbind(result_th3,result3)
      
    }
    result<-cbind(result_th1,result_th2,result_th3)
    
    return(result)
    
  }
  
  registerDoParallel(n_cores)
  theta_estimated<-foreach(r=1:n_cores,.combine = 'rbind',.packages = c('CVXR','dplyr','MASS','caret','mgcv'))%dopar%{
    
   covariate_missing_function(seed=core_seeds[r],n_cores,n_rep=1000,n=500,PS=1,MAR=1)
   
  }
  
  stopImplicitCluster()
  
  #saveRDS(theta_estimated,"Theta_MAR_n500_lowselectionbias.RDS")
  
    
    
  
  
  
  
  
