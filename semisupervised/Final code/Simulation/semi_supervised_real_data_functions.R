
#setwd("C:/Users/moushumi/Desktop/Codes_Reproducibility/Final semi-supervised code")
library(MASS)
library(mgcv)
library(CVXR)
library(caret)
library(dplyr)

######################################################################################################################################
####################################################PSSE function#####################################################################################################################
###############(1.1) polynomial generation function 

polynomial<- function(data,order) 
{
  #---------------------------Arguments----------------------------------------#
  # Purpose: This function is to produce the polynomial basis of X including 
  #           the intercept vector. 
  #
  # Input: 
  #       data: A matrix, whose each row is an observation of predictor vector.
  #       order: The polynomial order.
  #----------------------------------------------------------------------------#
  polynomial_Z<- rep(1,nrow(data))
  for (i in 1:order)
  {
    polynomial_i <- data^i
    polynomial_Z <- cbind(polynomial_Z,polynomial_i)
  }
  return(polynomial_Z)
}


###############(1.2) the data driven selector GBIC_ppo  
GBIC_ppo <- function (Y_new,covariates_matrix,order_poly){
  #-----------------------Arguments--------------------------------------------------# 
  # Purpose: GBIC_ppo criterion function 
  #
  # Input:
  #       Y_new: The first derivatives of the loss function L.
  #       covariates_matrix: A matrix, each row is an observation of predictor vector.
  #       order_poly: The polynomial order.
  #
  #---------------------------------------------------------------------------------#
  d <- ncol(Y_new)
  n <- nrow(Y_new)
  p <- ncol(covariates_matrix)
  
  polynomial_matrix <- polynomial(covariates_matrix, order_poly)
  
  gammahat <- lm(Y_new~polynomial_matrix-1)$coefficients
  residual_square <- (Y_new - polynomial_matrix %*% gammahat)^2
  sigmahat_square <- colSums(residual_square)/(n-p*order_poly-1)
  design_matrix <- t(polynomial_matrix)%*%polynomial_matrix/n
  
  
  if(class(try(solve(design_matrix),silent=T))[1]=="try-error"){
    
    GBIC<- 9999999
    
  }else {
    trace_det_AB = numeric()
    for (j in 1:d){
      Bhat_j<- t(polynomial_matrix*residual_square[,j])%*%polynomial_matrix/n
      cova_constrast <- 1/sigmahat_square[j]*solve(design_matrix)%*%Bhat_j
      trace_AB <- sum(diag(cova_constrast))
      det_AB <- det(cova_constrast)
      trace_det_AB=c(trace_det_AB,trace_AB-log(det_AB))
    }
    GBIC = 1/n*(d*(n-p*order_poly-1)+ n*sum(log(sigmahat_square)) + d*log(n)*p*order_poly+sum(trace_det_AB))
  }
  
  return(GBIC)
  
}


#################(1.3) the supervised estimate based on the given samples
supervised_one<- function(labelled_data){
  #---------------------Arguments-------------------------------------------#
  # Purpose: According the type of loss function, this function is used to  
  #         get the estimate for the target parameter only based on the given 
  #         labelled samples.
  #------------------------------------------------------------------------#
  
  hattheta_supervised <- lm(labelled_data[,1]~labelled_data[,-1])$coefficients
  
  return(hattheta_supervised)
}




#################(1.4)the first derivative of loss function L
L_first_derivative<- function(labelled_data,hattheta_supervised)
{
  #---------------------Arguments-------------------------------------------#
  # Purpose: According the type of loss function, this function is used to  
  #         determine the first-order derivatives of the loss function based on 
  #         an estimate of the target parameter and a labelled data set.
  #------------------------------------------------------------------------#
  n=nrow(labelled_data)
  
  
  residuals<- as.vector(labelled_data[,1]-cbind(rep(1,n),labelled_data[,-1])%*%hattheta_supervised)
  L_first_derivative <- residuals*cbind(rep(1,n),labelled_data[,-1])
  
  
  return(L_first_derivative)
}



PSSE1 <- function(labelled_data,unlabelled_data,c1=NULL,
                  type="linear",tau=0.5,alpha=NULL,gamma=10,sd=FALSE,Kfolds=5)
{
  #-----------------------------------------Arguments-----------------------------------------------------#
  # Purpose: This function is to implement our proposed semi-supervised method. 
  #           This function is the main function, and it also involves several sub-functions defined in
  #           (1.1) - (1.6) below.
  #
  # Input: 
  #      labelled_data: A matrix, whose each row is an observation of predictor vector and the response
  #                      variable and its first column is the observation vector of the response variable.  
  #      unlabelled_data: A matrix, whose each row is an observation of predictor vector. 
  #      c1: Weight. It is the weight to balance the contribution of labelled data and unlabelled data.
  #           Default value is nrow(labelled_data)/(nrow(labelled_data)+nrow(unlabelled_data)).
  #      type: A specification for the type of loss function $L$. We offer three choices: 
  #             "linear", "logistic", "quantile". Default value is "linear". 
  #      tau: The quantile level, which is only useful for "type=quantile". Default value is 0.5.
  #      alpha: The polynomial order. If it is NULL, then we use the data-driven selector GBIC_ppo to 
  #              determine the polynomial order. Default is NULL. 
  #      gamma: The maximum of possible polynomial order. Default is 10. Only if alpha is NULL, gamma will 
  #              be usful for the selection of polynomial order in the construction of Z. 
  #      sd: Return the pointwise standard deviation estimate based on the available samples? Default 
  #           value is FALSE.  
  #      Kfolds: K-folds cross validation method to avoiding the over-fitting during estimating variance, 
  #          which is only useful When "sd=TRUE".
  #
  # Output: 
  #       Hattheta: The estimate for the target parameter. 
  #       sd.of.hattheta: If "sd=TRUE", the elementwise standard deviation estimate for the estimator; 
  #                       otherwise, NULL.
  #       alpha: The selected polynomial order. 
  #------------------------------------------------------------------------------------------------------#
  
  # some basic parameters 
  n=nrow(labelled_data)
  N=nrow(unlabelled_data)
  p=ncol(labelled_data)-1
  
  
  if (p!=ncol(unlabelled_data)){
    print("the dimension of the labelled and unlabelled data is unmatched")
    
    return(NULL)
  }
  
  if(is.null(c1)==TRUE){
    c1=n/(n+N)
  }
  
  # whether the polynomial order needs to be selected or not 
  if (is.null(alpha)==TRUE){ 
    
    # data splitting 
    labelled_data_part1<- labelled_data[1:round(n/2),]
    labelled_data_part2<- labelled_data[(round(n/2)+1):n,]
    
    
    # the supervised estimate
    hattheta_supervised_part1<- supervised_one(labelled_data_part1)
    hattheta_supervised_part2<- supervised_one(labelled_data_part2)
    
    
    # the first derivative of the loss function L using the supervised estimate 
    L_first_derivative_part1<- L_first_derivative(labelled_data_part1,hattheta_supervised_part2)
    L_first_derivative_part2<- L_first_derivative(labelled_data_part2,hattheta_supervised_part1)
    
    
    # new data integrating the first derivatives of loss function
    Ynew<- rbind(L_first_derivative_part1,L_first_derivative_part2)
    covariates_labelled<- labelled_data[,-1]
    
    
    # GBIC_ppo to select the polynomial order 
    GBICscrores<-apply(as.matrix(1:gamma,1,gamma), 1, function(t) GBIC_ppo(Ynew,covariates_labelled,t))
    alpha<- which.min(GBICscrores)
    
  }
  
  # determine Z
  labelled_Z <- polynomial(labelled_data[,-1],alpha)
  unlabelled_Z <- polynomial(unlabelled_data,alpha)
  
  
  
  # weights (w_i) using in the our optimization problem to get the proposed semi-supervised estimate 
  unlabelled_Z_mean <- colMeans(unlabelled_Z)
  labelled_Z_seondmoment <- t(labelled_Z)%*%labelled_Z/n
  weights_loss <- as.vector(c1+(1-c1)*t(unlabelled_Z_mean)%*%solve(labelled_Z_seondmoment)%*%t(labelled_Z))
  
  
  
  # estimate of the target parameter
  if (type=="linear"){
    hattheta_PSSE <- as.vector(solve(t(cbind(rep(1,n),labelled_data[,-1]))%*%
                                       (weights_loss*cbind(rep(1,n),labelled_data[,-1])))%*%t(cbind(rep(1,n),labelled_data[,-1]))%*%
                                 (weights_loss*labelled_data[,1]))
  }
  
  if (type=="logistic"){
    hattheta_PSSE <- newton_iteration(labelled_data,tol=10^(-5),max_i=50,weights=weights_loss)
  }
  
  else if(type=="quantile"){
    
    beta_initial<- rq(as.vector(labelled_data[,1]) ~ labelled_data[,-1], tau=tau, weights=weights_loss*(weights_loss>=0), 
                      data = as.data.frame(labelled_data))$coefficients
    if(sum(weights_loss<0)>0)
    {
      tol=10^(-5)
      max_i=100
      rate=10^(-3)
      beta_initial<- gradient_descent(labelled_data,beta_initial=beta_initial,rate=rate,
                                      tol,max_i,weights=weights_loss,type="quantile",tau=tau)
      
    }
    hattheta_PSSE<- beta_initial
  }
  hattheta_PSSE = list("Hattheta"=hattheta_PSSE)
  
  ### estimating the sd
  if(sd==TRUE){
    
    if(type=="linear"){
      
      
      X_secondmoment_inverse=solve(t(cbind(rep(1,N),unlabelled_data))%*%cbind(rep(1,N),unlabelled_data)/N) 
      
      set.seed(20218080)
      index=createFolds(1:n, k = Kfolds) # data splitting
      W1_test = vector()
      
      
      for(k in 1:Kfolds){
        
        index_k=as.vector(index[[k]])
        Yt_train = labelled_data[-index_k,1]
        Xt_train = labelled_data[-index_k,-1]
        Zt_train = labelled_Z[-index_k,]
        Yt_test = labelled_data[index_k,1]
        Xt_test = labelled_data[index_k,-1]
        Zt_test = labelled_Z[index_k,]
        
        nrow_Xt_train = n-length(index_k)
        nrow_Xt_test = length(index_k)
        
        
        
        # projection matrix A(theta) estimation 
        
        hattheta_supervised_train_cv <- hattheta_PSSE[[1]]#lm(Yt_train~Xt_train)$coefficients
        L_firstder_train_cv <- as.vector(Yt_train-cbind(rep(1,nrow_Xt_train),Xt_train)%*%hattheta_supervised_train_cv)*
          cbind(rep(1,nrow_Xt_train),Xt_train)
        L_firstder_projection_cof <- solve(t(Zt_train)%*%Zt_train/nrow_Xt_train)%*%
          t(Zt_train)%*%L_firstder_train_cv/nrow_Xt_train
        
        
        # estimate variance by test data
        
        L_firstder_test_cv <- as.vector(Yt_test-cbind(rep(1,nrow_Xt_test),Xt_test)%*%hattheta_supervised_train_cv)*
          cbind(rep(1,nrow_Xt_test),Xt_test)
        
        W1_test_k <- L_firstder_test_cv+(c1-1)*Zt_test%*%L_firstder_projection_cof
        W1_test<- rbind(W1_test,W1_test_k)
        
        
        
      }
      
      
      
      L_firstder_total <- as.vector(labelled_data[,1]-(cbind(rep(1,n),labelled_data[,-1]))%*%hattheta_PSSE[[1]])*
        (cbind(rep(1,n),labelled_data[,-1]))
      L_firstder_projection_cof_total <- solve(t(labelled_Z)%*%labelled_Z/n)%*%t(labelled_Z)%*%L_firstder_total/n
      W2_total <- (1-c1)*unlabelled_Z%*%L_firstder_projection_cof_total
      
      W1_covariance <- t(W1_test)%*%W1_test/n
      W2_covariance <- t(W2_total)%*%W2_total/N
      Vc_hat_semi = W1_covariance+(n/N)*W2_covariance
      
      Var_matrix_hat = X_secondmoment_inverse%*%Vc_hat_semi%*%X_secondmoment_inverse 
      sd_hattheta_PSSE = sqrt(diag(Var_matrix_hat/n))
    }
    
    if(type=="logistic"){
      
      exp_linear_combined <- as.vector(cbind(rep(1,N),unlabelled_data)%*%hattheta_PSSE[[1]])
      weights_second_derivative <- 1/(1+exp(-exp_linear_combined))^2*exp(-exp_linear_combined)
      X_secondmoment_inverse=solve(t(weights_second_derivative*cbind(rep(1,N),unlabelled_data))%*%cbind(rep(1,N),unlabelled_data)/N) 
      
      set.seed(20218080)
      index=createFolds(1:n, k = Kfolds) # data splitting
      W1_test = vector()
      
      
      for(k in 1:Kfolds){
        
        index_k=as.vector(index[[k]])
        Yt_train = labelled_data[-index_k,1]
        Xt_train = labelled_data[-index_k,-1]
        Zt_train = labelled_Z[-index_k,]
        Yt_test = labelled_data[index_k,1]
        Xt_test = labelled_data[index_k,-1]
        Zt_test = labelled_Z[index_k,]
        
        nrow_Xt_train = n-length(index_k)
        nrow_Xt_test = length(index_k)
        
        
        
        # projection matrix A(theta) estimation 
        
        hattheta_supervised_train_cv <- hattheta_PSSE[[1]]#lm(Yt_train~Xt_train)$coefficients
        exp_linear_combined_train <- as.vector(cbind(rep(1,nrow_Xt_train),Xt_train)%*%hattheta_supervised_train_cv)
        L_firstder_train_cv <- as.vector(1/(1+exp(-exp_linear_combined_train))-Yt_train)*
          cbind(rep(1,nrow_Xt_train),Xt_train)
        L_firstder_projection_cof <- solve(t(Zt_train)%*%Zt_train/nrow_Xt_train)%*%
          t(Zt_train)%*%L_firstder_train_cv/nrow_Xt_train
        
        
        # estimate variance by test data
        
        exp_linear_combined_test <- as.vector(cbind(rep(1,nrow_Xt_test),Xt_test)%*%hattheta_supervised_train_cv)
        L_firstder_test_cv <- as.vector(1/(1+exp(-exp_linear_combined_test))-Yt_test)*
          cbind(rep(1,nrow_Xt_test),Xt_test)
        
        W1_test_k <- L_firstder_test_cv+(c1-1)*Zt_test%*%L_firstder_projection_cof
        W1_test<- rbind(W1_test,W1_test_k)
        
      }
      
      
      exp_linear_combined_total <- as.vector(cbind(rep(1,n),labelled_data[,-1])%*%hattheta_PSSE[[1]])
      L_firstder_total <- (1/(1+exp(-exp_linear_combined_total))-labelled_data[,1])*
        (cbind(rep(1,n),labelled_data[,-1]))
      L_firstder_projection_cof_total <- solve(t(labelled_Z)%*%labelled_Z/n)%*%t(labelled_Z)%*%L_firstder_total/n
      W2_total <- (1-c1)*unlabelled_Z%*%L_firstder_projection_cof_total
      
      W1_covariance <- t(W1_test)%*%W1_test/n
      W2_covariance <- t(W2_total)%*%W2_total/N
      Vc_hat_semi = W1_covariance+(n/N)*W2_covariance
      
      Var_matrix_hat = X_secondmoment_inverse%*%Vc_hat_semi%*%X_secondmoment_inverse 
      sd_hattheta_PSSE = sqrt(diag(Var_matrix_hat/n))
      
    }
    
    if(type=="quantile"){
      
      set.seed(20218080)
      index=createFolds(1:n, k = Kfolds) # data splitting
      W1_test = vector()
      
      
      for(k in 1:Kfolds){
        
        index_k=as.vector(index[[k]])
        Yt_train = labelled_data[-index_k,1]
        Xt_train = labelled_data[-index_k,-1]
        Zt_train = labelled_Z[-index_k,]
        Yt_test = labelled_data[index_k,1]
        Xt_test = labelled_data[index_k,-1]
        Zt_test = labelled_Z[index_k,]
        
        nrow_Xt_train = n-length(index_k)
        nrow_Xt_test = length(index_k)
        
        
        residual_train_indicator <- as.vector(Yt_train-cbind(rep(1,nrow_Xt_train),Xt_train)%*%hattheta_PSSE[[1]])<=0
        L_firstder_train_cv <- (residual_train_indicator-tau)*cbind(rep(1,nrow_Xt_train),Xt_train)
        L_firstder_projection_cof <- solve(t(Zt_train)%*%Zt_train)%*%
          t(Zt_train)%*%L_firstder_train_cv
        
        
        # estimate variance by test data
        
        residual_test_indicator <- as.vector(Yt_test-cbind(rep(1,nrow_Xt_test),Xt_test)%*%hattheta_PSSE[[1]])<=0
        L_firstder_test_cv <- (residual_test_indicator-tau)*cbind(rep(1,nrow_Xt_test),Xt_test)
        
        W1_test_k <- L_firstder_test_cv+(c1-1)*Zt_test%*%L_firstder_projection_cof
        W1_test<- rbind(W1_test,W1_test_k)
        
      }
      
      
      residual_total_indicator <- as.vector(labelled_data[,1]-cbind(rep(1,n),labelled_data[,-1])%*%hattheta_PSSE[[1]])<=0
      L_firstder_total <- (residual_total_indicator-tau)*(cbind(rep(1,n),labelled_data[,-1]))
      L_firstder_projection_cof_total <- solve(t(labelled_Z)%*%labelled_Z)%*%t(labelled_Z)%*%L_firstder_total
      W2_total <- (1-c1)*unlabelled_Z%*%L_firstder_projection_cof_total
      
      W1_covariance <- t(W1_test)%*%W1_test/n
      W2_covariance <- t(W2_total)%*%W2_total/N
      Vc_hat_semi = W1_covariance+(n/N)*W2_covariance
      
      
      
      ##estimating the second derivatives (M)
      B=2000
      G=mvrnorm(B,rep(0,(p+1)),diag(rep(1,(p+1))))
      theta_check_semi <- apply(1/sqrt(n)*G,1,function(t) hattheta_PSSE[[1]]+t)
      residual_semi_indicator= apply(cbind(rep(1,n),labelled_data[,-1])%*%theta_check_semi,2, function(t) as.vector(labelled_data[,1]-t)<=0)
      U_check_semi <- t(residual_semi_indicator-tau)%*%cbind(rep(1,n),labelled_data[,-1])/sqrt(n)
      
      hat_M_semi <- t(apply(U_check_semi,2, function(t) lm(t~G-1)$coefficients))
      
      
      if(class(try(solve(hat_M_semi),silent=T))[1]=="try-error"){
        X_secondmoment_inverse <- solve((hat_M_semi+t(hat_M_semi))/2)
      } else{
        X_secondmoment_inverse<- solve(hat_M_semi)
      }
      
      
      
      
      Var_matrix_hat = X_secondmoment_inverse%*%Vc_hat_semi%*%t(X_secondmoment_inverse)
      sd_hattheta_PSSE = sqrt(diag(Var_matrix_hat)/n)
    }
    
    
    hattheta_PSSE = append(hattheta_PSSE,list("sd.of.hattheta"=sd_hattheta_PSSE))
  }
  
  hattheta_PSSE = append(hattheta_PSSE,list("alpha"=alpha))
  
  return(hattheta_PSSE)
}





#####################DRESS###########################
DRESS<- function(labelled_data,unlabelled_data,L,Kfolds)
{
  #-------------------------Arguments-------------------------------------------#
  # Purpose: This function is to implement the method proposed by
  #         Kawakita and Kanamori (2013) for M-estimation,
  #         which uses the desity-ratio to improve the estimation
  #         efficiency.
  #
  # Input:
  #       labelled_data: Same as in the function "PSSE". 
  #       unlabelled_data: Same as in the function "PSSE".     
  #       type: Same as in the function "PSSE".
  #       tau: Same as in the function "PSSE".
  #        L: The polynomial order when we construct the polynomial function as
  #            base function for density ratio estimation.
  #       sd: Return the pointwise standard deviation estimate based on the available
  #           samples? Default value is FALSE.
  #       Kfolds: Same as in the function "PSSE". 
  #
  # Output: 
  #       Hattheta: The estimate for the target parameter. 
  #       sd.of.hattheta: If "sd=TRUE", the elementwise standard deviation 
  #                        estimate for the estimator; otherwise, NULL. 
  #       error: indicator of whether the parameters in density-ratio is properly
  #              estimated. "TRUE" represents "NOT"; "FALSE" represents "YES".
  #-----------------------------------------------------------------------------#
  
  n=nrow(labelled_data)
  p=ncol(labelled_data)-1
  N=nrow(unlabelled_data)
  
  
  base_labelled<- polynomial(labelled_data[,-1],L)
  base_unlabelled <- polynomial(unlabelled_data,L)
  alpha_first_derivative_unlabelled <- colMeans(base_unlabelled)
  
  ## estimating alpha (the parameter in density-ratio)
  tol=10^(-5)
  max_i = 50
  
  alpha_initial <- rep(0,L*p+1)
  alpha_distance<-10
  i=1
  error_svd=FALSE
  while(alpha_distance >tol& i<=max_i){
    exponential_phi_labelled <- as.vector(exp(base_labelled%*%alpha_initial))
    exponential_phi_labelled[which(exponential_phi_labelled==Inf)]=10^(-7)
    
    
    alpha_first_derivative <- colMeans(exponential_phi_labelled*base_labelled)-
      alpha_first_derivative_unlabelled
    
    alpha_second_derivative <- t(exponential_phi_labelled*base_labelled)%*%base_labelled/n
    
    
    if (max(abs(svd(alpha_second_derivative)$d))>99999|min(abs(svd(alpha_second_derivative)$d))<tol)
    {
      error_svd = TRUE
      alpha_initial<- rep(0,L*p+1)
      break
    }
    
    alpha_new <- as.vector(alpha_initial-solve(alpha_second_derivative)%*%alpha_first_derivative)
    
    
    alpha_distance <- sqrt(sum((alpha_new-alpha_initial)^2))
    alpha_initial<- alpha_new
    i=i+1
    #print(i)
  }
  
  
  exponential_phi_labelled <- as.vector(exp(base_labelled%*%alpha_initial))
  
  
  ## estimating theta*(the target parameter)

    hattheta_DRESS<- lm(labelled_data[,1]~labelled_data[,-1],weights=exponential_phi_labelled)$coefficients
    
  
  
  hattheta_DRESS = append(list("Hattheta"=hattheta_DRESS),list("error"=error_svd))
  
  
  ### estimating the sd

    
   
      
      
      X_secondmoment_inverse=solve(t(cbind(rep(1,N),unlabelled_data))%*%cbind(rep(1,N),unlabelled_data)/N) 
      
      set.seed(20218080)
      index=createFolds(1:n, k = Kfolds) # data splitting
      W1_test = vector()
      
      
      for(k in 1:Kfolds){
        
        index_k=as.vector(index[[k]])
        Yt_train = labelled_data[-index_k,1]
        Xt_train = labelled_data[-index_k,-1]
        Zt_train = base_labelled[-index_k,]
        Yt_test = labelled_data[index_k,1]
        Xt_test = labelled_data[index_k,-1]
        Zt_test = base_labelled[index_k,]
        
        nrow_Xt_train = n-length(index_k)
        nrow_Xt_test = length(index_k)
        
        
        
        # projection matrix A(theta) estimation 
        
        hattheta_supervised_train_cv <-  hattheta_DRESS[[1]]#lm(Yt_train~Xt_train)$coefficients
        L_firstder_train_cv <- as.vector(Yt_train-cbind(rep(1,nrow_Xt_train),Xt_train)%*%hattheta_supervised_train_cv)*
          cbind(rep(1,nrow_Xt_train),Xt_train)
        L_firstder_projection_cof <- solve(t(Zt_train)%*%Zt_train/nrow_Xt_train)%*%
          t(Zt_train)%*%L_firstder_train_cv/nrow_Xt_train
        
        
        # estimate variance by test data
        
        L_firstder_test_cv <- as.vector(Yt_test-cbind(rep(1,nrow_Xt_test),Xt_test)%*%hattheta_supervised_train_cv)*
          cbind(rep(1,nrow_Xt_test),Xt_test)
        
        W1_test_k <- L_firstder_test_cv-Zt_test%*%L_firstder_projection_cof
        W1_test<- rbind(W1_test,W1_test_k)
        
        
        
      }
      
      
      
      L_firstder_total <- as.vector(labelled_data[,1]-(cbind(rep(1,n),labelled_data[,-1]))%*%hattheta_DRESS[[1]])*
        (cbind(rep(1,n),labelled_data[,-1]))
      L_firstder_projection_cof_total <- solve(t(base_labelled)%*%base_labelled/n)%*%t(base_labelled)%*%L_firstder_total/n
      W2_total <- base_unlabelled%*%L_firstder_projection_cof_total
      
      W1_covariance <- t(W1_test)%*%W1_test/n
      W2_covariance <- t(W2_total)%*%W2_total/N
      Vc_hat_semi = W1_covariance+(n/N)*W2_covariance
      
      Var_matrix_hat = X_secondmoment_inverse%*%Vc_hat_semi%*%X_secondmoment_inverse 
      sd_hattheta_DRESS = sqrt(diag(Var_matrix_hat/n))
    
    
   
    
    hattheta_DRESS <- append(hattheta_DRESS,list("sd.of.hattheta"=sd_hattheta_DRESS))
  
  
  
  
  return(hattheta_DRESS)
}



###########################(2) PI proposed by Azriel et al. (2021) #################################
PI <- function(labelled_data,unlabelled_data)
{
  #---------------------------------Arguments------------------------------------------#
  # Purpose: This function is to implement the method proposed by Azriel et al.(2021)
  #          for linear working model. 
  #
  # Input:
  #       labelled_data: Same as in the function "PSSE". 
  #       unlabelled_data: Same as in the function "PSSE".     
  #
  # Output: 
  #        Hattheta: The estimate for the target parameter. 
  #------------------------------------------------------------------------------------#
  n=nrow(labelled_data)
  p=ncol(labelled_data)-1
  N=nrow(unlabelled_data)
  
  # combine all the covairates 
  X_combined <- rbind(labelled_data[,-1],unlabelled_data)
  X_labelled <- labelled_data[,-1]
  
  hat_beta_initial <- numeric()
  X_dot <- matrix(rep(0,n*p),n,p) 
  delta_tilde <- matrix(rep(0,n*p),n,p)
  
  for (j in 1:p)
  {
    ## first step 
    coefficients_negtive_j = lm(X_combined[,j]~X_combined[,-j])$coefficients #unlabeled data
    X_j_dot = X_labelled[,j] - cbind(rep(1,n),X_labelled[,-j]) %*% coefficients_negtive_j #labeled data projection error
    X_j_dot_total = X_combined[,j] - cbind(rep(1,n+N),X_combined[,-j]) %*% coefficients_negtive_j
    X_j_dot_square_sampleaverage=mean(X_j_dot_total^2)
    
    
    ## step step 
    W_j = as.vector(labelled_data[,1])*X_j_dot/X_j_dot_square_sampleaverage
    U1 = X_j_dot/X_j_dot_square_sampleaverage
    X_dot[,j] = as.vector(U1) 
    U = matrix(rep(0,n*p),n,p) 
    for (jj in 1:p)
    {
      if (jj==j)
      {
        U[,jj] = X_labelled[,jj]*X_j_dot/X_j_dot_square_sampleaverage-1
      }
      
      else
      {
        U[,jj] = X_labelled[,jj]*X_j_dot/X_j_dot_square_sampleaverage
      }
      
    }
    
    delta_tilde[,j] = W_j-as.matrix(cbind(rep(1,n),U1,U)) %*% (lm(W_j~U1+U)$coefficients)
    hat_beta_j = lm(W_j~U1+U)$coefficients[1]
    hat_beta_initial = c(hat_beta_initial,hat_beta_j)
  }
  hat_alpha = mean(labelled_data[,1])-t(hat_beta_initial)%*%colMeans(X_labelled)
  hat_theta_Azriel = as.vector(c(hat_alpha,hat_beta_initial))
  
  hat_theta_Azriel = list("Hattheta"=hat_theta_Azriel)
  
  
  
  return(hat_theta_Azriel)
}
#PI(data_labelled,data_unlabelled)








#############################################################################################################################
##############################################ET functions###############################################################
variance_theta_diag <- function(y, X, w_hat, theta_hat, ci_level = 0.95) {
  N <- length(y)
  p <- length(theta_hat)
  
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
    D[i, ] <- w_hat[i] * resid[i] * x_i
  }
  d_bar <- colMeans(D)
  
  # τ_hat
  tau_hat <- matrix(0, p, p)
  for (i in 1:N) {
    x_i <- as.numeric(X[i, ])
    tau_hat <- tau_hat + w_hat[i] * (x_i %o% x_i)
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


variance_theta_diag_scaled <- function(y, X, w_hat, theta_hat, ci_level = 0.95,
                                       scale_continuous = TRUE) {
  N <- length(y)
  p <- length(theta_hat)
  X <- as.matrix(X)
  
  # Ensure column names
  if (is.null(colnames(X))) colnames(X) <- paste0("X", 1:p)
  
  # Detect binary (0/1) columns (don't scale these)
  is_binary <- apply(X, 2, function(col) {
    vals <- unique(col)
    all(vals %in% c(0,1)) && length(vals) <= 2
  })
  
  # Choose columns to scale = non-binary (continuous) columns
  scale_idx <- if (scale_continuous) which(!is_binary) else integer(0)
  
  # SDs for columns we will scale (avoid divide-by-zero)
  sd_vec <- rep(1, p)
  if (length(scale_idx) > 0) {
    sd_tmp <- apply(X[, scale_idx, drop = FALSE], 2, sd, na.rm = TRUE)
    sd_tmp[!is.finite(sd_tmp) | sd_tmp == 0] <- 1
    sd_vec[scale_idx] <- sd_tmp
  }
  
  # Scale design WITHOUT centering: x* = x / sd
  Xs <- X
  if (length(scale_idx) > 0) {
    Xs[, scale_idx] <- sweep(X[, scale_idx, drop = FALSE], 2, sd_vec[scale_idx], "/")
  }
  
  # Transform coefficients to the scaled parameterization:
  # beta*_j = beta_j * sd_j  (intercept and dummies unchanged since sd=1)
  thetas <- theta_hat * sd_vec
  
  # Residuals computed on the scaled system give the same y - X beta
  resid <- as.numeric(y - Xs %*% thetas)
  
  # d_i rows (vectorized)
  D <- resid * w_hat * Xs
  d_bar <- colMeans(D)
  
  # τ_hat = (1/N) Σ w_i x_i x_i^T
  Xw <- sqrt(w_hat) * Xs
  tau_hat <- crossprod(Xw) / N
  
  # Middle term: (1/[N(N-1)]) Σ (d_i - d̄)(d_i - d̄)^T
  Dc <- sweep(D, 2, d_bar, "-")
  middle <- crossprod(Dc) / (N * (N - 1))
  
  # Sandwich variance on the scaled system
  tau_inv <- solve(tau_hat)
  var_theta_scaled <- tau_inv %*% middle %*% t(tau_inv)
  
  # Backscale variances to original units:
  # Var(beta_j) = Var(beta*_j) / sd_j^2
  variances <- diag(var_theta_scaled) / (sd_vec^2)
  names(variances) <- colnames(X)
  
  # SE and CI in original units
  se <- sqrt(variances)
  alpha <- 1 - ci_level
  z_val <- qnorm(1 - alpha / 2)
  lower <- theta_hat - z_val * se
  upper <- theta_hat + z_val * se
  
  # Result table in original units
  result <- data.frame(
    Variable = colnames(X),
    Estimate = theta_hat,
    Variance = variances,
    SE = se,
    CI_Lower = lower,
    CI_Upper = upper,
    row.names = NULL
  )
  
  return(result)
}


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
k_fold_function<-function(df,K,seed){
  res <- build_SU_folds(df, K , id_col = "ID", delta_col = "D", seed = seed)
  
  lab_idx_all <- which(df$D == 1)
  
  data_unlabeled<-list()
  for(k in 1:K){
    train_idx_lab <- setdiff(lab_idx_all, res$folds[[k]]$S_k_idx)
    validation_data_labeled <- df[res$folds[[k]]$S_k_idx,]
    validation_data_unlabeled<-df[res$folds[[k]]$U_k_idx,]
    
    train_data_labeled <-df[train_idx_lab,] ###get the data frame for trained data
    
    
    model.formula <- as.formula(paste("Y ~", paste(sprintf("`%s`", colnames( subset(train_data_labeled, select = -c(D, Y, ID,pi.hat)))), collapse = " + ")))
    
    lm_model<-lm(model.formula,data=as.data.frame( subset(train_data_labeled, select = -c(D, ID,pi.hat))))
    # gam_model <- gam(Y ~ s(age) + s(BMI) + s(SBP) + s(DBP) + sex + race,
    #             data=as.data.frame( subset(train_data_labeled, select = -c(D, ID,pi.hat))), method = "REML")
    
    y_hat<-predict( lm_model, newdata=as.data.frame(validation_data_unlabeled))
    #y_hat<-predict( gam_model, newdata=as.data.frame(validation_data_unlabeled))
    #rf_model <- randomForest(model.formula, data = train_data_labeled)
    #y_hat <- predict(rf_model, newdata = validation_data_unlabeled)
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

estimate_theta_EM_kfold_CVXR_ET <- function(th, data_full, K, seed, max.iter,eps) {
  
  # --- prepare cross-fitted data, stacked in fold order ---
  fold_hat <- k_fold_function(df = data_full, K = K, seed = seed)
  data_all <- do.call(rbind, fold_hat)
  n_k      <- sapply(fold_hat, nrow)
  ofs      <- c(0, cumsum(n_k))
  
  # features
  # define p properly:
  x_cols <- subset(data_all,select=-c(Y,D,pi.hat,ID,y.hat))
  p      <- length(x_cols)
  model.formula<-as.formula(paste("~",paste(colnames(x_cols),collapse="+")))
  new.data  <- as.data.frame(data_all[, colnames(x_cols), drop=FALSE])
  new.data1 <- (model.matrix(  model.formula,
                               data = as.data.frame(new.data)))
  
  
  y.hat <- data_all$y.hat
  D     <- data_all$D
  I1    <- which(D == 1)
  
  theta <- as.matrix(th)
  iter  <- 0
  
  repeat {
    iter <- iter + 1
    
    
    
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
    sol  <- CVXR::solve(prob,solver = "ECOS")
    
    
    if (sol$status != "optimal") {
      return(matrix(NA_real_, nrow=1, ncol=ncol(new.data1)))
    }
    
    W  <- as.numeric(sol$getValue(w))
    Q  <- new.data1[I1, , drop=FALSE]
    Y1 <- data_all$Y[I1]
    A<-t(Q) %*% diag(W) %*% Q
    B<-t(Q) %*% diag(W) %*% Y1
    th.new <- solve(A, B)
    #print(th.new)
    #print(max(abs(th.new - as.numeric(theta))))
    if (max(abs(th.new - as.numeric(theta))) < eps) {
      
      
      return(list(theta=matrix(th.new, ncol=1),w=W,estimate.table=variance_theta_diag(
        y = y,
        X = new.data1[I1, ],
        w_hat = W,
        theta_hat = th.new
      )))
    }
    if (iter >= max.iter) {
     # message("Maximum iterations reached.")
      return(list(theta=matrix(th.new, ncol=1),w=W,estimate.table=variance_theta_diag(
        y = y,
        X = new.data1[I1, ],
        w_hat = W,
        theta_hat = th.new
      )))
    }
    theta <- as.matrix(th.new)
  }
}



estimate_theta_EM_kfold_CVXR_HD <- function(th, data_full, K, seed, max.iter,eps) {
  
  # --- prepare cross-fitted data, stacked in fold order ---
  fold_hat <- k_fold_function(df = data_full, K = K, seed = seed)
  data_all <- do.call(rbind, fold_hat)
  n_k      <- sapply(fold_hat, nrow)
  ofs      <- c(0, cumsum(n_k))
  
  # features
  # define p properly:
  x_cols <- subset(data_all,select=-c(Y,D,pi.hat,ID,y.hat))
  p      <- length(x_cols)
  model.formula<-as.formula(paste("~",paste(colnames(x_cols),collapse="+")))
  new.data  <- as.data.frame(data_all[, colnames(x_cols), drop=FALSE])
  new.data1 <- (model.matrix(  model.formula,
                               data = as.data.frame(new.data)))
  
  
  y.hat <- data_all$y.hat
  D     <- data_all$D
  I1    <- which(D == 1)
  
  theta <- as.matrix(th)
  iter  <- 0
  
  repeat {
    iter <- iter + 1
    
    
    
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
    g<--sqrt(data_all$pi.hat)/2
    
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
    
    objective <- CVXR::Minimize(-sum(sqrt(w)))
    
    prob <- CVXR::Problem(objective, constr)
    sol  <- CVXR::solve(prob,solver = "ECOS")
    
    
    if (sol$status != "optimal") {
      return(matrix(NA_real_, nrow=1, ncol=ncol(new.data1)))
    }
    
    W  <- as.numeric(sol$getValue(w))
    Q  <- new.data1[I1, , drop=FALSE]
    Y1 <- data_all$Y[I1]
    A<-t(Q) %*% diag(W) %*% Q
    B<-t(Q) %*% diag(W) %*% Y1
    th.new <- solve(A, B)
    #print(th.new)
    #print(max(abs(th.new - as.numeric(theta))))
    if (max(abs(th.new - as.numeric(theta))) < eps) {
      
      
      return(list(theta=matrix(th.new, ncol=1),w=W,estimate.table=variance_theta_diag(
        y = y,
        X = new.data1[I1, ],
        w_hat = W,
        theta_hat = th.new
      )))
    }
    if (iter >= max.iter) {
      #message("Maximum iterations reached.")
      return(list(theta=matrix(th.new, ncol=1),w=W,estimate.table=variance_theta_diag(
        y = y,
        X = new.data1[I1, ],
        w_hat = W,
        theta_hat = th.new
      )))
    }
    theta <- as.matrix(th.new)
  }
}




estimate_theta_nested_kfold_CVXR_optim <- function(theta_init, data_full, K , seed ) {
  # ------------------------------
  # Define the objective function L(theta)
  # ------------------------------
  # --- prepare cross-fitted data, stacked in fold order ---
  fold_hat <- k_fold_function(df = data_full, K = K, seed = seed)
  data_all <- do.call(rbind, fold_hat)
  n_k      <- sapply(fold_hat, nrow)
  ofs      <- c(0, cumsum(n_k))
  
  # features
  # define p properly:
  x_cols <- subset(data_all,select=-c(Y,D,pi.hat,ID,y.hat))
  p      <- length(x_cols)
  model.formula<-as.formula(paste("~",paste(colnames(x_cols),collapse="+")))
  new.data  <- as.data.frame(data_all[, colnames(x_cols), drop=FALSE])
  new.data1 <- (model.matrix(  model.formula,
                               data = as.data.frame(new.data)))
  
  
  y.hat <- data_all$y.hat
  D     <- data_all$D
  I1    <- which(D == 1)
  
  g <- log(data_all$pi.hat)
  y <- data_all[I1, "Y"]
  
  L_fn <- function(theta) {
    
    
    # residual-based moments
    error <- y.hat - as.numeric(new.data1 %*% theta)
    H <- new.data1 * error
    A <- H
    r <- ncol(A)
    
    # --- CVXR setup ---
    w <- CVXR::Variable(length(I1), pos = TRUE)
    
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
    sol  <- CVXR::solve(prob,solver = "ECOS")
    
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
    control = list(reltol = 1e-3, maxit = 300)
  )
  
  theta_hat = optim_result$par
  
  # residual-based moments
  error <- y.hat - as.numeric(new.data1 %*%   theta_hat)
  H <- new.data1 * error
  A <- H
  r <- ncol(A)
  
  # --- CVXR setup ---
  w <- CVXR::Variable(length(I1), pos = TRUE)
  
  resid <- y - as.matrix(new.data1[I1, ]) %*%  theta_hat
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
  sol  <- CVXR::solve(prob,solver = "ECOS")
  
  if (sol$status != "optimal") {
    return(NA_real_)
  }
  
  
  
  W <- as.numeric(sol$getValue(w))
  
  # Return estimated theta and value
  
  list(
    theta_hat = optim_result$par,
    value = optim_result$value,
    convergence = optim_result$convergence,
    message = optim_result$message,
    estimate.table=variance_theta_diag(
      y = y,
      X = new.data1[I1, ],
      w_hat = W,
      theta_hat = optim_result$par
    )
  )
}











estimate_theta_nested_kfold_CVXR_nlm <- function(theta_init, data_full, K, seed) {
  
  # --- 1. Prepare data (unchanged) ---
  fold_hat <- k_fold_function(df = data_full, K = K, seed = seed)
  data_all <- do.call(rbind, fold_hat)
  n_k      <- sapply(fold_hat, nrow)
  ofs      <- c(0, cumsum(n_k))
  
  x_cols <- subset(data_all, select = -c(Y, D, pi.hat, ID, y.hat))
  p      <- length(x_cols)
  model.formula <- as.formula(paste("~", paste(colnames(x_cols), collapse = "+")))
  new.data  <- as.data.frame(data_all[, colnames(x_cols), drop = FALSE])
  new.data1 <- model.matrix(model.formula, data = new.data)
  
  y.hat <- data_all$y.hat
  D     <- data_all$D
  I1    <- which(D == 1)
  
  g <- log(data_all$pi.hat)
  y <- data_all[I1, "Y"]
  
  # --- 2. Define objective L(theta) (unchanged) ---
  L_fn <- function(theta) {
    
    error <- y.hat - as.numeric(new.data1 %*% theta)
    H <- new.data1 * error
    A <- H
    r <- ncol(A)
    
    w <- CVXR::Variable(length(I1), pos = TRUE)
    
    resid <- y - as.matrix(new.data1[I1, ]) %*% theta
    U_mat1 <- as.matrix(new.data1[I1, ]) * as.vector(resid)
    
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
    
    lamda <- lapply(constr, sol$getDualValue)
    lamda <- lamda[-3]
    lambda_hat <- -unlist(lamda)
    
    W <- as.numeric(sol$getValue(w))
    val <- -(sum(W)) + as.vector(U1 %*% lambda_hat)
    return(val)
  }
  
  # --- 3. Outer optimization replaced with nlminb() ---
  nlm_result <- tryCatch({
    nlminb(
      start = theta_init,
      objective = L_fn,
      control = list(rel.tol = 1e-2, x.tol = 1e-2, eval.max = 300, iter.max = 300)
    )
  }, error = function(e) {
    message("[nlminb] error: ", e$message)
    list(par = theta_init, objective = NA, convergence = 1, message = e$message)
  })
  
  theta_hat <- nlm_result$par
  
  # --- 4. Final CVXR solve for w_hat using estimated theta_hat ---
  error <- y.hat - as.numeric(new.data1 %*% theta_hat)
  H <- new.data1 * error
  A <- H
  r <- ncol(A)
  
  w <- CVXR::Variable(length(I1), pos = TRUE)
  resid <- y - as.matrix(new.data1[I1, ]) %*% theta_hat
  U_mat1 <- as.matrix(new.data1[I1, ]) * as.vector(resid)
  
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
    message("[Final CVXR] solver did not converge, returning NA for w_hat.")
    W <- rep(NA_real_, length(I1))
  } else {
    W <- as.numeric(sol$getValue(w))
  }
  
  # --- 5. Return ---
  list(
    theta_hat = theta_hat,
    value = nlm_result$objective,
    convergence = nlm_result$convergence,
    message = nlm_result$message,
    w_hat = W,
    estimate.table = variance_theta_diag(
      y = y,
      X = new.data1[I1, ],
      w_hat = W,
      theta_hat = theta_hat
    )
  )
}
