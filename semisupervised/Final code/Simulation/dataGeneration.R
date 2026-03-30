#---------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------#
#
#
#                                   Data generation functions 
#          (mainly for linear working model;logistic working model;quantile working model)
#
#---------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------#
GenerateData<- function(n,N,p,MAR){
sigma= 2  # the sd of the error term eta       
alpha0=1
alpha1=rep(1,p)
alpha2=rep(1,p)

X=mvrnorm(n+N,rep(0,p),diag(rep(1,p)))
#Y=alpha0+X%*%alpha1+cos(2*X)%*%alpha1+sin(X)%*%alpha1+(X^3-X^2+exp(X))%*%alpha2+rnorm(n+N,0,sigma) ###OM2
beta0=1
#Y=alpha0+X%*%alpha1+(X^3-X^2+exp(X))%*%alpha2+rnorm(n+N,0,sigma) ###OM2 final
#Y=beta0+X%*%alpha1+rnorm(n+N,0,sigma) ###OM1 final

Y=beta0+X%*%alpha1+rnorm(n+N,0,1) ###OM1 final


if(MAR==1){
  beta1<-rep(0.25,p+1)
#px=1/(1+exp(-(cbind(1,X)%*%beta1)))
  
  #px<-1/(1+exp(-(.25+X[,1]+.5*X[,2]-.5*X[,3]-0.1*X[,4])))##Final PS
  px <- 1 / (1 + exp(1 + X[, 1] + 0.5 * X[, 2] - 0.5 * X[, 3] - 0.1 * X[, 4]))
  D=rbinom(n+N,1,px)
  mean(D)
}else{
  D=rbinom(n+N,1,n/(n+N)) ###MCAR
}

data_full<-cbind(D,Y,X)


colnames(data_full)<-c("D","Y",paste0("X",1:p))
data_full<-as.data.frame(data_full)
data_full$Y[data_full$D == 0] <- NA
data_full$ID <- seq_len(nrow(data_full))
pi.hat1<-fitted(glm(D~.,data=subset(data_full, select = -c(Y, ID)),family=binomial()))

data_full$pi.hat<-  pi.hat1



data_labelled<-as.matrix(data_full[D==1,c("Y",paste0("X",1:p))])
data_unlabelled<-as.matrix(data_full[D==0,paste0("X",1:p)])

return(list(Data.labelled =data_labelled,Data.unlabelled=data_unlabelled,data_full=data_full))
}






GenerateData1<- function(n,N,p){
  #-----------------------------------------Arguments--------------------------------------------------------#
  # Purpose: This function is to generate the simulated data under different settings both in Section 4.1  
  #           of the paper and Section 3.2 of the supplementary material. 
  #
  # Input: 
  #      n: The labelled sample size.
  #      N: The unlabelled sample size. The default is 0. When "N=0", no unlabelled data are generated. 
  #      p: The dimension of the predictor vector. 
  #      option: The type of data settings. We offer nine choices: "i","ii","iii","W1","W2","W3",S1","S2","S3".
  #               These nine choices correspond to settings (i) - (iii) in Section 4.1 of the paper, 
  #               settings (W1) - (W2), settings (S1) - (S3) in Section 3.2 of the supplementary material 
  #               in sequence. 
  #
  # Output: 
  #       Data.labelled: A matrix, whose each row is an observation of predictor vector and the response
  #                      variable and its first column is the observation vector of the response variable.
  #       Data.unlabelled: A matrix, whose each row is an observation of predictor vector. 
  #----------------------------------------------------------------------------------------------------------#
  #print(paste("NOTE: The data are generated based on setting (",option,").",sep=""))
 
    ###parameters in setting (i)
    sigma= 2  # the sd of the error term eta       
    alpha0=1
    alpha1=rep(1,p)
    alpha2=rep(1,p)
    ### data generation including the labelled data set and the unlabelled data set 
    #X_labelled<-mvrnorm(n,rep(0,p),diag(rep(1,p)))
    
    X=mvrnorm(n+N,rep(0,p),diag(rep(1,p)))
    #Y=alpha0+X%*%alpha1+(X^3-X^2+exp(X))%*%alpha2+rnorm(n+N,0,sigma)
    #px=1/(1+exp(-(-1+X%*%beta1))) ###MAR
    
    if(MAR==1){
      #beta1<-rep(0.25,p)
      #px=1/(1+exp(-(-.5+X%*%beta1)))
      beta1<-rep(.5,p)
      #px=1/(1+exp(-(X%*%beta1)))
      #px<-1/(1+exp(-(.5+X[,1]+.5*X[,2]-.5*X[,3]-0.1*X[,4])))
      px <- 1 / (1 + exp(-( -0.25+0.8*X[,1] + 0.3*X[,2] - 0.3*X[,3] - 0.05*X[,4])))
      D=rbinom(n+N,1,px)
      mean(D)
    }else{
      D=rbinom(n+N,1,n/(n+N)) ###MCAR
    }
    
    
    X_labelled = X[D==1,]
    eta=rnorm(sum(D),0,sigma)
    #eta=rnorm(n,0,sigma)
   
    #Y_labelled<-Y[D==1]
    #Y_labelled=alpha0+X_labelled%*%alpha1+rnorm(sum(D),0,1)
    #Y=alpha0+X%*%alpha1+(X^3-X^2+exp(X))%*%alpha2+rnorm(n+N,0,sigma) ###OM2 final
    #Y_labelled=Y[D==1] ###OM2 final
    
    Y_labelled=alpha0+X_labelled%*%alpha1+(X_labelled^3-X_labelled^2+exp(X_labelled))%*%alpha2+eta ###OM2 final
    
    #Y_labelled=alpha0+X_labelled%*%alpha1+cos(2*X_labelled)%*%alpha1+sin(X_labelled)%*%alpha1+(X_labelled^3-X_labelled^2+exp(X_labelled))%*%alpha2+eta ###OM2
    
    #Y_labelled=alpha0+X_labelled%*%alpha1+eta ###OM1
    
    
    
    #Y=alpha0+X%*%alpha1+(X^3-X^2+exp(X))%*%alpha2+rnorm(N+n,0,sigma) ###OM2
    #Y_labelled<-Y[D==1]
    
    
    data_labelled=cbind(Y_labelled,X_labelled)  # labelled data set 
    
    # give the colnames of labelled data set and unlabelled data set 
    dimU<- ncol(data_labelled)-1
    colnames(data_labelled)<- c("Y",paste("X",1:dimU,sep=""))
    res_list<- list("Data.labelled"=data_labelled)
    if(N!=0){
      
      X_unlabelled=X[D==0,]
      
      
      data_unlabelled=cbind(X_unlabelled)
  
    colnames(data_unlabelled)<- paste("X",1:dimU,sep="")
    #colnames(data_unlabelled)<- c("Y",paste("U",1:dimU,sep=""))
    res_list<- append(res_list,list("Data.unlabelled"=data_unlabelled))
    }
    return(res_list)
}
    
    
    

GenerateData_new<- function(n,N,p){
  #-----------------------------------------Arguments--------------------------------------------------------#
  # Purpose: This function is to generate the simulated data under different settings both in Section 4.1  
  #           of the paper and Section 3.2 of the supplementary material. 
  #
  # Input: 
  #      n: The labelled sample size.
  #      N: The unlabelled sample size. The default is 0. When "N=0", no unlabelled data are generated. 
  #      p: The dimension of the predictor vector. 
  #      option: The type of data settings. We offer nine choices: "i","ii","iii","W1","W2","W3",S1","S2","S3".
  #               These nine choices correspond to settings (i) - (iii) in Section 4.1 of the paper, 
  #               settings (W1) - (W2), settings (S1) - (S3) in Section 3.2 of the supplementary material 
  #               in sequence. 
  #
  # Output: 
  #       Data.labelled: A matrix, whose each row is an observation of predictor vector and the response
  #                      variable and its first column is the observation vector of the response variable.
  #       Data.unlabelled: A matrix, whose each row is an observation of predictor vector. 
  #----------------------------------------------------------------------------------------------------------#
  #print(paste("NOTE: The data are generated based on setting (",option,").",sep=""))
  
  ###parameters in setting (i)
  sigma= 2  # the sd of the error term eta       
  alpha0=1
  alpha1=rep(1,p)
  alpha2=rep(1,p)
 
  X=mvrnorm(n+N,rep(0,p),diag(rep(1,p)))
  
  Y=alpha0+X%*%alpha1+(X^3-X^2+exp(X))%*%alpha2+rnorm(n+N,0,sigma) ###OM2 final
 
}


    
###############The covariance matrix of the predictor vector ###################
X_AR1_covmatrix<- function(rho,p){
  #---------------------------Arguments----------------------------------------#
  # Purpose: This function is to produce the covariance matrix of the predictors.
  #
  # Input:
  #      rho: A parameter
  #      p: The dimension of the predictor vector. 
  #----------------------------------------------------------------------------#
  x_covmatrix<- matrix(rep(0,p^2),p,p)
  for (ii in 1:p)
  {
    for (jj in 1:p){
      x_covmatrix[ii,jj]<- rho^(abs(ii-jj))
    }
  }
  return(x_covmatrix)
}

#####################A function related to setting (S1) #########################
g_S1<- function(u,eta)
{
  exp(u+2*cos(u)-0.3*eta)+0.1*exp(u+sin(u))+2.5*u^3 # more complicated with eta 
  
}
#####################A function related to setting (S2) #########################
g_S2<-function(u){
  0.9*u^4*(sin(u))^6-2*u^2+2.2    
}
#####################A function related to setting (S3) #########################
g_S3<- function(u,eta)
{
  0.5*exp(u+2*cos(u)-0.3*eta)+0.6*exp(u+sin(u))+u^3
}

