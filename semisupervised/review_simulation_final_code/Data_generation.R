GenerateData<- function(n,N,p,OR, MAR){
  sigma= 2  # the sd of the error term eta   
  
  alpha0=1
  alpha1=rep(1,p)
  alpha2=rep(1,p)
  
  X=mvrnorm(n+N,rep(0,p),diag(rep(1,p)))
  
  if(OR==1){
    
    Y=alpha0+X%*%alpha1+rnorm(n+N,0,1) ####OR1
    
  }else{
    Y=alpha0+X%*%alpha1+(X^3-X^2+exp(X))%*%alpha2+rnorm(n+N,0,2) ###OR2 
  }
  
  
  if(MAR==1){

    px <- 1 / (1 + exp(1 + X[, 1] + 0.5 * X[, 2] - 0.5 * X[, 3] - 0.1 * X[, 4])) ###MAR mechanism
    D=rbinom(n+N,1,px)

  }else{
    D=rbinom(n+N,1,n/(n+N)) ###MCAR mechanism
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



