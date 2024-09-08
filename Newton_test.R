      
      library(mgcv)
      library(caret)
      
      generate_data<-function(n,PS,MAR){
        
         x<-rnorm(n,0,1)
         z<-rexp(n,1)
        
        if(MAR==1){
          beta1<-1;beta2<-1;beta3<-2 
          y<-beta1+beta2*x+beta3*z+rnorm(n,0,1)
        }else{
          beta1<-1;beta2<-1;beta3<-2 
          y<-beta1+beta2*x+beta3*z+rnorm(n,0,1)
        }
        
        
        if(MAR==1){
          #eta1<--1;eta2<-.5;eta3<-0.5##PS1 with this choice 40% data is missing.final choice with exponential
          
          eta1<--1;eta2<-.25;eta3<-.25 ##PS1 with this choice 65% data is missing.final choice with exponential
          
        }else{
          eta1<--1;eta2<-0;eta3<-0 ##with this choice 70% data is missing. 
          
          }
        
        if(PS==1){
          
          px<-exp(eta1+eta2*y+eta3*x)/(1+exp(eta1+eta2*y+eta3*x))
          D<-rbinom(n,1,px)
          mean(D)
        } else {
          eta1<-0.5;eta2<-1;eta3<--.25 ###about 60% is missing
          #eta1<--.5;eta2<-.5;eta3<-1###20% IS MISSING FInal choice
         # eta1<--1.5;eta2<-.5;eta3<-.5
          px<-exp(eta1+eta2*x+eta3*(y-1)^2)/(1+exp(eta1+eta2*x+eta3*(y-1)^2))
          D<-rbinom(n,1,px)
          mean(D)
          }
        
        
        
        D.hat<-fitted(glm(D~y+x,family=binomial))
        
       
        dat<-data.frame(x=x,y=y,z=z,D=D,D.hat=D.hat)
        return(dat)
      }
      
      
    z.hat.function<-function(fold,q,new.dat0,new.dat1){
        
        pred0<-matrix(nrow = nrow(new.dat0), ncol =q)
        
        pred1<-matrix(nrow = nrow(new.dat1), ncol =1)
        
      
        for(k in 1:q){
          
          test_data1 <- new.dat1[fold[[k]],]
          train_data1 <- new.dat1[-fold[[k]],]
         
          model1<-gam(z~s(y,bs="tp")+s(x,bs="tp"),data=train_data1,method="REML")
          
          
          pred <- predict(model1, newdata=as.data.frame(test_data1))
          
          pred1[fold[[k]], ] <- as.numeric(pred)
          
          
          predU0<-as.numeric(predict(model1, newdata=as.data.frame(new.dat0)))
          pred0[,k]<-predU0
        }
        
        z.hat0<-as.numeric(apply(pred0,1,function(a)mean(a)))
        z.hat<-c(pred1,z.hat0)
        
        return(z.hat)
      }
      
      
  newton_EL_primal <- function(th,x,y,z,D.hat,z.hat,D, max.iter=50,eps=1e-6)
      {
        iter <- 0
        
        while (TRUE) {
          iter <- iter + 1
          
          theta<-as.matrix(th)
          b_theta1<-(y-theta[1]-x*theta[2]-theta[3]*z.hat)
          b_theta2<-x*(y-theta[1]-x*theta[2]-theta[3]*z.hat)
          b_theta3<-z.hat*(y-theta[1]-x*theta[2]-theta[3]*z.hat)
          
          
          
          H<-cbind(b_theta1,b_theta2,b_theta3,-(D.hat))
          
          
          
          
          ###Primal
          R=H[D==1,]
          p<-sum(D)
          
          w <- Variable(p,pos=TRUE)
          
          objective <- Minimize(-sum(log(w)))
          
          A<-cbind(rep(1,nrow(R)),R)
          b<-c(1,colMeans(H))
          
          constraint<-list(t(A)%*%w==b)
          problem <- Problem(objective, constraints = constraint)
          
          
          result3 <- solve(problem)
          if(result3$status=="optimal"){
            w2<-as.numeric(result3$getValue(w))
            
            new.dat<-cbind(data.frame(y=y,x=x,z=z)[D==1,],w2=w2)
            th.new<-as.matrix(as.numeric(coefficients(lm(y~x+z,data=new.dat,weights = w2))))
            
            print(max(abs(theta - th.new)))
            if (max(abs(theta - th.new)) < eps) {
              message("Convergence achieved.")
              return(th.new)
            }
            
            if (iter >= max.iter) {
              message("Maximum iterations reached. Returning current estimate.")
              return(as.matrix(rep(NA,length(th.new))))
            }
            print(iter)
            th <- th.new
            
          }else{
            message("Status is not optimal. Therefore returining NA values.")
            return(as.matrix(rep(NA,length(th.new))))
          }
        }
      }
      
      
      
      ###Hellinger distance (HD)
      
    newton_HD_primal <- function(th,x,y,z,D.hat,z.hat,D, max.iter=50,eps=1e-6)
      {
        iter <- 0
        
        while (TRUE) {
          iter <- iter + 1
          
          theta<-as.matrix(th)
          b_theta1<-(y-theta[1]-x*theta[2]-theta[3]*z.hat)
          b_theta2<-x*(y-theta[1]-x*theta[2]-theta[3]*z.hat)
          b_theta3<-z.hat*(y-theta[1]-x*theta[2]-theta[3]*z.hat)
          
          
          
          H<-cbind(b_theta1,b_theta2,b_theta3,2*sqrt(D.hat))
          
          
          
          
          ###Primal
          R=H[D==1,]
          p<-sum(D)
          
          w <- Variable(p,pos=TRUE)
          
          objective <- Minimize(-sum(4*sqrt(w)))
          
          A<-cbind(rep(1,nrow(R)),R)
          b<-c(1,colMeans(H))
          
          constraint<-list(t(A)%*%w==b)
          problem <- Problem(objective, constraints = constraint)
          
          
          result3 <- solve(problem)
          if(result3$status=="optimal"){
            w2<-as.numeric(result3$getValue(w))
            
            new.dat<-cbind(data.frame(y=y,x=x,z=z)[D==1,],w2=w2)
            th.new<-as.matrix(as.numeric(coefficients(lm(y~x+z,data=new.dat,weights = w2))))
            
            #print(max(abs(theta - th.new)))
            if (max(abs(theta - th.new)) < eps) {
              message("Convergence achieved.")
              return(th.new)
            }
            
            if (iter >= max.iter) {
              message("Maximum iterations reached. Returning current estimate.")
              return(as.matrix(rep(NA,length(th.new))))
            }
           # print(iter)
            th <- th.new
            
          }else{
            message("Status is not optimal. Therefore returining NA values.")
            return(as.matrix(rep(NA,length(th.new))))
          }
        }
      }
      
      
      
      
    newton_ET_primal <- function(th,x,y,z,D.hat,z.hat,D, max.iter=50,eps=1e-6)
      {
        iter <- 0
        
        while (TRUE) {
          iter <- iter + 1
          
          theta<-as.matrix(th)
          b_theta1<-(y-theta[1]-x*theta[2]-theta[3]*z.hat)
          b_theta2<-x*(y-theta[1]-x*theta[2]-theta[3]*z.hat)
          b_theta3<-z.hat*(y-theta[1]-x*theta[2]-theta[3]*z.hat)
          
          
          
          H<-cbind(b_theta1,b_theta2,b_theta3,-log(D.hat))
          
          
          
          
          ###Primal
          R=H[D==1,]
          p<-sum(D)
          
          w <- Variable(p,pos=TRUE)
          
          ones <- rep(1, p)
          
          # Define the objective function using kl_div
          objective <- Minimize(sum(kl_div(w, ones)-ones)) 
          
          
          
          A<-cbind(rep(1,nrow(R)),R)
          b<-c(1,colMeans(H))
          
          constraint<-list(t(A)%*%w==b)
          problem <- Problem(objective, constraints = constraint)
          
         
          result3 <- solve(problem)
          if(result3$status=="optimal"){
            w2<-as.numeric(result3$getValue(w))
           
            new.dat<-cbind(data.frame(y=y,x=x,z=z)[D==1,],w2=w2)
            th.new<-as.matrix(as.numeric(coefficients(lm(y~x+z,data=new.dat,weights = w2))))
            
            #print(max(abs(theta - th.new)))
            if (max(abs(theta - th.new)) < eps) {
              message("Convergence achieved.")
              return(th.new)
            }
            
            if (iter >= max.iter) {
              message("Maximum iterations reached. Returning current estimate.")
              return(as.matrix(rep(NA,length(th.new))))
            }
            #print(iter)
            th <- th.new
          
          }else{
            message("Status is not optimal. Therefore returining NA values.")
            return(as.matrix(rep(NA,length(th.new))))
          }
        }
      }
        
      
      