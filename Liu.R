#========Quantile-based estimation of Liu estimator in the Linear regression model========
#==========Authors: Muhammad Khan
#===========================Ph.D scholars in Statistics===============================
#set.seed(1234)
sigma1 = c(0.1,1,2,5)
ssize=c(25,50,100)
pred=c(4,10)
corr=c(0.85,0.95,0.99,0.999)
N=5000
#-----------------------------------------------------------------------------------------------------
MSE.r=matrix(0,16,length(corr))
MSE.s=array(0,dim = c(16,length(corr),length(sigma1)))
MSE.n=array(0,dim = c(16,length(corr),length(sigma1),length(ssize)))
MSE.p=array(0,dim = c(16,length(corr),length(sigma1),length(ssize),length(pred)))

MAE.r=matrix(0,16,length(corr))
MAE.s=array(0,dim = c(16,length(corr),length(sigma1)))
MAE.n=array(0,dim = c(16,length(corr),length(sigma1),length(ssize)))
MAE.p=array(0,dim = c(16,length(corr),length(sigma1),length(ssize),length(pred)))
#-----------------------------------------------------------------------------------------------------
# predictors loop starts here
for(b in 1:length(pred)){
  p=pred[b]
  # Sample size loop
  I=diag(p)
  for(a in 1:length(ssize)){
    n=ssize[a]
    # Generating standard normal random numbers
    z=matrix(0,n,p)
    for(i in 1:p){            #Generating stanndard normal random variable Z
      z[,i]=rnorm(n,0,1)
    }
    x=matrix(0,n,p)
    det.corr.x=rep(0,length(corr))
    # error variance loop
    for(s in 1:length(sigma1)){
      sigma=sigma1[s]
      # Correlation loop
      for(q in 1:length(corr)){
        r=corr[q]
        # Generating matrix of independent variables
        for(i in 1:p){            #Generating X matrix
          x[,i]=sqrt(1-(r^2))*z[,i]+r*z[,p]
        }
        x=scale(x,center = TRUE,scale = TRUE)   #Scaling x,,, standardizing x and y
        det.corr.x[q]=det(cor(x))
        eigen(cor(x))
        ev=eigen(cor(x))$values
        vec=eigen(cor(x))$vectors
        which.max(ev)
        beta=vec[,which.max(ev)]
        C=t(x)%*%(x)
        D=eigen(C)$vectors
        lamb=diag(t(D)%*%C%*%D)
        xstar=x%*%D
        alpha=t(D)%*%beta
        #-----------------------------------------------------------------------------------------------------------   
        betahat_ols=matrix(ncol=N, nrow=p)
        alphahatD1_M = matrix(ncol=N, nrow=p)
        alphahatD2_M = matrix(ncol=N, nrow=p)
        alphahatD3_M = matrix(ncol=N, nrow=p)
        alphahatD4_M = matrix(ncol=N, nrow=p)
        alphahatD5_M = matrix(ncol=N, nrow=p)
        alphahatD6_M = matrix(ncol=N, nrow=p)
        alphahatD7_M = matrix(ncol=N, nrow=p)
        alphahatD8_M = matrix(ncol=N, nrow=p)
        alphahatD9_M = matrix(ncol=N, nrow=p)
        alphahatD10_M = matrix(ncol=N, nrow=p)
        alphahatD11_M = matrix(ncol=N, nrow=p)
        alphahatD12_M = matrix(ncol=N, nrow=p)
        alphahatD13_M = matrix(ncol=N, nrow=p)
        alphahatD14_M = matrix(ncol=N, nrow=p)
        alphahatD15_M = matrix(ncol=N, nrow=p)
        
        #---------------------------------------------------------------------------------------------------------------  
        #Simulation loop starts here
        for(i in 1:N){
          e=rnorm(n,0,sigma)
          y=xstar%*%alpha+e
          betahat=solve(t(x)%*%x)%*%t(x)%*%y
          yhat=x%*%betahat
          sigmahat=(sum((y-yhat)^2))/(n-p-1)
          alphahat=solve(diag(lamb))%*%t(xstar)%*%y
          
          #OLS method
          betahat_ols[,i]=c(betahat)
          
          #Liu estimator 1993
          D.liu=abs(max(alphahat^2)-sigmahat)/(max(alphahat^2)+(1/max(ev)))
          D1=max(0,D.liu)
          #Shukur, mansson and sjolander 2015a
          D.sms=abs((alphahat^2)-sigmahat)/((alphahat^2)+(1/ev))
          D2=max(0,median(D.sms))
          D3=max(0,sum(D.sms)/p)
          #Khalaf and Shukur 2005
          D4=max(0,max(D.sms))
          D5=max(0,min(D.sms))
          
          #--------------------------------------------------------------------------------------
          #Proposed estimators
          #Based on D.sms inner part
          D6=max(0,abs(quantile(alphahat^2, probs = 0.00)-sigmahat)/(quantile(alphahat^2, probs = 0.00)+(1/quantile(ev, probs = 0.00))))
          D7=max(0,abs(quantile(alphahat^2, probs = 0.25)-sigmahat)/(quantile(alphahat^2, probs = 0.25)+(1/quantile(ev, probs = 0.25))))
          D8=max(0,abs(quantile(alphahat^2, probs = 0.50)-sigmahat)/(quantile(alphahat^2, probs = 0.50)+(1/quantile(ev, probs = 0.50))))
          D9=max(0,abs(quantile(alphahat^2, probs = 0.75)-sigmahat)/(quantile(alphahat^2, probs = 0.75)+(1/quantile(ev, probs = 0.75))))
          D10=max(0,abs(quantile(alphahat^2, probs = 1)-sigmahat)/(quantile(alphahat^2, probs = 1)+(1/quantile(ev, probs = 1))))
          
          wt=((pred+1/pred)^(1/pred))*(min(ev))
          MUS=D.sms*wt
          
          D11=max(0,mean(MUS))
          D12=max(0,median(MUS))
          D13=max(0,max(MUS)-min(MUS))
          D14=max(0,max(MUS))
          D15=max(0,min(MUS))
          #D.vector=c(D1,D2,D3,D4,D5,D6,D7,D8,D9,D10,D11,D12,D13,D14,D15)
          #--------------------------------------------------------------------------------------
          #Liu regression estimators
          alphahatD1=solve(t(xstar)%*%(xstar)+I)%*%(t(xstar)%*%(xstar)+D1*I)%*%alphahat
          alphahatD1_M[,i]=c(alphahatD1)
          
          alphahatD2=solve(t(xstar)%*%(xstar)+I)%*%(t(xstar)%*%(xstar)+D2*I)%*%alphahat
          alphahatD2_M[,i]=c(alphahatD2)
          
          alphahatD3=solve(t(xstar)%*%(xstar)+I)%*%(t(xstar)%*%(xstar)+D3*I)%*%alphahat
          alphahatD3_M[,i]=c(alphahatD3)
          
          alphahatD4=solve(t(xstar)%*%(xstar)+I)%*%(t(xstar)%*%(xstar)+D4*I)%*%alphahat
          alphahatD4_M[,i]=c(alphahatD4)
          
          alphahatD5=solve(t(xstar)%*%(xstar)+I)%*%(t(xstar)%*%(xstar)+D5*I)%*%alphahat
          alphahatD5_M[,i]=c(alphahatD5)
          
          alphahatD6=solve(t(xstar)%*%(xstar)+I)%*%(t(xstar)%*%(xstar)+D6*I)%*%alphahat
          alphahatD6_M[,i]=c(alphahatD6)
          
          alphahatD7=solve(t(xstar)%*%(xstar)+I)%*%(t(xstar)%*%(xstar)+D7*I)%*%alphahat
          alphahatD7_M[,i]=c(alphahatD7)
          
          alphahatD8=solve(t(xstar)%*%(xstar)+I)%*%(t(xstar)%*%(xstar)+D8*I)%*%alphahat
          alphahatD8_M[,i]=c(alphahatD8)
          
          alphahatD9=solve(t(xstar)%*%(xstar)+I)%*%(t(xstar)%*%(xstar)+D9*I)%*%alphahat
          alphahatD9_M[,i]=c(alphahatD9)
          
          alphahatD10=solve(t(xstar)%*%(xstar)+I)%*%(t(xstar)%*%(xstar)+D10*I)%*%alphahat
          alphahatD10_M[,i]=c(alphahatD10)
          alphahatD11=solve(t(xstar)%*%(xstar)+I)%*%(t(xstar)%*%(xstar)+D11*I)%*%alphahat
          alphahatD11_M[,i]=c(alphahatD11)
           alphahatD12=solve(t(xstar)%*%(xstar)+I)%*%(t(xstar)%*%(xstar)+D12*I)%*%alphahat
          alphahatD12_M[,i]=c(alphahatD12)
           alphahatD13=solve(t(xstar)%*%(xstar)+I)%*%(t(xstar)%*%(xstar)+D13*I)%*%alphahat
          alphahatD13_M[,i]=c(alphahatD13)
           alphahatD14=solve(t(xstar)%*%(xstar)+I)%*%(t(xstar)%*%(xstar)+D14*I)%*%alphahat
          alphahatD14_M[,i]=c(alphahatD14)
            alphahatD15=solve(t(xstar)%*%(xstar)+I)%*%(t(xstar)%*%(xstar)+D15*I)%*%alphahat
          alphahatD15_M[,i]=c(alphahatD15)
          
         
          #------------------------------------------------------------------------------------------------------          
          
        }     #Loop for N=5000 Runs ends
        
        #Mean square error estimation
        MSE_ols=sum((betahat_ols-c(beta))^2)/N
        
        MSE1<-sum((alphahatD1_M-c(alpha))^2)/N
        MSE2<-sum((alphahatD2_M-c(alpha))^2)/N
        MSE3<-sum((alphahatD3_M-c(alpha))^2)/N
        MSE4<-sum((alphahatD4_M-c(alpha))^2)/N
        MSE5<-sum((alphahatD5_M-c(alpha))^2)/N
        MSE6<-sum((alphahatD6_M-c(alpha))^2)/N
        MSE7<-sum((alphahatD7_M-c(alpha))^2)/N
        MSE8<-sum((alphahatD8_M-c(alpha))^2)/N
        MSE9<-sum((alphahatD9_M-c(alpha))^2)/N
        
        MSE10<-sum((alphahatD10_M-c(alpha))^2)/N
        MSE11<-sum((alphahatD11_M-c(alpha))^2)/N
        MSE12<-sum((alphahatD12_M-c(alpha))^2)/N
        MSE13<-sum((alphahatD13_M-c(alpha))^2)/N
        MSE14<-sum((alphahatD14_M-c(alpha))^2)/N
        
        MSE15<-sum((alphahatD15_M-c(alpha))^2)/N
        
        
        
        MSE=c(MSE_ols,MSE1,MSE2,MSE3,MSE4,MSE5,MSE6,MSE7,MSE8,MSE9,MSE10,MSE11,MSE12,MSE13,MSE14,MSE15)
        
        MSE=round(MSE,10)
        as.matrix(MSE)
        MSE.r[,q]=MSE
        #----------------------------------------------------------------------------------------------------------------------
        #Mean absolute error estimation
        MAE_ols=sum(abs(betahat_ols-c(beta)))/N
        
        MAE1<-sum(abs(alphahatD1_M-c(alpha)))/N
        MAE2<-sum(abs(alphahatD2_M-c(alpha)))/N
        MAE3<-sum(abs(alphahatD3_M-c(alpha)))/N
        MAE4<-sum(abs(alphahatD4_M-c(alpha)))/N
        MAE5<-sum(abs(alphahatD5_M-c(alpha)))/N
        MAE6<-sum(abs(alphahatD6_M-c(alpha)))/N
        MAE7<-sum(abs(alphahatD7_M-c(alpha)))/N
        MAE8<-sum(abs(alphahatD8_M-c(alpha)))/N
        MAE9<-sum(abs(alphahatD9_M-c(alpha)))/N
        
        MAE10<-sum(abs(alphahatD10_M-c(alpha)))/N
        MAE11<-sum(abs(alphahatD11_M-c(alpha)))/N
        MAE12<-sum(abs(alphahatD12_M-c(alpha)))/N
        MAE13<-sum(abs(alphahatD13_M-c(alpha)))/N
        MAE14<-sum(abs(alphahatD14_M-c(alpha)))/N
        
        MAE15<-sum(abs(alphahatD15_M-c(alpha)))/N
        
        MAE=c(MAE_ols,MAE1,MAE2,MAE3,MAE4,MAE5,MAE6,MAE7,MAE8,MAE9,MAE10,MAE11,MAE12,MAE13,MAE14,MAE15)
        
        MAE=round(MAE,4)
        as.matrix(MAE)
        MAE.r[,q]=MAE
        
      }   #End of correlation loop
      MSE.s[,,s]=MSE.r
      MAE.s[,,s]=MAE.r
    }     #End of error variance loop
    MSE.n[,,,a]=MSE.s
    MAE.n[,,,a]=MAE.s
  }       #End of sample size loop
  MSE.p[,,,,b]=MSE.n
  MAE.p[,,,,b]=MAE.n
}         #End of predictors loop

col.names=c("0.85","0.95","0.99","0.999")
row.names=c("OLS","D1","D2","D3","D4","D5","D6","D7","D8","D9","D10","D11","D12","D13","D14","D15")
matrix.names=c("0.1","1","2","5")
matrix.names2=c("25","50","100")
matrix.names3=c("4","10")
dimnames(MSE.p)=list(row.names,col.names,matrix.names,matrix.names2,matrix.names3)
dimnames(MAE.p)=list(row.names,col.names,matrix.names,matrix.names2,matrix.names3)

write.csv(MSE.p,file = "D:/PHD/MSE.quantile1.csv")
write.csv(MAE.p,file = "D:/PHD/MAE.quantile1.csv")



