rm(list=ls())
set.seed(1990)
library(MASS)
rm(list=ls())
set.seed(1990)
library(MASS)
p=8                                     #No of explanatory variables
I=diag(p)
setwd("E:/PHD/Real data")      #To change working directory
data=read.csv("US crime data Table10.18.csv",header=TRUE)
x1=data$X2
x2=data$X3
x3=data$X4
x4=data$X6
x5=data$X7
x6=data$X10
x7=data$X11
x8=data$X12

n<- length(x1)
x<-cbind(x1,x2,x3,x4,x5,x6,x7,x8)
x=scale(x,center = TRUE,scale = TRUE)   #Scaling x,,, standardizing x and y
y=data$Y
y=(y-mean(y))/sd(y)
#------------------------------------------
evec=eigen(cor(x))
ev=evec$values
e.vec=evec$vectors
beta=e.vec[,which.max(ev)]
c=t(x)%*%x
D=eigen(c)$vectors
xstar=x%*%D
#Z=x%*%D
lam=t(xstar)%*%xstar
lam=round(lam,4)
lamda=diag(lam)
alpha=t(D)%*%beta

beta.hat=solve(t(x)%*%x)%*%t(x)%*%y
y.hat=x%*%beta.hat
sigmahat=sum((y-y.hat)^2)/(n-p-1)
alphahat=solve(lam)%*%t(xstar)%*%y

#OLS estimator
alphahat.ols=alphahat
alphahat.ols=c(alphahat.ols)
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
          wt=((p+1/p)^(1/p))*(min(ev))
          MUS=D.sms*wt
          D11=max(0,mean(MUS))
          D12=max(0,median(MUS))
          D13=max(0,max(MUS)-min(MUS))
          D14=max(0,max(MUS))
          D15=max(0,min(MUS))
D.vector=c(D1,D2,D3,D4,D5,D6,D7,D8,D9,D10,D11,D12,D13,D14,D15)
#Liu regression estimators
alphahatD1=solve(t(xstar)%*%(xstar)+I)%*%(t(xstar)%*%(xstar)+D1*I)%*%alphahat
alphahatD1=c(alphahatD1)

alphahatD2=solve(t(xstar)%*%(xstar)+I)%*%(t(xstar)%*%(xstar)+D2*I)%*%alphahat
alphahatD2=c(alphahatD2)

alphahatD3=solve(t(xstar)%*%(xstar)+I)%*%(t(xstar)%*%(xstar)+D3*I)%*%alphahat
alphahatD3=c(alphahatD3)

alphahatD4=solve(t(xstar)%*%(xstar)+I)%*%(t(xstar)%*%(xstar)+D4*I)%*%alphahat
alphahatD4=c(alphahatD4)

alphahatD5=solve(t(xstar)%*%(xstar)+I)%*%(t(xstar)%*%(xstar)+D5*I)%*%alphahat
alphahatD5=c(alphahatD5)

alphahatD6=solve(t(xstar)%*%(xstar)+I)%*%(t(xstar)%*%(xstar)+D6*I)%*%alphahat
alphahatD6=c(alphahatD6)

alphahatD7=solve(t(xstar)%*%(xstar)+I)%*%(t(xstar)%*%(xstar)+D7*I)%*%alphahat
alphahatD7=c(alphahatD7)

alphahatD8=solve(t(xstar)%*%(xstar)+I)%*%(t(xstar)%*%(xstar)+D8*I)%*%alphahat
alphahatD8=c(alphahatD8)

alphahatD9=solve(t(xstar)%*%(xstar)+I)%*%(t(xstar)%*%(xstar)+D9*I)%*%alphahat
alphahatD9=c(alphahatD9)

alphahatD10=solve(t(xstar)%*%(xstar)+I)%*%(t(xstar)%*%(xstar)+D10*I)%*%alphahat
alphahatD10=c(alphahatD10)
####################################
alphahatD11=solve(t(xstar)%*%(xstar)+I)%*%(t(xstar)%*%(xstar)+D11*I)%*%alphahat
alphahatD11=c(alphahatD11)

alphahatD12=solve(t(xstar)%*%(xstar)+I)%*%(t(xstar)%*%(xstar)+D12*I)%*%alphahat
alphahatD12=c(alphahatD12)

alphahatD13=solve(t(xstar)%*%(xstar)+I)%*%(t(xstar)%*%(xstar)+D13*I)%*%alphahat
alphahatD13=c(alphahatD13)

alphahatD14=solve(t(xstar)%*%(xstar)+I)%*%(t(xstar)%*%(xstar)+D14*I)%*%alphahat
alphahatD14=c(alphahatD14)

alphahatD15=solve(t(xstar)%*%(xstar)+I)%*%(t(xstar)%*%(xstar)+D15*I)%*%alphahat
alphahatD15=c(alphahatD15)
#Computing regression coefficients
comb.coeff=rbind(alphahat.ols,alphahatD1,alphahatD2,alphahatD3,alphahatD4,alphahatD5,alphahatD6,
                 alphahatD7,alphahatD8,alphahatD9,alphahatD10,alphahatD11,alphahatD12,alphahatD13,alphahatD14,alphahatD15)
comb.coeff=round(comb.coeff,5)
write.csv(comb.coeff,file = "E:/PHD/coefficents.chemical.csv")
write.csv(D.vector,file = "E:/PHD/D-vector.csv")

#------------------------------------------------------------------------------------------------------          
#MSE computation
MSE.OLS=sigmahat*sum(1/ev)
#MSE of all Liu estimators

a=rep(0,length(D.vector)); b=rep(0,length(D.vector))                            #to obtain MSE of ridge estimators
for(i in 1:length(D.vector)){
  
  a[i]=sum(((ev+D.vector[i])^2)/(ev*((ev+1)^2)))
  b[i]=((D.vector[i]-1)^2)*sum((alphahat^2)/((ev+1)^2))
  
}
MSE.Liu=sigmahat*a+b
MSE=c(MSE.OLS,MSE.Liu)
MSE=as.matrix(MSE)
MSE=round(MSE,4)
row.names(MSE)=c("OLS","D1","D2","D3","D4","D5","D6","D7","D8","D9","D10","D11","D12","D13","D14","D15")
colnames(MSE)=c("MSE")
write.csv(MSE,file ="E:/MSE.tabacco.csv")
#Prediction interval
x0 <- c(122,170,160,1080,170,60,695,28)
as.vector(x0)
y0.ols=t(x0)%*%alphahat.ols
y0.D1=t(x0)%*%alphahatD1
y0.D2=t(x0)%*%alphahatD2
y0.D3=t(x0)%*%alphahatD3
y0.D4=t(x0)%*%alphahatD4
y0.D5=t(x0)%*%alphahatD5
y0.D6=t(x0)%*%alphahatD6
y0.D7=t(x0)%*%alphahatD7
y0.D8=t(x0)%*%alphahatD8
y0.D9=t(x0)%*%alphahatD9
y0.D10=t(x0)%*%alphahatD10
#############################
y0.D11=t(x0)%*%alphahatD11
y0.D12=t(x0)%*%alphahatD12
y0.D13=t(x0)%*%alphahatD13
y0.D14=t(x0)%*%alphahatD14
y0.D15=t(x0)%*%alphahatD15
y0.vector=c(y0.ols,y0.D1,y0.D2,y0.D3,y0.D4,y0.D5,y0.D6,y0.D7,y0.D8,y0.D9,y0.D10,y0.D11,y0.D12,y0.D13,y0.D14,y0.D15)
y0.vector=as.vector(y0.vector)
var.y0=sigmahat*(t(x0)%*%solve(t(x)%*%x)%*%x0)
se.y0.ols=sqrt(var.y0)

QD1=solve(t(x)%*%x+I)%*%(t(x)%*%x+D1*I)
QD2=solve(t(x)%*%x+I)%*%(t(x)%*%x+D2*I)
QD3=solve(t(x)%*%x+I)%*%(t(x)%*%x+D3*I)
QD4=solve(t(x)%*%x+I)%*%(t(x)%*%x+D4*I)
QD5=solve(t(x)%*%x+I)%*%(t(x)%*%x+D5*I)
QD6=solve(t(x)%*%x+I)%*%(t(x)%*%x+D6*I)
QD7=solve(t(x)%*%x+I)%*%(t(x)%*%x+D7*I)
QD8=solve(t(x)%*%x+I)%*%(t(x)%*%x+D8*I)
QD9=solve(t(x)%*%x+I)%*%(t(x)%*%x+D9*I)
QD10=solve(t(x)%*%x+I)%*%(t(x)%*%x+D10*I)
#############
QD11=solve(t(x)%*%x+I)%*%(t(x)%*%x+D11*I)
QD12=solve(t(x)%*%x+I)%*%(t(x)%*%x+D12*I)
QD13=solve(t(x)%*%x+I)%*%(t(x)%*%x+D13*I)
QD14=solve(t(x)%*%x+I)%*%(t(x)%*%x+D14*I)
QD15=solve(t(x)%*%x+I)%*%(t(x)%*%x+D15*I)

se.y0.D1=sqrt(sigmahat*(t(x0)%*%QD1%*%solve(t(x)%*%x)%*%t(QD1)%*%x0))
se.y0.D2=sqrt(sigmahat*(t(x0)%*%QD2%*%solve(t(x)%*%x)%*%t(QD2)%*%x0))
se.y0.D3=sqrt(sigmahat*(t(x0)%*%QD3%*%solve(t(x)%*%x)%*%t(QD3)%*%x0))
se.y0.D4=sqrt(sigmahat*(t(x0)%*%QD4%*%solve(t(x)%*%x)%*%t(QD4)%*%x0))
se.y0.D5=sqrt(sigmahat*(t(x0)%*%QD5%*%solve(t(x)%*%x)%*%t(QD5)%*%x0))
se.y0.D6=sqrt(sigmahat*(t(x0)%*%QD6%*%solve(t(x)%*%x)%*%t(QD6)%*%x0))
se.y0.D7=sqrt(sigmahat*(t(x0)%*%QD7%*%solve(t(x)%*%x)%*%t(QD7)%*%x0))
se.y0.D8=sqrt(sigmahat*(t(x0)%*%QD8%*%solve(t(x)%*%x)%*%t(QD8)%*%x0))
se.y0.D9=sqrt(sigmahat*(t(x0)%*%QD9%*%solve(t(x)%*%x)%*%t(QD9)%*%x0))
se.y0.D10=sqrt(sigmahat*(t(x0)%*%QD10%*%solve(t(x)%*%x)%*%t(QD10)%*%x0))
####################
se.y0.D11=sqrt(sigmahat*(t(x0)%*%QD11%*%solve(t(x)%*%x)%*%t(QD11)%*%x0))
se.y0.D12=sqrt(sigmahat*(t(x0)%*%QD12%*%solve(t(x)%*%x)%*%t(QD12)%*%x0))
se.y0.D13=sqrt(sigmahat*(t(x0)%*%QD13%*%solve(t(x)%*%x)%*%t(QD13)%*%x0))
se.y0.D14=sqrt(sigmahat*(t(x0)%*%QD14%*%solve(t(x)%*%x)%*%t(QD14)%*%x0))
se.y0.D15=sqrt(sigmahat*(t(x0)%*%QD15%*%solve(t(x)%*%x)%*%t(QD15)%*%x0))
se.y0.vector=c(se.y0.ols,se.y0.D1,se.y0.D2,se.y0.D3,se.y0.D4,se.y0.D5,se.y0.D6,se.y0.D7,se.y0.D8,se.y0.D9,se.y0.D10,se.y0.D11,se.y0.D12,se.y0.D13,se.y0.D14,se.y0.D15)
se.y0.vector=as.vector(se.y0.vector)
v.df=n-p
t.table=qt(0.975, df=v.df, lower.tail = TRUE, log.p = FALSE)
LCL_y0=rep(0,length(y0.vector))
UCL_y0=rep(0,length(y0.vector))
for (i in 1:length(y0.vector)) {
  LCL_y0[i]=y0.vector[i]-(t.table*se.y0.vector[i])
  UCL_y0[i]=y0.vector[i]+(t.table*se.y0.vector[i])
}
Pred.int=round(cbind(LCL_y0,UCL_y0),5)
row.names(Pred.int)=c("OLS","D1","D2","D3","D4","D5","D6","D7","D8","D9","D10","D11","D12","D13","D14","D15")
write.csv(Pred.int,file = "D:/PHD/Pred.int.csv")
# Difference (width) of prediction intervals
PI_width = UCL_y0 - LCL_y0
PI_width = round(PI_width, 5)

# Combine with existing results (optional)
PI_results = cbind(LCL_y0, UCL_y0, PI_width)
PI_results = round(PI_results, 10)

row.names(PI_results) = c("OLS","D1","D2","D3","D4","D5","D6","D7","D8","D9","D10","D11","D12","D13","D14","D15")
colnames(PI_results) = c("LCL","UCL","Width")

# Save to CSV
write.csv(PI_results, file = "D:/PHD/PI_width.csv")


#OLS
xnam <- paste0("x", 1:p) 
fmla<- as.formula(paste("y ~ ", paste(xnam, collapse= "+"))) 
model<-lm(fmla)
sigma2<-sum(model$residuals^2)/(n-p-1)
sigmahat.e=round(sqrt(sigma2),4)
#eigen values
ev1=eigen(cor(x))$values
ev1=round(ev1,4)
max.per.ev1=(max(ev1)/8)*100
max.per.ev1=round(max.per.ev1,2)

cor.x=round(cor(x),4)
cor.y=round(cor(x,y),4)
cor.xy=cbind(cor.x,cor.y)
write.csv(cor.xy,file = "E:/Corr.US.Crime.csv")
cor.xy
round(det(cor(x)),6)
MSE
ev1
max.per.ev1
sigma2
sigmahat.e
CN=max(ev1)/min(ev1)        #Computing condition number
CN1=sqrt(CN)
CN
CN1
#Checking normality of errors
res.ols=model$residuals
hist(res.ols)
#lines(density(res.ols),col="red")
ks.test(res.ols,"pnorm",3)
#library(nortest)
shapiro.test(res.ols)
library(boot)
round(k3.linear(res.ols),4)

#Checking outliers
#Influence measures
inflm.model <- influence.measures(model)
which(apply(inflm.model$is.inf, 1, any))      #outliers
summary(inflm.model)                          #Summary of outlying observations
#----------------------------------------------------------
# Variance Inflation Factor (VIF)
#----------------------------------------------------------

# Install car package if not already installed
# install.packages("car")

library(car)

# Calculate VIF values
VIF <- vif(model)

# Round VIF values
VIF <- round(VIF, 4)
VIF
