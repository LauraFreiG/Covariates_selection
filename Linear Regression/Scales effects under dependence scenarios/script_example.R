############################################
# Example script
############################################


# Data generation
#---------------------------
# Scenario 1
stand="yes" # applying univariate standardization
source("scenario_1.R")
p=100 # nº covariates
n=300 # sample size (n=25,50,100,150,300)
s=10 # relevant terms
scen="b" # scenario (scen="a","b","c")
data=data_gener(p,n,s,scen)

# # Scenario 2
# stand="yes" # applying univariate standardization
# source("scenario_2.R")
# p=100 # nº covariates
# n=300 # sample size (n=25,50,100,150,300)
# rho=0.5 # rho value (rho=0.5,0.9)
# scen="a" # scenario (scen="a","b")
# data=data_gener(p,n,rho,scen)

# # Scenario 3
# stand="yes" # applying univariate standardization
# source("scenario_3.R")
# p=100 # nº covariates
# n=300 # sample size (n=25,50,100,150,300)
# rho=0.5 # rho value (rho=0.5,0.9)
# scen="a" # scenario (scen="a","b")
# data=data_gener(p,n,rho,scen)

x=data$x
if(stand=="yes"){
  x=apply(x,2,function(z){z/sd(z)})
}
y=data$y



# Covariates selection techniques
#-----------------------------------
sel=rep(0,11) #nº of selected covariates
MSE=rep(0,11) #estimated MSE
pdev=rep(0,11) #estimated %Dev
nam_cov=paste("X",1:p,sep="")
cov_sel=list() #names of selected covariates


#LASSO.min
library(glmnet)
nfold=10
cv_fit_LR <- cv.glmnet(x, y,family="gaussian", alpha = 1, standardize=FALSE, nfolds=nfold,
                       intercept=FALSE, standardize.response=FALSE)
lambda=cv_fit_LR$lambda.min
fit_LR.min=glmnet(x,y, family="gaussian", alpha = 1, standardize=FALSE, lambda=lambda,
                  intercept=FALSE, standardize.response=FALSE)
sel[1]=sum(fit_LR.min$beta!=0)
ind=which(fit_LR.min$beta!=0)
cov_sel[[1]]=nam_cov[ind]
if(length(ind)==0){
  MSE[1]=sum(y^2)/n
  pdev[1]=0
}else{
  RSS=sum((y-predict(lm(y~x[,ind])))^2)
  MSE[1]=RSS/n
  pdev[1]=(sum(y^2)-RSS)/sum(y^2)
}



#LASSO.1se
library(glmnet)
nfold=10
cv_fit_LR <- cv.glmnet(x, y,family="gaussian", alpha = 1, standardize=FALSE, nfolds=nfold,
                       intercept=FALSE, standardize.response=FALSE)
lambda=cv_fit_LR$lambda.1se
fit_LR.1se=glmnet(x,y, family="gaussian", alpha = 1, standardize=FALSE, lambda=lambda,
                  intercept=FALSE, standardize.response=FALSE)
sel[2]=sum(fit_LR.1se$beta!=0)
ind=which(fit_LR.1se$beta!=0)
cov_sel[[2]]=nam_cov[ind]
if(length(ind)==0){
  MSE[2]=sum(y^2)/n
  pdev[2]=0
}else{
  RSS=sum((y-predict(lm(y~x[,ind])))^2)
  MSE[2]=RSS/n
  pdev[2]=(sum(y^2)-RSS)/sum(y^2)
}



#LASSO.BIC
library(glmnet)
nfold=10 
cv_fit_LR <- cv.glmnet(x, y,family="gaussian", alpha = 1, standardize=FALSE, nfolds=nfold,
                       intercept=FALSE, standardize.response=FALSE)
#Computing BIC values for lambda grid of glmnet---
lambda_seq=cv_fit_LR$lambda
BIC=rep(0,length(lambda_seq))
for(i in  1:length(BIC)){
  lambda_aux=cv_fit_LR$glmnet.fit$lambda[i]
  indices=which(cv_fit_LR$glmnet.fit$beta[,i]!=0)
  if(length(indices)!=0){ 
    MSE_aux=sum((y-predict(lm(y~x[,indices])))^2)/n
    BIC[i]=n*log(MSE_aux)+length(indices)*log(n)
  }
}
BIC[BIC==0]=max(BIC[BIC!=0])+1 #do not select empty models
BIC_opt=min(BIC); if(length(BIC_opt)>1){BIC_opt=BIC_opt[1]}
#In case of a tie, simpler model
ind_aux=which(BIC==BIC_opt); if(length(ind_aux)>1){ind_aux=ind_aux[1]}
lambda_opt=cv_fit_LR$lambda[ind_aux]
coeficientes=which(cv_fit_LR$glmnet.fit$beta[,ind_aux]!=0)
sel[3]=length(coeficientes)
ind=coeficientes
cov_sel[[3]]=nam_cov[ind]
if(length(ind)==0){
  MSE[3]=sum(y^2)/n
  pdev[3]=0
}else{
  RSS=sum((y-predict(lm(y~x[,ind])))^2)
  MSE[3]=RSS/n
  pdev[3]=(sum(y^2)-RSS)/sum(y^2)
}



#AdapL.min
library(glmnet)
nfold=10
cv.ridge <- cv.glmnet(x,y,family="gaussian", alpha=0, parallel=F, standardize=FALSE,
                      standardize.response=FALSE, intercept=FALSE, nfold=nfold)
w3 <- 1/abs(matrix(coef(cv.ridge, s=cv.ridge$lambda.min)[, 1][2:(p+1)] ))^1
w3[w3[,1] == Inf] <- 999999999 ## Replacing values estimated as Infinite for 999999999
beta=rep(0,p)
iter=TRUE
while(iter==TRUE){
  cv_fit_LR <- cv.glmnet(x, y,family="gaussian", alpha = 1, standardize=FALSE, 
                         standardize.response=FALSE, intercept=FALSE, nfold=nfold,
                         penalty.factor=w3)
  w <- cv_fit_LR$lambda.min
  fit_LR=glmnet(x,y, family="gaussian", alpha = 1, standardize=FALSE, lambda=w,
                standardize.response=FALSE, intercept=FALSE,penalty.factor=w3)
  norma=norm(beta,"2")
  if(norma<1e-12){norma=1e-6}
  dif=norm(beta-as.vector(fit_LR$beta),"2")/norma
  if(dif<1e-2){iter=FALSE}
  beta=as.vector(fit_LR$beta)
  w3<-1/abs(beta)
  w3[w3 == Inf] <- 999999999 ## Replacing values estimated as Infinite for 999999999
}
fit_AdapL.min=fit_LR
sel[4]=sum(fit_AdapL.min$beta!=0)
ind=which(fit_AdapL.min$beta!=0)
cov_sel[[4]]=nam_cov[ind]
if(length(ind)==0){
  MSE[4]=sum(y^2)/n
  pdev[4]=0
}else{
  RSS=sum((y-predict(lm(y~x[,ind])))^2)
  MSE[4]=RSS/n
  pdev[4]=(sum(y^2)-RSS)/sum(y^2)
}



#AdapL.1se
library(glmnet)
nfold=10
cv.ridge <- cv.glmnet(x,y,family="gaussian", alpha=0, parallel=F, standardize=FALSE,
                      standardize.response=FALSE, intercept=FALSE, nfold=nfold)
w3 <- 1/abs(matrix(coef(cv.ridge, s=cv.ridge$lambda.1se)[, 1][2:(p+1)] ))^1
w3[w3[,1] == Inf] <- 999999999 ## Replacing values estimated as Infinite for 999999999
beta=rep(0,p)
iter=TRUE
while(iter==TRUE){
  cv_fit_LR <- cv.glmnet(x, y,family="gaussian", alpha = 1, standardize=FALSE, 
                         standardize.response=FALSE, intercept=FALSE, nfold=nfold,
                         penalty.factor=w3)
  w <- cv_fit_LR$lambda.1se
  fit_LR=glmnet(x,y, family="gaussian", alpha = 1, standardize=FALSE, lambda=w,
                standardize.response=FALSE, intercept=FALSE,penalty.factor=w3)
  norma=norm(beta,"2")
  if(norma<1e-12){norma=1e-6}
  dif=norm(beta-as.vector(fit_LR$beta),"2")/norma
  if(dif<1e-2){iter=FALSE}
  beta=as.vector(fit_LR$beta)
  w3<-1/abs(beta)
  w3[w3 == Inf] <- 999999999 ## Replacing values estimated as Infinite for 999999999
}
fit_AdapL.1se=fit_LR
sel[5]=sum(fit_AdapL.1se$beta!=0)
ind=which(fit_AdapL.1se$beta!=0)
cov_sel[[5]]=nam_cov[ind]
if(length(ind)==0){
  MSE[5]=sum(y^2)/n
  pdev[5]=0
}else{
  RSS=sum((y-predict(lm(y~x[,ind])))^2)
  MSE[5]=RSS/n
  pdev[5]=(sum(y^2)-RSS)/sum(y^2)
}



#SCAD
library(ncvreg)
nfold=10
cv_fit_SCAD <- cv.ncvreg(x, y,family="gaussian", penalty="SCAD", nfolds=nfold)
lambda=cv_fit_SCAD$lambda.min
fit_SCAD=ncvreg(x, y,family="gaussian", penalty="SCAD", lambda=lambda)
sel[6]=sum(fit_SCAD$beta[-1]!=0)
ind=which(fit_SCAD$beta[-1]!=0)
cov_sel[[6]]=nam_cov[ind]
if(length(ind)==0){
  MSE[6]=sum(y^2)/n
  pdev[6]=0
}else{
  RSS=sum((y-predict(lm(y~x[,ind])))^2)
  MSE[6]=RSS/n
  pdev[6]=(sum(y^2)-RSS)/sum(y^2)
}



#Dant
library(flare)
k=10
seq_lambda=slim(x,y,nlambda=100,method="dantzig",verbose=FALSE)$lambda
#Function to calculate the MSE value for each lambda value of the grid
val<-function(l){
  lambda_aux=l
  proof2<-function(kk){
    ind_muest=((n/k)*(kk-1)+1):((n/k)*(kk))
    x_valid=x[ind_muest,]; y_valid=y[ind_muest]
    x_entre=x[-ind_muest,]; y_entre=y[-ind_muest]
    dantzig=slim(x_entre,y_entre,lambda=lambda_aux,method="dantzig",verbose=FALSE)
    MSE_aux=sum( (y_valid-x_valid%*%dantzig$beta)^2 )/length(y_valid)
    return(MSE_aux)
  }
  MSE_aux=apply(as.matrix(1:k),1,proof2)
}
#Estimation of the lambda min value
grid_MSE=colMeans(apply(as.matrix(seq_lambda),1,val))
#lambda minimizing MSE estimated by k-fold cross validation
ind_lambda=which(grid_MSE==min(grid_MSE))
lambda_min=seq_lambda[ind_lambda]
#Adjustment with optimal parameters
fit_Dantzig=slim(x,y,method="dantzig",lambda=lambda_min,verbose=FALSE)
sel[7]=sum(fit_Dantzig$beta!=0)
ind=which(fit_Dantzig$beta!=0)
cov_sel[[7]]=nam_cov[ind]
if(length(ind)==0){
  MSE[7]=sum(y^2)/n
  pdev[7]=0
}else{
  RSS=sum((y-predict(lm(y~x[,ind])))^2)
  MSE[7]=RSS/n
  pdev[7]=(sum(y^2)-RSS)/sum(y^2)
}



#RelaxL
library(relaxo)
nfold=10
error=TRUE
while(error){
  cv_fit_RelaxL <- tryCatch(cvrelaxo(x, y, phi=seq(1e-4, 1, length = 10), k=nfold, 
                                     warn=FALSE),
                            error= function(z){TRUE})
  if(length(cv_fit_RelaxL)!=1){error=FALSE }
}
sel[8]=sum(cv_fit_RelaxL$beta!=0)
ind=which(cv_fit_RelaxL$beta!=0)
cov_sel[[8]]=nam_cov[ind]
if(length(ind)==0){
  MSE[8]=sum(y^2)/n
  pdev[8]=0
}else{
  RSS=sum((y-predict(lm(y~x[,ind])))^2)
  MSE[8]=RSS/n
  pdev[8]=(sum(y^2)-RSS)/sum(y^2)
}



#SqrtL
library(flare)
k=10
seq_lambda=slim(x,y,nlambda=100,verbose=FALSE)$lambda
#Function to calculate the MSE value for each lambda value of the grid
val<-function(l){
  lambda_aux=l
  proof2<-function(kk){
    ind_muest=((n/k)*(kk-1)+1):((n/k)*(kk))
    x_valid=x[ind_muest,]; y_valid=y[ind_muest]
    x_entre=x[-ind_muest,]; y_entre=y[-ind_muest]
    dantzig=slim(x_entre,y_entre,lambda=lambda_aux,verbose=FALSE)
    MSE_aux=sum( (y_valid-x_valid%*%dantzig$beta)^2 )/length(y_valid)
    return(MSE_aux)
  }
  MSE_aux=apply(as.matrix(1:k),1,proof2)
}
#Estimation of the lambda min value
grid_MSE=colMeans(apply(as.matrix(seq_lambda),1,val))
#lambda minimizing MSE estimated by k-fold cross validation
ind_lambda=which(grid_MSE==min(grid_MSE))
lambda_min=seq_lambda[ind_lambda]
#Adjustment with optimal parameters
fit_square=slim(x,y,lambda=lambda_min,verbose=FALSE)
sel[9]=sum(fit_square$beta!=0)
ind=which(fit_square$beta!=0)
cov_sel[[9]]=nam_cov[ind]
if(length(ind)==0){
  MSE[9]=sum(y^2)/n
  pdev[9]=0
}else{
  RSS=sum((y-predict(lm(y~x[,ind])))^2)
  MSE[9]=RSS/n
  pdev[9]=(sum(y^2)-RSS)/sum(y^2)
}



#ScalL
library(scalreg)
fit_ScalL=scalreg(x,y)
sel[10]=sum(fit_ScalL$coefficients!=0)
ind=which(fit_ScalL$coefficients!=0)
cov_sel[[10]]=nam_cov[ind]
if(length(ind)==0){
  MSE[10]=sum(y^2)/n
  pdev[10]=0
}else{
  RSS=sum((y-predict(lm(y~x[,ind])))^2)
  MSE[10]=RSS/n
  pdev[10]=(sum(y^2)-RSS)/sum(y^2)
}



#DC.VS
source("gam.vs.R")
dat=data.frame("x"=x,"y"=y)
sal=gam.vs(resp="y",dat=dat)
fit_gam=sal[[3]]
i_pred=sal[[5]]
sel[11]=sum(i_pred!=0)
ind=which(i_pred!=0)
cov_sel[[11]]=nam_cov[ind]
if(length(ind)==0){
  MSE[11]=sum(y^2)/n
  pdev[11]=0
}else{
  RSS=sum((y-predict(lm(y~x[,ind])))^2)
  MSE[11]=RSS/n
  pdev[11]=(sum(y^2)-RSS)/sum(y^2)
}



#nº selected covariates:
sel

#names of selected covariates:
names(cov_sel)=c("LASSO.min","LASSO.1se","LASSO.BIC","#AdapL.min","#AdapL.1se",
                 "SCAD","Dant","RelaxL","SqrtL","ScalL","DC.VS")
cov_sel

#estimated MSE:
MSE

#estimated % of deviance:
pdev





