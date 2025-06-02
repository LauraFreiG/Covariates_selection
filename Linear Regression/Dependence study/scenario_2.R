############################################
# Code for scenario 2: dependence by blocks
############################################

data_gener<-function(p,n,s,rho){

# p=100 # nº covariates
# n=25 # sample size (n=25,50,100,200,400)
# s=10 # relevant terms (s=10,15,20)
# rho=0.5 # rho value (rho=0.5,0.9)

#beta vector
nbeta=s #first s non-zero coefficients
beta=c(rep(1,nbeta),rep(0,p-nbeta))

#explanation of the 90% of deviance
per=0.9
suma=ifelse(s>10,2*rho*(s-10),0)
dev=sqrt( ((1-per)/per)*(rho*s+suma) ) # error deviation


#data generation
library(MASS)
Sig=c()
for(i in 1:p){
  Sig_aux=rep(0,p)
  for(j in 1:p){
    Sig_aux[j]=ifelse(i%%10 == j%%10, rho, 0)    
  }
  Sig=cbind(Sig,Sig_aux)
}
x=mvrnorm(n=n,mu=rep(0,p),Sigma=Sig)
x=x-matrix(colMeans(x),n,p,byrow=TRUE)
x=apply(x,2,function(z){z/sd(z)})
eps=rnorm(n=n, sd=dev)
y=x%*%beta+eps; y=as.numeric(y)
y=y-mean(y)


return(list("x"=x,"y"=y))}