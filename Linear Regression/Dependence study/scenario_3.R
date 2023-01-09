############################################
# Code for scenario 3: Toeplitz covariance
############################################

data_gener<-function(p,n,scen,rho){

# p=100 # nº covariates
# n=25 # sample size (n=25,50,100,200,400)
# scen="a" # scenario (scen="a","b")
if(scen=="a"){
  s=15
}else if(scen=="b"){
  s=10
}else{
  stop("Wrong scenario")
}
# rho=0.5 # rho value (rho=0.5,0.9)


#beta vector
beta=rep(0,p)
pos=ifelse(rep(s==15,s), 1:s , seq(1,100,by=10))
beta[pos]=0.5


#explanation of the 90% of deviance
per=0.9
if(scen=="a"){
  suma1=0
  for(jjj in 1:(s-1)){
    kkk=(jjj+1):15
    suma1=suma1+sum(rho^abs(jjj-kkk))
  }
  suma=suma1
}

if(scen=="b"){
  suma2=0
  secuencia=seq(1,100,by=10)
  for(jjj in 1:(s-1)){
    jjj=secuencia[jjj]
    kkk=secuencia[which(secuencia>jjj)]
    suma2=suma2+sum(rho^abs(jjj-kkk))
  }
  suma=suma2
}
dev=sqrt( ((1-per)/per)*(s*(0.5^2)+2*(0.5^2)*suma) )


#data generation
library(MASS)
Sig=matrix(0,p,p)
valores=rho^{0:(p-1)}
for(i in 1:p){
  Sig[i,1:p>=i]=Sig[1:p>=i,i]=valores[1:sum(1:p>=i)]
}
x=mvrnorm(n=n,mu=rep(0,p),Sigma=Sig)
x=x-matrix(colMeans(x),n,p,byrow=TRUE)
x=apply(x,2,function(z){z/sd(z)})
eps=rnorm(n=n, sd=dev)
y=x%*%beta+eps; y=as.numeric(y)
y=y-mean(y)


return(list("x"=x,"y"=y))}
