############################################################
# Code for scenario 2: Toeplitz covariance with unit scales
############################################################

data_gener<-function(p,n,rho,scen){

# p=100 # nº covariates
# n=400 # sample size (n=25,50,100,200,400)
# rho=0.5 # rho value (rho=0.5,0.9)
# scen="a" # scenario (scen="a","b")
  
  

beta=rep(0,p)
Sig=diag(1,p,p) #Sigma matrix
if(scen=="a"){
  s=15
  beta[1:s]=0.5
  sum1=0
  for(jjj in 1:(s-1)){
    kkk=(jjj+1):15
    sum1=sum1+sum(rho^abs(jjj-kkk))
  }
  sum_total=sum1
}else if(scen=="b"){
  s=10
  beta[seq(3,30,by=3)]=0.5
  sum2=0
  seq_vec=seq(3,30,by=3)
  seq_vec=seq_vec[1:s]
  for(jjj in 1:(s-1)){  
    jjj=seq_vec[jjj]
    kkk=seq_vec[which(seq_vec>jjj)]
    sum2=sum2+sum(rho^abs(jjj-kkk))
  }
  sum_total=sum2
}else{
  stop("Wrong scenario")
}  


#explanation of the 90% of deviance
per=0.9 #explanation of the 90% of deviance
dev=sqrt( ((1-per)/per)*(s*(0.5^2)+2*(0.5^2)*sum_total) ) # error deviation


#data generation
library(MASS)
Sig=matrix(0,p,p)
values=rho^{0:(p-1)}
for(i in 1:p){
  Sig[i,1:p>=i]=Sig[1:p>=i,i]=values[1:sum(1:p>=i)]
}
x=mvrnorm(n=n,mu=rep(0,p),Sigma=Sig)
x=x-matrix(colMeans(x),n,p,byrow=TRUE)
eps=rnorm(n=n, sd=dev)
y=x%*%beta+eps; y=as.numeric(y)
y=y-mean(y)


return(list("x"=x,"y"=y))}