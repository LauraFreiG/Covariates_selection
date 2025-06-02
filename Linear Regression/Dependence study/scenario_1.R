############################################
# Code for scenario 1: orthogonal design
############################################

data_gener<-function(p,n,s){
  
# p=100 # nº covariates
# n=25 # sample size (n=25,50,100,200,400)
# s=10 # relevant terms (s=10,15,20)

#beta vector
nbeta=s #first s non-zero coefficients
beta=c(rep(1.25,nbeta),rep(0,p-nbeta))

#explanation of the 90% of deviance
per=0.9
dev=sqrt((1-per)/per*sum(beta^2)) # error deviation


#data generation
x=matrix(rnorm(p * n), n)
x=x-matrix(colMeans(x),n,p,byrow=TRUE)
x=apply(x,2,function(z){z/sd(z)})
eps=rnorm(n=n, sd=dev)
y=x%*%beta+eps; y=as.numeric(y)
y=y-mean(y)


return(list("x"=x,"y"=y))}