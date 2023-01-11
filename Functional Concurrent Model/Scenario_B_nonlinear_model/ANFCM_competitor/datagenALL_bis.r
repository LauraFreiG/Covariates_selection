######################################################
# datagenALL.r - April/25/2014			#
# Script created by Janet Kim			#
# Contents:					#
# (i)	X, W					#
# (ii)	F.nlin, F.nlin2, F.lin, F.lin2, F.anova1, F.anova2	#
# (iii)	build.AR, Error_resp			#
# (iv)	datagen, anova.datagen			#
#						#
# @ Sample size in training data:			#
#   n subjects, m observation points per subject		#
# @ Sample size in test data:				#
#   n.test subjects, m observation points per subject	#
######################################################

# Functional Covariate without errors
X = function(n,t,Q=1){ 
  # dim(t) = m-by-1
  # Q: number of multiple functional predictors
  m = length(t)
  a0 = rnorm(n)
  a1 = rnorm(n, 0, sd=0.85)
  a2 = rnorm(n, 0, sd=0.70)
  ones = rep(1,m)
  phi1 = sqrt(2)*sin(pi*t)				# m-by-1
  phi2 = sqrt(2)*cos(pi*t)				# m-by-1
  X.list = as.list(NULL)
  for (q in 1:Q){
    # X.list[[q]] = 2^(1-q)*(a0%*%t(ones)+a1%*%t(phi1)+a2%*%t(phi2))
    X.list[[q]] = (a0%*%t(ones)+a1%*%t(phi1)+a2%*%t(phi2))/sqrt(2^(q-1))
  }
  X.list
}

# Functional Covariate with errors
W = function(x,tau){ # dim(x) = n-by-m
  n = nrow(x)
  m = ncol(x)
  error.x = mvrnorm(n, mu=rep(0,m), Sigma=tau*diag(m))	# n-by-m
  obs = x + error.x					# n-by-m
  obs									
}

# True Nonlinear Function of x and t
F.nlin = function(x,t){1 + x + t + 2*x^2*t}
F.nlin2 = function(x,t){0.75*exp(x*t)}
# True Linear Function f x and t
F.lin = function(x,t){1 + x + t}
F.lin2 = function(x,t){x*t}
# Null model used in hypothesis testing
F.anova1 = function(x,t,d){1+2*t+t^2+d*(x*t/8)}
# F.anova2 = function(x1,x2,t,d){2*t+t^2+x1*sin(pi*t)/4+d*2*cos(x2*t)}
F.anova2 = function(x2,t,d){2*t+t^2+d*2*cos(x2*t)}

# Build AR(1) m-by-m covariance matrix for the random errors
build.AR = function(m, sigma2, rho){ 
  power = 0:(m-1)
  H= abs(outer(power, power, "-"))			# m-by-m
  V = sigma2*(rho^H)				# m-by-m
  V
}

# Create Random Errors
Error_resp = function(n, t, error_type, sig2noise, rho){ # dim(t) = m-by-1
  m = length(t)
  Error0 = mvrnorm(n, mu=rep(0,m), Sigma=sig2noise*diag(m))
  if (error_type==1) {# Independent
    Error = Error0
  }
  if (error_type==2) {# AR(1) + WN
    Error.temp = mvrnorm(n, mu=rep(0,m), Sigma=build.AR(m, sig2noise, rho))
    Error = Error.temp+ Error0
  }
  if (error_type==3) {# KL Expansion + WN
    x1i = sqrt(2)*rnorm(n)
    x2i = rnorm(n,0,sd=0.75)
    Error = x1i%*%t(sqrt(2)*cos(pi*t))+x2i%*%t(sqrt(2)*sin(pi*t))+ Error0
  }  		
  Error						# n-by-m
}


# Generate training data set 
datagen = function(n, m, mi_low=NULL, mi_up=NULL, 
                   F=list(F.lin), Q=1, error_type, tau, sig2noise, rho=NULL, seed){
  # F: a list of true functions, e.g., Y(t) = F.nlin(x1,t) + F.nlin2(x2,t)
  # Dense design: "m" time points per curve
  # Sparse Design: "mi_low" ~ "mi_up" time points per curve

  set.seed(seed*10)
  t=seq(0, 1, len=m)					# m-by-1
  trep = matrix(rep(t,n), nrow=n, byrow=TRUE)       	# n-by-m
  
  # a list of functional covariates
  Xfull = X(n,t,Q)
  Wfull = NULL
  for (q in 1:Q){ Wfull[[q]] = W(Xfull[[q]], tau)}  		# Noisy Covariate; n-by-m
  
  if (is.null(mi_low)&is.null(mi_up)){
    Xeval <- lapply(c(1:Q), function(q) Xfull[[q]])
    Weval <- lapply(c(1:Q), function(q) Wfull[[q]])
  } else{
    Xeval = lapply(c(1:Q), function(q) NA*Xfull[[q]])
    Weval = lapply(c(1:Q), function(q) NA*Wfull[[q]])
    mi = sample(mi_low:mi_up, n, replace=TRUE)
    for (i in 1:n){
      ind = sample(1:m, mi[i], replace=FALSE)
      for (q in 1:Q){
        Xeval[[q]][i, ind] = Xfull[[q]][i, ind]
        Weval[[q]][i, ind] = Wfull[[q]][i, ind]
      }
    }
  }
  
  # True Function F evaluated at "Xfull" and "trep"; n-by-m
  Ffull = lapply(c(1:Q), function(q) F[[q]](Xfull[[q]], trep))
  # True Function F evaluated at "Xeval" and "trep"; n-by-m
  Feval = lapply(c(1:Q), function(q) F[[q]](Xeval[[q]], trep))
  
  # Random Error; n-by-m
  Error = Error_resp(n,t,error_type,sig2noise,rho)
  
  # Response based on Ffull; n-by-m
  Yfull = Reduce('+', lapply(c(1:Q), function(q) Ffull[[q]])) + Error
  # Response based on Feval; n-by-m 
  Yeval = Reduce('+', lapply(c(1:Q), function(q) Feval[[q]])) + Error

  result = list(Xfull = Xfull, Wfull = Wfull, Xeval = Xeval, Weval = Weval, 
                Error = Error, Ffull = Ffull, Feval = Feval, Yfull = Yfull, Yeval = Yeval) 
}

# Generate data set based on null model of hypothesis testing
anova.datagen = function(n, m, mi_low=NULL, mi_up=NULL, F, Q=2, d, error_type, tau, sig2noise, rho, seed){ 
  # Dense design: "m" time points per curve
  # Sparse Design: "mi_low" ~ "mi_up" time points per curve
  
  set.seed(seed)
  t=seq(0, 1, len=m)                                  		# m-by-1
  trep = matrix(rep(t,n), nrow=n, byrow=TRUE)         	# n-by-m
  
  # a list of functional covariates
  Xfull = X(n,t,Q)
  Wfull = NULL
  for (q in 1:Q){ Wfull[[q]] = W(Xfull[[q]], tau)}  		# Noisy Covariate; n-by-m
  
  if (is.null(mi_low)&is.null(mi_up)){
    Xeval <- lapply(c(1:Q), function(q) Xfull[[q]])
    Weval <- lapply(c(1:Q), function(q) Wfull[[q]])
  } else{
    Xeval = lapply(c(1:Q), function(q) NA*Xfull[[q]])
    Weval = lapply(c(1:Q), function(q) NA*Wfull[[q]])
    mi = sample(mi_low:mi_up, n, replace=TRUE)
    for (i in 1:n){
      ind = sample(1:m, mi[i], replace=FALSE)
      for (q in 1:Q){
        Xeval[[q]][i, ind] = Xfull[[q]][i, ind]
        Weval[[q]][i, ind] = Wfull[[q]][i, ind]
      }
    }
  }
  
  if (length(Xfull)==1){ 
    Fanova = F(Xeval[[1]],trep,d) 
  } else{
    # Fanova = F(Xeval[[1]],Xeval[[2]],trep,d)
    Fanova = F(Xeval[[2]],trep,d)
  }
   
  # Random Error; n-by-m
  Error = Error_resp(n,t,error_type,sig2noise,rho)
  # Response based on Feval; n-by-m 
  Yanova = Fanova + Error
  
  result = list(Xfull = Xfull, Wfull = Wfull, Xeval = Xeval, Weval = Weval, Error = Error, Fanova = Fanova, Yanova = Yanova) 
}
