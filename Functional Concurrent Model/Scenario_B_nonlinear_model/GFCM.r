######################################################
# GFCM.r - May/27/2016				#
# Script created by Janet Kim              		#
# Contents:                                   			#
# (i) G                                       			#
# (ii) trf                                        			#
# (iii) GFCM                                   			#
# (iv) Std.est                               			#
# (v) CP						#
# @ Sample size in training data:			#
#   n subjects, m observation points per subject		#
# @ Sample size in test data:				#
#   n.test subjects, m observation points per subject	#
######################################################

# Estimate covariance of random error 
G = function(r, n, m, pve, iind.na=NULL){# "r" is residuals obtained from gam fit; dim(r) = n*m-by-1
  res.mat = matrix(r, ncol=m, byrow=T)
  res.sm = fpca.sc(res.mat, pve=pve, var=TRUE)
  eval = res.sm$evalues/m			# eigenvalues
  efunc = (res.sm$efunctions)*sqrt(m)		# eigenfunctions
  sig2 = res.sm$sigma2			# variance of white noise error
  
  if (!any(is.na(r))){ 
    if (length(eval) == 1){ G.tilde = eval*efunc%*%t(efunc)+sig2*diag(m) }
    if (length(eval) > 1) { G.tilde = efunc%*%diag(eval)%*%t(efunc)+sig2*diag(m) }
  } else{
    G.tilde = NULL
    for (i in 1:n){
      efunc0 = efunc[-iind.na[[i]],]
      mi.len = nrow(as.matrix(efunc0))
      if (length(eval) == 1){ G.tilde[[i]] = eval*efunc0%*%t(efunc0)+sig2*diag(mi.len) }
      if (length(eval) > 1) { G.tilde[[i]] = efunc0%*%diag(eval)%*%t(efunc0)+sig2*diag(mi.len) }
    }
  }
  G.tilde
}

trf = function(w, wtest=NULL ,pve=pve){# Smooth noisy covariate w by applying FPCA
  w.sm = fpca.sc(w, pve=pve, var=TRUE)
  xhat = w.sm$Yhat                                                                # n-by-m
  m=ncol(xhat)
  w.eval = w.sm$evalues/m
  w.efunc = (w.sm$efunctions)*sqrt(m)
  
  # Centering/Scaling transformation of xhat 
  mu.hat = w.sm$mu                                                                # m-by-1
  if (length(w.eval) == 1){ sig2.hat = diag(w.eval*w.efunc%*%t(w.efunc)) }
  if (length(w.eval) > 1){ sig2.hat = diag(w.efunc%*%diag(w.eval)%*%t(w.efunc)) }
  xhat.tr = scale(xhat, center = mu.hat, scale = sqrt(sig2.hat))  # n-by-m
  if (is.null(wtest)){xtest.tr = NULL}
  if (!is.null(wtest)){
    xtest = fpca.sc(Y=w, Y.pred = wtest, pve=pve, var=TRUE)$Yhat
    xtest.tr = scale(xtest, center = mu.hat, scale = sqrt(sig2.hat))
  }
  
  result = list(xhat.tr = xhat.tr, xtest.tr = xtest.tr)
}

# Fit GFCM or linear FCM
GFCM = function(y, w, nbasis=list(7, 7), pve=0.99, fit_opt=list("nonlinear"), opt="gam", newdata){
  # dim(y) = dim(w) = n-by-m; dim(t) = m-by-1
  # "nbasis": a list of number of basis functions used in the analysis
  #  e.g., if Q=2 and fit_opt=list("nonlienar", "nonlinear"), nbasis=list(c(7,7), c(7,7)).
  #          if Q=2 and fit_opt=list("nonlinear", "linear"), nbasis=list(c(7,7), 7).
  # xtest = newdata$X, ytest = newdata$Y  
  
  n = nrow(y); m = ncol(y); Q = length(w)
  t = seq(0, 1, len=m)					# m-by-1
  trep = matrix(rep(t, n), nrow=n, byrow=TRUE)			# n*m-by-1
  xtest = newdata$X					# a list of covariates in test data; n.test-by-m
  ytest = newdata$Y					# Response in test data; n.test-by-m
  
  # Smooth and transform the noisy covariate
  w.sm = lapply(c(1:Q), function(q) trf(w[[q]], pve=pve, xtest[[q]]))
  xhat.tr = lapply(c(1:Q), function(q) w.sm[[q]]$xhat.tr) # n-by-m
  
  # Vectorize the y's and xhat.tr's
  iind.na = NULL
  if (!any(is.na(y))){
    xhat.tr.vec = lapply(c(1:Q), function(q) as.vector(t(xhat.tr[[q]])))	# n*m-by-1
    t.vec = as.vector(t(trep))					# n*m-by-1
    y.vec = as.vector(t(y))					# n*m-by-1    
  } else{
    # Save index of "NA" for each subject
    iind.na = lapply(c(1:n), function(k) na.action(na.omit(y[k,])))
    # Save index of "NA" for all subjects
    ind.na = which(is.na(as.vector(t(y))))			# N=n*m-length(ind.na)
    xhat.tr.vec = lapply(c(1:Q), function(q) as.vector(t(xhat.tr[[q]]))[-ind.na])# N-by-1
    t.vec = as.vector(t(trep))[-ind.na]				# N-by-1
    y.vec = as.vector(t(y))[-ind.na]				# N-by-1
  }
  
  x1 = xhat.tr.vec[[1]]                     				# N-by-1
  formula = as.list(NULL)
  if (fit_opt[[1]] == "nonlinear"){# Fit GFCM
    formula[[1]] = paste("y.vec~","te(x1, t.vec, bs=\"ps\", k=nbasis[[1]])")
  }
  if (fit_opt[[1]] == "linear"){# Fit linear FCM
    formula[[1]] = paste("y.vec~","s(t.vec, bs=\"ps\", k=nbasis[[1]])+s(t.vec, by=x1, bs=\"ps\", k=nbasis[[1]])")
  }
  if (Q == 1){ formula.fit = as.formula(formula[[1]]) }
  if (Q == 2){
    x2 = xhat.tr.vec[[2]]
    if (fit_opt[[2]] == "nonlinear"){ formula[[2]] = paste("+te(x2, t.vec, bs=\"ps\", k=nbasis[[2]])") }
    if (fit_opt[[2]] == "linear"){ formula[[2]] = paste("+s(t.vec, by=x2, bs=\"ps\", k=nbasis[[2]])") }
    formula.fit = as.formula(paste(formula[[1]],formula[[2]]))
  }
  if (opt=="gam"){ fit = gam(formula.fit, method="REML") }
  if (opt=="bam"){ fit = bam(formula.fit, method="REML") }
  
  # Obtain fitted values
  Fstar.fit = fit$fitted.values
  
  # Do prediction using the transformed test data
  xtest.tr.vec = lapply(1:Q, function(q) as.vector(t(w.sm[[q]]$xtest.tr)))
  if (Q == 1){newd = list(x1 = xtest.tr.vec[[1]], t.vec=rep(t, nrow(ytest)))}
  if (Q == 2){newd = list(x1 = xtest.tr.vec[[1]], x2 = xtest.tr.vec[[2]], t.vec=rep(t, nrow(ytest)))}
  Fstar.pred = predict.gam(fit, newdata=newd)			# n.test*m-by-1
  
  # Obtain residuals
  Res = fit$residuals					# Residuals from training data
  Resp = as.vector(t(ytest))-Fstar.pred				# Residuals from test data
  
  # Estimate covariance of random errors
  if (!any(is.na(y))){ 
    r <- Res 
  } else{
    r = rep(NA, n*m)
    ind = (1:(n*m))[-ind.na]					# Index non-NA elements
    r[ind] <- Res
  }
  
  # Estimate RMSPE (prediction error)
  RMSPE_in = sqrt(mean(Res^2))
  RMSPE_out = sqrt(mean(Resp^2))
  
  # Estimate var(ytest) and var(Fstar.pred)
  if (!any(is.na(y))){
    CovRes = G(r, n, m, pve=pve)                  			# m-by-m
    VarFp = Std.est(fit, newd, CovRes, CovResp=NULL, n, n.test=nrow(ytest), m)
  } else{
    CovRes = G(r, n, m, pve=pve, iind.na=iind.na)			# mi-by-mi
    CovResp = G(Res, n, m, pve=pve)         			# m-by-m
    VarFp = Std.est(fit, newd, CovRes, CovResp=CovResp, n, n.test=nrow(ytest), m, iind.na=iind.na)
  }
  
  # Estimate ICP for new observation
  ICP85 = CP(ytest, Fstar.pred, VarFp, 0.15);  ICP90 = CP(ytest, Fstar.pred, VarFp, 0.1);  ICP95 = CP(ytest, Fstar.pred, VarFp, 0.05)
  ICPall = c(ICP85$ICP, ICP90$ICP, ICP95$ICP)
  
  # Integrated length (IL) and minimum and maximum standard errors of prediction bands
  IL= c(mean(rowMeans(2*ICP85$MOE)), mean(rowMeans(2*ICP90$MOE)), mean(rowMeans(2*ICP95$MOE)))
  minSE = c(mean(apply(2*ICP85$MOE,1,min)), mean(apply(2*ICP90$MOE,1,min)), mean(apply(2*ICP95$MOE,1,min)))
  maxSE = c(mean(apply(2*ICP85$MOE,1,max)), mean(apply(2*ICP90$MOE,1,max)), mean(apply(2*ICP95$MOE,1,max)))
  
  # range of SE; refer to main paper
  R85 = c(minSE[1], maxSE[1]); R90 = c(minSE[2], maxSE[2]); R95 = c(minSE[3], maxSE[3])
  R = rbind(R85, R90, R95)
  
  result = list(Fstar.fit = Fstar.fit, Fstar.pred = Fstar.pred, RMSPE_in = RMSPE_in, RMSPE_out = RMSPE_out,
                CovRes = CovRes, VarFp = VarFp, ICPall = ICPall, IL = IL, R = R)  
}

# Estimate variance of parameters, variance of predicted response
Std.est = function(fit, newd, CovRes, CovResp=NULL, n, n.test, m, iind.na=NULL){
  Z = model.matrix(fit)
  Xp = predict.gam(fit, newdata=newd, type="lpmatrix")	# nrow(Xp) = ntest*m
  H = vcov(fit, dispersion=1)
  
  sandwich = 0
  if (length(fit$residuals) == n*m){
    for (i in 1:n){ sandwich = sandwich + t(Z[(m*(i-1)+1):(m*i),])%*%CovRes%*%Z[(m*(i-1)+1):(m*i),] }
    G.diag = rep(diag(CovRes), n.test)			# n.test*m-by-1
  } else{
    for (i in 1:n){
      mi.len = m -length(iind.na[[i]])
      sandwich = sandwich + t(Z[(1:mi.len),])%*%CovRes[[i]]%*%Z[(1:mi.len),]
      Z = Z[-(1:mi.len),]
    }
    G.diag = rep(diag(CovResp), n.test)			# n.test*m-by-1
  }
  CovTheta = H%*%sandwich%*%t(H)			# Est. of Cov(Theta); r-by-r
  VarFp = G.diag + diag(Xp%*%CovTheta%*%t(Xp))
  return(VarFp)
}

CP = function(ytest, Fstar.pred, Var.est, alpha){
  # dim(ytest) = n.test-by-m; dim(Fstar.pred) = n.test*m-by-1
  # VarFp: varaince for prediction; dim(VarFp) = n.test*m-by-1
  pred = matrix(Fstar.pred, nrow = nrow(ytest), ncol=ncol(ytest), byrow=TRUE)
  StdF.mat = matrix(sqrt(Var.est), nrow = nrow(ytest), ncol=ncol(ytest), byrow=TRUE)
  MOE = qnorm(1-alpha/2)*StdF.mat		# margin of error at a prescribed nominal level; n.test-by-m
  F.llim = pred-MOE                     		# n.test-by-m
  F.ulim = pred+MOE			# n.test-by-m
  Logic = (F.llim<ytest) & (ytest<F.ulim)	# n.test-by-m
  ICP = mean(rowMeans(Logic))		# Integrated coverage prob; n.test-by-m
  result = list(MOE = MOE, ICP = ICP)
}