#########################################################
# test.r - May/20/2016					#
# Script created by Janet Kim				#
# Contents: Hypothesis testing				#
# (i)  test.anova					#
# (ii) rejectH0prob					#
# @ Sample size:					#
#   n subjects, m observation points/subject		#
# @ Number of bootstrap = B = 200			#
#########################################################

test.anova = function(y, w, test.type, nbasis.null, nbasis.full, pve=pve, B, seed.B, l){
  # test.type = 1: Algorithm 1 in the main paper
  # test.type = 2: Algorithm 2 in the main paper
  # "seed.B": B-by-1 vector of seed numbers used to generate bootstrap data
  source("GFCM.r")
  n = nrow(y);  m = ncol(y); Q = length(w)
  t = seq(0, 1, len=m)
  trep = matrix(rep(t, n), nrow=n, byrow=TRUE)

  if(m<35){
    y.sm = fpca.face(y, pve=pve, knots=floor(m/2))
  }else{
    y.sm = fpca.face(y, pve=pve)
  }
  ytilde = y - matrix(y.sm$mu, nrow=n, ncol=m, byrow=TRUE)
  
  # Smooth and transform the noisy covariate:
  xhat.tr = lapply(c(1:Q), function(q) trf(w[[q]], pve=pve)$xhat.tr) # n-by-m
  
  # Vectorize the y's and xhat.tr's
  if (!any(is.na(y))){# dense design
    xhat.tr.vec = lapply(c(1:Q), function(q) as.vector(t(xhat.tr[[q]])))	# n*m-by-1
    t.vec = as.vector(t(trep))							      # n*m-by-1
    y.vec = as.vector(t(ytilde))							          # n*m-by-1    
  } else{
    # Save index of "NA" for all subjects
    ind.na = which(is.na(as.vector(t(y))))				# N=n*m-length(ind.na)
    xhat.tr.vec = lapply(c(1:Q), function(q) as.vector(t(xhat.tr[[q]]))[-ind.na])# N-by-1
    t.vec = as.vector(t(trep))[-ind.na]						# N-by-1
    y.vec = as.vector(t(ytilde))[-ind.na]						  # N-by-1
  }
  
  # Fit the null and the full model
  if (test.type == 1){
    x1 = xhat.tr.vec[[1]]
    fit.null = gam(y.vec ~ s(t.vec, bs="ps", k=nbasis.null), method="REML")
    fit.full = gam(y.vec ~ te(x1, t.vec, bs="ps", k=nbasis.full), method="REML")
  }
  if (test.type == 2){
    x1 = xhat.tr.vec[[1]]; x2 = xhat.tr.vec[[2]]
    fit.null = gam(y.vec ~ te(x1, t.vec, bs="ps", k=nbasis.full), method="REML")
    fit.full = gam(y.vec ~ te(x1, t.vec, bs="ps", k=nbasis.full) + te(x2, t.vec, bs="ps", k=nbasis.full), 
                   method="REML")
  }
  
  # Anova test
  test = anova.gam(fit.null, fit.full, freq=TRUE, test="F", p.type=0)
  if (test$Deviance[2] <= 0){Fst = 0}
  if (test$Deviance[2] > 0){Fst = test$F[2]}
  
  # Residuals from the full model fit
  if (!any(is.na(y))){# dense design 
    Res = fit.full$residuals 
  } else{# sparse design
    Res = rep(NA, n*m); Res[(1:(n*m))[-ind.na]] <- fit.full$residuals
  }
  
  # Bootstrap:
  FstB  = NULL; count.zero = 0
  for (b in 1:B){ 
    set.seed(seed.B[b]*100) 
    ind.B = sample(c(1:n), n, replace = TRUE)				            # Resample subject index
    ResB = (matrix(Res, nrow=n, ncol=m, byrow=TRUE))[ind.B,]		# Boobstrap the residuals

        
    # Reconstruct the bootstrap data:
    if (!any(is.na(y))){# dense design 
      tB <- t.vec							                            # n*m-by-1
      xhat.tr.vecB <- xhat.tr.vec
      yB <- fit.null$fitted.values + as.vector(t(ResB))   # n*m-by-1
    } else{# sparse design
      ind.naB = which(is.na(as.vector(t(ResB))))          # index of "NA" in bootstrap data
      tB <- as.vector(t(trep))[-ind.naB]
      xhat.tr.vecB <- lapply(c(1:Q), function(q) as.vector(t(xhat.tr[[q]]))[-ind.naB])
      if (test.type == 1){
        yB <- predict(fit.null, newdata=list(t.vec = tB)) + as.vector(t(ResB))[-ind.naB]
      }
      if (test.type ==2){
        yB <- predict(fit.null, newdata=list(t.vec = tB, x1 = xhat.tr.vecB[[1]])) + as.vector(t(ResB))[-ind.naB]
      }
    }
    # fit the null and the full model using the bootstrap data
    if (test.type == 1){
      x1B <- xhat.tr.vecB[[1]]
      fit.nullB = gam(yB ~ s(tB, bs="ps", k=nbasis.null), method="REML")
      fit.fullB = gam(yB ~ te(x1B, tB, bs="ps", k=nbasis.full), method="REML")
    }
    if (test.type == 2){
      x1B <- xhat.tr.vecB[[1]]; x2B <- xhat.tr.vecB[[2]]
      fit.nullB = gam(yB ~ te(x1B, tB, bs="ps", k=nbasis.null), method="REML")
      fit.fullB = gam(yB ~ te(x1B, tB, bs="ps", k=nbasis.full) + te(x2B, tB, bs="ps", k=nbasis.full), method="REML")
    }
    
    # Anova test
    testB = anova.gam(fit.nullB, fit.fullB, freq=TRUE, test="F", p.type=0)
    if (testB$Deviance[2] <= 0){ FstB = c(FstB, 0); count.zero = count.zero+1 }
    if (testB$Deviance[2] > 0){ FstB = c(FstB, testB$F[2]) }
    cat("Bootstrap iter =", b, "", "FstB =", testB$F[2], "\n")
  } 
  
  # Bootstrap p-value
  pvalue = mean(na.omit(Fst < FstB))
  cat("Iter =", l, "", "pvalue =", pvalue, "", "count.zero =", count.zero, "", "\n")
  
  Result = list(Fst = Fst, FstB = FstB, pvalue = pvalue, count.zero = count.zero)
}


# Calculate probability of rejecting the null hypothesis
rejectH0prob = function(p, alpha){ 
  mean(p < alpha) 
}
