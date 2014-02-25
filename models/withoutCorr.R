# Validate teh full model tau and phi
rm(list = ls())

library(stats4)
library(plyr)
library(doMC)

registerDoMC(cores = 5)
#load('./../computeResids/fullResid.Rdata')
load('./../computeResids/briansResid.Rdata')
bobData = read.csv('./../data/tauPhiBob.csv')

cyResids$M = cyResids$Mag

fullLogLik = function(df, t1, t2, s1, s2){
  logLik = 0
  eqids = unique(df$EQID)
  
  for(i in 1:length(eqids)){
    eqid = eqids[i]
    idx = which(df$EQID == eqid)
    if(length(idx) > 4){
      M = unique(df$M[idx])[1]
      tau = t1 + (t2 - t1)/1.5 * (min(max(M,5),6.5) - 5)
      phi = s1 + (s2 - s1)/1.5 * (min(max(M,5),6.5) - 5)
      N = length(idx)
      C = diag(N)*phi^2 + matrix(1, N, N) * tau^2
      
      detC = det(C)
      if(detC > 0 ){
        Cinv = chol2inv(chol(C))
        logLik = logLik - 0.5*log(detC) - 0.5 * t(df$resids[idx]) %*% Cinv %*% df$resids[idx]
      }else{
        logLik = logLik - 10000
      }
    }
  }
  return(-1 * logLik)
}

computeSigma = function(df){
  # starting values
  d = data.frame(t1 = 0.4, t2 = 0.3, s1 = 0.7, s2 = 0.6)
  mlePhiTau = mle(fullLogLik, start = list(t1 = d$t1, t2 = d$t2, s1 = d$s1, s2 = d$s2), fixed = list(df = df))
  d$t1 = abs(mlePhiTau@coef[[1]])
  d$t2 = abs(mlePhiTau@coef[[2]])
  d$s1 = abs(mlePhiTau@coef[[3]])
  d$s2 = abs(mlePhiTau@coef[[4]])
  return(d)
}

compFor = c("T0.010S","T0.200S","T0.500S","T1.000S","T2.000S")
data = subset(cyResids, variable %in% compFor)
sigmas = ddply(data, "variable", computeSigma, .parallel = TRUE)

extractPeriod = function(per){
  per = sub("T","",per)
  return(as.numeric(sub("S","",per)))
}

sigmas$periods = sapply(sigmas$variable, extractPeriod)

save(sigmas, file = "withoutCorr.Rdata")