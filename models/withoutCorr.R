# Validate teh full model tau and phi
rm(list = ls())

library(stats4)
library(plyr)
library(doMC)

registerDoMC(cores = 6)
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
      tau = t1 + (t2 - t1)/2.25 * (min(max(M,5),7.25) - 5)
      phi = s1 + (s2 - s1)/2.25 * (min(max(M,5),7.25) - 5)
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
  mlePhiTau = mle(fullLogLik, start = list(t1 = d$t1, t2 = d$t2, s1 = d$s1, s2 = d$s2), 
                  fixed = list(df = df))
  d$t1 = abs(mlePhiTau@coef[[1]])
  d$t2 = abs(mlePhiTau@coef[[2]])
  d$s1 = abs(mlePhiTau@coef[[3]])
  d$s2 = abs(mlePhiTau@coef[[4]])
  return(d)
}

# Remove the Tottori event from computation
cyResids = subset(cyResids, EQID != 176)

# Only use the selected periods
compFor = c("T0.010S", "T0.020S", "T0.030S", "T0.040S", "T0.050S", "T0.075S",
            "T0.100S", "T0.120S", "T0.150S", "T0.170S", "T0.200S", "T0.250S",
            "T0.300S", "T0.400S", "T0.500S", "T0.750S", "T1.000S", "T1.500S",
            "T2.000S", "T3.000S", "T4.000S", "T5.000S", "T7.500S", "T10.000S")
data = subset(cyResids, variable %in% compFor)

# Compute the phi and taus
sigmas = ddply(data, "variable", computeSigma, .parallel = TRUE)

# Add numeric periods to the dataframe
extractPeriod = function(per){
  per = sub("T","",per)
  return(as.numeric(sub("S","",per)))
}
sigmas$periods = sapply(sigmas$variable, extractPeriod)

save(sigmas, file = "withoutCorr.Rdata")