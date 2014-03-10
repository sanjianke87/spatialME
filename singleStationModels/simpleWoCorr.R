# Validate teh full model tau and phi
rm(list = ls())

library(stats4)
library(plyr)
library(doMC)

registerDoMC(cores = 6)
load('./../computeResids/briansResid.Rdata')

cyResids$M = cyResids$Mag

fullLogLik = function(df, tau, phi){
  logLik = 0
  eqids = unique(df$EQID)
  
  for(i in 1:length(eqids)){
    eqid = eqids[i]
    idx = which(df$EQID == eqid)
    if(length(idx) > 4){
      M = unique(df$M[idx])[1]
      N = length(idx)
      C = diag(N)*phi^2 + matrix(1, N, N) * tau^2
      detC = det(C)
      if(detC > 0 & is.finite(detC)){
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
  d = data.frame(tau = 0.4, phi = 0.6)
  mlePhiTau = mle(fullLogLik, start = list(tau = d$tau, phi = d$phi), fixed = list(df = df))
  d$tau = abs(mlePhiTau@coef[[1]])
  d$phi = abs(mlePhiTau@coef[[2]])
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

save(sigmas, file = "simpleWoCorr.Rdata")