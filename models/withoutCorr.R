# Validate teh full model tau and phi
rm(list = ls())

library(stats4)
library(plyr)
library(ggplot2)
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
  #print(paste(t1,t2,s1,s2,logLik))
  return(-1 * logLik)
}

computeSigma = function(df){
  print(paste("Computing sigmas for", unique(df$variable)))
  # starting values
  d = data.frame(t1 = 0.4, t2 = 0.3, s1 = 0.7, s2 = 0.6)
  mlePhiTau = mle(fullLogLik, start = list(t1 = d$t1, t2 = d$t2, s1 = d$s1, s2 = d$s2), fixed = list(df = df))
  d$t1 = abs(mlePhiTau@coef[[1]])
  d$t2 = abs(mlePhiTau@coef[[2]])
  d$s1 = abs(mlePhiTau@coef[[3]])
  d$s2 = abs(mlePhiTau@coef[[4]])
  print(paste("Computation done for", unique(df$variable)))
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

# plot Tau
p = ggplot()
p = p + geom_point(data = sigmas, aes(x  = periods, y = t1), color = "red")
p = p + geom_point(data = sigmas, aes(x  = periods, y = t2), color = "blue")
p = p + geom_line(data = bobData, aes(x  = per, y = tau1), color = "red")
p = p + geom_line(data = bobData, aes(x  = per, y = tau2), color = "blue")
p = p + scale_x_log10() + scale_y_continuous(limits = c(0,1))
p + theme_bw(base_size = 28)

# plot phi
p = ggplot()
p = p + geom_point(data = sigmas, aes(x  = periods, y = s1), color = "red")
p = p + geom_point(data = sigmas, aes(x  = periods, y = s2), color = "blue")
p = p + geom_line(data = bobData, aes(x  = per, y = phi1), color = "red")
p = p + geom_line(data = bobData, aes(x  = per, y = phi2), color = "blue")
p = p + scale_x_log10() + scale_y_continuous(limits = c(0,1))
p + theme_bw(base_size = 28)

save(sigmas, file = "fullLogLik.Rdata")