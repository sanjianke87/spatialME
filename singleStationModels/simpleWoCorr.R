# Validate teh full model tau and phi
rm(list = ls())

library(stats4)
library(plyr)
library(doMC)

registerDoMC(cores = 6)
load('./../computeResids/briansResid.Rdata')

cyResids$M = cyResids$Mag

fullLogLik = function(df, tau, phiSS, phiS2S, mat1, mat2, mat3){
  logLik = 0  
  C = mat1*phiSS^2 + mat2*tau^2 + mat3 * phiS2S^2
  detC = det(C)
  if(detC > 0 & is.finite(detC)){
    Cinv = chol2inv(chol(C))
    logLik = logLik - 0.5*log(detC) - 0.5 * t(df$resids) %*% Cinv %*% df$resids
  }else{
    logLik = logLik - 10000
  }
  #print(paste(tau, phiSS, phiS2S, logLik))
  return(-1 * logLik)
}

getMat1 = function(N){
  i = seq(1,N)
  j = seq(1,N)
  mat1 = sparseMatrix(i = i, j =j, x = 1, dims = c(N,N))
  return(mat1)
}

getMat2 = function(df){
  N = length(df[,1])
  eqids = unique(df$EQID)
  counter = 0
  Is = c()
  Js = c()
  for(i in 1:length(eqids)){
    eqid = eqids[i]
    idx = which(df$EQID == eqid)
    idxStart = counter + 1
    idxEnd = counter + length(idx)
    for(j in idxStart:idxEnd){
      Is = c(Is, rep(j , length(idx)))
      Js = c(Js, idxStart:idxEnd)
    }
    counter = counter + length(idx)
  }
  mat2 = sparseMatrix(i = Is, j = Js, x = 1, dims = c(N,N))
}

getLL = function(df){
  return(paste(df$lat,"_",df$lon, sep = ""))
}

getMat3 = function(df){
  N = length(df[,1])
  ll = daply(df, c(), getLL)
  Is = c()
  Js = c()
  for(i in 1:length(ll)){
    idx = which(ll == ll[i])
    Is = c(Is, rep(i, length(idx)))
    Js = c(Js, idx)
  }
  mat3 = sparseMatrix(i = Is, j = Js, x = 1, dims = c(N,N))
}

computeSigma = function(df){
  # starting values
  d = data.frame(tau = 0.4, phiSS = 0.4, phiS2S = 0.4)
  
  print("Preparing mat1")
  mat1 = getMat1(length(df[,1]))
  print("Preparing mat2")
  mat2 = getMat2(df)
  print("Preparing mat3")
  mat3 = getMat3(df)
  
  per = unique(df$variable)[1]
  print(paste("Starting MLE for",per))
  mlePhiTau = mle(fullLogLik, start = list(tau = d$tau, phiSS = d$phiSS, phiS2S = d$phiS2S), 
                  fixed = list(df = df, mat1 = mat1, mat2 = mat2, mat3 = mat3))
  d$tau = abs(mlePhiTau@coef[[1]])
  d$phiSS = abs(mlePhiTau@coef[[2]])
  d$phiS2S = abs(mlePhiTau@coef[[2]])
  return(d)
}

# Remove the Tottori event from computation
cyResids = subset(cyResids, EQID != 176)

# Only use the selected periods

#compFor = c("T0.010S", "T0.020S", "T0.030S", "T0.040S", "T0.050S", "T0.075S",
#            "T0.100S", "T0.120S", "T0.150S", "T0.170S", "T0.200S", "T0.250S",
#            "T0.300S", "T0.400S", "T0.500S", "T0.750S", "T1.000S", "T1.500S",
#            "T2.000S", "T3.000S", "T4.000S", "T5.000S", "T7.500S", "T10.000S")

compFor = c("T0.010S")

data = subset(cyResids, variable %in% compFor)

# Compute the phi and taus
system.time(sigmas <- ddply(data, "variable", computeSigma))#, .parallel = TRUE)

# Add numeric periods to the dataframe
extractPeriod = function(per){
  per = sub("T","",per)
  return(as.numeric(sub("S","",per)))
}
sigmas$periods = sapply(sigmas$variable, extractPeriod)

save(sigmas, file = "simpleWoCorr.Rdata")