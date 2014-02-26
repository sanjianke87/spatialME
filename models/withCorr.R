# Validate the full model tau and phi
rm(list = ls())

library(stats4)
library(plyr)
library(doMC)

registerDoMC(cores = 6)
load('./../computeResids/briansResid.Rdata')
bobData = read.csv('./../data/tauPhiBob.csv')

cyResids$M = cyResids$Mag

negLogLik = function(df, distMat, t1, t2, s1, s2, range){
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
      corrMat = exp(-3*distMat[[eqid]]/range)
      C = corrMat*phi^2 + matrix(1, N, N) * tau^2
      
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

computeSigma = function(df, distMat){
  per = unique(df$variable)[1]
  distMat = distMat[[per]]
  period = sub("T","",per)
  period = as.numeric(sub("S","",period))
  startRange = 1
  if(period < 1){
    startRange = 40.7 - 15.0 * period  
  }else{
    startRange = 22.0 + 3.7*period
  }
  # starting values
  d = data.frame(t1 = 0.4, t2 = 0.3, s1 = 0.7, s2 = 0.6, range = startRange)
  mlePhiTau = mle(negLogLik, start = list(t1 = d$t1, t2 = d$t2, s1 = d$s1, s2 = d$s2, range = d$range), 
                  fixed = list(df = df, distMat = distMat))
  d$t1 = abs(mlePhiTau@coef[[1]])
  d$t2 = abs(mlePhiTau@coef[[2]])
  d$s1 = abs(mlePhiTau@coef[[3]])
  d$s2 = abs(mlePhiTau@coef[[4]])
  d$range = abs(mlePhiTau@coef[[5]])
  return(d)
}

getDistance = function(latLon1, latLon2){
  degreesToRadian = pi/180.0
  phi1 = (90 - latLon1$lat) * degreesToRadian
  phi2 = (90 - latLon2$lat) * degreesToRadian
  theta1 = latLon1$lon*degreesToRadian
  theta2 = latLon2$lon*degreesToRadian
  cosValue = sin(phi1) * sin(phi2) * cos(theta1 - theta2) + cos(phi1)*cos(phi2)
  arc = acos(cosValue)
  dist = 6373*arc
  return(dist)
}

computeDistanceMat = function(df, latLonData){
  eqids = unique(df$EQID)
  distMats = list()
  for(eqid in eqids){
    idx = which(df$EQID == eqid)
    distMat = matrix(0,length(idx),length(idx))
    for(i in 1:length(idx)){
      for(j in 1:i){
        if(i == j){
          distMat[i,j] = 0
        }else{
          idxI = which(latLonData$seq == df$SeqNo[idx[i]])
          idxJ = which(latLonData$seq == df$SeqNo[idx[j]])
          distMat[i,j] = getDistance(latLonData[idxI,], latLonData[idxJ,])
          distMat[j,i] = distMat[i,j]
          if(is.nan(distMat[i,j])){
            print(paste("NAN",latLonData[idxI,],latLonData[idxJ,]))
          }
        }
      }
    }
    distMats[[eqid]] = distMat
  }
  return(distMats)
}

print("Step 1: Preparing Data")

# Remove the Tottori event from computation
cyResids = subset(cyResids, EQID != 176)

# Only use the selected periods
compFor = c("T0.010S")
data = subset(cyResids, variable %in% compFor)

print("Step 2: Preparing Distance Matrix")

latLonData = read.csv("./../data/NGAW2_latLon.csv")
distanceMat = dlply(data, "variable", computeDistanceMat, latLonData, .parallel = TRUE)

print("Step 3: Maximum Likelihood")

# Compute the phi and taus
sigmas = ddply(.data = data, .variables = c("variable"), .fun = computeSigma, distanceMat, .parallel = TRUE)

# Add numeric periods to the dataframe
extractPeriod = function(per){
  per = sub("T","",per)
  return(as.numeric(sub("S","",per)))
}
sigmas$periods = sapply(sigmas$variable, extractPeriod)

save(sigmas, file = "withCorr.Rdata")