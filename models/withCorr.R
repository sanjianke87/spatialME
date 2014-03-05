# Validate the full model tau and phi
rm(list = ls())

library(stats4)
library(plyr)
library(doMC)

registerDoMC(cores = 5)
load('./../computeResids/briansResid.Rdata')
bobData = read.csv('./../data/tauPhiBob.csv')

load("./withoutCorr.Rdata")

sigmas_old = sigmas

cyResids$M = cyResids$Mag

negLogLik = function(t1, t2, s1, s2, range, df, distMat){
  logLik = 0
  t1 = abs(t1)
  t2 = abs(t2)
  s1 = abs(s1)
  s2 = abs(s2)
  range = abs(range)
  
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
      #print(eqid)
      C = corrMat*phi^2 + matrix(1, N, N) * tau^2 + matrix(runif(N*N, 1e-6, 5e-6),N,N)
      detC = det(C)
      if(detC > 0 & is.finite(detC)){
        update = tryCatch({
          Cinv = chol2inv(chol(C))
          0.5*log(detC) + 0.5 * t(df$resids[idx]) %*% Cinv %*% df$resids[idx]
        }, error = function(e){
          #print("HERE")
          return(10000)
        })
        logLik = logLik - update
      }else{
        logLik = logLik - 10000
      }
    }
  }
  #print(paste(t1,t2,s1,s2,range,logLik))
  return(-1 * logLik)
}

computeSigma = function(df, distMats){
  per = unique(df$variable)[1]
  period = sub("T","",per)
  period = as.numeric(sub("S","",period))
  startRange = 1
  if(period < 1){
    startRange = 8.5 + 17.2 * period  
  }else{
    startRange = 22.0 + 3.7 * period
  }
  idx = which(sigmas_old$periods == period)
  print(per)
  print(names(distMats))
  distMat = distMats[[as.character(per)]]
  # starting values
  d = data.frame(t1 = sigmas_old$t1[idx], t2 = sigmas_old$t2[idx], s1 = sigmas_old$s1[idx], s2 = sigmas_old$s2[idx], range = startRange)
  mlePhiTau = mle(negLogLik, start = list(t1 = d$t1, t2 = d$t2, s1 = d$s1, s2 = d$s2), lower = c(0.05,0.05,0.05,0.05), upper = c(0.6,0.6,1.4,1.4),
                  method = "L-BFGS-B",fixed = list(df = df, distMat = distMat, range = startRange))
  d$t1 = abs(mlePhiTau@coef[[1]])
  d$t2 = abs(mlePhiTau@coef[[2]])
  d$s1 = abs(mlePhiTau@coef[[3]])
  d$s2 = abs(mlePhiTau@coef[[4]])
  #d$range = abs(mlePhiTau@coef[[5]])
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
#compFor = c("T0.500S")
#compFor = c("T0.010S", "T0.020S", "T0.030S", "T0.040S", "T0.050S", "T0.075S",
            #"T0.100S", "T0.120S", "T0.150S", "T0.170S", "T0.200S", "T0.250S",
            #"T0.300S", "T0.400S", "T0.500S", "T0.750S", "T1.000S", "T1.500S",
            #"T2.000S", "T3.000S", "T4.000S", "T5.000S", "T7.500S", "T10.000S")

compFor = c("T0.010S", "T0.200S", "T0.500S", "T1.000S", "T2.000S")


data = subset(cyResids, variable %in% compFor)

print("Step 2: Preparing Distance Matrix")

#latLonData = read.csv("./../data/NGAW2_latLon.csv")
#distanceMat = dlply(data, "variable", computeDistanceMat, latLonData, .parallel = TRUE)
#save(distanceMat, file = "distanceMat.Rdata")
load("distanceMat.Rdata")
print("Step 3: Maximum Likelihood")

# Compute the phi and taus
sigmas = list()
for(i in 1:100){
  sigmas[[i]] = ddply(.data = data, .variables = c("variable"), .fun = computeSigma, distanceMat)#, .parallel = TRUE)
}

print(sigmas)
# Add numeric periods to the dataframe
#extractPeriod = function(per){
#  per = sub("T","",per)
#  return(as.numeric(sub("S","",per)))
#}
#sigmas$periods = sapply(sigmas$variable, extractPeriod)

save(sigmas, file = "withCorr50.Rdata")
