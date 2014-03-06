# Plot the sigmas and taus from the homo and heteroskedastic models

rm(list = ls())
library(plyr)
library(ggplot2)

load('./../models//withoutCorr.Rdata')
withoutHetero = sigmas
load('./../models//withCorr.Rdata')
withHetero = sigmas

load('./../simpleModels//simpleWoCorr.Rdata')
withoutHomo = sigmas
load('./../simpleModels//simpleCorr.Rdata')
withHomo = sigmas

withCorr = withHetero
withCorr$tau = withHomo$tau
withCorr$phi = withHomo$phi
withCorr$type = "withCorrelation"

withoutCorr = withoutHetero
withoutCorr$tau = withoutHomo$tau
withoutCorr$phi = withoutHomo$phi
withoutCorr$type = "withoutCorrelation"

data = rbind(withCorr, withoutCorr)

heteroSigma = function(M, s1, s2){
  return(s1 + (s2 - s1)/2.25 * (min(max(M,5),7.25) - 5))
}

makePlot = function(df){
  idxWith = which(df$type == "withCorrelation")
  idxWithout = which(df$type == "withoutCorrelation")
  
  mags = seq(4,8,0.1)
  N = length(mags)
  
  tau_with_homo = rep(df$tau[idxWith], N)
  phi_with_homo = rep(df$phi[idxWith], N)
  tau_without_homo = rep(df$tau[idxWithout], N)
  phi_without_homo = rep(df$phi[idxWithout], N)
  
  tau_with_hetero = sapply(mags, heteroSigma, df$t1[idxWith], df$t2[idxWith])
  phi_with_hetero = sapply(mags, heteroSigma, df$s1[idxWith], df$s2[idxWith])
  tau_without_hetero = sapply(mags, heteroSigma, df$t1[idxWithout], df$t2[idxWithout])
  phi_without_hetero = sapply(mags, heteroSigma, df$s1[idxWithout], df$s2[idxWithout])
  
  dphi = data.frame(M = mags, phi = c(phi_with_homo, phi_without_homo, phi_with_hetero, phi_without_hetero),
                    corr = c(rep("withCorr",N), rep("withoutCorr",N), rep("withCorr",N), rep("withoutCorr",N)),
                    model = c(rep("homo", 2*N), rep("hetero", 2*N)))
  dtau = data.frame(M = mags, tau = c(tau_with_homo, tau_without_homo, tau_with_hetero, tau_without_hetero),
                    corr = c(rep("withCorr",N), rep("withoutCorr",N), rep("withCorr",N), rep("withoutCorr",N)),
                    model = c(rep("homo", 2*N), rep("hetero", 2*N)))
  
  fnamePhi = paste("./plots/phi_", per, ".jpg", sep = "")
  
  p = ggplot(dphi, aes(x = M, y = phi, color = corr, linetype = model))
  p = p + geom_line() + them_bw(base_size = 28)
  ggsave(plot = p, file = fnamePhi)
}