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

withCorr = withCorr[,-6]

data = rbind(withCorr, withoutCorr)

heteroSigma = function(M, s1, s2){
  return(s1 + (s2 - s1)/2.25 * (min(max(M,5),7.25) - 5))
}

makePlot = function(df){
  idxWith = which(df$type == "withCorrelation")
  idxWithout = which(df$type == "withoutCorrelation")
  
  per = unique(df$variable)[1]
  
  mags = seq(4,8,0.25)
  N = length(mags)
  
  tau_with_homo = rep(df$tau[idxWith], N)
  phi_with_homo = rep(df$phi[idxWith], N)
  tau_without_homo = rep(df$tau[idxWithout], N)
  phi_without_homo = rep(df$phi[idxWithout], N)
  
  tau_with_hetero = sapply(mags, heteroSigma, df$t1[idxWith], df$t2[idxWith])
  phi_with_hetero = sapply(mags, heteroSigma, df$s1[idxWith], df$s2[idxWith])
  tau_without_hetero = sapply(mags, heteroSigma, df$t1[idxWithout], df$t2[idxWithout])
  phi_without_hetero = sapply(mags, heteroSigma, df$s1[idxWithout], df$s2[idxWithout])
  
  tau_homo = tau_with_homo/tau_without_homo
  phi_homo = phi_with_homo/phi_without_homo
  sig_homo = sqrt(tau_with_homo^2 + phi_with_homo^2)/sqrt(tau_without_homo^2 + phi_without_homo^2)
  
  tau_hetero = tau_with_hetero/tau_without_hetero
  phi_hetero = phi_with_hetero/phi_without_hetero
  sig_hetero = sqrt(tau_with_hetero^2 + phi_with_hetero^2)/sqrt(tau_without_hetero^2 + phi_without_hetero^2)
  
  dphi = data.frame(M = mags, phi = c(phi_homo, phi_hetero),
                    model = c(rep("homoSk", N), rep("heteroSk", N)))
  
  dtau = data.frame(M = mags, tau = c(tau_homo, tau_hetero),
                    model = c(rep("homoSk", N), rep("heteroSk", N)))
  
  dsig = data.frame(M = mags, sig = c(sig_homo, sig_hetero),
                    model = c(rep("homoSk", N), rep("heteroSk", N)))
  
  fnamePhi = paste("./ratioPhi/phi_", per, ".jpg", sep = "")
  p = ggplot(dphi, aes(x = M, y = phi, linetype = model))
  p = p + geom_line() + theme_bw(base_size = 28) + ylab(expression(phi)) + scale_y_continuous(lim = c(0.5,1.5))
  ggsave(plot = p, file = fnamePhi)
  
  fnameTau = paste("./ratioTau/tau_", per, ".jpg", sep = "")
  p = ggplot(dtau, aes(x = M, y = tau, linetype = model))
  p = p + geom_line() + theme_bw(base_size = 28) + ylab(expression(tau))+ scale_y_continuous(lim = c(0.5,1.5))
  ggsave(plot = p, file = fnameTau)
  
  fnameSig = paste("./ratioSig/sig_", per, ".jpg", sep = "")
  p = ggplot(dsig, aes(x = M, y = sig, linetype = model))
  p = p + geom_line() + theme_bw(base_size = 28) + ylab(expression(sigma))+ scale_y_continuous(lim = c(0.5,1.5))
  ggsave(plot = p, file = fnameSig)
}

d_ply(data, "variable", makePlot)