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
withCorr$p1 = withCorr$s1
withCorr$p2 = withCorr$s2
withCorr$s1 = sqrt(withCorr$t1^2 + withCorr$p1^2)
withCorr$s2 = sqrt(withCorr$t2^2 + withCorr$p2^2)
withCorr$sig = sqrt(withCorr$tau^2 + withCorr$phi^2)

withoutCorr = withoutHetero
withoutCorr$tau = withoutHomo$tau
withoutCorr$phi = withoutHomo$phi
withoutCorr$p1 = withoutCorr$s1
withoutCorr$p2 = withoutCorr$s2
withoutCorr$s1 = sqrt(withoutCorr$t1^2 + withoutCorr$p1^2)
withoutCorr$s2 = sqrt(withoutCorr$t2^2 + withoutCorr$p2^2)
withoutCorr$sig = sqrt(withoutCorr$tau^2 + withoutCorr$phi^2)

ratioT1 = withCorr$t1/withoutCorr$t1
ratioT2 = withCorr$t2/withoutCorr$t2

ratioP1 = withCorr$p1/withoutCorr$p1
ratioP2 = withCorr$p2/withoutCorr$p2

ratioS1 = withCorr$s1/withoutCorr$s1
ratioS2 = withCorr$s2/withoutCorr$s2

ratioTAU = withCorr$tau/withoutCorr$tau
ratioPHI = withCorr$phi/withoutCorr$phi
ratioSIG = withCorr$sig/withoutCorr$sig

ratios = data.frame(periods = withCorr$periods, t1 = ratioT1, t2 = ratioT2, tau = ratioTAU, p1 = ratioP1, p2 = ratioP2,
                    phi = ratioPHI, s1 = ratioS1, s2 = ratioS2, sig = ratioSIG)

write.csv(ratios, file = "ratios.csv")