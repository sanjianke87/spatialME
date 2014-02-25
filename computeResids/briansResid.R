# Save the given resids
rm(list = ls())

library(reshape)

briansResids = read.csv('./../data/U08.residT.csv')

cyData = melt(briansResids, c("X","SeqNo","EQID","EQNAME","Mag","TopOfRupt","dTopOfRupt","EQRegion","HypLat","HypLon",
                        "Mech", "RV", "NM", "Dip", "AS", "FS", "Width", "JBD", "ClstD", "Rt",
                        "SrcSiteA", "Vs30", "Z1", "dZ1", "cIEP"))

cyData = subset(cyData, !is.na(value))

extractPeriod = function(per){
  per = sub("T","",per)
  return(as.numeric(sub("S","",per)))
}

cyData$periods = lapply(cyData$variable, extractPeriod)
cyData$resids = cyData$value

cyResids = cyData
save(cyResids, file = "briansResid.Rdata")
