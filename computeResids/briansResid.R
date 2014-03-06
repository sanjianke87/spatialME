# Save the given resids
rm(list = ls())

library(reshape)

briansResids = read.csv('./../data/U08.residT.csv')
latLonData = read.csv("./../data/NGAW2_latLon.csv")

idx = rep(0,length(briansResids[,1]))
for(i in 1:length(idx)){
  idx[i] = which(latLonData$seq == briansResids$SeqNo[i])
}
briansResids$lat = latLonData$lat[idx]
briansResids$lon = latLonData$lon[idx]

x = data.frame(eqid = briansResids$EQID, lat = briansResids$lat, lon = briansResids$lon)

briansResids = briansResids[!duplicated(x),]

cyData = melt(briansResids, c("X","SeqNo","EQID","EQNAME","Mag","TopOfRupt","dTopOfRupt","EQRegion","HypLat","HypLon",
                        "Mech", "RV", "NM", "Dip", "AS", "FS", "Width", "JBD", "ClstD", "Rt",
                        "SrcSiteA", "Vs30", "Z1", "dZ1", "cIEP", "lat", "lon"))

cyData = subset(cyData, !is.na(value))

extractPeriod = function(per){
  per = sub("T","",per)
  return(as.numeric(sub("S","",per)))
}

cyData$periods = lapply(cyData$variable, extractPeriod)
cyData$resids = cyData$value

cyResids = cyData
save(cyResids, file = "briansResid.Rdata")
