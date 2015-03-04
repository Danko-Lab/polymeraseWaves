##
## callWaves.R -- identify the position of the leading edge of the wave of Pol II as it moves 
##                across gene bodies.
##

require(groHMM)
require(bigWig)
source("~/src/polymeraseWave.bw.R")

## BigWig files.
pth <- "/home/cgd24/storage/home/work/jay/data/bigWigs/"
wtNHSpl <- paste(pth,"WT_NHS_BRs_pl.bigWig",sep="")
wtNHSmn <- paste(pth,"WT_NHS_BRs_mn.bigWig",sep="")
wt12pl <- paste(pth,"WT_12HS_BRs_pl.bigWig",sep="")
wt12mn <- paste(pth,"WT_12HS_BRs_mn.bigWig",sep="")
wt60pl <- paste(pth,"WT_60HS_BRs_pl.bigWig",sep="")
wt60mn <- paste(pth,"WT_60HS_BRs_mn.bigWig",sep="")

## Bed files.
down <- read.table("../data/beds/DownRegGenes_WT_12minHS.bed")
up   <- read.table("../data/beds/UpRegGenes_WT_12minHS.bed")

minSize=60000
down <- down[(down[,3]-down[,2]) > minSize,]
up   <- up[(up[,3]-up[,2]) > minSize,]

## Try to run it ... 
# 12m
approx=24000
r12m.1.d <- polymeraseWaveBW(wt12pl, wt12mn, wtNHSpl, wtNHSmn, up[,c(1:3,6,4:5)], TSmooth= 20, approxDist=approx, returnVal="alldata", prefix="IMG/up.12m.")

# 60m
approx=120000
up   <- up[(up[,3]-up[,2]) > 150000,]
r60m.1.d <- polymeraseWaveBW(wt60pl, wt60mn, wtNHSpl, wtNHSmn, up[c(1:3),c(1:3,6,4:5)], approxDist=approx, TSmooth= 20, returnVal="alldata", prefix="IMG/up.60m.")


#r60m.1.d <- polymeraseWave(wtNr1, wt60r1, down[,c(1:3,6,4:5)], approxDist=approx, returnVal="alldata")
#r60m.2.d <- polymeraseWave(wtNr2, wt60r2, down[,c(1:3,6,4:5)], approxDist=approx, returnVal="alldata")

r12m.1.u <- polymeraseWave(wt12r1, wtNr1, up[,c(1:3,6,4:5)], approxDist=approx, returnVal="alldata")
r12m.2.u <- polymeraseWave(wt12r2, wtNr2, up[,c(1:3,6,4:5)], approxDist=approx, returnVal="alldata")

#r60m.1.u <- polymeraseWave(wt60r1, wtNr1, up[,c(1:3,6,4:5)], approxDist=approx, returnVal="alldata")
#r60m.2.u <- polymeraseWave(wt60r2, wtNr2, up[,c(1:3,6,4:5)], approxDist=approx, returnVal="alldata")

save.image("waves.RData")

write.table(cbind(r12m.1.d[[1]], down[,c(1:3,6,4:5)]), "Rate.12m.DN.r1.Rflat")
write.table(cbind(r12m.2.d[[1]], down[,c(1:3,6,4:5)]), "Rate.12m.DN.r2.Rflat")
write.table(cbind(r12m.1.u[[1]], up[,c(1:3,6,4:5)]), "Rate.12m.UP.r1.Rflat")
write.table(cbind(r12m.2.u[[1]], up[,c(1:3,6,4:5)]), "Rate.12m.UP.r2.Rflat")

