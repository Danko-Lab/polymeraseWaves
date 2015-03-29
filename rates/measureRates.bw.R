##
## callWaves.R -- identify the position of the leading edge of the wave of Pol II as it moves 
##                across gene bodies.
##

require(groHMM)
require(bigWig)
source("../polymeraseWave.bw.R")

## BigWig files.
pth <- "/local/storage/projects/polymeraseWaves/data/bigWigs/"

wtNHSpl <- paste(pth,"WT_NHS_BRs_pl.bigWig",sep="")
wtNHSmn <- paste(pth,"WT_NHS_BRs_mn.bigWig",sep="")
wt12pl  <- paste(pth,"WT_12HS_BRs_pl.bigWig",sep="")
wt12mn  <- paste(pth,"WT_12HS_BRs_mn.bigWig",sep="")
wt60pl  <- paste(pth,"WT_60HS_BRs_pl.bigWig",sep="")
wt60mn  <- paste(pth,"WT_60HS_BRs_mn.bigWig",sep="")

dkoNHSpl <- paste(pth,"HsfDKO_NHS_BRs_pl.bigWig",sep="")
dkoNHSmn <- paste(pth,"HsfDKO_NHS_BRs_mn.bigWig",sep="")
dko12pl  <- paste(pth,"HsfDKO_12HS_BRs_pl.bigWig",sep="")
dko12mn  <- paste(pth,"HsfDKO_12HS_BRs_mn.bigWig",sep="")
dko60pl  <- paste(pth,"HsfDKO_60HS_BRs_pl.bigWig",sep="")
dko60mn  <- paste(pth,"HsfDKO_60HS_BRs_mn.bigWig",sep="")

## Bed files.
readBed <- function(filename, minSize) {
  dataf <- read.table(filename)
  dataf <- dataf[(dataf[,3]-dataf[,2]) > minSize,]
  dataf
}

minSize=120000
upwt <- readBed("../data/beds/UpRegGenes_WT_12minHS.bed", minSize)
updko<- readBed("../data/beds/UpRegGenes_HsfDKO_12minHS.bed", minSize)

## Try to run it ... 
# 12m
approx=24000
dko12m.d.list <- polymeraseWaveBW(dko12pl, dko12mn, dkoNHSpl, dkoNHSmn, updko[,c(1:3,6,4:5)], TSmooth= 20, approxDist=approx, returnVal="alldata", prefix="IMG/dko12.")
wt12m.d.list <- polymeraseWaveBW(wt12pl, wt12mn, wtNHSpl, wtNHSmn, upwt[,c(1:3,6,4:5)], TSmooth= 20, approxDist=approx, returnVal="alldata", prefix="IMG/wt12.")

approx=100000
dko60m.d.list <- polymeraseWaveBW(dko60pl, dko60mn, dkoNHSpl, dkoNHSmn, updko[,c(1:3,6,4:5)], TSmooth= 20, approxDist=approx, returnVal="alldata", prefix="IMG/dko60.")
wt60m.d.list <- polymeraseWaveBW(wt60pl, wt60mn, wtNHSpl, wtNHSmn, upwt[,c(1:3,6,4:5)], TSmooth= 20, approxDist=approx, returnVal="alldata", prefix="IMG/wt60.")

save.image("waves.12m.60m.wt.dko.RData")

q("no")

dko12m.d.list <- dko12m.d.list[[NROW(dko12m.d.list)]]
wt12m.d.list <- wt12m.d.list[[NROW(wt12m.d.list)]]
dko60m.d.list <- dko60m.d.list[[NROW(dko60m.d.list)]]
wt60m.d.list <- wt60m.d.list[[NROW(wt60m.d.list)]]

which.cleanup <- function(f.d) {
  return(f.d$minOfMax & f.d$minOfAvg & f.d$KLdivParametric > 1)
}

dko.d <- which.cleanup(dko12m.d.list) & which.cleanup(dko60m.d.list)
wt.d <- which.cleanup(wt12m.d.list) & which.cleanup(wt60m.d.list)

rates.wt <- wt60m.d.list$EndWave - wt12m.d.list$EndWave
rates.dko <- dko60m.d.list$EndWave - dko12m.d.list$EndWave

## Compare distances in WT and DKO.
pdf("WTvDKO.rates.pdf")

par(mfrow=c(2,1))
hist(dko.d$Rate, breaks=seq(0,75000,5000))
hist(wt.d$Rate, breaks=seq(0,75000,5000))

par(mfrow=c(1,1))
dko.cdf <- ecdf(dko.d$Rate)
wt.cdf  <- ecdf(wt.d$Rate)

plot(dko.cdf, col="dark red", xlim=c(0, 75000), ylim=c(0,1))
par(new=TRUE)
plot(wt.cdf, col="black", xlim=c(0, 75000), ylim=c(0,1))
ks.test(dko.d$Rate, wt.d$Rate) ## NOT significant.

boxplot(dko.d$Rate, wt.d$Rate, names=c("DKO", "WT"))
wilcox.test(dko.d$Rate, wt.d$Rate) ## Means *ARE* significantly different.

dev.off()

## Let's define genes that are direct HSF targets (in their promoter).

