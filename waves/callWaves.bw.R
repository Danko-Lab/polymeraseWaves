##
## callWaves.R -- identify the position of the leading edge of the wave of Pol II as it moves 
##                across gene bodies.
##

require(groHMM)
require(bigWig)
source("../polymeraseWave.bw.R")

## BigWig files.
pth <- "/home/cgd24/storage/home/work/polymeraseWaves/data/bigWigs/"

wtNHSpl <- paste(pth,"WT_NHS_BRs_pl.bigWig",sep="")
wtNHSmn <- paste(pth,"WT_NHS_BRs_mn.bigWig",sep="")
wt12pl  <- paste(pth,"WT_12HS_BRs_pl.bigWig",sep="")
wt12mn  <- paste(pth,"WT_12HS_BRs_mn.bigWig",sep="")

dkoNHSpl <- paste(pth,"HsfDKO_NHS_BRs_pl.bigWig",sep="")
dkoNHSmn <- paste(pth,"HsfDKO_NHS_BRs_mn.bigWig",sep="")
dko12pl  <- paste(pth,"HsfDKO_12HS_BRs_pl.bigWig",sep="")
dko12mn  <- paste(pth,"HsfDKO_12HS_BRs_mn.bigWig",sep="")

## Bed files.
readBed <- function(filename, minSize) {
  dataf <- read.table(filename)
  dataf <- dataf[(dataf[,3]-dataf[,2]) > minSize,]
  dataf
}

minSize=60000
upwt <- readBed("../data/beds/UpRegGenes_WT_12minHS.bed", minSize)
updko<- readBed("../data/beds/UpRegGenes_HsfDKO_12minHS.bed", minSize)

## Try to run it ... 
# 12m
approx=24000
dko12m.d.list <- polymeraseWaveBW(dko12pl, dko12mn, dkoNHSpl, dkoNHSmn, updko[,c(1:3,6,4:5)], TSmooth= 20, approxDist=approx, returnVal="alldata", prefix="IMG/up.dko.12m.")
wt12m.d.list <- polymeraseWaveBW(wt12pl, wt12mn, wtNHSpl, wtNHSmn, upwt[,c(1:3,6,4:5)], TSmooth= 20, approxDist=approx, returnVal="alldata", prefix="IMG/up.wt.12m.")

save.image("waves.12m.wt.dko.RData")

cleanup <- function(f.d.list) {
  f.d <- f.d.list[[NROW(f.d.list)]]
  f.d[f.d$minOfMax & f.d$minOfAvg & f.d$KLdivParametric > 1,]
}

dko.d <- cleanup(dko12m.d.list)
wt.d <- cleanup(wt12m.d.list)

## Compare distances in WT and DKO.
pdf("WTvDKO.pdf")

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

