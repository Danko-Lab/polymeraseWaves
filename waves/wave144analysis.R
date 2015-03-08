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
wt144pl <- paste(pth,"WT_144sHS_BRs_pl.bigWig",sep="")
wt144mn <- paste(pth,"WT_144sHS_BRs_mn.bigWig",sep="")

koNHSpl <- paste(pth,"Hsf1KO_NHS_BRs_pl.bigWig",sep="")
koNHSmn <- paste(pth,"Hsf1KO_NHS_BRs_mn.bigWig",sep="")
ko144pl <- paste(pth,"Hsf1KO_144sHS_BRs_pl.bigWig",sep="")
ko144mn <- paste(pth,"Hsf1KO_144sHS_BRs_mn.bigWig",sep="")

## Bed files.
bedPth <- "/home/cgd24/storage/home/work/polymeraseWaves/data/beds/"

readBed <- function(filename, minSize) {
  dataf <- read.table(filename)
  dataf <- dataf[(dataf[,3]-dataf[,2]) > minSize,]
  dataf
}

cleanup <- function(f.d.list) {
#  f.d <- f.d.list[[NROW(f.d.list)]] ## ONLY used for alldata
  f.d[f.d$minOfMax & f.d$minOfAvg & f.d$KLdivParametric > 1,]
}

#################################
## Compare wildtype and HSF1ko
minSize=20000
f144hs_WT_HSF1ko   <- readBed(paste(bedPth, "144sHS_UpregGenes_WT-Hsf1KO.txt", sep=""), minSize=minSize)[,c(1:3,6,4:5)]
f144hs_WT_only     <- readBed(paste(bedPth, "144sHS_UpregGenes_WTonly.txt", sep=""), minSize=minSize)[,c(1:3,6,4:5)]
f144hs_HSF1ko_only <- readBed(paste(bedPth, "144sHS_UpregGenes_hsf1KOonly.txt", sep=""), minSize=minSize)[,c(1:3,6,4:5)]

approx=5000
wt144_both <- polymeraseWaveBW(wt144pl, wt144mn, wtNHSpl, wtNHSmn, f144hs_WT_HSF1ko, TSmooth= 20, approxDist=approx, returnVal="simple", prefix=NULL)
ko144_both <- polymeraseWaveBW(ko144pl, ko144mn, wtNHSpl, wtNHSmn, f144hs_WT_HSF1ko, TSmooth= 20, approxDist=approx, returnVal="alldata", prefix=NULL)

wt144_only <- polymeraseWaveBW(wt144pl, wt144mn, wtNHSpl, wtNHSmn, f144hs_WT_only, TSmooth= 20, approxDist=approx, returnVal="alldata", prefix=NULL)
ko144_only <- polymeraseWaveBW(ko144pl, ko144mn, wtNHSpl, wtNHSmn, f144hs_HSF1ko_only, TSmooth= 20, approxDist=approx, returnVal="alldata", prefix=NULL)

save.image("wt.ko.144s.RData")

pdf("WTvHSF1KO.144s.pdf")
  hist(cleanup(wt144_both$Rate), breaks=seq(0,40000,3000))
  hist(cleanup(ko144_both$Rate), breaks=seq(0,40000,3000))

  hist(cleanup(wt144_only$Rate), breaks=seq(0,40000,3000))
  hist(cleanup(ko144_only$Rate), breaks=seq(0,40000,3000))

  indxBoth <- wt144_both$minOfMax & wt144_both$minOfAvg & wt144_both$KLdivParametric > 1 & ko144_both$minOfMax & ko144_both$minOfAvg & ko144_both$KLdivParametric > 1
  wt.cdf <- ecdf(wt144_both$Rate[indxBoth])
  ko.cdf <- ecdf(ko144_both$Rate[indxBoth])
  plot(wt.cdf, col="dark red", xlim=c(0, 40000), ylim=c(0,1))
  par(new=TRUE)
  plot(ko.cdf, col="black", xlim=c(0, 40000), ylim=c(0,1))
  ks.test(ko144_both$Rate[indxBoth], wt144_both$Rate[indxBoth]) ## NOT significant.

  boxplot(ko144_both$Rate[indxBoth], wt144_both$Rate[indxBoth], names=c("HSF1KO", "WT"))
  wilcox.test(ko144_both$Rate[indxBoth], wt144_both$Rate[indxBoth])
dev.off()

