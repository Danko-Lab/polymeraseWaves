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
wt12pl <- paste(pth,"WT_12HS_BRs_pl.bigWig",sep="")
wt12mn <- paste(pth,"WT_12HS_BRs_mn.bigWig",sep="")

koNHSpl <- paste(pth,"Hsf1KO_NHS_BRs_pl.bigWig",sep="")
koNHSmn <- paste(pth,"Hsf1KO_NHS_BRs_mn.bigWig",sep="")
ko144pl <- paste(pth,"Hsf1KO_144sHS_BRs_pl.bigWig",sep="")
ko144mn <- paste(pth,"Hsf1KO_144sHS_BRs_mn.bigWig",sep="")
ko12pl <- paste(pth,"Hsf1KO_12HS_BRs_pl.bigWig",sep="")
ko12mn <- paste(pth,"Hsf1KO_12HS_BRs_mn.bigWig",sep="")

## Bed files.
bedPth <- "/home/cgd24/storage/home/work/polymeraseWaves/data/beds/"

readBed <- function(filename, minSize) {
  dataf <- read.table(filename)
  dataf <- dataf[(dataf[,3]-dataf[,2]) > minSize,]
  dataf
}

cleanup <- function(f.d) {
#  f.d <- f.d.list[[NROW(f.d.list)]] ## ONLY used for alldata
  f.d[f.d$minOfMax & f.d$minOfAvg & f.d$KLdivParametric > 1,]
}

#################################
## Compare wildtype and HSF1ko
minSize=20000
f144hs_WT_HSF1ko   <- readBed(paste(bedPth, "144sHS_UpregGenes_WT-Hsf1KO.txt", sep=""), minSize=minSize)[,c(1:3,6,4:5)]
f144hs_WT_only     <- readBed(paste(bedPth, "144sHS_UpregGenes_WTonly.txt", sep=""), minSize=minSize)[,c(1:3,6,4:5)]
f144hs_HSF1ko_only <- readBed(paste(bedPth, "144sHS_UpregGenes_Hsf1KOonly.txt", sep=""), minSize=minSize)[,c(1:3,6,4:5)]

approx=5000
wt144_both <- polymeraseWaveBW(wt144pl, wt144mn, wtNHSpl, wtNHSmn, f144hs_WT_HSF1ko, TSmooth= 20, approxDist=approx, returnVal="simple", prefix="IMG/wt144_both.")
ko144_both <- polymeraseWaveBW(ko144pl, ko144mn, wtNHSpl, wtNHSmn, f144hs_WT_HSF1ko, TSmooth= 20, approxDist=approx, returnVal="simple", prefix="IMG/ko144_both.")

wt144_only <- polymeraseWaveBW(wt144pl, wt144mn, wtNHSpl, wtNHSmn, f144hs_WT_only, TSmooth= 20, approxDist=approx, returnVal="simple", prefix="IMG/wt144_only.")
ko144_only <- polymeraseWaveBW(ko144pl, ko144mn, wtNHSpl, wtNHSmn, f144hs_HSF1ko_only, TSmooth= 20, approxDist=approx, returnVal="simple", prefix="IMG/ko144_only.")

pdf("WTvHSF1KO.144s.pdf")
  hist(cleanup(wt144_both)$Rate, breaks=seq(0,40000,3000))
  hist(cleanup(ko144_both)$Rate, breaks=seq(0,40000,3000))

  hist(cleanup(wt144_only)$Rate, breaks=seq(0,40000,3000))
  hist(cleanup(ko144_only)$Rate, breaks=seq(0,40000,3000))

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

#################################
## Compare wildtype and HSF1ko
minSize=20000
f12_AND_144.wt <- readBed(paste(bedPth, "12HS_and_144sHS_UpregGenes_WT.txt", sep=""), minSize=minSize)[,c(1:3,6,4:5)]
f12_NOT_144.wt <- readBed(paste(bedPth, "12HS_NOT_144sHS_UpregGenes_WT.txt", sep=""), minSize=minSize)[,c(1:3,6,4:5)]
f12_AND_144.ko <- readBed(paste(bedPth, "12HS_and_144sHS_UpregGenes_Hsf1KO.txt", sep=""), minSize=minSize)[,c(1:3,6,4:5)]
f12_NOT_144.ko <- readBed(paste(bedPth, "12HS_NOT_144sHS_UpregGenes_Hsf1KO.txt", sep=""), minSize=minSize)[,c(1:3,6,4:5)]

##
r12_AND_144.wt <- polymeraseWaveBW(wt12pl, wt12mn, wtNHSpl, wtNHSmn, f12_AND_144.wt, TSmooth= 20, approxDist=approx, returnVal="simple", prefix="IMG/WT_12_AND_144.")
r12_NOT_144.wt <- polymeraseWaveBW(wt12pl, wt12mn, wtNHSpl, wtNHSmn, f12_NOT_144.wt, TSmooth= 20, approxDist=approx, returnVal="simple", prefix="IMG/WT_12_NOT_144.")
r12_AND_144.ko <- polymeraseWaveBW(ko12pl, ko12mn, koNHSpl, koNHSmn, f12_AND_144.ko, TSmooth= 20, approxDist=approx, returnVal="simple", prefix="IMG/KO_12_AND_144.")
r12_NOT_144.ko <- polymeraseWaveBW(ko12pl, ko12mn, koNHSpl, koNHSmn, f12_NOT_144.ko, TSmooth= 20, approxDist=approx, returnVal="simple", prefix="IMG/KO_12_NOT_144.")

##
pdf("144s_V_12m.pdf")
  wt.and.cdf <- ecdf(cleanup(r12_AND_144.wt)$Rate)
  wt.not.cdf <- ecdf(cleanup(r12_NOT_144.wt)$Rate)
  plot(wt.and.cdf, col="dark red", xlim=c(0, 40000), ylim=c(0,1))
  par(new=TRUE)
  plot(wt.not.cdf, col="black", xlim=c(0, 40000), ylim=c(0,1))
  ks.test(cleanup(r12_AND_144.wt)$Rate, cleanup(r12_NOT_144.wt)$Rate) ## NOT significant.

  ko.and.cdf <- ecdf(cleanup(r12_AND_144.ko)$Rate)
  ko.not.cdf <- ecdf(cleanup(r12_NOT_144.ko)$Rate)
  plot(ko.and.cdf, col="dark red", xlim=c(0, 40000), ylim=c(0,1))
  par(new=TRUE)
  plot(ko.not.cdf, col="black", xlim=c(0, 40000), ylim=c(0,1))
  ks.test(cleanup(r12_AND_144.ko)$Rate, cleanup(r12_NOT_144.ko)$Rate) ## NOT significant.
dev.off()
#########################

## Write out data...
save.image("wave144data.RData")

write.table(cbind(f144hs_WT_HSF1ko, wt144_both), "WT.144s.both.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(cbind(f144hs_WT_HSF1ko, ko144_both), "KO.144s.both.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(cbind(f144hs_WT_only, wt144_only), "WT.144s.only.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(cbind(f144hs_HSF1ko_only, ko144_only), "KO.144s.only.tsv", sep="\t", quote=FALSE, row.names=FALSE)

write.table(cbind(f12_AND_144.wt, r12_AND_144.wt), "WT.12_and_144.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(cbind(f12_NOT_144.wt, r12_NOT_144.wt), "WT.12_not_144.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(cbind(f12_AND_144.ko, r12_AND_144.ko), "KO.12_and_144.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(cbind(f12_NOT_144.ko, r12_NOT_144.ko), "KO.12_not_144.tsv", sep="\t", quote=FALSE, row.names=FALSE)


