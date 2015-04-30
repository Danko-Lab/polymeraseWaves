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
Hsp_bed  <- readBed(paste(bedPth, "HspGenes_mm10.bed", sep=""), minSize=minSize)[,c(1:3,6,4:5)]
Cyt_bed  <- readBed(paste(bedPth, "CytoskeletonGenes_UpHC_WT-144sHS_unique.bed", sep=""), minSize=minSize)[,c(1:3,6,4:5)]

approx=5000
hsp_144 <- polymeraseWaveBW(wt144pl, wt144mn, wtNHSpl, wtNHSmn, Hsp_bed, TSmooth= 20, approxDist=approx, returnVal="simple", prefix="IMG/Hsp.")
cyt_144 <- polymeraseWaveBW(wt144pl, wt144mn, wtNHSpl, wtNHSmn, Cyt_bed, TSmooth= 20, approxDist=approx, returnVal="simple", prefix="IMG/Cyt.")

pdf("Hsf_Cyto.144s.pdf")
  hist(cleanup(hsp_144)$Rate, breaks=seq(0,40000,3000))
  hist(cleanup(cyt_144)$Rate, breaks=seq(0,40000,3000))

  indx_hsp <- hsp_144$minOfMax & hsp_144$minOfAvg & hsp_144$KLdivParametric > 1 
  indx_cyt <- cyt_144$minOfMax & cyt_144$minOfAvg & cyt_144$KLdivParametric > 1
  hsp.cdf <- ecdf(hsp_144$Rate[indx_hsp])
  cyt.cdf <- ecdf(cyt_144$Rate[indx_cyt])
  plot(hsp.cdf, col="dark red", xlim=c(0, 40000), ylim=c(0,1))
  par(new=TRUE)
  plot(cyt.cdf, col="black", xlim=c(0, 40000), ylim=c(0,1))
  ks.test(hsp_144$Rate[indx_hsp], cyt_144$Rate[indx_cyt]) ## NOT significant.

  boxplot(hsp_144$Rate[indx_hsp], cyt_144$Rate[indx_cyt], names=c("HSF1KO", "WT"))
  wilcox.test(hsp_144$Rate[indx_hsp], cyt_144$Rate[indx_cyt])
dev.off()

## Write out data...
save.image("cyto_hsf.RData")

write.table(cbind(Hsp_bed, hsp_144), "HSP.144s.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(cbind(Cyt_bed, cyt_144), "CYT.144s.tsv", sep="\t", quote=FALSE, row.names=FALSE)

