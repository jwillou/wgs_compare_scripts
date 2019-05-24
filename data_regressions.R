setwd("/Users/jannawilloughby/GDrive/WGS_divergence/data/")
library(scales)
library(vioplot)

#read in data
data = read.table("data_may23.csv", header=T, sep=",")

#class
classes = unique(data$class)
OUT = NULL
for(c in 1:length(classes)){
  t = data[data$class==as.character(classes[c]),,drop=F]
  d = c(mean(t$co1_K80,  na.rm=T), sd(t$co1_K80,  na.rm=T), nrow(t[!is.na(t$co1_K80),, drop=F]), (sd(t$co1_K80,  na.rm=T)/(sqrt(nrow(t[!is.na(t$co1_K80),, drop=F])))), 
        mean(t$cytb_K80, na.rm=T), sd(t$cytb_K80, na.rm=T), nrow(t[!is.na(t$cytb_K80),,drop=F]), (sd(t$cytb_K80, na.rm=T)/(sqrt(nrow(t[!is.na(t$cytb_K80),,drop=F])))), 
        mean(t$mtW_K80, na.rm=T),  sd(t$mtW_K80,  na.rm=T), nrow(t[!is.na(t$mtW_K80),, drop=F]), (sd(t$mtW_K80,  na.rm=T)/(sqrt(nrow(t[!is.na(t$mtW_K80),, drop=F])))), 
        mean(t$gen_K80, na.rm=T),  sd(t$gen_K80,  na.rm=T), nrow(t[!is.na(t$gen_K80),, drop=F]), (sd(t$gen_K80,  na.rm=T)/(sqrt(nrow(t[!is.na(t$gen_K80),, drop=F])))), 
        as.character(classes[c]))
  OUT = rbind(OUT, d)
}
class.sum = as.data.frame(OUT)
colnames(class.sum) = c("co1M", "co1SD", "co1N", "co1SE", "cytbM", "cytbSD", "cytbN", "cytbSE,", "mtwM", "mtwSD", "mtwN", "mtwSE", "genM", "genSD", "genN", "genSE")
rownames(class.sum) = seq(1, nrow(class.sum), 1)

#fix factors
for(c in 1:16){
  class.sum[,c] = as.numeric(as.character(class.sum[,c]))
}

Mcols = c(1, 5, 9, 13)
labs  = c("co1", "cytb", "mt", "genome")
l = 0
for(c in Mcols){
  l = l + 1
  plot(-100,-100, xlim=c(1.65,4.40), ylim=c(0,0.15), xlab="Class", ylab="divergence rate (Kimura 1980)", axes=T, main=labs[l])
  colors5 = c("saddlebrown", "goldenrod3", "dodgerblue3", "firebrick3", "darkorchid3")
  b = 0.35
  o = 0.1
  for(r in 2:4){
    rect(xleft=c((r-b),(r-b)), xright=c((r+b),(r+b)), ybottom=c(0,0), ytop=c(class.sum[r,c],class.sum[r,c]), col=alpha(colors5[r], 0.5))
    segments(x0=r,   y0=(class.sum[r,c]-class.sum[r,(c+3)]), x1=r,   y1=(class.sum[r,c]+class.sum[r,(c+3)]), lwd=2) #vertical
    segments(x0=r-o, y0=(class.sum[r,c]-class.sum[r,(c+3)]), x1=r+o, y1=(class.sum[r,c]-class.sum[r,(c+3)]), lwd=2) #bottom
    segments(x0=r-o, y0=(class.sum[r,c]+class.sum[r,(c+3)]), x1=r+o, y1=(class.sum[r,c]+class.sum[r,(c+3)]), lwd=2) #top
  }
}

#individuals
Mcols = c(10, 16, 22, 34)
labs    = c("co1", "cytb", "mt", "genome")
colors6 = c("saddlebrown", "goldenrod3", "dodgerblue3",    "firebrick3", "chartreuse3",    "darkorchid3")
classes = c("Reptilia",    "Aves",       "Actinopterygii", "Mammalia",   "Chondrichthyes", "Amphibia")
l = 0
for(c in Mcols){
  l = l + 1
  t = data[!is.na(data[,Mcols]),]
  t = t[!is.na(t[,1]),]
  plot(-100,-100, xlim=c(0, (nrow(t)+1)), ylim=c(0,0.3), xlab="pairs", ylab="divergence rate (Kimura 1980)", axes=T, main=labs[l])
  p  = 1
  c6 = 0
  for(g in classes){
    c6 = c6 + 1
    tt = t[t$class==as.character(g),,drop=F]
    points(x=seq(p, (p+nrow(tt)-1), 1), y=tt[,c], col=alpha(colors6[c6], 0.5), pch=19)
    p = p + nrow(tt)
  }
}

Mcols1 = c(10, 22)
Mcols2 = c(16, 34)
nlabs = 
l = 0
for(c in 1:length(Mcols1)){
  l = l + 1
  t = data[!is.na(data[,Mcols2[c]]),]
  t = t[!is.na(t[,1]),]
  p  = 1
  c6 = 0
  for(g in classes){
    c6 = c6 + 1
    tt = t[t$class==as.character(g),,drop=F]
    vioplot(x=tt[,Mcols1[c]], col = alpha(colors6[c6], 0.5), plotCentre = "line", side = "right", ylim=c(0,0.3), main=labs[l], xlab=classes[c6])
    vioplot(x=tt[,Mcols2[c]], col = alpha(colors6[c6], 0.5), plotCentre = "line", side = "left", ylim=c(0,0.3), add=T)
    p = p + nrow(tt)
  }
}



