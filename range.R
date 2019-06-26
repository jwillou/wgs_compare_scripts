setwd("/Users/jannawilloughby/GDrive/WGS_divergence/data/range")
library(scales)
library(vioplot)
library(phytools)

#read in data
data = read.table("../data_june19.csv", header=T, sep=",")

pdf("biospp.pdf", height=5, width=5, onefile=T)
sink("biospp_lm.txt")
#how does species range location/overlap influence divergence
colors6 = c("saddlebrown", "goldenrod3", "dodgerblue2",    "firebrick2", "chartreuse3",    "darkorchid3") #
classes = c("Reptilia",    "Aves",       "Actinopterygii", "Mammalia",   "Chondrichthyes", "Amphibia") #    

#1. does a larger range overlap mean that species are less diverged? - NO
data$perrange = (data$gArea_Int / (data$gArea_sp1 + data$gArea_sp2))
plot(-100,-100, xlim=c(0,0.35), ylim=c(0,0.15), xlab="prop. range overlap", ylab="nuclear divergence rate")
#add regression line
x0 = 0
x1 = max(data$perrange, na.rm=T)
y0 = (coef(lm(genS_K80~perrange, data=data[data$perrange>0,]))[2]*x0)+coef(lm(genS_K80~perrange, data=data[data$perrange>0,]))[1]
y1 = (coef(lm(genS_K80~perrange, data=data[data$perrange>0,]))[2]*x1)+coef(lm(genS_K80~perrange, data=data[data$perrange>0,]))[1]
segments(x0=x0, x1=x1, y1=y1, y0=y0, lty=1, col="grey20", lwd=3)
for(c in 1:length(classes)){
  t = data[data$class==as.character(classes[c]),,drop=F]
  t$perrange[t$perrange==0] = NA
  t = t[!is.na(t$perrange),,drop=FALSE]
  points(t$perrange, t$genS_K80, col=alpha(colors6[c], 0.5), pch=19, cex=1)
}
print(summary(lm(genS_K80~perrange, data=data)))

plot(-100,-100, xlim=c(0,0.35), ylim=c(0,0.3), xlab="prop. range overlap", ylab="mitochondrial divergence rate")
#add regression line
x0 = 0
x1 = max(data$perrange, na.rm=T)
y0 = (coef(lm(mtW_K80~perrange, data=data[data$perrange>0,]))[2]*x0)+coef(lm(mtW_K80~perrange, data=data[data$perrange>0,]))[1]
y1 = (coef(lm(mtW_K80~perrange, data=data[data$perrange>0,]))[2]*x1)+coef(lm(mtW_K80~perrange, data=data[data$perrange>0,]))[1]
segments(x0=x0, x1=x1, y1=y1, y0=y0, lty=1, col="grey20", lwd=3)
for(c in 1:length(classes)){
  t = data[data$class==as.character(classes[c]),,drop=F]
  t$perrange[t$perrange==0] = NA
  t = t[!is.na(t$perrange),,drop=FALSE]
  points(t$perrange, t$mtW_K80, col=alpha(colors6[c], 0.5), pch=19, cex=1)
}
print(summary(lm(mtW_K80~perrange, data=data)))

#2. if speices ranges are farther apart, estimated by the shortest distance between range edges, are they more diverged? - NO nuc, YES mt
data$rdist_km = data$gDist/1000
plot(-100,-100, xlim=c(0,13000), ylim=c(0,0.15), xlab="smallest distance between range edges (km)", ylab="divergence rate (Kimura 1980)")
#add regression line
x0 = 0
x1 = max(data$rdist_km, na.rm=T)
y0 = (coef(lm(genS_K80~rdist_km, data=data[data$rdist_km>0,]))[2]*x0)+coef(lm(genS_K80~rdist_km, data=data[data$rdist_km>0,]))[1]
y1 = (coef(lm(genS_K80~rdist_km, data=data[data$rdist_km>0,]))[2]*x1)+coef(lm(genS_K80~rdist_km, data=data[data$rdist_km>0,]))[1]
segments(x0=x0, x1=x1, y1=y1, y0=y0, lty=1, col="grey20", lwd=3)
for(c in 1:length(classes)){
  t = data[data$class==as.character(classes[c]),,drop=F]
  t$rdist_km[t$rdist_km==0] = NA
  t = t[!is.na(t$rdist_km),,drop=FALSE]
  points(t$rdist_km, t$genS_K80, col=alpha(colors6[c], 0.5), pch=19, cex=1)
}
print(summary(lm(genS_K80~log(rdist_km+1), data=data)))

plot(-100,-100, xlim=c(0,13000), ylim=c(0,0.3), xlab="smallest distance between range edges (km)", ylab="mitochondrial divergence rate")
#add regression line
x0 = 0
x1 = max(data$rdist_km, na.rm=T)
y0 = (coef(lm(mtW_K80~rdist_km, data=data[data$rdist_km>0,]))[2]*x0)+coef(lm(mtW_K80~rdist_km, data=data[data$rdist_km>0,]))[1]
y1 = (coef(lm(mtW_K80~rdist_km, data=data[data$rdist_km>0,]))[2]*x1)+coef(lm(mtW_K80~rdist_km, data=data[data$rdist_km>0,]))[1]
segments(x0=x0, x1=x1, y1=y1, y0=y0, lty=1, col="grey20", lwd=3)
for(c in 1:length(classes)){
  t = data[data$class==as.character(classes[c]),,drop=F]
  t$rdist_km[t$rdist_km==0] = NA
  t = t[!is.na(t$rdist_km),,drop=FALSE]
  points(t$rdist_km, t$mtW_K80, col=alpha(colors6[c], 0.5), pch=19, cex=1)
}
print(summary(lm(mtW_K80~log(rdist_km+1), data=data)))

#3a. when species ranges overlap more/less, are hybrids more/less likely? -NO
data$hyb = rep(0, nrow(data))
data$hyb[data$hybridize=="Y"] = 1

plot(-100, -100, xlim=c(0,0.5), ylim=c(0,1), xlab="prop. range overlap", ylab="hybridize")
for(c in 1:length(classes)){
  t = data[data$class==as.character(classes[c]),,drop=F]
  t = t[!is.na(t$perrange),,drop=FALSE]
  points(t$perrange, t$hyb, col=alpha(colors6[c], 0.5), pch=19, cex=1)
}
print(summary(glm(hyb~perrange, data=data, family="binomial")))

#3b. does range distance change probability of hybridization? - YES
plot(-100, -100, xlim=c(0,13000), ylim=c(0,1), xlab="smallest distance between range edges (km)", ylab="hybridize")
for(c in 1:length(classes)){
  t = data[data$class==as.character(classes[c]),,drop=F]
  t = t[!is.na(t$rdist_km),,drop=FALSE]
  points(t$rdist_km, t$hyb, col=alpha(colors6[c], 0.5), pch=19, cex=1)
}
print(summary(glm(hyb~log(rdist_km+1), data=data, family="binomial")))
dev.off()
sink()







####classes
for(c in c("Mammalia", "Aves")){
  t = data[!is.na(data$rdist_km),]
}
plot(-100,-100, xlim=c(0,13000), ylim=c(0,0.15), xlab="smallest distance between range edges (km)", ylab="divergence rate (Kimura 1980)")
#add regression line
x0 = 0
x1 = max(data$rdist_km, na.rm=T)
y0 = (coef(lm(genS_K80~rdist_km, data=data[data$rdist_km>0,]))[2]*x0)+coef(lm(genS_K80~rdist_km, data=data[data$rdist_km>0,]))[1]
y1 = (coef(lm(genS_K80~rdist_km, data=data[data$rdist_km>0,]))[2]*x1)+coef(lm(genS_K80~rdist_km, data=data[data$rdist_km>0,]))[1]
segments(x0=x0, x1=x1, y1=y1, y0=y0, lty=1, col="grey20", lwd=3)
for(c in 1:length(classes)){
  t = data[data$class==as.character(classes[c]),,drop=F]
  t$rdist_km[t$rdist_km==0] = NA
  t = t[!is.na(t$rdist_km),,drop=FALSE]
  points(t$rdist_km, t$genS_K80, col=alpha(colors6[c], 0.5), pch=19, cex=1)
}
print(summary(lm(genS_K80~log(rdist_km+1), data=data)))

plot(-100,-100, xlim=c(0,13000), ylim=c(0,0.3), xlab="smallest distance between range edges (km)", ylab="mitochondrial divergence rate")
#add regression line
x0 = 0
x1 = max(data$rdist_km, na.rm=T)
y0 = (coef(lm(mtW_K80~rdist_km, data=data[data$rdist_km>0,]))[2]*x0)+coef(lm(mtW_K80~rdist_km, data=data[data$rdist_km>0,]))[1]
y1 = (coef(lm(mtW_K80~rdist_km, data=data[data$rdist_km>0,]))[2]*x1)+coef(lm(mtW_K80~rdist_km, data=data[data$rdist_km>0,]))[1]
segments(x0=x0, x1=x1, y1=y1, y0=y0, lty=1, col="grey20", lwd=3)
for(c in 1:length(classes)){
  t = data[data$class==as.character(classes[c]),,drop=F]
  t$rdist_km[t$rdist_km==0] = NA
  t = t[!is.na(t$rdist_km),,drop=FALSE]
  points(t$rdist_km, t$mtW_K80, col=alpha(colors6[c], 0.5), pch=19, cex=1)
}
print(summary(lm(mtW_K80~log(rdist_km+1), data=data)))
