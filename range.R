setwd("/Users/jannawilloughby/GDrive/WGS_divergence/data/range")
library(scales)
library(vioplot)
library(phytools)
colors  = c("saddlebrown", "goldenrod3", "dodgerblue2",    "firebrick2", "chartreuse3",    "darkorchid3") #
classes = c("Reptilia",    "Aves",       "Actinopterygii", "Mammalia",   "Chondrichthyes", "Amphibia") #    

#read in data
data = read.table("../data_june19.csv", header=T, sep=",")

#1. does a larger range overlap mean that species are less diverged? - NO
data$perrange = (data$gArea_Int / (data$gArea_sp1 + data$gArea_sp2))

#nuclear
pdf("rangeoverlap_genomic.pdf", height=5, width=5, onefile=T)
sink("rangeoverlap_genomic.txt")
plot(-100,-100, xlim=c(0,0.35), ylim=c(0,0.10), xlab="prop. range overlap", ylab="nuclear divergence rate")
#mammal
i = 4
t = data[data$class==classes[i],]
t$gen_K80[t$gen_K80==0] = 0.00001
t = data.frame(gen_K80=log(t$gen_K80), perrange=t$perrange)
t = t[t$perrange > 0,]
t = t[complete.cases(t),]
tlm = lm(gen_K80~perrange, data=t)
shapiro.test(tlm$residuals)
print(classes[i])
print(summary(tlm))
plotcurve = data.frame(x=t$perrange, y=predict(tlm))
plotcurve$y = exp(plotcurve$y)
plotcurve = plotcurve[order(plotcurve$x),]
#lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
points(x=t$perrange, y=exp(t$gen_K80), col=alpha(colors[i], 0.7), pch=16, cex=1)

#bird
i = 2
t = data[data$class==classes[i],]
t$gen_K80[t$gen_K80==0] = 0.00001
t = data.frame(gen_K80=log(t$gen_K80), perrange=t$perrange)
t = t[t$perrange > 0,]
t = t[complete.cases(t),]
tlm = lm(gen_K80~perrange, data=t)
shapiro.test(tlm$residuals)
print(classes[i])
print(summary(tlm))
plotcurve = data.frame(x=t$perrange, y=predict(tlm))
plotcurve$y = exp(plotcurve$y)
plotcurve = plotcurve[order(plotcurve$x),]
#lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
points(x=t$perrange, y=exp(t$gen_K80), col=alpha(colors[i], 0.7), pch=16, cex=1)

dev.off()
sink()

#mt
pdf("rangeoverlap_mt.pdf", height=5, width=5, onefile=T)
sink("rangeoverlap_mt.txt")
plot(-100,-100, xlim=c(0,0.35), ylim=c(0,0.20), xlab="prop. range overlap", ylab="mt divergence")

#mammal
i = 4
t = data[data$class==classes[i],]
t = data.frame(mtW_K80=log(t$mtW_K80), perrange=t$perrange)
t = t[t$perrange > 0,]
t = t[complete.cases(t),]
tlm = lm(mtW_K80~perrange, data=t)
shapiro.test(tlm$residuals)
print(classes[i])
print(summary(tlm))
plotcurve = data.frame(x=t$perrange, y=predict(tlm))
plotcurve$y = exp(plotcurve$y)
plotcurve = plotcurve[order(plotcurve$x),]
#lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
points(x=t$perrange, y=exp(t$mtW_K80), col=alpha(colors[i], 0.7), pch=16, cex=1)

#bird
i = 2
t = data[data$class==classes[i],]
t = data.frame(mtW_K80=log(t$mtW_K80), perrange=t$perrange)
t = t[t$perrange > 0,]
t = t[complete.cases(t),]
tlm = lm(mtW_K80~perrange, data=t)
shapiro.test(tlm$residuals)
print(classes[i])
print(summary(tlm))
plotcurve = data.frame(x=t$perrange, y=predict(tlm))
plotcurve$y = exp(plotcurve$y)
plotcurve = plotcurve[order(plotcurve$x),]
#lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
points(x=t$perrange, y=exp(t$mtW_K80), col=alpha(colors[i], 0.7), pch=16, cex=1)
dev.off()
sink()

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
  t = data[data$class==c,]
  print(c)
  print(summary(lm(genS_K80~log(rdist_km+1), data=t)))
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
