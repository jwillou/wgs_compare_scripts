setwd("/Users/jannawilloughby/GDrive/WGS_divergence/data/range")
library(scales)
library(vioplot)
library(phytools)
colors  = c("saddlebrown", "goldenrod3", "dodgerblue2",    "firebrick2", "chartreuse3",    "darkorchid3") #
classes = c("Reptilia",    "Aves",       "Actinopterygii", "Mammalia",   "Chondrichthyes", "Amphibia") #    

#read in data
data = read.table("../data_july02.csv", header=T, sep=",")

#1. does a larger range overlap mean that species are less diverged? - NO
data$perrange = (data$gArea_Int / (data$gArea_sp1 + data$gArea_sp2))

#nuclear
pdf("rangeoverlap_genomic.pdf", height=5, width=5, onefile=T)
sink("rangeoverlap_genomic.txt")
plot(-100,-100, xlim=c(0,0.35), ylim=c(0,0.10), xlab="prop. range overlap", ylab="nuclear divergence rate")
#mammal
i = 4
t = data[data$class==classes[i],]
t = data.frame(gen_K80=sqrt(t$gen_K80), perrange=t$perrange)
t = t[t$perrange > 0,]
t = t[complete.cases(t),]
tlm = lm(gen_K80~perrange, data=t)
shapiro.test(tlm$residuals)
print(classes[i])
print(summary(tlm))
plotcurve = data.frame(x=t$perrange, y=predict(tlm))
plotcurve$y = (plotcurve$y)^2
plotcurve = plotcurve[order(plotcurve$x),]
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
points(x=t$perrange, y=(t$gen_K80)^2, col=alpha(colors[i], 0.7), pch=16, cex=1)

#bird
i = 2
t = data[data$class==classes[i],]
t = data.frame(gen_K80=sqrt(t$gen_K80), perrange=t$perrange)
t = t[t$perrange > 0,]
t = t[complete.cases(t),]
tlm = lm(gen_K80~perrange, data=t)
shapiro.test(tlm$residuals)
print(classes[i])
print(summary(tlm))
plotcurve = data.frame(x=t$perrange, y=predict(tlm))
plotcurve$y = (plotcurve$y)^2
plotcurve = plotcurve[order(plotcurve$x),]
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
points(x=t$perrange, y=(t$gen_K80)^2, col=alpha(colors[i], 0.7), pch=16, cex=1)

dev.off()
sink()

#mt
pdf("rangeoverlap_mt.pdf", height=5, width=5, onefile=T)
sink("rangeoverlap_mt.txt")
plot(-100,-100, xlim=c(0,0.35), ylim=c(0,0.20), xlab="prop. range overlap", ylab="mt divergence")

#mammal
i = 4
t = data[data$class==classes[i],]
t = data.frame(mtW_K80=sqrt(t$mtW_K80), perrange=t$perrange)
t = t[t$perrange > 0,]
t = t[complete.cases(t),]
tlm = lm(mtW_K80~perrange, data=t)
shapiro.test(tlm$residuals)
print(classes[i])
print(summary(tlm))
plotcurve = data.frame(x=t$perrange, y=predict(tlm))
plotcurve$y = (plotcurve$y)^2
plotcurve = plotcurve[order(plotcurve$x),]
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
points(x=t$perrange, y=(t$mtW_K80)^2, col=alpha(colors[i], 0.7), pch=16, cex=1)

#bird
i = 2
t = data[data$class==classes[i],]
t = data.frame(mtW_K80=sqrt(t$mtW_K80), perrange=t$perrange)
t = t[t$perrange > 0,]
t = t[complete.cases(t),]
tlm = lm(mtW_K80~perrange, data=t)
shapiro.test(tlm$residuals)
print(classes[i])
print(summary(tlm))
plotcurve = data.frame(x=t$perrange, y=predict(tlm))
plotcurve$y = (plotcurve$y)^2
plotcurve = plotcurve[order(plotcurve$x),]
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
points(x=t$perrange, y=(t$mtW_K80)^2, col=alpha(colors[i], 0.7), pch=16, cex=1)
dev.off()
sink()

#2. if speices ranges are farther apart, estimated by the shortest distance between range edges, are they more diverged? - yes bird nuc, YES mt
data$rdist_km = data$gDist/1000

#nuclear
pdf("rangedist_genomic.pdf", height=5, width=5, onefile=T)
sink("rangedist_genomic.txt")
plot(-100,-100, xlim=c(0,13000), ylim=c(0,0.1), xlab="smallest distance between range edges (km)", ylab="nuclear divergence")
#mammal
i = 4
t = data[data$class==classes[i],]
t$rdist_km[t$rdist_km==0] = NA
t = data.frame(gen_K80=sqrt(t$gen_K80), rdist_km=log(t$rdist_km))
t = t[complete.cases(t),]
tlm = lm(gen_K80~rdist_km, data=t)
shapiro.test(tlm$residuals)
print(classes[i])
print(summary(tlm))
plotcurve = data.frame(x=t$rdist_km, y=predict(tlm))
plotcurve$y = (plotcurve$y)^2
plotcurve$x = exp(plotcurve$x)
plotcurve = plotcurve[order(plotcurve$x),]
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
points(x=exp(t$rdist_km), y=(t$gen_K80)^2, col=alpha(colors[i], 0.7), pch=16, cex=1)

#bird
i = 2
t = data[data$class==classes[i],]
t$rdist_km[t$rdist_km==0] = NA
t = data.frame(gen_K80=sqrt(t$gen_K80), rdist_km=log(t$rdist_km))
t = t[complete.cases(t),]
tlm = lm(gen_K80~rdist_km, data=t)
shapiro.test(tlm$residuals)
print(classes[i])
print(summary(tlm))
plotcurve = data.frame(x=t$rdist_km, y=predict(tlm))
plotcurve$y = (plotcurve$y)^2
plotcurve$x = exp(plotcurve$x)
plotcurve = plotcurve[order(plotcurve$x),]
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
points(x=exp(t$rdist_km), y=(t$gen_K80)^2, col=alpha(colors[i], 0.7), pch=16, cex=1)

dev.off()
sink()

pdf("rangedist_mt.pdf", height=5, width=5, onefile=T)
sink("rangedist_mt.txt")
plot(-100,-100, xlim=c(0,13000), ylim=c(0,0.2), xlab="smallest distance between range edges (km)", ylab="mitochondrial divergence")
#mammal
i = 4
t = data[data$class==classes[i],]
t$rdist_km[t$rdist_km==0] = NA
t = data.frame(mtW_K80=sqrt(t$mtW_K80), rdist_km=log(t$rdist_km))
t = t[complete.cases(t),]
tlm = lm(mtW_K80~rdist_km, data=t)
shapiro.test(tlm$residuals)
print(classes[i])
print(summary(tlm))
plotcurve = data.frame(x=t$rdist_km, y=predict(tlm))
plotcurve$y = (plotcurve$y)^2
plotcurve$x = exp(plotcurve$x)
plotcurve = plotcurve[order(plotcurve$x),]
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
points(x=exp(t$rdist_km), y=(t$mtW_K80)^2, col=alpha(colors[i], 0.7), pch=16, cex=1)

#bird
i = 2
t = data[data$class==classes[i],]
t$rdist_km[t$rdist_km==0] = NA
t = data.frame(mtW_K80=sqrt(t$mtW_K80), rdist_km=log(t$rdist_km))
t = t[complete.cases(t),]
tlm = lm(mtW_K80~rdist_km, data=t)
shapiro.test(tlm$residuals)
print(classes[i])
print(summary(tlm))
plotcurve = data.frame(x=t$rdist_km, y=predict(tlm))
plotcurve$y = (plotcurve$y)^2
plotcurve$x = exp(plotcurve$x)
plotcurve = plotcurve[order(plotcurve$x),]
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
points(x=exp(t$rdist_km), y=(t$mtW_K80)^2, col=alpha(colors[i], 0.7), pch=16, cex=1)

dev.off()
sink()

