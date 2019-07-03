setwd("/Users/jannawilloughby/GDrive/WGS_divergence/data/nuc_mt/")
library(scales)

#read in data
data = read.table("../data_july02.csv", header=T, sep=",")

#all classes - mt vs. genomic
colors6 = c("saddlebrown", "goldenrod3", "dodgerblue2",    "firebrick2", "chartreuse3",    "darkorchid3") #
classes = c("Reptilia",    "Aves",       "Actinopterygii", "Mammalia",   "Chondrichthyes", "Amphibia") #    

sink("nuc_mt_div_regression.txt")
pdf("nuc_mt_div.pdf", height=5, width=5)
plot(-100, -100, ylim=c(-0.005,0.180), xlim=c(0,0.37), xlab="mitochondrial divergence rate", ylab="nuclear divergence rate", axes=F)
segments(x0=-0.005, x1=0.33, y0=-0.005, y1=-0.005)
segments(x0=-0.005, x1=-0.005, y0=-0.005, y1=0.16)
segments(x0=-0.005, x1=0.33, y0=0.16, y1=0.16)
segments(x0= 0.33,  x1=0.33, y0=-0.005, y1=0.16)
axis(side=2, labels=T, tick=T, pos=-0.005, at=seq(0,0.15,0.05))
axis(side=1, labels=T, tick=T, pos=-0.005, at=seq(0,0.30,0.1))
segments(x0=0, x1=0.15, y0=0, y1=0.15, col="grey50", lty=2)
d = data.frame(y=sqrt(data$gen_K80), x = sqrt(data$mtW_K80))
d = d[complete.cases(d),]
lmr = lm(y~x, data=d) #regression
print("all")
print(summary(lmr))
shapiro.test(lmr$residuals)
plotcurve = data.frame(x=d$x, y=predict(lmr))
plotcurve$x = plotcurve$x^2
plotcurve$y = plotcurve$y^2
plotcurve = plotcurve[order(plotcurve$x),]
lines(plotcurve$x, plotcurve$y, lwd=3, col="black")
ci = 0
for(c in classes){
  ci = ci + 1
  t = data[data$class==as.character(c),,drop=F]
  points(y=t$gen_K80, x=t$mtW_K80, pch=19, col=alpha(colors6[ci], 0.5))
  print(c)
  print(mean(c(t$gen_K80/t$mtW_K80), na.rm=T))
  print(length(t$mtW_K80[!is.na(t$mtW_K80)]))
  print(length(t$gen_K80[!is.na(t$gen_K80)]))
}
dev.off()

#%diff
mean(c(data$gen_K80/data$mtW_K80), na.rm=T)

#mammals
pdf("nuc_mt_div_mamm.pdf", height=5, width=5)
tdata = data[data$class=="Mammalia",]
plot(-100, -100, ylim=c(-0.005,0.180), xlim=c(0,0.37), xlab="mitochondrial divergence rate", ylab="nuclear divergence rate", axes=F)
segments(x0=-0.005, x1=0.33, y0=-0.005, y1=-0.005)
segments(x0=-0.005, x1=-0.005, y0=-0.005, y1=0.16)
segments(x0=-0.005, x1=0.33, y0=0.16, y1=0.16)
segments(x0= 0.33,  x1=0.33, y0=-0.005, y1=0.16)
axis(side=2, labels=T, tick=T, pos=-0.005, at=seq(0,0.15,0.05))
axis(side=1, labels=T, tick=T, pos=-0.005, at=seq(0,0.30,0.1))
segments(x0=0, x1=0.15, y0=0, y1=0.15, col="grey50", lty=2)
d = data.frame(y=sqrt(tdata$gen_K80), x = sqrt(tdata$mtW_K80))
d = d[complete.cases(d),]
lmr = lm(y~x, data=d) #regression
print("Mammalia")
print(summary(lmr))
shapiro.test(lmr$residuals)
plotcurve = data.frame(x=d$x, y=predict(lmr))
plotcurve$x = plotcurve$x^2
plotcurve$y = plotcurve$y^2
plotcurve = plotcurve[order(plotcurve$x),]
ci = 4
c = "Mammalia"
t = data[data$class==as.character(c),,drop=F]
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors6[ci], 1))
points(y=t$gen_K80, x=t$mtW_K80, pch=19, col=alpha(colors6[ci], 0.5))
print(mean(c(t$gen_K80/t$mtW_K80), na.rm=T))
print(length(t$mtW_K80[!is.na(t$mtW_K80)]))
print(length(t$gen_K80[!is.na(t$gen_K80)]))
dev.off()

#fish
pdf("nuc_mt_div_fish.pdf", height=5, width=5)
tdata = data[data$class=="Actinopterygii",]
plot(-100, -100, ylim=c(-0.005,0.180), xlim=c(0,0.37), xlab="mitochondrial divergence rate", ylab="nuclear divergence rate", axes=F)
segments(x0=-0.005, x1=0.33, y0=-0.005, y1=-0.005)
segments(x0=-0.005, x1=-0.005, y0=-0.005, y1=0.16)
segments(x0=-0.005, x1=0.33, y0=0.16, y1=0.16)
segments(x0= 0.33,  x1=0.33, y0=-0.005, y1=0.16)
axis(side=2, labels=T, tick=T, pos=-0.005, at=seq(0,0.15,0.05))
axis(side=1, labels=T, tick=T, pos=-0.005, at=seq(0,0.30,0.1))
segments(x0=0, x1=0.15, y0=0, y1=0.15, col="grey50", lty=2)
lmr = lm(gen_K80~mtW_K80, data=tdata) #regression
d = data.frame(y=sqrt(tdata$gen_K80), x = sqrt(tdata$mtW_K80))
d = d[complete.cases(d),]
lmr = lm(y~x, data=d) #regression
print("Actinopterygii")
print(summary(lmr))
shapiro.test(lmr$residuals)
plotcurve = data.frame(x=d$x, y=predict(lmr))
plotcurve$x = plotcurve$x^2
plotcurve$y = plotcurve$y^2
plotcurve = plotcurve[order(plotcurve$x),]
ci = 3
c = "Actinopterygii"
t = data[data$class==as.character(c),,drop=F]
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors6[ci], 1))
points(y=t$gen_K80, x=t$mtW_K80, pch=19, col=alpha(colors6[ci], 0.5))
print(mean(c(t$gen_K80/t$mtW_K80), na.rm=T))
print(length(t$mtW_K80[!is.na(t$mtW_K80)]))
print(length(t$gen_K80[!is.na(t$gen_K80)]))
dev.off()

#mammals
pdf("nuc_mt_div_bird.pdf", height=5, width=5)
tdata = data[data$class=="Aves",]
plot(-100, -100, ylim=c(-0.005,0.180), xlim=c(0,0.37), xlab="mitochondrial divergence rate", ylab="nuclear divergence rate", axes=F)
segments(x0=-0.005, x1=0.33, y0=-0.005, y1=-0.005)
segments(x0=-0.005, x1=-0.005, y0=-0.005, y1=0.16)
segments(x0=-0.005, x1=0.33, y0=0.16, y1=0.16)
segments(x0= 0.33,  x1=0.33, y0=-0.005, y1=0.16)
axis(side=2, labels=T, tick=T, pos=-0.005, at=seq(0,0.15,0.05))
axis(side=1, labels=T, tick=T, pos=-0.005, at=seq(0,0.30,0.1))
segments(x0=0, x1=0.15, y0=0, y1=0.15, col="grey50", lty=2)
lmr = lm(gen_K80~mtW_K80, data=tdata) #regression
d = data.frame(y=sqrt(tdata$gen_K80), x = sqrt(tdata$mtW_K80))
d = d[complete.cases(d),]
lmr = lm(y~x, data=d) #regression
print("Aves")
print(summary(lmr))
shapiro.test(lmr$residuals)
plotcurve = data.frame(x=d$x, y=predict(lmr))
plotcurve$x = plotcurve$x^2
plotcurve$y = plotcurve$y^2
plotcurve = plotcurve[order(plotcurve$x),]
ci = 2
c = "Aves"
t = data[data$class==as.character(c),,drop=F]
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors6[ci], 1))
points(y=t$gen_K80, x=t$mtW_K80, pch=19, col=alpha(colors6[ci], 0.5))
print(mean(c(t$gen_K80/t$mtW_K80), na.rm=T))
print(length(t$mtW_K80[!is.na(t$mtW_K80)]))
print(length(t$gen_K80[!is.na(t$gen_K80)]))
dev.off()
