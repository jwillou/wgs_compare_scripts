setwd("/Users/jannawilloughby/GDrive/WGS_divergence/data/lifehistory/")
library(scales)
library(vioplot)
library(phytools)

#read in data
data = read.table("../data_june19.csv", header=T, sep=",")

#body size
data$meanbd = apply(cbind(data$sp1_bodysize, data$sp2_bodysize), 1, mean)
data$meanbd = data$meanbd/1000
colors  = c("firebrick2") #
classes = c("Mammalia") #  
pdf("bodysize_genomic.pdf", height=4.25, width=4)
sink(file="bodysize_genomic.txt")
plot(-100, -100, xlim=c(0, 300), ylim=c(-0.01, 0.12),xlab=NA, ylab="nuclear divergence rate", axes=F)
axis(side=1, at=seq(0,300,100), labels=T, tick=T, pos=-0.01)
axis(side=2, at=seq(0,0.10,0.05), labels=T, tick=T, pos=-10)
segments(x0=-10,x1=310, y0=-0.01, y1=-0.01)
segments(x0=-10,x1=-10, y0=-0.01, y1=0.11)
segments(x0=-10,x1=310, y0=0.11,  y1=0.11)
segments(x0=310,x1=310, y0=-0.01, y1=0.11)
#mammals
i = 1
t = data[data$class==classes[i],]
t = data.frame(gen_K80=t$gen_K80, meanbd=t$meanbd)
t = t[complete.cases(t),]
tlm = lm(gen_K80~meanbd, data=t)
x0=min(t$meanbd, na.rm=T)
x1=max(t$meanbd, na.rm=T)
y0=(tlm$coefficients[2]*min(t$meanbd, na.rm=T)) + tlm$coefficients[1]
y1=(tlm$coefficients[2]*max(t$meanbd, na.rm=T)) + tlm$coefficients[1]
segments(x0=x0, x1=x1, y0=y0, y1=y1, col=colors[i], lwd=3)
points(x=t$meanbd, y=t$gen_K80, col=alpha(colors[i], 0.8), pch=16, cex=1)
print(classes[i])
print(summary(tlm))
sink()  
dev.off()  

#mt
pdf("gentime_mt.pdf", height=4.25, width=4)
sink(file="gentime_mt.txt")
plot(-100, -100, xlim=c(0,17), ylim=c(-0.01, 0.32),xlab=NA, ylab="nuclear divergence rate", axes=F)
axis(side=1, at=seq(0,15,5), labels=T, tick=T, pos=-0.01)
axis(side=2, at=seq(0,0.30,0.1), labels=T, tick=T, pos=-0.40)
segments(x0=-0.40,x1=16,   y0=-0.01, y1=-0.01)
segments(x0=-0.40,x1=-0.40,y0=-0.01, y1=0.32)
segments(x0=-0.40,x1=16,   y0=0.32,  y1=0.32)
segments(x0=16,   x1=16,   y0=-0.01, y1=0.32)
#fish
i = 2
t = data[data$class==classes[i],]
t = data.frame(mt_K80=t$mt_K80, meangen=t$meangen)
t = t[complete.cases(t),]
tlm = lm(mt_K80~meangen, data=t)
x0=min(t$meangen, na.rm=T)
x1=max(t$meangen, na.rm=T)
y0=(tlm$coefficients[2]*min(t$meangen, na.rm=T)) + tlm$coefficients[1]
y1=(tlm$coefficients[2]*max(t$meangen, na.rm=T)) + tlm$coefficients[1]
segments(x0=x0, x1=x1, y0=y0, y1=y1, col=colors[i], lwd=1.5)
points(x=t$meangen, y=t$mt_K80, col=alpha(colors[i], 0.8), pch=16, cex=1)
print(classes[i])
print(summary(tlm))
#birds
i = 3
t = data[data$class==classes[i],]
t = data.frame(mt_K80=t$mt_K80, meangen=t$meangen)
t$mt_K80[t$meangen>15] = NA
t$meangen[t$meangen>15] = NA
t = t[complete.cases(t),]
tlm = lm(mt_K80~meangen, data=t)
x0=min(t$meangen, na.rm=T)
x1=max(t$meangen, na.rm=T)
y0=(tlm$coefficients[2]*min(t$meangen, na.rm=T)) + tlm$coefficients[1]
y1=(tlm$coefficients[2]*max(t$meangen, na.rm=T)) + tlm$coefficients[1]
segments(x0=x0, x1=x1, y0=y0, y1=y1, col=colors[i], lwd=1.5)
points(x=t$meangen, y=t$mt_K80, col=alpha(colors[i], 0.8), pch=16, cex=1)
print(classes[i])
print(summary(tlm))
#mammals
i = 1
t = data[data$class==classes[i],]
t = data.frame(mt_K80=t$mt_K80, meangen=t$meangen)
t$mt_K80[t$meangen>15] = NA
t$meangen[t$meangen>15] = NA
t = t[complete.cases(t),]
tlm = lm(mt_K80~meangen, data=t)
x0=min(t$meangen, na.rm=T)
x1=max(t$meangen, na.rm=T)
y0=(tlm$coefficients[2]*min(t$meangen, na.rm=T)) + tlm$coefficients[1]
y1=(tlm$coefficients[2]*max(t$meangen, na.rm=T)) + tlm$coefficients[1]
segments(x0=x0, x1=x1, y0=y0, y1=y1, col=colors[i], lwd=1.5)
points(x=t$meangen, y=t$mt_K80, col=alpha(colors[i], 0.8), pch=16, cex=1)
print(classes[i])
print(summary(tlm))
sink()  
dev.off()  



#### generation time ####
data$meangen = apply(cbind(data$sp1_gen_iucn, data$sp2_gen_iucn), 1, mean)

colors  = c("firebrick2", "dodgerblue2",     "goldenrod3") #
classes = c("Mammalia",   "Actinopterygii",  "Aves") #  
pdf("gentime_genomic.pdf", height=4.25, width=4)
sink(file="gentime_genomic.txt")
plot(-100, -100, xlim=c(0,17), ylim=c(-0.01, 0.12),xlab=NA, ylab="nuclear divergence rate", axes=F)
axis(side=1, at=seq(0,15,5), labels=T, tick=T, pos=-0.01)
axis(side=2, at=seq(0,0.10,0.05), labels=T, tick=T, pos=-0.40)
segments(x0=-0.40,x1=16,   y0=-0.01, y1=-0.01)
segments(x0=-0.40,x1=-0.40,y0=-0.01, y1=0.11)
segments(x0=-0.40,x1=16,   y0=0.11,  y1=0.11)
segments(x0=16,   x1=16,   y0=-0.01, y1=0.11)
#fish
i = 2
t = data[data$class==classes[i],]
t = data.frame(gen_K80=t$gen_K80, meangen=t$meangen)
t = t[complete.cases(t),]
tlm = lm(gen_K80~meangen, data=t)
x0=min(t$meangen, na.rm=T)
x1=max(t$meangen, na.rm=T)
y0=(tlm$coefficients[2]*min(t$meangen, na.rm=T)) + tlm$coefficients[1]
y1=(tlm$coefficients[2]*max(t$meangen, na.rm=T)) + tlm$coefficients[1]
segments(x0=x0, x1=x1, y0=y0, y1=y1, col=colors[i], lwd=1.5)
points(x=t$meangen, y=t$gen_K80, col=alpha(colors[i], 0.8), pch=16, cex=1)
print(classes[i])
print(summary(tlm))
#birds
i = 3
t = data[data$class==classes[i],]
t = data.frame(gen_K80=t$gen_K80, meangen=t$meangen)
t = t[complete.cases(t),]
t$gen_K80[t$gen_K80>0.03 | t$meangen>17] = NA
t$meangen[t$gen_K80>0.03 | t$meangen>17] = NA
tlm = lm(gen_K80~meangen, data=t)
x0=min(t$meangen, na.rm=T)
x1=max(t$meangen, na.rm=T)
y0=(tlm$coefficients[2]*min(t$meangen, na.rm=T)) + tlm$coefficients[1]
y1=(tlm$coefficients[2]*max(t$meangen, na.rm=T)) + tlm$coefficients[1]
segments(x0=x0, x1=x1, y0=y0, y1=y1, col=colors[i], lwd=1.5)
points(x=t$meangen, y=t$gen_K80, col=alpha(colors[i], 0.8), pch=16, cex=1)
print(classes[i])
print(summary(tlm))
#mammals
i = 1
t = data[data$class==classes[i],]
t = data.frame(gen_K80=t$gen_K80, meangen=t$meangen)
t = t[complete.cases(t),]
t$gen_K80[t$meangen>15] = NA
t$meangen[t$meangen>15] = NA
tlm = lm(gen_K80~meangen, data=t)
x0=min(t$meangen, na.rm=T)
x1=max(t$meangen, na.rm=T)
y0=(tlm$coefficients[2]*min(t$meangen, na.rm=T)) + tlm$coefficients[1]
y1=(tlm$coefficients[2]*max(t$meangen, na.rm=T)) + tlm$coefficients[1]
segments(x0=x0, x1=x1, y0=y0, y1=y1, col=colors[i], lwd=3)
points(x=t$meangen, y=t$gen_K80, col=alpha(colors[i], 0.8), pch=16, cex=1)
print(classes[i])
print(summary(tlm))
sink()  
dev.off()  

#mt
pdf("gentime_mt.pdf", height=4.25, width=4)
sink(file="gentime_mt.txt")
plot(-100, -100, xlim=c(0,17), ylim=c(-0.01, 0.32),xlab=NA, ylab="nuclear divergence rate", axes=F)
axis(side=1, at=seq(0,15,5), labels=T, tick=T, pos=-0.01)
axis(side=2, at=seq(0,0.30,0.1), labels=T, tick=T, pos=-0.40)
segments(x0=-0.40,x1=16,   y0=-0.01, y1=-0.01)
segments(x0=-0.40,x1=-0.40,y0=-0.01, y1=0.32)
segments(x0=-0.40,x1=16,   y0=0.32,  y1=0.32)
segments(x0=16,   x1=16,   y0=-0.01, y1=0.32)
#fish
i = 2
t = data[data$class==classes[i],]
t = data.frame(mt_K80=t$mt_K80, meangen=t$meangen)
t = t[complete.cases(t),]
tlm = lm(mt_K80~meangen, data=t)
x0=min(t$meangen, na.rm=T)
x1=max(t$meangen, na.rm=T)
y0=(tlm$coefficients[2]*min(t$meangen, na.rm=T)) + tlm$coefficients[1]
y1=(tlm$coefficients[2]*max(t$meangen, na.rm=T)) + tlm$coefficients[1]
segments(x0=x0, x1=x1, y0=y0, y1=y1, col=colors[i], lwd=1.5)
points(x=t$meangen, y=t$mt_K80, col=alpha(colors[i], 0.8), pch=16, cex=1)
print(classes[i])
print(summary(tlm))
#birds
i = 3
t = data[data$class==classes[i],]
t = data.frame(mt_K80=t$mt_K80, meangen=t$meangen)
t$mt_K80[t$meangen>15] = NA
t$meangen[t$meangen>15] = NA
t = t[complete.cases(t),]
tlm = lm(mt_K80~meangen, data=t)
x0=min(t$meangen, na.rm=T)
x1=max(t$meangen, na.rm=T)
y0=(tlm$coefficients[2]*min(t$meangen, na.rm=T)) + tlm$coefficients[1]
y1=(tlm$coefficients[2]*max(t$meangen, na.rm=T)) + tlm$coefficients[1]
segments(x0=x0, x1=x1, y0=y0, y1=y1, col=colors[i], lwd=1.5)
points(x=t$meangen, y=t$mt_K80, col=alpha(colors[i], 0.8), pch=16, cex=1)
print(classes[i])
print(summary(tlm))
#mammals
i = 1
t = data[data$class==classes[i],]
t = data.frame(mt_K80=t$mt_K80, meangen=t$meangen)
t$mt_K80[t$meangen>15] = NA
t$meangen[t$meangen>15] = NA
t = t[complete.cases(t),]
tlm = lm(mt_K80~meangen, data=t)
x0=min(t$meangen, na.rm=T)
x1=max(t$meangen, na.rm=T)
y0=(tlm$coefficients[2]*min(t$meangen, na.rm=T)) + tlm$coefficients[1]
y1=(tlm$coefficients[2]*max(t$meangen, na.rm=T)) + tlm$coefficients[1]
segments(x0=x0, x1=x1, y0=y0, y1=y1, col=colors[i], lwd=1.5)
points(x=t$meangen, y=t$mt_K80, col=alpha(colors[i], 0.8), pch=16, cex=1)
print(classes[i])
print(summary(tlm))
sink()  
dev.off()  





