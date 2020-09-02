setwd("/Users/jannawilloughby/GDrive/WGS_divergence/data/lifehistory/")
library(scales)
library(vioplot)
library(phytools)

#read in data
data = read.table("../data_feb7.csv", header=T, sep=",")

####body size####
sink(file="bodysize_genomic.txt")

data$meanbd = apply(cbind(data$sp1_bodysize, data$sp2_bodysize), 1, mean)
for(r in 1:nrow(data)){
  if(is.na(data$meanbd[r])){
    if(!is.na(data$sp1_bodysize[r])){
      data$meanbd[r] = data$sp1_bodysize[r]
      next
    }
    if(!is.na(data$sp2_bodysize[r])){
      data$meanbd[r] = data$sp2_bodysize[r]
      next
    }
  }
}

colors  = c("firebrick2", "goldenrod3", "dodgerblue2") #
classes = c("Mammalia",   "Aves",       "Actinopterygii") #  

#nuclear
#mammals
pdf("bodysize_genomic_mammal.pdf", height=5.25, width=5)
par(bg=NA)
plot(-100, -100, xlim=c(0, 501000), ylim=c(-0.01, 0.12),xlab=NA, ylab="nuclear divergence", axes=F)
axis(side=1, at=seq(0,500000,100000), labels=T, tick=T, pos=-0.01)
axis(side=2, at=seq(0,0.10,0.05), labels=T, tick=T, pos=-10000)
segments(x0=-10000,x1=505000, y0=-0.01, y1=-0.01)
segments(x0=-10000,x1=-10000, y0=-0.01, y1=0.11)
segments(x0=-10000,x1=505000, y0=0.11,  y1=0.11)
segments(x0=505000,x1=505000, y0=-0.01, y1=0.11)

i = 1
t = data[data$class==classes[i],]
t = t[t$meanbd < 800000, ]
t = data.frame(meanbd=log(t$meanbd), gen_K80=sqrt(t$gen_K80))
t = t[complete.cases(t),]
tlm = lm(gen_K80~meanbd, data=t)
print("mammal body size")
print(summary(tlm))
shapiro.test(tlm$residuals)
plotcurve = data.frame(x=t$meanbd, y=predict(tlm))
plotcurve$y = (plotcurve$y)^2
plotcurve$x = exp(plotcurve$x)
plotcurve = plotcurve[order(plotcurve$x),]
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
points(x=exp(t$meanbd), y=(t$gen_K80)^2, col=alpha(colors[i], 0.7), pch=16, cex=1)
dev.off()

pdf("bodysize_genomic_mammal_inset.pdf", height=4, width=4)
par(bg="white")
plot(-100, -100, xlim=c(0, 15), ylim=c(-0.3, 0.3),xlab=NA, ylab="nuclear divergence")
points(x=(t$meanbd), y=(t$gen_K80), col=alpha(colors[i], 0.7), pch=16, cex=1)
segments(x0=min(t$meanbd, na.rm=T), x1=max(t$meanbd, na.rm=T), y0=(coefficients(tlm)[2]*min(t$meanbd, na.rm=T)+coefficients(tlm)[1]), y1=(coefficients(tlm)[2]*max(t$meanbd, na.rm=T)+coefficients(tlm)[1]), lwd=3, col=alpha(colors[i]))
dev.off()

#birds
pdf("bodysize_genomic_bird.pdf", height=5.25, width=5)
par(bg=NA)
plot(-100, -100, xlim=c(0, 6050), ylim=c(-0.01, 0.12),xlab=NA, ylab="nuclear divergence", axes=F)
axis(side=1, at=seq(0,6000,1000), labels=T, tick=T, pos=-0.01)
axis(side=2, at=seq(0,0.10,0.05), labels=T, tick=T, pos=-100)
segments(x0=-100,x1=6040, y0=-0.01, y1=-0.01)
segments(x0=-100,x1=-100, y0=-0.01, y1=0.11)
segments(x0=-100,x1=6040, y0=0.11,  y1=0.11)
segments(x0=6040,x1=6040, y0=-0.01, y1=0.11)

i = 2
t = data[data$class==classes[i],]
t = data.frame(meanbd=log(t$meanbd), gen_K80=sqrt(t$gen_K80))
t = t[complete.cases(t),]
tlm = lm(gen_K80~meanbd, data=t)
print("bird body size")
print(summary(tlm))
shapiro.test(tlm$residuals)
plotcurve = data.frame(x=t$meanbd, y=predict(tlm))
plotcurve$y = (plotcurve$y)^2
plotcurve$x = exp(plotcurve$x)
plotcurve = plotcurve[order(plotcurve$x),]
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
points(x=exp(t$meanbd), y=(t$gen_K80)^2, col=alpha(colors[i], 0.7), pch=16, cex=1)
dev.off()

pdf("bodysize_genomic_bird_inset.pdf", height=4, width=4)
par(bg="white")
plot(-100, -100, xlim=c(0, 10), ylim=c(-0.3, 0.3),xlab=NA, ylab="nuclear divergence")
points(x=(t$meanbd), y=(t$gen_K80), col=alpha(colors[i], 0.7), pch=16, cex=1)
segments(x0=min(t$meanbd, na.rm=T), x1=max(t$meanbd, na.rm=T), y0=(coefficients(tlm)[2]*min(t$meanbd, na.rm=T)+coefficients(tlm)[1]), y1=(coefficients(tlm)[2]*max(t$meanbd, na.rm=T)+coefficients(tlm)[1]), lwd=3, col=alpha(colors[i]))
dev.off()

#fish
pdf("bodysize_genomic_fish.pdf", height=5.25, width=5)
par(bg=NA)
plot(-100, -100, xlim=c(0, 410), ylim=c(-0.01, 0.12),xlab=NA, ylab="nuclear divergence", axes=F)
axis(side=1, at=seq(0,400,100), labels=T, tick=T, pos=-0.01)
axis(side=2, at=seq(0,0.10,0.05), labels=T, tick=T, pos=-5)
segments(x0=-5,x1=405, y0=-0.01, y1=-0.01)
segments(x0=-5,x1=-5, y0=-0.01, y1=0.11)
segments(x0=-5,x1=405, y0=0.11,  y1=0.11)
segments(x0=405,x1=405, y0=-0.01, y1=0.11)
i = 3
t = data[data$class==classes[i],]
t = data.frame(meanbd=log(t$meanbd), gen_K80=sqrt(t$gen_K80))
t = t[complete.cases(t),]
tlm = lm(gen_K80~meanbd, data=t)
print("fish body size")
print(summary(tlm))
shapiro.test(tlm$residuals)
plotcurve = data.frame(x=t$meanbd, y=predict(tlm))
plotcurve$y = (plotcurve$y)^2
plotcurve$x = exp(plotcurve$x)
plotcurve = plotcurve[order(plotcurve$x),]
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
t = t[t$gen_K80 < sqrt(0.1), ]
points(x=exp(t$meanbd), y=(t$gen_K80)^2, col=alpha(colors[i], 0.7), pch=16, cex=1)
dev.off()

pdf("bodysize_genomic_fish_inset.pdf", height=4, width=4)
par(bg="white")
plot(-100, -100, xlim=c(0, 6), ylim=c(-0.01, 0.3),xlab=NA, ylab="nuclear divergence")
points(x=(t$meanbd), y=(t$gen_K80), col=alpha(colors[i], 0.7), pch=16, cex=1)
segments(x0=min(t$meanbd, na.rm=T), x1=max(t$meanbd, na.rm=T), y0=(coefficients(tlm)[2]*min(t$meanbd, na.rm=T)+coefficients(tlm)[1]), y1=(coefficients(tlm)[2]*max(t$meanbd, na.rm=T)+coefficients(tlm)[1]), lwd=3, col=alpha(colors[i]))
dev.off()

sink()  

#mt
sink(file="bodysize_mt.txt")

#mammal
pdf("bodysize_mt_mammal.pdf", height=5.25, width=5)
par(bg=NA)
plot(-100, -100, xlim=c(0, 501000), ylim=c(-0.01, 0.32),xlab=NA, ylab="mt divergence", axes=F)
axis(side=1, at=seq(0,500000,100000), labels=T, tick=T, pos=-0.01)
axis(side=2, at=seq(0,0.30,0.1), labels=T, tick=T, pos=-10000)
segments(x0=-10000,x1=505000, y0=-0.01, y1=-0.01)
segments(x0=-10000,x1=-10000, y0=-0.01, y1=0.31)
segments(x0=-10000,x1=505000, y0=0.31,  y1=0.31)
segments(x0=505000,x1=505000, y0=-0.01, y1=0.31)
i = 1
t = data[data$class==classes[i],]
t = t[t$meanbd < 800000, ]
t = data.frame(meanbd=log(t$meanbd), mtW_K80=sqrt(t$mtW_K80))
t = t[complete.cases(t),]
tlm = lm(mtW_K80~meanbd, data=t)
print("mammal body size")
print(summary(tlm))
shapiro.test(tlm$residuals)
plotcurve = data.frame(x=t$meanbd, y=predict(tlm))
plotcurve$y = (plotcurve$y)^2
plotcurve$x = exp(plotcurve$x)
plotcurve = plotcurve[order(plotcurve$x),]
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
points(x=exp(t$meanbd), y=(t$mtW_K80)^2, col=alpha(colors[i], 0.7), pch=16, cex=1)
dev.off()

pdf("bodysize_mt_mammal_inset.pdf", height=4, width=4)
par(bg="white")
plot(-100, -100, xlim=c(0, 15), ylim=c(0, 0.55),xlab=NA, ylab="mt divergence")
points(x=(t$meanbd), y=(t$mtW_K80), col=alpha(colors[i], 0.7), pch=16, cex=1)
segments(x0=min(t$meanbd, na.rm=T), x1=max(t$meanbd, na.rm=T), y0=(coefficients(tlm)[2]*min(t$meanbd, na.rm=T)+coefficients(tlm)[1]), y1=(coefficients(tlm)[2]*max(t$meanbd, na.rm=T)+coefficients(tlm)[1]), lwd=3, col=alpha(colors[i]))
dev.off()

#bird
pdf("bodysize_mt_bird.pdf", height=5.25, width=5)
par(bg=NA)
plot(-100, -100, xlim=c(0, 6050), ylim=c(-0.01, 0.32),xlab=NA, ylab="mt divergence", axes=F)
axis(side=1, at=seq(0,6000,1000), labels=T, tick=T, pos=-0.01)
axis(side=2, at=seq(0,0.30,0.1), labels=T, tick=T, pos=-100)
segments(x0=-100,x1=6040, y0=-0.01, y1=-0.01)
segments(x0=-100,x1=-100, y0=-0.01, y1=0.31)
segments(x0=-100,x1=6040, y0=0.31,  y1=0.31)
segments(x0=6040,x1=6040, y0=-0.01, y1=0.31)
i = 2
t = data[data$class==classes[i],]
t = data.frame(meanbd=log(t$meanbd), mtW_K80=sqrt(t$mtW_K80))
t = t[complete.cases(t),]
tlm = lm(mtW_K80~meanbd, data=t)
print("bird body size")
print(summary(tlm))
shapiro.test(tlm$residuals)
plotcurve = data.frame(x=t$meanbd, y=predict(tlm))
plotcurve$y = (plotcurve$y)^2
plotcurve$x = exp(plotcurve$x)
plotcurve = plotcurve[order(plotcurve$x),]
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
points(x=exp(t$meanbd), y=(t$mtW_K80)^2, col=alpha(colors[i], 0.7), pch=16, cex=1)
dev.off()

pdf("bodysize_mt_bird_inset.pdf", height=4, width=4)
par(bg="white")
plot(-100, -100, xlim=c(0, 10), ylim=c(0, 0.55),xlab=NA, ylab="mt divergence")
points(x=(t$meanbd), y=(t$mtW_K80), col=alpha(colors[i], 0.7), pch=16, cex=1)
segments(x0=min(t$meanbd, na.rm=T), x1=max(t$meanbd, na.rm=T), y0=(coefficients(tlm)[2]*min(t$meanbd, na.rm=T)+coefficients(tlm)[1]), y1=(coefficients(tlm)[2]*max(t$meanbd, na.rm=T)+coefficients(tlm)[1]), lwd=3, col=alpha(colors[i]))
dev.off()

#fish
pdf("bodysize_mt_fish.pdf", height=5.25, width=5)
par(bg=NA)
plot(-100, -100, xlim=c(0, 410), ylim=c(-0.01, 0.32),xlab=NA, ylab="mt divergence", axes=F)
axis(side=1, at=seq(0,400,100), labels=T, tick=T, pos=-0.01)
axis(side=2, at=seq(0,0.30,0.1), labels=T, tick=T, pos=-5)
segments(x0=-5,x1=405, y0=-0.01, y1=-0.01)
segments(x0=-5,x1=-5, y0=-0.01, y1=0.31)
segments(x0=-5,x1=405, y0=0.31,  y1=0.31)
segments(x0=405,x1=405, y0=-0.01, y1=0.31)
i = 3
t = data[data$class==classes[i],]
t = data.frame(meanbd=log(t$meanbd), mtW_K80=sqrt(t$mtW_K80))
t = t[complete.cases(t),]
tlm = lm(mtW_K80~meanbd, data=t)
print("fish body size")
print(summary(tlm))
shapiro.test(tlm$residuals)
plotcurve = data.frame(x=t$meanbd, y=predict(tlm))
plotcurve$y = (plotcurve$y)^2
plotcurve$x = exp(plotcurve$x)
plotcurve = plotcurve[order(plotcurve$x),]
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
points(x=exp(t$meanbd), y=(t$mtW_K80)^2, col=alpha(colors[i], 0.7), pch=16, cex=1)
dev.off()

pdf("bodysize_mt_fish_inset.pdf", height=4, width=4)
par(bg="white")
plot(-100, -100, xlim=c(0, 10), ylim=c(0, 0.55),xlab=NA, ylab="mt divergence")
points(x=(t$meanbd), y=(t$mtW_K80), col=alpha(colors[i], 0.7), pch=16, cex=1)
segments(x0=min(t$meanbd, na.rm=T), x1=max(t$meanbd, na.rm=T), y0=(coefficients(tlm)[2]*min(t$meanbd, na.rm=T)+coefficients(tlm)[1]), y1=(coefficients(tlm)[2]*max(t$meanbd, na.rm=T)+coefficients(tlm)[1]), lwd=3, col=alpha(colors[i]))
dev.off()

sink()  



#### generation time ####
data$meangen = apply(cbind(data$sp1_gen_iucn, data$sp2_gen_iucn), 1, mean)
for(r in 1:nrow(data)){
  if(is.na(data$meangen[r])){
    if(!is.na(data$sp1_gen_iucn[r])){
      data$meangen[r] = data$sp1_gen_iucn[r]
      next
    }
    if(!is.na(data$sp2_gen_iucn[r])){
      data$meangen[r] = data$sp2_gen_iucn[r]
      next
    }
  }
}
colors  = c("firebrick2", "dodgerblue2",     "goldenrod3") #
classes = c("Mammalia",   "Actinopterygii",  "Aves") #  

pdf("gentime_genomic.pdf", height=5.25, width=5)
sink(file="gentime_genomic.txt")

plot(-100, -100, xlim=c(0,22), ylim=c(-0.01, 0.12),xlab=NA, ylab="nuclear divergence", axes=F)
axis(side=1, at=seq(0,20,5), labels=T, tick=T, pos=-0.01)
axis(side=2, at=seq(0,0.10,0.05), labels=T, tick=T, pos=-0.40)
segments(x0=-0.40,x1=21,   y0=-0.01, y1=-0.01)
segments(x0=-0.40,x1=-0.40,y0=-0.01, y1=0.11)
segments(x0=-0.40,x1=21,   y0=0.11,  y1=0.11)
segments(x0=21,   x1=21,   y0=-0.01, y1=0.11)

#mammal
i = 1
t = data[data$class==classes[i],]
t = t[t$meangen < 20, ]
t = data.frame(gen_K80=sqrt(t$gen_K80), meangen=t$meangen)
t = t[complete.cases(t),]
tlm = lm(gen_K80~meangen, data=t)
shapiro.test(tlm$residuals)
print(classes[i])
print(summary(tlm))
plotcurve = data.frame(x=t$meangen, y=predict(tlm))
plotcurve$y = (plotcurve$y)^2
plotcurve = plotcurve[order(plotcurve$x),]
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
points(x=t$meangen, y=(t$gen_K80)^2, col=alpha(colors[i], 0.7), pch=16, cex=1)

#fish
i = 2
t = data[data$class==classes[i],]
t = data.frame(gen_K80=sqrt(t$gen_K80), meangen=t$meangen)
t = t[complete.cases(t),]
tlm = lm(gen_K80~meangen, data=t)
shapiro.test(tlm$residuals)
print(classes[i])
print(summary(tlm))
plotcurve = data.frame(x=t$meangen, y=predict(tlm))
plotcurve$y = (plotcurve$y)^2
plotcurve = plotcurve[order(plotcurve$x),]
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
points(x=t$meangen, y=(t$gen_K80)^2, col=alpha(colors[i], 0.7), pch=16, cex=1)

#birds
i = 3
t = data[data$class==classes[i],]
t = t[t$gen_K80 < 0.025, ]
t = data.frame(gen_K80=sqrt(t$gen_K80), meangen=t$meangen)
t = t[complete.cases(t),]
tlm = lm(gen_K80~meangen, data=t)
shapiro.test(tlm$residuals)
print(classes[i])
print(summary(tlm))
plotcurve = data.frame(x=t$meangen, y=predict(tlm))
plotcurve$y = (plotcurve$y)^2
plotcurve = plotcurve[order(plotcurve$x),]
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
points(x=t$meangen, y=(t$gen_K80)^2, col=alpha(colors[i], 0.7), pch=16, cex=1)

sink()  
dev.off()  

#mt
pdf("gentime_mt.pdf", height=5.25, width=5)
sink(file="gentime_mt.txt")

plot(-100, -100, xlim=c(0,22), ylim=c(-0.01, 0.22),xlab=NA, ylab="mt divergence", axes=F)
axis(side=1, at=seq(0,20,5), labels=T, tick=T, pos=-0.01)
axis(side=2, at=seq(0,0.20,0.1), labels=T, tick=T, pos=-0.40)
segments(x0=-0.40,x1=21,   y0=-0.01, y1=-0.01)
segments(x0=-0.40,x1=-0.40,y0=-0.01, y1=0.22)
segments(x0=-0.40,x1=21,   y0=0.22,  y1=0.22)
segments(x0=21,   x1=21,   y0=-0.01, y1=0.22)

#fish
i = 2
t = data[data$class==classes[i],]
t = data.frame(mtW_K80=sqrt(t$mtW_K80), meangen=t$meangen)
t = t[complete.cases(t),]
tlm = lm(mtW_K80~meangen, data=t)
shapiro.test(tlm$residuals)
print(classes[i])
print(summary(tlm))
plotcurve = data.frame(x=t$meangen, y=predict(tlm))
plotcurve$y = (plotcurve$y)^2
plotcurve = plotcurve[order(plotcurve$x),]
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
points(x=t$meangen, y=(t$mtW_K80)^2, col=alpha(colors[i], 0.7), pch=16, cex=1)

#birds
i = 3
t = data[data$class==classes[i],]
t = data.frame(mtW_K80=sqrt(t$mtW_K80), meangen=t$meangen)
t = t[complete.cases(t),]
tlm = lm(mtW_K80~meangen, data=t)
shapiro.test(tlm$residuals)
print(classes[i])
print(summary(tlm))
plotcurve = data.frame(x=t$meangen, y=predict(tlm))
plotcurve$y = (plotcurve$y)^2
plotcurve = plotcurve[order(plotcurve$x),]
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
points(x=t$meangen, y=(t$mtW_K80)^2, col=alpha(colors[i], 0.7), pch=16, cex=1)

#mammals
i = 1
t = data[data$class==classes[i],]
t = data.frame(mtW_K80=sqrt(t$mtW_K80), meangen=t$meangen)
t = t[t$meangen < 20,]
t = t[complete.cases(t),]
tlm = lm(mtW_K80~meangen, data=t)
shapiro.test(tlm$residuals)
print(classes[i])
print(summary(tlm))
plotcurve = data.frame(x=t$meangen, y=predict(tlm))
plotcurve$y = (plotcurve$y)^2
plotcurve = plotcurve[order(plotcurve$x),]
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
points(x=t$meangen, y=(t$mtW_K80)^2, col=alpha(colors[i], 0.7), pch=16, cex=1)
sink()  
dev.off()  

