setwd("/Users/jannawilloughby/GDrive/WGS_divergence/data/lifehistory/")
library(scales)
library(vioplot)
library(phytools)

#read in data
data = read.table("../data_july02.csv", header=T, sep=",")

####body size####
#mammals
pdf("bodysize_genomic.pdf", height=5.25, width=5)
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
data$meanbd = data$meanbd/1000
t = data[!is.na(data$meanbd),]

colors  = c("firebrick2", "goldenrod3") #
classes = c("Mammalia",   "Aves") #  

#nuclear
plot(-100, -100, xlim=c(0, 260), ylim=c(-0.01, 0.12),xlab=NA, ylab="nuclear divergence rate", axes=F)
axis(side=1, at=seq(0,250,50), labels=T, tick=T, pos=-0.01)
axis(side=2, at=seq(0,0.10,0.05), labels=T, tick=T, pos=-10)
segments(x0=-10,x1=255, y0=-0.01, y1=-0.01)
segments(x0=-10,x1=-10, y0=-0.01, y1=0.11)
segments(x0=-10,x1=255, y0=0.11,  y1=0.11)
segments(x0=255,x1=255, y0=-0.01, y1=0.11)

#mammal
i = 1
t = data[data$class==classes[i],]
t = t[t$meanbd < 250, ]
t = data.frame(meanbd=log(t$meanbd), gen_K80=sqrt(t$gen_K80))
#t = t[t$meanbd > 0.1, ]
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

#bird
i = 2
t = data[data$class==classes[i],]
t = data.frame(meanbd=log(t$meanbd), gen_K80=sqrt(t$gen_K80))
#t = t[t$meanbd > 0.1, ]
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

sink()  
dev.off()  

#mt
pdf("bodysize_mt.pdf", height=5.25, width=5)
sink(file="bodysize_mt.txt")
plot(-100, -100, xlim=c(0, 260), ylim=c(-0.01, 0.32),xlab=NA, ylab="mitochondrial divergence", axes=F)
axis(side=1, at=seq(0,250,50), labels=T, tick=T, pos=-0.01)
axis(side=2, at=seq(0,0.30,0.1), labels=T, tick=T, pos=-10)
segments(x0=-10,x1=255,y0=-0.01, y1=-0.01)
segments(x0=-10,x1=-10,y0=-0.01, y1=0.32)
segments(x0=-10,x1=255,y0=0.32,  y1=0.32)
segments(x0=255,x1=255,y0=-0.01, y1=0.32)

#mammal
i = 1
t = data[data$class==classes[i],]
t = t[t$meanbd < 250, ]
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

#bird
i = 2
t = data[data$class==classes[i],]
t = t[t$meanbd < 250, ]
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

sink()  
dev.off()  

#NEED TO ADD BIRDS


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
t = data.frame(gen_K80=log(t$gen_K80), meangen=t$meangen)
t = t[complete.cases(t),]
tlm = lm(gen_K80~meangen, data=t)
print(classes[i])
print(summary(tlm))
plotcurve = data.frame(x=t$meangen, y=predict(tlm))
plotcurve$y = exp(plotcurve$y)
plotcurve = plotcurve[order(plotcurve$x),]
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
points(x=t$meangen, y=exp(t$gen_K80), col=alpha(colors[i], 0.7), pch=16, cex=1)

#fish
i = 2
t = data[data$class==classes[i],]
t = data.frame(gen_K80=log(t$gen_K80), meangen=t$meangen)
t = t[complete.cases(t),]
tlm = lm(gen_K80~meangen, data=t)
print(classes[i])
print(summary(tlm))
plotcurve = data.frame(x=t$meangen, y=predict(tlm))
plotcurve$y = exp(plotcurve$y)
plotcurve = plotcurve[order(plotcurve$x),]
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
points(x=t$meangen, y=exp(t$gen_K80), col=alpha(colors[i], 0.7), pch=16, cex=1)

#birds
i = 3
t = data[data$class==classes[i],]
t = t[t$gen_K80 < 0.025, ]
t = data.frame(gen_K80=log(t$gen_K80), meangen=t$meangen)
t = t[complete.cases(t),]
tlm = lm(gen_K80~meangen, data=t)
print(classes[i])
print(summary(tlm))
plotcurve = data.frame(x=t$meangen, y=predict(tlm))
plotcurve$y = exp(plotcurve$y)
plotcurve = plotcurve[order(plotcurve$x),]
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
points(x=t$meangen, y=exp(t$gen_K80), col=alpha(colors[i], 0.7), pch=16, cex=1)

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
t = data.frame(mtW_K80=log(t$mtW_K80), meangen=t$meangen)
t = t[complete.cases(t),]
tlm = lm(mtW_K80~meangen, data=t)
print(classes[i])
print(summary(tlm))
plotcurve = data.frame(x=t$meangen, y=predict(tlm))
plotcurve$y = exp(plotcurve$y)
plotcurve = plotcurve[order(plotcurve$x),]
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
points(x=t$meangen, y=exp(t$mtW_K80), col=alpha(colors[i], 0.7), pch=16, cex=1)

#birds
i = 3
t = data[data$class==classes[i],]
t = data.frame(mtW_K80=log(t$mtW_K80), meangen=t$meangen)
t = t[complete.cases(t),]
tlm = lm(mtW_K80~meangen, data=t)
print(classes[i])
print(summary(tlm))
plotcurve = data.frame(x=t$meangen, y=predict(tlm))
plotcurve$y = exp(plotcurve$y)
plotcurve = plotcurve[order(plotcurve$x),]
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
points(x=t$meangen, y=exp(t$mtW_K80), col=alpha(colors[i], 0.7), pch=16, cex=1)

#mammals
i = 1
t = data[data$class==classes[i],]
t = data.frame(mtW_K80=log(t$mtW_K80), meangen=t$meangen)
t = t[t$meangen < 20,]
t = t[complete.cases(t),]
tlm = lm(mtW_K80~meangen, data=t)
print(classes[i])
print(summary(tlm))
plotcurve = data.frame(x=t$meangen, y=predict(tlm))
plotcurve$y = exp(plotcurve$y)
plotcurve = plotcurve[order(plotcurve$x),]
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
points(x=t$meangen, y=exp(t$mtW_K80), col=alpha(colors[i], 0.7), pch=16, cex=1)
sink()  
dev.off()  

