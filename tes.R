setwd("/Users/jannawilloughby/GDrive/WGS_divergence/data/tes/")
library(scales)
library(vioplot)
library(phytools)

#read in data
data = read.table("../data_sept10.csv", header=T, sep=",")

####total repetitive elements####
data$meante = apply(cbind(data$sp1_prop, data$sp2_prop), 1, mean, na.rm=T)
for(r in 1:nrow(data)){
  if(is.na(data$meante[r])){
    if(!is.na(data$sp1_prop[r])){
      data$meante[r] = data$sp1_prop[r]
      next
    }
    if(!is.na(data$sp2_prop[r])){
      data$meante[r] = data$sp2_prop[r]
      next
    }
  }
}

colors  = c("firebrick2", "dodgerblue2",     "goldenrod3") #
classes = c("Mammalia",   "Actinopterygii",  "Aves") #  

sink(file="te_genomic.txt")
pdf("te_genomic.pdf", height=5.25, width=5)

plot(-100, -100, xlim=c(-0.1,0.66), ylim=c(-0.01, 0.12),xlab=NA, ylab="nuclear divergence rate", axes=F)
axis(side=1, at=seq(0,0.6,0.1), labels=T, tick=T, pos=-0.01)
axis(side=2, at=seq(0,0.1,0.05), labels=T, tick=T, pos=-0.05)
segments(x0=-0.05, x1=0.63,  y0=-0.01, y1=-0.01)
segments(x0=-0.05, x1=-0.05, y0=-0.01, y1=0.11)
segments(x0=-0.05, x1=0.63,  y0=0.11, y1=0.11)
segments(x0=0.63,  x1=0.63,  y0=-0.01, y1=0.11)

#mammal
i = 1
t = data[data$class==classes[i],]
t = t[t$gen_K80 > 0 & t$gen_K80 < 0.06,]
t = data.frame(gen_K80=sqrt(t$gen_K80), meante=t$meante)
t = t[complete.cases(t),]
tlm = lm(gen_K80~meante, data=t)
print(classes[i])
print(summary(tlm))
plotcurve = data.frame(x=t$meante, y=predict(tlm))
plotcurve$y = (plotcurve$y)^2
plotcurve = plotcurve[order(plotcurve$x),]
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
points(x=t$meante, y=(t$gen_K80)^2, col=alpha(colors[i], 0.6), pch=16, cex=1)

#fish
i = 2
t = data[data$class==classes[i],]
t = t[t$gen_K80 < 0.08,]
t = data.frame(gen_K80=sqrt(t$gen_K80), meante=(t$meante))
t = t[complete.cases(t),]
tlm = lm(gen_K80~meante, data=t)
print(classes[i])
print(summary(tlm))
plotcurve = data.frame(x=t$meante, y=predict(tlm))
plotcurve$y = (plotcurve$y)^2
plotcurve = plotcurve[order(plotcurve$x),]
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
points(x=t$meante, y=(t$gen_K80)^2, col=alpha(colors[i], 0.6), pch=16, cex=1)

#bird
i = 3
t = data[data$class==classes[i],]
t = t[t$gen_K80 < 0.025,]
t = data.frame(gen_K80=sqrt(t$gen_K80), meante=t$meante)
t = t[complete.cases(t),]
tlm = lm(gen_K80~meante, data=t)
print(classes[i])
print(summary(tlm))
plotcurve = data.frame(x=t$meante, y=predict(tlm))
plotcurve$y = (plotcurve$y)^2
plotcurve = plotcurve[order(plotcurve$x),]
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
points(x=t$meante, y=(t$gen_K80)^2, col=alpha(colors[i], 0.6), pch=16, cex=1)

sink()
dev.off()

####non ltr####
data$meannl = apply(cbind((data$sp1_sines+data$sp1_lines)/data$sp1_tlen, (data$sp1_sines+data$sp1_lines)/data$sp2_tlen), 1, mean, na.rm=T)

colors  = c("firebrick2", "dodgerblue2",     "goldenrod3") #
classes = c("Mammalia",   "Actinopterygii",  "Aves") #  

sink(file="nonltr_genomic.txt")
pdf("nonltr_genomic.pdf", height=5.25, width=5)

plot(-100, -100, xlim=c(-0.1,0.58), ylim=c(-0.01, 0.12),xlab=NA, ylab="nuclear divergence rate", axes=F)
axis(side=1, at=seq(0,0.4,0.1), labels=T, tick=T, pos=-0.01)
axis(side=2, at=seq(0,0.1,0.05), labels=T, tick=T, pos=-0.05)
segments(x0=-0.05, x1=0.45,  y0=-0.01, y1=-0.01)
segments(x0=-0.05, x1=-0.05, y0=-0.01, y1=0.11)
segments(x0=-0.05, x1=0.45,  y0=0.11, y1=0.11)
segments(x0=0.45,  x1=0.45,  y0=-0.01, y1=0.11)

#mammal
i = 1
t = data[data$class==classes[i],]
t = t[t$gen_K80 > 0 & t$gen_K80 < 0.06,]
t = data.frame(gen_K80=sqrt(t$gen_K80), meante=(t$meannl))
t = t[complete.cases(t),]
tlm = lm(gen_K80~meante, data=t)
print(classes[i])
print(summary(tlm))
plotcurve = data.frame(x=t$meante, y=predict(tlm))
plotcurve$y = (plotcurve$y)^2
plotcurve = plotcurve[order(plotcurve$x),]
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
points(x=(t$meante), y=(t$gen_K80)^2, col=alpha(colors[i], 0.6), pch=16, cex=1)

#fish
i = 2
t = data[data$class==classes[i],]
t = t[t$gen_K80 < 0.08,]
t = data.frame(gen_K80=sqrt(t$gen_K80), meante=sqrt(t$meannl))
t = t[complete.cases(t),]
tlm = lm(gen_K80~meante, data=t)
print(classes[i])
print(summary(tlm))
plotcurve = data.frame(x=t$meante, y=predict(tlm))
plotcurve$y = (plotcurve$y)^2
plotcurve = plotcurve[order(plotcurve$x),]
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
points(x=t$meante, y=(t$gen_K80)^2, col=alpha(colors[i], 0.6), pch=16, cex=1)

#bird
i = 3
t = data[data$class==classes[i],]
t = t[t$gen_K80 < 0.025,]
t = data.frame(gen_K80=sqrt(t$gen_K80), meante=sqrt(t$meannl))
t = t[complete.cases(t),]
tlm = lm(gen_K80~meante, data=t)
print(classes[i])
print(summary(tlm))
plotcurve = data.frame(x=t$meante, y=predict(tlm))
plotcurve$y = (plotcurve$y)^2
plotcurve = plotcurve[order(plotcurve$x),]
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
points(x=t$meante, y=(t$gen_K80)^2, col=alpha(colors[i], 0.6), pch=16, cex=1)

sink()
dev.off()

####ltrs####
data$meanlt = apply(cbind(data$sp1_ltrs/data$sp1_tlen, data$sp2_ltrs/data$sp2_tlen), 1, mean, na.rm=T)
for(r in 1:nrow(data)){
  if(is.na(data$meanlt[r])){
    if(!is.na(data$sp1_ltrs[r])){
      data$meanlt[r] = data$sp1_ltrs[r]/data$sp1_tlen[r]
      next
    }
    if(!is.na(data$sp2_ltrs[r])){
      data$meanlt[r] = data$sp2_ltrs[r]/data$sp2_tlen[r]
      next
    }
  }
}

colors  = c("firebrick2", "dodgerblue2",     "goldenrod3") #
classes = c("Mammalia",   "Actinopterygii",  "Aves") #  

sink(file="ltrs_genomic.txt")
pdf("ltrs_genomic.pdf", height=5.25, width=5)

plot(-100, -100, xlim=c(-0.1,0.58), ylim=c(-0.01, 0.12),xlab=NA, ylab="nuclear divergence rate", axes=F)
axis(side=1, at=seq(0,0.4,0.1), labels=T, tick=T, pos=-0.01)
axis(side=2, at=seq(0,0.1,0.05), labels=T, tick=T, pos=-0.05)
segments(x0=-0.05, x1=0.45,  y0=-0.01, y1=-0.01)
segments(x0=-0.05, x1=-0.05, y0=-0.01, y1=0.11)
segments(x0=-0.05, x1=0.45,  y0=0.11, y1=0.11)
segments(x0=0.45,  x1=0.45,  y0=-0.01, y1=0.11)

#mammal
i = 1
t = data[data$class==classes[i],]
t = t[t$gen_K80 > 0 & t$gen_K80 < 0.06,]
t = data.frame(gen_K80=sqrt(t$gen_K80), meante=sqrt(t$meanlt))
t = t[complete.cases(t),]
tlm = lm(gen_K80~meante, data=t)
print(classes[i])
print(summary(tlm))
plotcurve = data.frame(x=t$meante, y=predict(tlm))
plotcurve$y = (plotcurve$y)^2
plotcurve = plotcurve[order(plotcurve$x),]
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
points(x=t$meante, y=(t$gen_K80)^2, col=alpha(colors[i], 0.6), pch=16, cex=1)

#fish
i = 2
t = data[data$class==classes[i],]
t = t[t$gen_K80 < 0.08,]
t = data.frame(gen_K80=sqrt(t$gen_K80), meante=sqrt(t$meanlt))
t = t[complete.cases(t),]
tlm = lm(gen_K80~meante, data=t)
print(classes[i])
print(summary(tlm))
plotcurve = data.frame(x=t$meante, y=predict(tlm))
plotcurve$y = (plotcurve$y)^2
plotcurve = plotcurve[order(plotcurve$x),]
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
points(x=t$meante, y=(t$gen_K80)^2, col=alpha(colors[i], 0.6), pch=16, cex=1)

#bird
i = 3
t = data[data$class==classes[i],]
t = t[t$gen_K80 < 0.025,]
t = data.frame(gen_K80=sqrt(t$gen_K80), meante=sqrt(t$meanlt))
t = t[complete.cases(t),]
tlm = lm(gen_K80~meante, data=t)
print(classes[i])
print(summary(tlm))
plotcurve = data.frame(x=t$meante, y=predict(tlm))
plotcurve$y = (plotcurve$y)^2
plotcurve = plotcurve[order(plotcurve$x),]
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
points(x=t$meante, y=(t$gen_K80)^2, col=alpha(colors[i], 0.6), pch=16, cex=1)

sink()
dev.off()



