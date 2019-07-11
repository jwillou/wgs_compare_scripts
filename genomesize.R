setwd("/Users/jannawilloughby/GDrive/WGS_divergence/data/genomesize")
library(scales)
library(vioplot)
library(phytools)

#read in data
data = read.table("../data_july02.csv", header=T, sep=",")

####genome size####
data$meanc = apply(cbind(data$sp1_cvalue, data$sp2_cvalue), 1, mean, na.rm=T)
for(r in 1:nrow(data)){
  if(is.na(data$meanc[r])){
    if(!is.na(data$sp1_cvalue[r])){
      data$meanc[r] = data$sp1_cvalue[r]
      next
    }
    if(!is.na(data$sp2_cvalue[r])){
      data$meanc[r] = data$sp2_cvalue[r]
      next
    }
  }
}


colors  = c("firebrick2", "dodgerblue2",     "goldenrod3") #
classes = c("Mammalia",   "Actinopterygii",  "Aves") #  

sink(file="genomesize_genomic.txt")
pdf("genomesize_genomic_blk.pdf", height=5.25, width=5)

plot(-100, -100, xlim=c(-0.1, 5.3), ylim=c(-0.01, 0.12),xlab=NA, ylab="nuclear divergence rate", axes=F)
axis(side=1, at=seq(0,5,1), labels=T, tick=T, pos=-0.01)
axis(side=2, at=seq(0,0.1,0.05), labels=T, tick=T, pos=-0.1)
segments(x0=-0.1, x1=5.2,  y0=-0.01, y1=-0.01)
segments(x0=-0.1, x1=-0.1, y0=-0.01, y1=0.11)
segments(x0=-0.1, x1=5.2,  y0=0.11, y1=0.11)
segments(x0=5.2,  x1=5.2,  y0=-0.01, y1=0.11)

#mammal
i = 1
t = data[data$class==classes[i],]
t = t[t$gen_K80 > 0,]
t = data.frame(gen_K80=sqrt(t$gen_K80), meanc=t$meanc)
t = t[complete.cases(t),]
tlm = lm(gen_K80~meanc, data=t)
print(classes[i])
print(summary(tlm))
plotcurve = data.frame(x=t$meanc, y=predict(tlm))
plotcurve$y = (plotcurve$y)^2
plotcurve = plotcurve[order(plotcurve$x),]
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
points(x=t$meanc, y=(t$gen_K80)^2, col=alpha(colors[i], 0.6), pch=16, cex=1)

#fish
i = 2
t = data[data$class==classes[i],]
t = t[t$gen_K80 < 0.08,]
t = data.frame(gen_K80=sqrt(t$gen_K80), meanc=t$meanc)
t = t[complete.cases(t),]
tlm = lm(gen_K80~meanc, data=t)
print(classes[i])
print(summary(tlm))
plotcurve = data.frame(x=t$meanc, y=predict(tlm))
plotcurve$y = (plotcurve$y)^2
plotcurve = plotcurve[order(plotcurve$x),]
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
points(x=t$meanc, y=(t$gen_K80)^2, col=alpha(colors[i], 0.6), pch=16, cex=1)

#bird
i = 3
t = data[data$class==classes[i],]
t = t[t$gen_K80 < 0.025,]
t = data.frame(gen_K80=sqrt(t$gen_K80), meanc=t$meanc)
t = t[complete.cases(t),]
tlm = lm(gen_K80~meanc, data=t)
print(classes[i])
print(summary(tlm))
plotcurve = data.frame(x=t$meanc, y=predict(tlm))
plotcurve$y = (plotcurve$y)^2
plotcurve = plotcurve[order(plotcurve$x),]
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha(colors[i], 1))
points(x=t$meanc, y=(t$gen_K80)^2, col=alpha(colors[i], 0.6), pch=16, cex=1)

sink()
dev.off()
