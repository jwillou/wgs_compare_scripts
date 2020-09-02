setwd("/Users/jannawilloughby/GDrive/WGS_divergence/data/")
library(scales)
library(vioplot)
library(phytools)

#read in data
data = read.table("../data_feb7.csv", header=T, sep=",")

#nuclear data
data$meanGAL = apply(cbind(data$sp1_assembly_length, data$sp2_assembly_length), 1, mean, na.rm=T)
data$propGAL = data$gen_length/data$meanGAL

t = data.frame(propGAL=log(data$propGAL), gen_K80=sqrt(data$gen_K80))
tlm = lm(gen_K80~propGAL, data=t)
summary(tlm)

pdf("bp_nucdiv.pdf", height=5.25, width=5)
plotcurve = data.frame(x=t$propGAL, y=predict(tlm))
plotcurve$y = (plotcurve$y)^2
plotcurve$x = exp(plotcurve$x)
plotcurve = plotcurve[order(plotcurve$x),]
plot(-100, -100, xlim=c(0,1), ylim=c(0,0.15), xlab="prop. aligned base pairs", ylab="nuclear divergence")
lines(plotcurve$x, plotcurve$y, lwd=3, col=alpha("black", 1))
points(x=exp(t$propGAL), y=(t$gen_K80)^2, col=alpha("black", 0.7), pch=16, cex=1)
dev.off()
