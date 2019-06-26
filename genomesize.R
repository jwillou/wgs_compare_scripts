setwd("/Users/jannawilloughby/GDrive/WGS_divergence/data/genomesize")
library(scales)
library(vioplot)
library(phytools)

#read in data
data = read.table("../data_june19.csv", header=T, sep=",")

####genome size####
data$meanc = apply(cbind(data$sp1_cvalue, data$sp2_cvalue), 1, mean, na.rm=T)
colors  = c("firebrick2", "dodgerblue2") #
classes = c("Mammalia",   "Actinopterygii") #  
sink(file="genomesize_genomic.txt")
print(summary(lm(gen_K80~meanc, data=data[data$class=="Mammalia",])))
print(summary(lm(gen_K80~meanc, data=data[data$class=="Actinopterygii" & data$gen_K80 < 0.08,])))
sink()
pdf("genomesize_genomic.pdf", height=5.25, width=5)
plot(-100, -100, xlim=c(-0.1, 5.3), ylim=c(-0.01, 0.17),xlab=NA, ylab="nuclear divergence rate", axes=F)
axis(side=1, at=seq(0,5,1), labels=T, tick=T, pos=-0.01)
axis(side=2, at=seq(0,0.15,0.05), labels=T, tick=T, pos=-0.1)
segments(x0=-0.1, x1=5.2,  y0=-0.01, y1=-0.01)
segments(x0=-0.1, x1=-0.1, y0=-0.01, y1=0.16)
segments(x0=-0.1, x1=5.2,  y0=0.16, y1=0.16)
segments(x0=5.2,  x1=5.2,  y0=-0.01, y1=0.16)
tlm = lm(gen_K80~meanc, data=data[data$class=="Mammalia",])
x0=min(data$meanc[data$class=="Mammalia"], na.rm=T)
x1=max(data$meanc[data$class=="Mammalia"], na.rm=T)
y0=(tlm$coefficients[2]*min(data$meanc[data$class=="Mammalia"], na.rm=T)) + tlm$coefficients[1]
y1=(tlm$coefficients[2]*max(data$meanc[data$class=="Mammalia"], na.rm=T)) + tlm$coefficients[1]
segments(x0=x0, x1=x1, y0=y0, y1=y1, col=colors[1], lwd=1.5)
points(x=data$meanc[data$class=="Mammalia"], y=data$gen_K80[data$class=="Mammalia"], col=alpha(colors[1], 0.8), pch=19)
tlm = lm(gen_K80~meanc, data=data[data$class=="Actinopterygii" & data$gen_K80 < 0.08,])
x0=min(data$meanc[data$class=="Actinopterygii" & data$gen_K80 < 0.08], na.rm=T)
x1=max(data$meanc[data$class=="Actinopterygii" & data$gen_K80 < 0.08], na.rm=T)
y0=(tlm$coefficients[2]*min(data$meanc[data$class=="Actinopterygii" & data$gen_K80 < 0.08], na.rm=T)) + tlm$coefficients[1]
y1=(tlm$coefficients[2]*max(data$meanc[data$class=="Actinopterygii" & data$gen_K80 < 0.08], na.rm=T)) + tlm$coefficients[1]
segments(x0=x0, x1=x1, y0=y0, y1=y1, col=colors[2], lwd=3)
points(x=data$meanc[data$class=="Actinopterygii" & data$gen_K80 < 0.08], y=data$gen_K80[data$class=="Actinopterygii" & data$gen_K80 < 0.08], col=alpha(colors[2], 0.8), pch=19)
dev.off()
