setwd("/Users/jannawilloughby/GDrive/WGS_divergence/data/sisterspp/")
library(scales)
library(vioplot)
library(phytools)

#read in data
data = read.table("../data_june19.csv", header=T, sep=",")

classes = row.names(table(data$class, data$sisterspp)[table(data$class, data$sisterspp)[,4]>0,])
OUT = NULL
for(c in classes){
  t = data[data$class==c,]
  ts = t$gen_K80[t$sisterspp=="yes"]
  tg = t$gen_K80[t$sisterspp=="no"]
  OUT = rbind(OUT, c(c, length(ts), round(mean(ts, na.rm=T), 4), round(quantile(ts, probs=0.025, na.rm=T), 4), round(quantile(ts, probs=0.975, na.rm=T), 4), length(tg), round(mean(tg, na.rm=T), 4), round(quantile(tg, probs=0.025, na.rm=T), 4), round(quantile(tg, probs=0.975, na.rm=T), 4)))
}
colnames(OUT) = c("class", "n_ss", "mean_ss", "lowerq_ss", "upperq_ss", "n_gs", "mean_gs", "lowerq_gs", "upperq_gs")
sumgen = OUT

OUT = NULL
tdata = data[!is.na(data$mt_K80),]
classes = row.names(table(tdata$class, tdata$sisterspp)[table(tdata$class, tdata$sisterspp)[,4]>0,])
for(c in classes){
  t = tdata[tdata$class==c,]
  ts = t$mt_K80[t$sisterspp=="yes"]
  tg = t$mt_K80[t$sisterspp=="no"]
  OUT = rbind(OUT, c(c, length(ts), round(mean(ts, na.rm=T), 4), round(quantile(ts, probs=0.025, na.rm=T), 4), round(quantile(ts, probs=0.975, na.rm=T), 4), length(tg), round(mean(tg, na.rm=T), 4), round(quantile(tg, probs=0.025, na.rm=T), 4), round(quantile(tg, probs=0.975, na.rm=T), 4)))
}
colnames(OUT) = c("class", "n_ss", "mean_ss", "lowerq_ss", "upperq_ss", "n_gs", "mean_gs", "lowerq_gs", "upperq_gs")
summt = OUT

#make plots
colors  = c("firebrick2", "dodgerblue2",     "goldenrod3") #
classes = c("Mammalia",   "Actinopterygii",  "Aves") #  

#genomic
sumgen = as.data.frame(sumgen)
sumgen$lowerq_gs = as.numeric(as.character(sumgen$lowerq_gs))
sumgen$upperq_gs = as.numeric(as.character(sumgen$upperq_gs))
sumgen$mean_gs   = as.numeric(as.character(sumgen$mean_gs))
sumgen$mean_ss   = as.numeric(as.character(sumgen$mean_ss))
write.table(sumgen, "sissp_genomic.csv", sep=",", col.names=T, row.names=F)
pdf("sisp_genomic.pdf", height=5, width=4)
plot(-100, -100, xlim=c(0.5,3.5), ylim=c(-0.01, 0.17), xlab=NA, ylab="nuclear divergence rate", axes=F)
axis(side=1, at=seq(1,3,1), labels=c("Mammalia",   "Actinopterygii",  "Aves"), tick=T, pos=-0.01)
axis(side=2, at=seq(0,0.15,0.05), labels=T, tick=T, pos=0.40)
segments(x0=0.40,x1=3.5, y0=-0.01,y1=-0.01)
segments(x0=0.40,x1=0.40,y0=-0.01,y1=0.16)
segments(x0=0.40,x1=3.5, y0=0.16,y1=0.16)
segments(x0=3.5,x1=3.5,y0=-0.01,y1=0.16)
i = 0
for(c in classes){
  i = i + 1 
  t = data[data$class==as.character(c) & !is.na(data$gen_K80) & data$sisterspp=="no",]
  g = sumgen[sumgen$class==c,]
  polygon(x=c((i-0.1),(i+0.1),(i+0.1),(i-0.1)), y=c(g$lowerq_gs, g$lowerq_gs, g$upperq_gs, g$upperq_gs), col=alpha(colors[i], 0.4), border=F, lwd=2)
  segments(x0=(i-0.2), x1=(i+0.2), y0=g$mean_gs, y1=g$mean_gs, lty=1, lwd=3, col="grey50")
  points(x=i+sample(seq(-0.1,0.1,0.01), nrow(t), replace=T), y=t$gen_K80, col=alpha(colors[i], 0.8), pch=16, cex=0.75)
  t = data[data$class==as.character(c) & !is.na(data$gen_K80) & data$sisterspp=="yes",]
  points(x=i+sample(seq(-0.1,0.1,0.01), nrow(t), replace=T), y=t$gen_K80, col=alpha("black", 0.8), pch=16, cex=1)
  segments(x0=(i-0.2), x1=(i+0.2), y0=g$mean_ss, y1=g$mean_ss, lty=1, lwd=3, col="black")
}
dev.off()

#mt
summt = as.data.frame(summt)
summt$lowerq_gs = as.numeric(as.character(summt$lowerq_gs))
summt$upperq_gs = as.numeric(as.character(summt$upperq_gs))
summt$mean_gs   = as.numeric(as.character(summt$mean_gs))
summt$mean_ss   = as.numeric(as.character(summt$mean_ss))
write.table(summt, "sissp_mt.csv", sep=",", col.names=T, row.names=F)
pdf("sisp_mt.pdf", height=5, width=4)
plot(-100, -100, xlim=c(0.5,3.5), ylim=c(-0.01, 0.34), xlab=NA, ylab="mitochondrial divergence rate", axes=F)
axis(side=1, at=seq(1,3,1), labels=c("Mammalia",   "Actinopterygii",  "Aves"), tick=T, pos=-0.01)
axis(side=2, at=seq(0,0.3,0.1), labels=T, tick=T, pos=0.40)
segments(x0=0.40,x1=3.5, y0=-0.01,y1=-0.01)
segments(x0=0.40,x1=0.40,y0=-0.01,y1=0.33)
segments(x0=0.40,x1=3.5, y0=0.33,y1=0.33)
segments(x0=3.5,x1=3.5,y0=-0.01,y1=0.33)
i = 0
for(c in classes){
  i = i + 1 
  t = data[data$class==as.character(c) & !is.na(data$mt_K80) & data$sisterspp=="no",]
  g = summt[summt$class==c,]
  polygon(x=c((i-0.1),(i+0.1),(i+0.1),(i-0.1)), y=c(g$lowerq_gs, g$lowerq_gs, g$upperq_gs, g$upperq_gs), col=alpha(colors[i], 0.4), border=F, lwd=2)
  segments(x0=(i-0.2), x1=(i+0.2), y0=g$mean_gs, y1=g$mean_gs, lty=1, lwd=3, col="grey50")
  points(x=i+sample(seq(-0.1,0.1,0.01), nrow(t), replace=T), y=t$mt_K80, col=alpha(colors[i], 0.8), pch=16, cex=0.75)
  t = data[data$class==as.character(c) & !is.na(data$mt_K80) & data$sisterspp=="yes",]
  points(x=i+sample(seq(-0.1,0.1,0.01), nrow(t), replace=T), y=t$mt_K80, col=alpha("black", 0.8), pch=16, cex=1)
  segments(x0=(i-0.2), x1=(i+0.2), y0=g$mean_ss, y1=g$mean_ss, lty=1, lwd=3, col="black")
}
dev.off()


