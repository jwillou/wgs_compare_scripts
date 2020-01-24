setwd("/Users/jannawilloughby/GDrive/WGS_divergence/data/nuc_mt/")
library(scales)

#read in data
data = read.table("../data_sept10.csv", header=T, sep=",")

#all classes - mt vs. genomic
colors6 = c("saddlebrown", "goldenrod3", "dodgerblue2",    "firebrick2", "chartreuse3",    "darkorchid3") #
classes = c("Reptilia",    "Aves",       "Actinopterygii", "Mammalia",   "Chondrichthyes", "Amphibia") #    

####mt_div####
pdf("nuc_mt_div_expected20.pdf", height=5, width=5)
plot(-100, -100, ylim=c(-0.005,0.180), xlim=c(0,0.37), xlab="mitochondrial divergence rate", ylab="nuclear divergence rate", axes=F)
segments(x0=-0.005, x1=0.33, y0=-0.005, y1=-0.005)
segments(x0=-0.005, x1=-0.005, y0=-0.005, y1=0.16)
segments(x0=-0.005, x1=0.33, y0=0.16, y1=0.16)
segments(x0= 0.33,  x1=0.33, y0=-0.005, y1=0.16)
axis(side=2, labels=T, tick=T, pos=-0.005, at=seq(0,0.15,0.05))
axis(side=1, labels=T, tick=T, pos=-0.005, at=seq(0,0.30,0.1))
segments(x0=0, x1=0.15, y0=0, y1=0.15, col="grey50", lty=2)
d = data.frame(y=sqrt(data$gen_K80), x = sqrt(data$gen_K80*20))
d = d[complete.cases(d),]
lmr = lm(y~x, data=d) #regression
plotcurve = data.frame(x=d$x, y=predict(lmr))
plotcurve$x = plotcurve$x^2
plotcurve$y = plotcurve$y^2
plotcurve = plotcurve[order(plotcurve$x),]
plotcurve = plotcurve[plotcurve$x <=0.3,]
lines(plotcurve$x, plotcurve$y, lwd=3, col="black")
plotcurve.expected = plotcurve
ci = 0
for(c in classes){
  ci = ci + 1
  t = data[data$class==as.character(c),,drop=F]
  tt = data.frame(y=t$gen_K80, x=(t$gen_K80*sample(seq(10,30,1), length(t$gen_K80), replace=T)))
  tt = tt[tt$x<=0.3,]
  points(y=tt$y, x=tt$x, pch=19, col=alpha(colors6[ci], 0.5))
}
dev.off()

sink("nuc_mt_div_regression.txt")
#pdf("nuc_mt_div_ga.pdf", height=5, width=5)
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
#lines(plotcurve.expected$x, plotcurve.expected$y, lwd=3, col="black", lty=2)
ci = 0
for(c in classes){
  ci = ci + 1
  #if(ci==1 | ci==5 | ci==6){next}
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
#%diff
1-mean(c(t$gen_K80/t$mtW_K80), na.rm=T)

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
1-mean(c(tdata$gen_K80/tdata$mtW_K80), na.rm=T)

#birds
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
tt = cbind(tdata$gen_K80,tdata$mtW_K80)
tt = tt[complete.cases(tt),]
tt = as.data.frame(tt)
colnames(tt) = c("gen_K80", "mtW_K80")
1-mean(c(tt$gen_K80/tt$mtW_K80), na.rm=T)
range(c(tt$gen_K80), na.rm=T)
range(c(tt$mtW_K80), na.rm=T)


c = "Amphibia"
t = data[data$class==as.character(c),,drop=F]
1-mean(c(t$gen_K80/t$mtW_K80), na.rm=T)
sink()

####mt_genes_COI####
sink("mt_gene_COI_regression.txt")
pdf("mt_gene_COI.pdf", height=5, width=5)
plot(-100, -100, ylim=c(0,0.28), xlim=c(0,0.28), ylab="mitochondrial divergence", xlab="COI divergence", axes=F)
segments(x0=-0.005, x1=0.28, y0=-0.005, y1=-0.005)
segments(x0=-0.005, x1=-0.005, y0=-0.005, y1=0.28)
segments(x0=-0.005, x1=0.28, y0=0.28, y1=0.28)
segments(x0= 0.28,  x1=0.28, y0=-0.005, y1=0.28)
axis(side=2, labels=T, tick=T, pos=-0.005, at=seq(0,0.25,0.05))
axis(side=1, labels=T, tick=T, pos=-0.005, at=seq(0,0.25,0.05))
segments(x0=0, x1=0.25, y0=0, y1=0.25, col="grey50", lty=2)
d = data.frame(y=sqrt(data$mtW_K80), x = sqrt(data$co1_K80))
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
  points(y=t$mtW_K80, x=t$co1_K80, pch=19, col=alpha(colors6[ci], 0.5))
}
dev.off()

#%diff
t = cbind(data$mtW_K80, data$co1_K80)
t = t[t[,1]>0 & t[,2]>0,]
mean(c(t[,1]/t[,2]), na.rm=T)

#mammals
tdata = data[data$class=="Mammalia",]
d = data.frame(y=sqrt(tdata$mtW_K80), x = sqrt(tdata$co1_K80))
d = d[complete.cases(d),]
lmr = lm(y~x, data=d) #regression
print("Mammalia")
print(summary(lmr))
mean(c(tdata$mtW_K80/tdata$co1_K80), na.rm=T)

#fish
tdata = data[data$class=="Actinopterygii",]
d = data.frame(y=sqrt(tdata$mtW_K80), x = sqrt(tdata$co1_K80))
d = d[complete.cases(d),]
lmr = lm(y~x, data=d) #regression
print("Actinopterygii")
print(summary(lmr))
mean(c(tdata$mtW_K80/tdata$co1_K80), na.rm=T)

#birds
tdata = data[data$class=="Aves",]
d = data.frame(y=sqrt(tdata$mtW_K80), x = sqrt(tdata$co1_K80))
d = d[complete.cases(d),]
lmr = lm(y~x, data=d) #regression
print("Aves")
print(summary(lmr))
t = tdata$mtW_K80/tdata$co1_K80
t = t[!is.infinite(t)]
mean(c(t), na.rm=T)

dev.off()
sink()



####mt_genes_cytb####
sink("mt_gene_cytb_regression.txt")
pdf("mt_gene_cytb.pdf", height=5, width=5)
plot(-100, -100, ylim=c(0,0.28), xlim=c(0,0.28), ylab="mitochondrial divergence", xlab="cytb divergence", axes=F)
segments(x0=-0.005, x1=0.28, y0=-0.005, y1=-0.005)
segments(x0=-0.005, x1=-0.005, y0=-0.005, y1=0.28)
segments(x0=-0.005, x1=0.28, y0=0.28, y1=0.28)
segments(x0= 0.28,  x1=0.28, y0=-0.005, y1=0.28)
axis(side=2, labels=T, tick=T, pos=-0.005, at=seq(0,0.25,0.05))
axis(side=1, labels=T, tick=T, pos=-0.005, at=seq(0,0.25,0.05))
segments(x0=0, x1=0.25, y0=0, y1=0.25, col="grey50", lty=2)
d = data.frame(y=sqrt(data$mtW_K80), x = sqrt(data$cytb_K80))
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
  points(y=t$mtW_K80, x=t$cytb_K80, pch=19, col=alpha(colors6[ci], 0.5))
}
dev.off()

#%diff
t = cbind(data$mtW_K80, data$cytb_K80)
t = t[t[,1]>0 & t[,2]>0,]
mean(c(t[,1]/t[,2]), na.rm=T)

#mammals
tdata = data[data$class=="Mammalia",]
d = data.frame(y=sqrt(tdata$mtW_K80), x = sqrt(tdata$cytb_K80))
d = d[complete.cases(d),]
lmr = lm(y~x, data=d) #regression
print("Mammalia")
print(summary(lmr))
mean(c(tdata$mtW_K80/tdata$cytb_K80), na.rm=T)

#fish
tdata = data[data$class=="Actinopterygii",]
d = data.frame(y=sqrt(tdata$mtW_K80), x = sqrt(tdata$cytb_K80))
d = d[complete.cases(d),]
lmr = lm(y~x, data=d) #regression
print("Actinopterygii")
print(summary(lmr))
mean(c(tdata$mtW_K80/tdata$cytb_K80), na.rm=T)

#birds
tdata = data[data$class=="Aves",]
d = data.frame(y=sqrt(tdata$mtW_K80), x = sqrt(tdata$cytb_K80))
d = d[complete.cases(d),]
lmr = lm(y~x, data=d) #regression
print("Aves")
print(summary(lmr))
t = tdata$mtW_K80/tdata$cytb_K80
t = t[!is.infinite(t)]
mean(c(t), na.rm=T)

dev.off()
sink()

####end####