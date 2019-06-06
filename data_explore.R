setwd("/Users/jannawilloughby/GDrive/WGS_divergence/data/explore/")
library(scales)

#read in data
data = read.table("../data_may23.csv", header=T, sep=",")

#how correlated are various estimates within the same region comparisons?
write.table(c("correlation between divergence estimates"), "div_corr.txt", sep="\t", row.names=F, col.names=F, append=F, quote=F)
write.table(" ", "div_corr.txt", sep="\t", row.names=F, col.names=F, append=T, quote=F)
lowcol = c( 9, 15, 21, 27, 33, 39)
uppcol = c(13, 19, 25, 31, 37, 43)
ordcol = c("cox1", "cytb", "mt - whole", "mt - blast", "gen - blast", "gen - subsampled")
divsts = c("raw", "K80", "K81", "F84", "TN93")
for(g in 1:length(lowcol)){
  temp = data[,lowcol[g]:uppcol[g]]
  temp = temp[complete.cases(temp),]
  tempcor = cor(temp, method=c("spearman"))
  diag(tempcor) = NA
  tempcor = round(tempcor, 6)
  colnames(tempcor) = divsts
  rownames(tempcor) = divsts
  write.table(ordcol[g], "div_corr.txt", sep="\t", row.names=F, col.names=F, append=T, quote=F)
  write.table(paste("min", round(range(tempcor, na.rm=T)[1], 6)), "div_corr.txt", sep="\t", row.names=F, col.names=F, append=T, quote=F)
  write.table(paste("max", round(range(tempcor, na.rm=T)[2], 6)), "div_corr.txt", sep="\t", row.names=F, col.names=F, append=T, quote=F)
  write.table(tempcor, "div_corr.txt", sep="\t", row.names=T, col.names=T, append=T, quote=F)
  write.table(" ", "div_corr.txt", sep="\t", row.names=F, col.names=F, append=T, quote=F)
}

#nice figure for mt vs. genomic
colors6 = c("saddlebrown", "goldenrod3", "dodgerblue2",    "firebrick2", "chartreuse3",    "darkorchid3") #
classes = c("Reptilia",    "Aves",       "Actinopterygii", "Mammalia",   "Chondrichthyes", "Amphibia") #    

pdf("../nuc_mt_div.pdf", height=5, width=5)
plot(-100, -100, ylim=c(-0.005,0.180), xlim=c(0,0.37), xlab="mitochondrial divergence rate", ylab="nuclear divergence rate", axes=F)
segments(x0=-0.005, x1=0.33, y0=-0.005, y1=-0.005)
segments(x0=-0.005, x1=-0.005, y0=-0.005, y1=0.16)
segments(x0=-0.005, x1=0.33, y0=0.16, y1=0.16)
segments(x0= 0.33,  x1=0.33, y0=-0.005, y1=0.16)
axis(side=2, labels=T, tick=T, pos=-0.005, at=seq(0,0.15,0.05))
axis(side=1, labels=T, tick=T, pos=-0.005, at=seq(0,0.30,0.1))
segments(x0=0, x1=0.15, y0=0, y1=0.15, col="grey50", lty=2)
lmr = lm(gen_K80~mtW_K80, data=data) #regression
x0 = min(data$mtW_K80, na.rm=T)
x1 = max(data$mtW_K80, na.rm=T)
y0 = (coef(lmr)[2]*x0)+coef(lmr)[1]
y1 = (coef(lmr)[2]*x1)+coef(lmr)[1]
segments(x0=x0, x1=x1, y1=y1, y0=y0, lty=1, col="grey20", lwd=3)

ci = 0
for(c in classes){
  ci = ci + 1
  t = data[data$class==as.character(c),,drop=F]
  points(y=t$gen_K80, x=t$mtW_K80, pch=19, col=alpha(colors6[ci], 0.5))
  print(mean(c(t$gen_K80/t$mtW_K80), na.rm=T))
  print(length(t$mtW_K80[!is.na(t$mtW_K80)]))
  print(length(t$gen_K80[!is.na(t$gen_K80)]))
}
dev.off()
summary(lmr)

#%diff
mean(c(data$gen_K80/data$mtW_K80), na.rm=T)


#does blasting change estimates very much?
pdf("div_blasteffect.pdf", height=4, width=4, onefile=T)
plot(-100, -100, xlim=c(0,0.3), ylim=c(0,0.3), xlab="mt - whole", ylab="mt - blast", main="divergence - Kimura 2 parameter")
segments(0,0,1,1, lty=2, col="grey50")
points(x=data$mtW_K80, y=data$mt_K80, pch=20, col=alpha("firebrick3", 0.5), cex=1.25)
dev.off()

#does having a smaller genome assembly change estimates very much?
pdf("div_crapassemblyeffect.pdf", height=4, width=4, onefile=T)
plot(-100, -100, xlim=c(0,0.3), ylim=c(0,0.3), xlab="gen - blast", ylab="gen - subsampled", main="divergence - Kimura 2 parameter")
segments(0,0,1,1, lty=2, col="grey50")
points(x=data$gen_K80, y=data$genS_K80, pch=20, col=alpha("firebrick3", 0.5), cex=1.25)
dev.off()

#how similar are different divergence from different genetic data?

#prep comparisons data frame
K80cols = c(10, 16, 22, 34)
compare = expand.grid(K80cols, K80cols)
compare = compare[compare$Var1 != compare$Var2,]
compare$comb = paste(t(apply(compare, 1, sort))[,1], t(apply(compare, 1, sort))[,2])
ucombs = unique(compare$comb)
OUT = NULL
for(r in 1:nrow(compare)){
  t = compare[compare$comb==as.character(ucombs[r]),,drop=F]
  t = t[!is.na(t[,1])]
  OUT = rbind(OUT, t[1,])
}
compare = as.data.frame(OUT)

#plot comparisons
ordcol = data.frame(name = c("cox1", "cytb", "mt - whole", "gen - blast"), index = c(10, 16, 22, 34))
pdf("div_genetic.pdf", height=4, width=4, onefile=T)
for(p in 1:nrow(compare)){
  plot(-100, -100, xlim=c(0,0.2), ylim=c(0,0.2), xlab=as.character(ordcol$name[ordcol$index==compare$Var2[p]]), ylab=as.character(ordcol$name[ordcol$index==compare$Var1[p]]))
  segments(0,0,1,1, lty=2, col="grey50")
  points(x=data[,compare$Var2[p]], y=data[,compare$Var1[p]], pch=20, col=alpha("firebrick3", 0.5), cex=1.25)
  abline(summary(lm(data[,compare$Var1[p]]~data[,compare$Var2[p]])), col="grey50")
}
dev.off()
sink("div_genetic_lm.txt")
for(p in 1:nrow(compare)){
  print(paste("y=", as.character(ordcol$name[ordcol$index==compare$Var1[p]]), "; x=", as.character(ordcol$name[ordcol$index==compare$Var2[p]]), sep=""))
  print(summary(lm(data[,compare$Var1[p]]~data[,compare$Var2[p]])))
}
sink()

sink("div_genetic_glm.txt")
for(p in 1:nrow(compare)){
  print(paste("y=", as.character(ordcol$name[ordcol$index==compare$Var1[p]]), "; x=", as.character(ordcol$name[ordcol$index==compare$Var2[p]]), sep=""))
  print(summary(glm(data[,compare$Var1[p]]~data[,compare$Var2[p]])))
}
sink()
