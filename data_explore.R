setwd("/Users/jannawilloughby/GDrive/WGS_divergence/data/explore/")
library(scales)

#read in data
data = read.table("data_may23.csv", header=T, sep=",")

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
