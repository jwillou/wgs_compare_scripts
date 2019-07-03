setwd("/Users/jannawilloughby/GDrive/WGS_divergence/data/family/")
library(scales)
library(vioplot)
library(phytools)

#read in data
data = read.table("../data_july02.csv", header=T, sep=",")

####class####
pdf("family/classcompare.pdf", height=5, width=5, onefile=T)
classes = unique(data$class)
OUT = NULL
for(c in 1:length(classes)){
  t = data[data$class==as.character(classes[c]),,drop=F]
  d = c(mean(t$co1_K80,  na.rm=T), sd(t$co1_K80,  na.rm=T), nrow(t[!is.na(t$co1_K80),, drop=F]), (sd(t$co1_K80,  na.rm=T)/(sqrt(nrow(t[!is.na(t$co1_K80),, drop=F])))), 
        mean(t$cytb_K80, na.rm=T), sd(t$cytb_K80, na.rm=T), nrow(t[!is.na(t$cytb_K80),,drop=F]), (sd(t$cytb_K80, na.rm=T)/(sqrt(nrow(t[!is.na(t$cytb_K80),,drop=F])))), 
        mean(t$mtW_K80, na.rm=T),  sd(t$mtW_K80,  na.rm=T), nrow(t[!is.na(t$mtW_K80),, drop=F]), (sd(t$mtW_K80,  na.rm=T)/(sqrt(nrow(t[!is.na(t$mtW_K80),, drop=F])))), 
        mean(t$gen_K80, na.rm=T),  sd(t$gen_K80,  na.rm=T), nrow(t[!is.na(t$gen_K80),, drop=F]), (sd(t$gen_K80,  na.rm=T)/(sqrt(nrow(t[!is.na(t$gen_K80),, drop=F])))), 
        as.character(classes[c]))
  OUT = rbind(OUT, d)
}
class.sum = as.data.frame(OUT)
colnames(class.sum) = c("co1M", "co1SD", "co1N", "co1SE", "cytbM", "cytbSD", "cytbN", "cytbSE,", "mtwM", "mtwSD", "mtwN", "mtwSE", "genM", "genSD", "genN", "genSE")
rownames(class.sum) = seq(1, nrow(class.sum), 1)

#fix factors
for(c in 1:16){
  class.sum[,c] = as.numeric(as.character(class.sum[,c]))
}

Mcols = c(1, 5, 9, 13)
labs  = c("co1", "cytb", "mt", "genome")
l = 0

for(c in Mcols){
  l = l + 1
  plot(-100,-100, xlim=c(1.65,4.40), ylim=c(0,0.15), xlab="Class", ylab="divergence (Kimura 1980)", axes=T, main=labs[l])
  colors5 = c("saddlebrown", "goldenrod3", "dodgerblue3", "firebrick3", "darkorchid3")
  b = 0.35
  o = 0.1
  for(r in 2:4){
    rect(xleft=c((r-b),(r-b)), xright=c((r+b),(r+b)), ybottom=c(0,0), ytop=c(class.sum[r,c],class.sum[r,c]), col=alpha(colors5[r], 0.5))
    segments(x0=r,   y0=(class.sum[r,c]-class.sum[r,(c+3)]), x1=r,   y1=(class.sum[r,c]+class.sum[r,(c+3)]), lwd=2) #vertical
    segments(x0=r-o, y0=(class.sum[r,c]-class.sum[r,(c+3)]), x1=r+o, y1=(class.sum[r,c]-class.sum[r,(c+3)]), lwd=2) #bottom
    segments(x0=r-o, y0=(class.sum[r,c]+class.sum[r,(c+3)]), x1=r+o, y1=(class.sum[r,c]+class.sum[r,(c+3)]), lwd=2) #top
  }
}
write.table(class.sum, "class/classcompare.csv", sep=",", row.names=F, col.names=T)

#by class with raw/density data
Mcols = c(12, 18, 24, 36)
labs    = c("co1", "cytb", "mt", "genome")
colors6 = c("goldenrod3", "dodgerblue3",    "firebrick3") #"saddlebrown", , "chartreuse3",    "darkorchid3"
classes = c("Aves",       "Actinopterygii", "Mammalia") #"Reptilia",    ,   "Chondrichthyes", "Amphibia"
l = 0
for(c in Mcols){
  l = l + 1
  t = data[!is.na(data[,Mcols]),]
  t = t[!is.na(t[,1]),]
  if(labs[l]=="genome"){
    plot(-100,-100, xlim=c(0.65,(length(classes)+0.4)), ylim=c(0,0.15), xlab="Class", ylab="divergence rate (Kimura 1980)", axes=T, main=labs[l])
  }else{
    plot(-100,-100, xlim=c(0.65,(length(classes)+0.4)), ylim=c(0,0.3), xlab="Class", ylab="divergence rate (Kimura 1980)", axes=T, main=labs[l])
    
  }
  b = 0.35
  o = 0.05
  j = 0.05
  c6 = 0
  for(g in classes){
    c6 = c6 + 1
    tt = t[t$class==as.character(g),,drop=F]
    tt = tt[!is.na(tt[,Mcols]),,drop=F]
    points(x=(c6-sample(seq(j, 0.2, 0.01), nrow(tt), replace=T)), y=tt[,c], col=alpha(colors6[c6], 0.7), pch=19, cex=0.75)
    vioplot(x=tt[,c], col = alpha(colors6[c6], 0.9), plotCentre = "line", side = "right", ylim=c(0,0.3), add=T, at=(c6+j))
  }
}
#all classes combined
l = 0
for(c in Mcols){
  l = l + 1
  t = data[!is.na(data[,Mcols]),]
  t = t[!is.na(t[,1]),]
  if(labs[l]=="genome"){
    plot(-100,-100, xlim=c(0.65,(length(classes)+0.4)), ylim=c(0,0.15), xlab="Class", ylab="divergence rate (Kimura 1980)", axes=T, main=labs[l])
  }else{
    plot(-100,-100, xlim=c(0.65,(length(classes)+0.4)), ylim=c(0,0.3), xlab="Class", ylab="divergence rate (Kimura 1980)", axes=T, main=labs[l])
    
  }
  b = 0.35
  o = 0.05
  j = 0.05
  c6 = 0
  g  
  c6 = c6 + 1
  tt = t
  tt = tt[!is.na(tt[,Mcols]),,drop=F]
  points(x=(c6-sample(seq(j, 0.2, 0.01), nrow(tt), replace=T)), y=tt[,c], col=alpha("grey50", 0.7), pch=19, cex=0.75)
  vioplot(x=tt[,c], col = alpha("grey50", 0.9), plotCentre = "line", side = "right", ylim=c(0,0.3), add=T, at=(c6+j))
}
dev.off()

####family####
colors6 = c("saddlebrown", "goldenrod3", "dodgerblue3",    "firebrick3", "chartreuse3",    "darkorchid3") #
classes = c("Reptilia",    "Aves",       "Actinopterygii", "Mammalia",   "Chondrichthyes", "Amphibia") #    

#create tree for next plot
allfams = names(table(data$family))
pdata = NULL
for(f in allfams){
  t = data[data$family==as.character(f),,drop=F]
  t = t[!is.na(t$family),]
  pdata = rbind(pdata, t[1,])
}
as.phylo.formula2 = function (x, data = parent.frame(), ...){
  err <- "Formula must be of the kind \"~A1/A2/.../An\"."
  if (length(x) != 2) 
    stop(err)
  if (x[[1]] != "~") 
    stop(err)
  f <- x[[2]]
  taxo <- list()
  while (length(f) == 3) {
    if (f[[1]] != "/") 
      stop(err)
    if (!is.factor(data[[deparse(f[[3]])]])) 
      stop(paste("Variable", deparse(f[[3]]), "must be a factor."))
    taxo[[deparse(f[[3]])]] <- data[[deparse(f[[3]])]]
    if (length(f) > 1) 
      f <- f[[2]]
  }
  if (!is.factor(data[[deparse(f)]])) 
    stop(paste("Variable", deparse(f), "must be a factor."))
  taxo[[deparse(f)]] <- data[[deparse(f)]]
  taxo.data <- as.data.frame(taxo)
  leaves.names <- as.character(taxo.data[, 1])
  taxo.data[, 1] <- 1:nrow(taxo.data)
  f.rec <- function(subtaxo) {
    u <- ncol(subtaxo)
    levels <- unique(subtaxo[, u])
    if (u == 1) {
      if (length(levels) != nrow(subtaxo)) 
        warning("Error, leaves names are not unique.")
      return(as.character(subtaxo[, 1]))
    }
    t <- character(length(levels))
    for (l in 1:length(levels)) {
      x <- f.rec(subtaxo[subtaxo[, u] == levels[l], ][1:(u - 1)])
      t[l] <- paste("(", paste(x, collapse = ","), ")", sep = "")
    }
    return(t)
  }
  string <- paste("(", paste(f.rec(taxo.data), collapse = ","),");", sep = "")
  phy <- read.newick(text = string) ## so that singles will be read without error
  phy$edge.length <- rep(1,nrow(phy$edge))
  phy <- collapse.singles(phy)
  phy$tip.label <- leaves.names[as.numeric(phy$tip.label)]
  return(phy)
}
atree  = as.phylo.formula2(~class/order/family, data=pdata)
pdf("family_tree.pdf", height=10, width=5)
par(bg=NA)
plot(atree, cex=0.75)
dev.off()
famorder = rev(atree$tip.label)

#set up color order
classcol = data.frame(class=classes, color=colors6)
colorder = NULL
for(f in famorder){
  t = data[data$family==as.character(f),]
  t = t$class[1]
  colorder = c(colorder, as.character(classcol$color[classcol$class==as.character(t)]))
}

#genomic divergence
Mcols   = 36 #(12, 18, 24, 36)
#set up plot
pdf("family_divergence_g.pdf", height=6, width=8)
par(bg=NA)
plot(-100,-100, xlim=c(1,59), ylim=c(-0.10,0.155), xlab="Family", ylab="divergence (Kimura 1980)", axes=F)
axis(side=1, at=seq(1,59,1), pos=-0.005, labels=F, lwd=0.75)
segments(x0=0, x1=60,y0=-0.005, y1=-0.005, lwd=1)
axis(side=2, at=seq(0,0.15,0.05), pos=0, labels=T, lwd=0.75)
segments(x0=0, x1=0,y0=-0.005, y1=0.155, lwd=1)
segments(x0=0, x1=60,y0=0.155, y1=0.155, lwd=1)
segments(x0=60, x1=60,y0=0.155, y1=-0.005, lwd=1)

#find cutoff value, add to plot
all = mean(c(data[,Mcols]), na.rm=T) + sd(c(data[,Mcols]), na.rm=T)
segments(x0=1,x1=59,y0=all,y1=all, lty=2, col="black", lwd=2)
#mamm = mean(c(data[data$class=="Mammalia",Mcols]), na.rm=T) + sd(c(data[data$class=="Mammalia",Mcols]), na.rm=T)
#segments(x0=3,x1=28,y0=mamm,y1=mamm, lty=1, col="firebrick3", lwd=2)
#rept = mean(c(data[data$class=="Reptilia",Mcols]), na.rm=T) + sd(c(data[data$class=="Reptilia",Mcols]), na.rm=T)
#segments(x0=29,x1=32,y0=rept,y1=rept, lty=1, col="saddlebrown", lwd=2)
#fish = mean(c(data[data$class=="Actinopterygii",Mcols]), na.rm=T) + sd(c(data[data$class=="Actinopterygii",Mcols]), na.rm=T)
#segments(x0=33,x1=47,y0=fish,y1=fish, lty=1, col="dodgerblue3", lwd=2)
#bird = mean(c(data[data$class=="Aves",Mcols]), na.rm=T) + sd(c(data[data$class=="Aves",Mcols]), na.rm=T)
#segments(x0=48,x1=59,y0=bird,y1=bird, lty=1, col="goldenrod3", lwd=2)

#add polygon shading
for(i in seq(1, 60, 2)){
  polygon(x=c((i-0.5), (i+0.5), (i+0.5), (i-0.5)), y=c(-0.10,-0.10,0.155,0.155), col=alpha("grey50", 0.25), border=F)
}

#plot data by family/class
ci = 0
b = 0.35
o = 0.05
j = 0.005
c6 = 0
for(f in famorder){
  ci = ci + 1
  c6 = c6 + 1
  tt = data[data$family==as.character(f),,drop=F]
  tt = tt[!is.na(tt[,Mcols]),,drop=F]
  points(x=(c6-sample(seq(j, 0.2, 0.01), nrow(tt), replace=T)), y=tt[,Mcols], col=alpha(colorder[ci], 1), pch=19, cex=0.75) #
}
dev.off()

#summarize by class with half violin plots
pdf("family_violin_g.pdf", height=8, width=3)
par(bg=NA)
plot(-100, -100, xlim=c(0.5, 3), ylim=c(0,0.15), axes=T)
tt = data[data$class=="Mammalia",,drop=F]
vioplot(x=tt[,Mcols], col = alpha("firebrick3", 1), plotCentre = "line", side = "right", ylim=c(-0.10,0.155), add=T, at=0.5)
tt = data[data$class=="Reptilia",,drop=F]
vioplot(x=tt[,Mcols], col = alpha("saddlebrown", 1), plotCentre = "line", side = "right", ylim=c(-0.10,0.155), add=T, at=1)
tt = data[data$class=="Actinopterygii",,drop=F]
vioplot(x=tt[,Mcols], col = alpha("dodgerblue3", 1), plotCentre = "line", side = "right", ylim=c(-0.10,0.155), add=T, at=1.5)
tt = data[data$class=="Aves",,drop=F]
vioplot(x=tt[,Mcols], col = alpha("goldenrod3", 1), plotCentre = "line", side = "right", ylim=c(-0.10,0.155), add=T, at=2)
tt = data
vioplot(x=tt[,Mcols], col = alpha("grey50", 1), plotCentre = "line", side = "right", ylim=c(-0.10,0.155), add=T, at=2.5)
dev.off()

#mtdna divergence
Mcols   = 24
#set up plot
pdf("family_divergence_m.pdf", height=5, width=8)
par(bg=NA)
plot(-100,-100, xlim=c(1,59), ylim=c(-0.10,0.33), xlab="Family", ylab="divergence rate (Kimura 1980)", axes=F)
axis(side=1, at=seq(1,59,1), pos=-0.005, labels=F, lwd=0.75)
segments(x0=0, x1=60,y0=-0.005, y1=-0.005, lwd=1)
axis(side=2, at=seq(0,0.3,0.05), pos=0, labels=T, lwd=0.75)
segments(x0=0, x1=0,y0=-0.005, y1=0.33, lwd=1)
segments(x0=0, x1=60,y0=0.33, y1=0.33, lwd=1)
segments(x0=60, x1=60,y0=0.33, y1=-0.005, lwd=1)

#find cutoff value, add to plot
all = mean(c(data[,Mcols]), na.rm=T) + sd(c(data[,Mcols]), na.rm=T)
segments(x0=1,x1=59,y0=all,y1=all, lty=2, col="black", lwd=2)

#add polygon shading
for(i in seq(1, 60, 2)){
  polygon(x=c((i-0.5), (i+0.5), (i+0.5), (i-0.5)), y=c(-0.10,-0.10,0.33,0.33), col=alpha("grey50", 0.25), border=F)
}

#plot data by family/class
ci = 0
b = 0.35
o = 0.05
j = 0.005
c6 = 0
for(f in famorder){
  ci = ci + 1
  c6 = c6 + 1
  tt = data[data$family==as.character(f),,drop=F]
  tt = tt[!is.na(tt[,Mcols]),,drop=F]
  points(x=(c6-sample(seq(j, 0.2, 0.01), nrow(tt), replace=T)), y=tt[,Mcols], col=alpha(colorder[ci], 1), pch=19, cex=0.75) #
}
dev.off()

#summarize by class with half violin plots
pdf("family_violin_m.pdf", height=8, width=3)
par(bg=NA)
plot(-100, -100, xlim=c(0.5, 3), ylim=c(0,0.3), axes=T)
tt = data[data$class=="Mammalia",,drop=F]
vioplot(x=tt[,Mcols], col = alpha("firebrick3", 1), plotCentre = "line", side = "right", ylim=c(-0.10,0.155), add=T, at=0.5)
tt = data[data$class=="Reptilia",,drop=F]
vioplot(x=tt[,Mcols], col = alpha("saddlebrown", 1), plotCentre = "line", side = "right", ylim=c(-0.10,0.155), add=T, at=1)
tt = data[data$class=="Actinopterygii",,drop=F]
vioplot(x=tt[,Mcols], col = alpha("dodgerblue3", 1), plotCentre = "line", side = "right", ylim=c(-0.10,0.155), add=T, at=1.5)
tt = data[data$class=="Aves",,drop=F]
vioplot(x=tt[,Mcols], col = alpha("goldenrod3", 1), plotCentre = "line", side = "right", ylim=c(-0.10,0.155), add=T, at=2)
tt = data
vioplot(x=tt[,Mcols], col = alpha("grey50", 1), plotCentre = "line", side = "right", ylim=c(-0.10,0.155), add=T, at=2.5)
dev.off()

#coi divergence
Mcols   = 12
#set up plot
pdf("family_divergence_coi.pdf", height=5, width=8)
par(bg=NA)
plot(-100,-100, xlim=c(1,59), ylim=c(-0.10,0.33), xlab="Family", ylab="divergence rate (Kimura 1980)", axes=F)
axis(side=1, at=seq(1,59,1), pos=-0.005, labels=F, lwd=0.75)
segments(x0=0, x1=60,y0=-0.005, y1=-0.005, lwd=1)
axis(side=2, at=seq(0,0.3,0.05), pos=0, labels=T, lwd=0.75)
segments(x0=0, x1=0,y0=-0.005, y1=0.33, lwd=1)
segments(x0=0, x1=60,y0=0.33, y1=0.33, lwd=1)
segments(x0=60, x1=60,y0=0.33, y1=-0.005, lwd=1)

#find cutoff value, add to plot
all = mean(c(data[,Mcols]), na.rm=T) + sd(c(data[,Mcols]), na.rm=T)
segments(x0=1,x1=59,y0=all,y1=all, lty=2, col="black", lwd=2)

#add polygon shading
for(i in seq(1, 60, 2)){
  polygon(x=c((i-0.5), (i+0.5), (i+0.5), (i-0.5)), y=c(-0.10,-0.10,0.33,0.33), col=alpha("grey50", 0.25), border=F)
}

#plot data by family/class
ci = 0
b = 0.35
o = 0.05
j = 0.005
c6 = 0
for(f in famorder){
  ci = ci + 1
  c6 = c6 + 1
  tt = data[data$family==as.character(f),,drop=F]
  tt = tt[!is.na(tt[,Mcols]),,drop=F]
  points(x=(c6-sample(seq(j, 0.2, 0.01), nrow(tt), replace=T)), y=tt[,Mcols], col=alpha(colorder[ci], 1), pch=19, cex=0.75) #
}
dev.off()

#summarize by class with half violin plots
pdf("family_violin_coi.pdf", height=8, width=3)
par(bg=NA)
plot(-100, -100, xlim=c(0.5, 3), ylim=c(0,0.3), axes=T)
tt = data[data$class=="Mammalia",,drop=F]
vioplot(x=tt[,Mcols], col = alpha("firebrick3", 1), plotCentre = "line", side = "right", ylim=c(-0.10,0.155), add=T, at=0.5)
tt = data[data$class=="Reptilia",,drop=F]
vioplot(x=tt[,Mcols], col = alpha("saddlebrown", 1), plotCentre = "line", side = "right", ylim=c(-0.10,0.155), add=T, at=1)
tt = data[data$class=="Actinopterygii",,drop=F]
vioplot(x=tt[,Mcols], col = alpha("dodgerblue3", 1), plotCentre = "line", side = "right", ylim=c(-0.10,0.155), add=T, at=1.5)
tt = data[data$class=="Aves",,drop=F]
vioplot(x=tt[,Mcols], col = alpha("goldenrod3", 1), plotCentre = "line", side = "right", ylim=c(-0.10,0.155), add=T, at=2)
tt = data
vioplot(x=tt[,Mcols], col = alpha("grey50", 1), plotCentre = "line", side = "right", ylim=c(-0.10,0.155), add=T, at=2.5)
dev.off()

#cytb divergence
Mcols   = 18
#set up plot
pdf("family_divergence_cytb.pdf", height=5, width=8)
par(bg=NA)
plot(-100,-100, xlim=c(1,59), ylim=c(-0.10,0.33), xlab="Family", ylab="divergence rate (Kimura 1980)", axes=F)
axis(side=1, at=seq(1,59,1), pos=-0.005, labels=F, lwd=0.75)
segments(x0=0, x1=60,y0=-0.005, y1=-0.005, lwd=1)
axis(side=2, at=seq(0,0.3,0.05), pos=0, labels=T, lwd=0.75)
segments(x0=0, x1=0,y0=-0.005, y1=0.33, lwd=1)
segments(x0=0, x1=60,y0=0.33, y1=0.33, lwd=1)
segments(x0=60, x1=60,y0=0.33, y1=-0.005, lwd=1)

#find cutoff value, add to plot
all = mean(c(data[,Mcols]), na.rm=T) + sd(c(data[,Mcols]), na.rm=T)
segments(x0=1,x1=59,y0=all,y1=all, lty=2, col="black", lwd=2)

#add polygon shading
for(i in seq(1, 60, 2)){
  polygon(x=c((i-0.5), (i+0.5), (i+0.5), (i-0.5)), y=c(-0.10,-0.10,0.33,0.33), col=alpha("grey50", 0.25), border=F)
}

#plot data by family/class
ci = 0
b = 0.35
o = 0.05
j = 0.005
c6 = 0
for(f in famorder){
  ci = ci + 1
  c6 = c6 + 1
  tt = data[data$family==as.character(f),,drop=F]
  tt = tt[!is.na(tt[,Mcols]),,drop=F]
  points(x=(c6-sample(seq(j, 0.2, 0.01), nrow(tt), replace=T)), y=tt[,Mcols], col=alpha(colorder[ci], 1), pch=19, cex=0.75) #
}
dev.off()

#summarize by class with half violin plots
pdf("family_violin_cytb.pdf", height=8, width=3)
par(bg=NA)
plot(-100, -100, xlim=c(0.5, 3), ylim=c(0,0.3), axes=T)
tt = data[data$class=="Mammalia",,drop=F]
vioplot(x=tt[,Mcols], col = alpha("firebrick3", 1), plotCentre = "line", side = "right", ylim=c(-0.10,0.155), add=T, at=0.5)
tt = data[data$class=="Reptilia",,drop=F]
vioplot(x=tt[,Mcols], col = alpha("saddlebrown", 1), plotCentre = "line", side = "right", ylim=c(-0.10,0.155), add=T, at=1)
tt = data[data$class=="Actinopterygii",,drop=F]
vioplot(x=tt[,Mcols], col = alpha("dodgerblue3", 1), plotCentre = "line", side = "right", ylim=c(-0.10,0.155), add=T, at=1.5)
tt = data[data$class=="Aves",,drop=F]
vioplot(x=tt[,Mcols], col = alpha("goldenrod3", 1), plotCentre = "line", side = "right", ylim=c(-0.10,0.155), add=T, at=2)
tt = data
vioplot(x=tt[,Mcols], col = alpha("grey50", 1), plotCentre = "line", side = "right", ylim=c(-0.10,0.155), add=T, at=2.5)
dev.off()


