setwd("/Users/jannawilloughby/GDrive/WGS_divergence/data/nuisance")
library(scales)
library(vioplot)
library(phytools)

#read in data
data = read.table("../data_july02.csv", header=T, sep=",")####assembly size####

####assembly quality####
data$assembly = rep(0, nrow(data))
for(r in 1:nrow(data)){
  t1 = as.character(data$sp1_assembly[r])
  t2 = as.character(data$sp2_assembly[r])
  if(t1=="Chromosome" & t2=="Chromosome"){
    data$assembly[r] = "CC"
  }
  if(t1=="Chromosome" & t2=="Scaffold" | t1=="Scaffold" & t2=="Chromosome" ){
    data$assembly[r] = "CS"
  }
  if(t1=="Scaffold" & t2=="Scaffold"){
    data$assembly[r] = "SS"
  }
  if(t1=="Chromosome" & t2=="Contig" | t1=="Contig" & t2=="Chromosome" ){
    data$assembly[r] = "Cc"
  }
  if(t1=="Contig" & t2=="Scaffold" | t1=="Scaffold" & t2=="Contig" ){
    data$assembly[r] = "cS"
  }
  if(t1=="Contig" & t2=="Contig"){
    data$assembly[r] = "cc"
  }
}
sink("assemblyqual.txt")
print(table(data$assembly))

print("CC")
hist(data$gen_K80[data$assembly=="CC"], breaks=seq(0,0.15, 0.01), xlim=c(0,0.15))
print(mean(data$gen_K80[data$assembly=="CC"], na.rm=T))
print(quantile(data$gen_K80[data$assembly=="CC"], probs=c(0.085, 0.915), na.rm=T))

print("CS")
hist(data$gen_K80[data$assembly=="CS"], breaks=seq(0,0.15, 0.01), xlim=c(0,0.15))
print(mean(data$gen_K80[data$assembly=="CS"], na.rm=T))
print(quantile(data$gen_K80[data$assembly=="CS"], probs=c(0.085, 0.915), na.rm=T))

print("SS")
hist(data$gen_K80[data$assembly=="SS"], breaks=seq(0,0.15, 0.01), xlim=c(0,0.15))
print(mean(data$gen_K80[data$assembly=="SS"], na.rm=T))
print(quantile(data$gen_K80[data$assembly=="SS"], probs=c(0.085, 0.915), na.rm=T))

print("cS")
hist(data$gen_K80[data$assembly=="cS"], breaks=seq(0,0.15, 0.01), xlim=c(0,0.15))
print(mean(data$gen_K80[data$assembly=="cS"], na.rm=T))
print(quantile(data$gen_K80[data$assembly=="cS"], probs=c(0.085, 0.915), na.rm=T))
sink()

####number of chromosomes####
data$meannc = apply(cbind(data$sp1_chrom_num, data$sp2_chrom_num), 1, mean)


sink("numberchroms.txt")
print("all")
print(summary(lm(gen_K80~meannc, data=data)))

#classes
classes = c("Mammalia",   "Reptilia",    "Actinopterygii", "Aves")
for(c in classes){
  print(c)
  print(summary(lm(gen_K80~meannc, data=data[data$class==c,])))
}
sink()





