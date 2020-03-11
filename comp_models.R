setwd("/Users/jannawilloughby/GDrive/WGS_divergence/data/")
library(scales)
library(vioplot)
library(phytools)

#read in data
data = read.table("data_feb7.csv", header=T, sep=",")

#update variables
data$meanbd = apply(cbind(data$sp1_bodysize, data$sp2_bodysize), 1, mean)
for(r in 1:nrow(data)){
  if(is.na(data$meanbd[r])){
    if(!is.na(data$sp1_bodysize[r])){
      data$meanbd[r] = data$sp1_bodysize[r]
      next
    }
    if(!is.na(data$sp2_bodysize[r])){
      data$meanbd[r] = data$sp2_bodysize[r]
      next
    }
  }
}
data$meangen = apply(cbind(data$sp1_gen_iucn, data$sp2_gen_iucn), 1, mean)
for(r in 1:nrow(data)){
  if(is.na(data$meangen[r])){
    if(!is.na(data$sp1_gen_iucn[r])){
      data$meangen[r] = data$sp1_gen_iucn[r]
      next
    }
    if(!is.na(data$sp2_gen_iucn[r])){
      data$meangen[r] = data$sp2_gen_iucn[r]
      next
    }
  }
}
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
data$perrange = (data$gArea_Int / (data$gArea_sp1 + data$gArea_sp2))
data$rdist_km = data$gDist/1000

#remove outliers or data not included in regressions
data$perrange[data$perrange==0] = NA
data$rdist_km[data$rdist_km==0] = NA
data$meanbd[data$meanbd > 800000] = NA
data$meangen[data$meangen > 20 & data$class=="Mammalia"] = NA

OUT = NULL
#nuclear genome size - nuclear divergence
classes = c("Mammalia",   "Actinopterygii",  "Aves") #  mammal, bird, fish
for(c in classes){
  t = data[data$class==c,]
  if(c=="Mammalia"){
    t = t[t$gen_K80 > 0,]
  }
  if(c=="Actinopterygii"){
    t = t[t$gen_K80 < 0.08,]
  }
  if(c=="Aves"){
    t = t[t$gen_K80 < 0.025,]
  }
  t = data.frame(y=sqrt(t$gen_K80), x=t$meanc)
  t = t[complete.cases(t),]
  tlm = lm(y~x, data=t)
  tmp = c(c, "genome size", "raw", "nuc", summary(tlm)$coefficients[1,1], summary(tlm)$coefficients[1,2], summary(tlm)$coefficients[2,1], summary(tlm)$coefficients[2,2], summary(tlm)$fstatistic, summary(tlm)$r.squared)
  t$x = scale(t$x, scale=T, center=T)[,1]
  tlm = lm(y~x, data=t)
  tmp = c(tmp, c(c, "genome size", "std", "nuc", summary(tlm)$coefficients[1,1], summary(tlm)$coefficients[1,2], summary(tlm)$coefficients[2,1], summary(tlm)$coefficients[2,2], summary(tlm)$fstatistic, summary(tlm)$r.squared, summary(tlm)$r.squared/(1-summary(tlm)$r.squared)))
  OUT = rbind(OUT, tmp)
}

#mean generation time - nuclear and mt divergence
classes = c("Mammalia",   "Actinopterygii",  "Aves") #mammal, bird, fish
for(c in classes){
  t = data[data$class==c,]
  if(c=="Aves"){
    t = t[t$gen_K80 < 0.025, ]
  }
  t = data.frame(y=sqrt(t$gen_K80), x=t$meangen)
  t = t[complete.cases(t),]
  tlm = lm(y~x, data=t)
  tmp = c(c, "generationtime", "raw", "nuc", summary(tlm)$coefficients[1,1], summary(tlm)$coefficients[1,2], summary(tlm)$coefficients[2,1], summary(tlm)$coefficients[2,2], summary(tlm)$fstatistic, summary(tlm)$r.squared)
  t$x = scale(t$x, scale=T, center=T)[,1]
  tlm = lm(y~x, data=t)
  tmp = c(tmp, c(c, "generationtime", "std", "nuc", summary(tlm)$coefficients[1,1], summary(tlm)$coefficients[1,2], summary(tlm)$coefficients[2,1], summary(tlm)$coefficients[2,2], summary(tlm)$fstatistic, summary(tlm)$r.squared, summary(tlm)$r.squared/(1-summary(tlm)$r.squared)))
  OUT = rbind(OUT, tmp)
  
  t = data[data$class==c,]
  t = data.frame(y=sqrt(t$mtW_K80), x=t$meangen)
  t = t[complete.cases(t),]
  tlm = lm(y~x, data=t)
  tmp = c(c, "generationtime", "raw", "mt", summary(tlm)$coefficients[1,1], summary(tlm)$coefficients[1,2], summary(tlm)$coefficients[2,1], summary(tlm)$coefficients[2,2], summary(tlm)$fstatistic, summary(tlm)$r.squared)
  t$x = scale(t$x, scale=T, center=T)[,1]
  tlm = lm(y~x, data=t)
  tmp = c(tmp, c(c, "generationtime", "std", "mt", summary(tlm)$coefficients[1,1], summary(tlm)$coefficients[1,2], summary(tlm)$coefficients[2,1], summary(tlm)$coefficients[2,2], summary(tlm)$fstatistic, summary(tlm)$r.squared, summary(tlm)$r.squared/(1-summary(tlm)$r.squared)))
  OUT = rbind(OUT, tmp)
}

#mean body size - nuclear and mt divergence
classes = c("Mammalia", "Aves", "Actinopterygii") #mammal, bird
for(c in classes){
  t = data[data$class==c,]
  if(c=="Mammalia"){
    t$meanbd = t$meanbd/1000
  }
  t = data.frame(y=sqrt(t$gen_K80), x=sqrt(t$meanbd))
  t = t[complete.cases(t),]
  tlm = lm(y~x, data=t)
  tmp = c(c, "bodysize", "raw", "nuc", summary(tlm)$coefficients[1,1], summary(tlm)$coefficients[1,2], summary(tlm)$coefficients[2,1], summary(tlm)$coefficients[2,2], summary(tlm)$fstatistic, summary(tlm)$r.squared)
  t$x = scale(t$x, scale=T, center=T)[,1]
  tlm = lm(y~x, data=t)
  tmp = c(tmp, c(c, "bodysize", "std", "nuc", summary(tlm)$coefficients[1,1], summary(tlm)$coefficients[1,2], summary(tlm)$coefficients[2,1], summary(tlm)$coefficients[2,2], summary(tlm)$fstatistic, summary(tlm)$r.squared, summary(tlm)$r.squared/(1-summary(tlm)$r.squared)))
  OUT = rbind(OUT, tmp)
  
  t = data[data$class==c,]
  if(c=="Mammalia"){
    t$meanbd = t$meanbd/1000
  }
  t = data.frame(y=sqrt(t$mtW_K80), x=sqrt(t$meanbd))
  t = t[complete.cases(t),]
  tlm = lm(y~x, data=t)
  tmp = c(c, "bodysize", "raw", "mt", summary(tlm)$coefficients[1,1], summary(tlm)$coefficients[1,2], summary(tlm)$coefficients[2,1], summary(tlm)$coefficients[2,2], summary(tlm)$fstatistic, summary(tlm)$r.squared)
  t$x = scale(t$x, scale=T, center=T)[,1]
  tlm = lm(y~x, data=t)
  tmp = c(tmp, c(c, "bodysize", "std", "mt", summary(tlm)$coefficients[1,1], summary(tlm)$coefficients[1,2], summary(tlm)$coefficients[2,1], summary(tlm)$coefficients[2,2], summary(tlm)$fstatistic, summary(tlm)$r.squared, summary(tlm)$r.squared/(1-summary(tlm)$r.squared)))
  OUT = rbind(OUT, tmp)
}

#mean tes - nuclear 
classes = c("Mammalia", "Aves", "Actinopterygii") #mammal, bird
for(c in classes){
  t = data[data$class==c,]
  t = data.frame(y=sqrt(t$gen_K80), x=t$meante)
  t = t[complete.cases(t),]
  tlm = lm(y~x, data=t)
  tmp = c(c, "te", "raw", "nuc", summary(tlm)$coefficients[1,1], summary(tlm)$coefficients[1,2], summary(tlm)$coefficients[2,1], summary(tlm)$coefficients[2,2], summary(tlm)$fstatistic, summary(tlm)$r.squared)
  t$x = scale(t$x, scale=T, center=T)[,1]
  tlm = lm(y~x, data=t)
  tmp = c(tmp, c(c, "te", "std", "nuc", summary(tlm)$coefficients[1,1], summary(tlm)$coefficients[1,2], summary(tlm)$coefficients[2,1], summary(tlm)$coefficients[2,2], summary(tlm)$fstatistic, summary(tlm)$r.squared, summary(tlm)$r.squared/(1-summary(tlm)$r.squared)))
  OUT = rbind(OUT, tmp)
}

#range distance - nuclear and mt divergence
classes = c("Mammalia", "Aves") #mammal, bird
for(c in classes){
  t = data[data$class==c,]
  t = data.frame(y=sqrt(t$gen_K80), x=log(t$rdist_km))
  t = t[complete.cases(t),]
  tlm = lm(y~x, data=t)
  tmp = c(c, "mindist", "raw", "nuc", summary(tlm)$coefficients[1,1], summary(tlm)$coefficients[1,2], summary(tlm)$coefficients[2,1], summary(tlm)$coefficients[2,2], summary(tlm)$fstatistic, summary(tlm)$r.squared)
  t$x = scale(t$x, scale=T, center=T)[,1]
  tlm = lm(y~x, data=t)
  tmp = c(tmp, c(c, "mindist", "std", "nuc", summary(tlm)$coefficients[1,1], summary(tlm)$coefficients[1,2], summary(tlm)$coefficients[2,1], summary(tlm)$coefficients[2,2], summary(tlm)$fstatistic, summary(tlm)$r.squared, summary(tlm)$r.squared/(1-summary(tlm)$r.squared)))
  OUT = rbind(OUT, tmp)
  
  t = data[data$class==c,]
  t = data.frame(y=sqrt(t$mtW_K80), x=log(t$rdist_km))
  t = t[complete.cases(t),]
  tlm = lm(y~x, data=t)
  tmp = c(c, "mindist", "raw", "mt", summary(tlm)$coefficients[1,1], summary(tlm)$coefficients[1,2], summary(tlm)$coefficients[2,1], summary(tlm)$coefficients[2,2], summary(tlm)$fstatistic, summary(tlm)$r.squared)
  t$x = scale(t$x, scale=T, center=T)[,1]
  tlm = lm(y~x, data=t)
  tmp = c(tmp, c(c, "mindist", "std", "mt", summary(tlm)$coefficients[1,1], summary(tlm)$coefficients[1,2], summary(tlm)$coefficients[2,1], summary(tlm)$coefficients[2,2], summary(tlm)$fstatistic, summary(tlm)$r.squared, summary(tlm)$r.squared/(1-summary(tlm)$r.squared)))
  OUT = rbind(OUT, tmp)
}

#range overlap - nuclear and mt divergence
classes = c("Mammalia", "Aves") #mammal, bird
for(c in classes){
  t = data[data$class==c,]
  t = data.frame(y=sqrt(t$gen_K80), x=t$perrange)
  t = t[complete.cases(t),]
  tlm = lm(y~x, data=t)
  tmp = c(c, "rangeoverlap", "raw", "nuc", summary(tlm)$coefficients[1,1], summary(tlm)$coefficients[1,2], summary(tlm)$coefficients[2,1], summary(tlm)$coefficients[2,2], summary(tlm)$fstatistic, summary(tlm)$r.squared)
  t$x = scale(t$x, scale=T, center=T)[,1]
  tlm = lm(y~x, data=t)
  tmp = c(tmp, c(c, "rangeoverlap", "std", "nuc", summary(tlm)$coefficients[1,1], summary(tlm)$coefficients[1,2], summary(tlm)$coefficients[2,1], summary(tlm)$coefficients[2,2], summary(tlm)$fstatistic, summary(tlm)$r.squared, summary(tlm)$r.squared/(1-summary(tlm)$r.squared)))
  OUT = rbind(OUT, tmp)
  
  t = data[data$class==c,]
  t = data.frame(y=sqrt(t$mtW_K80), x=t$perrange)
  t = t[complete.cases(t),]
  tlm = lm(y~x, data=t)
  tmp = c(c, "rangeoverlap", "raw", "mt", summary(tlm)$coefficients[1,1], summary(tlm)$coefficients[1,2], summary(tlm)$coefficients[2,1], summary(tlm)$coefficients[2,2], summary(tlm)$fstatistic, summary(tlm)$r.squared)
  t$x = scale(t$x, scale=T, center=T)[,1]
  tlm = lm(y~x, data=t)
  tmp = c(tmp, c(c, "rangeoverlap", "std", "mt", summary(tlm)$coefficients[1,1], summary(tlm)$coefficients[1,2], summary(tlm)$coefficients[2,1], summary(tlm)$coefficients[2,2], summary(tlm)$fstatistic, summary(tlm)$r.squared, summary(tlm)$r.squared/(1-summary(tlm)$r.squared)))
  OUT = rbind(OUT, tmp)
}

colnames(OUT) = c("class", "x", "regtype", "y", "intM", "intSE", "slopeM", "slopeSE", "f", "df1", "df2", "r2","class", "x", "regtype", "y", "intM", "intSE", "slopeM", "slopeSE", "f", "df1", "df2", "r2", "cohens")
write.table(OUT, "statscompare.csv", sep=",", row.names=F, col.names=T)
