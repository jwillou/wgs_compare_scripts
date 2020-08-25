setwd("/Users/jannawilloughby/GDrive/WGS_divergence/data/")

#read in data
data = read.table("data_feb7.csv", header=T, sep=",")

####update variables####
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
data$perrange[data$perrange==0] = NA
data$rdist_km = data$gDist/1000
data$rdist_km[data$rdist_km==0] = NA

#remove outliers or data not included in regressions
data$perrange[data$perrange==0] = NA
data$rdist_km[data$rdist_km==0] = NA
data$meanbd[data$meanbd > 800000] = NA
data$meangen[data$meangen > 20 & data$class=="Mammalia"] = NA

####predictor correlations####
sdata = cbind(data[,1:7], data.frame(co1_K80=data$co1_K80, cytb_K80=data$cytb_K80, mtW_K80=data$mtW_K80, gen_K80=data$gen_K80, meanbd=data$meanbd, meangen=data$meangen, meanc=data$meanc, meante=data$meante, perrange=data$perrange, rdist_km=data$rdist_km))
sdata = as.data.frame(sdata)

####pairwise correlations between predictors####
sink(file="/Users/jannawilloughby/GDrive/WGS_divergence/data/correlations.txt")
preds = c("meanbd", "meangen", "meanc", "meante", "perrange", "rdist_km")
corrs = NULL
for(p in 1:length(preds)){
  if(p==length(preds)){next}
  t = sdata[,colnames(sdata)==preds[p]]
  for(pp in (p+1):length(preds)){
    q = sdata[,colnames(sdata)==preds[pp]]
    df = data.frame(t=t, q=q)
    df = df[complete.cases(df),]
    print(paste("correlation:", preds[p], ", ", preds[pp], sep=""))
    print(cor(df,method="spearman")[1,2])
    corrs = rbind(corrs, c(preds[p],  preds[pp], cor(df,method="spearman")[1,2]))
  }
}
sink()
colnames(corrs) = c("pred1", "pred2", "spearman")
corrs = as.data.frame(corrs)
corrs$spearman = as.numeric(as.character(corrs$spearman))
corrs
range(corrs$spearman, na.rm=T)

####sign flip####
lm(sdata$gen_K80~meanbd+)

