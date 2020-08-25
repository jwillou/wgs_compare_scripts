setwd("/Users/jannawilloughby/GDrive/WGS_divergence/data/")

#read in data
data = read.table("data_feb7.csv", header=T, sep=",")

ALLCOR = NULL
#nuclear
df = data[,35:39]
df = df[complete.cases(df),]
cor(df)
t = cor(df)
diag(t) = NA
ALLCOR = rbind(ALLCOR,t)

#mt
df = data[,29:33]
df = df[complete.cases(df),]
cor(df)
t = cor(df)
diag(t) = NA
ALLCOR = rbind(ALLCOR,t)

#cytb
df = data[,17:21]
df = df[complete.cases(df),]
cor(df)
t = cor(df)
diag(t) = NA
ALLCOR = rbind(ALLCOR,t)

#co1
df = data[,11:15]
df = df[complete.cases(df),]
cor(df)
t = cor(df)
diag(t) = NA
ALLCOR = rbind(ALLCOR,t)

range(ALLCOR, na.rm=T)
