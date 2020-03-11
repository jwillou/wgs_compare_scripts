setwd("/Users/jannawilloughby/GDrive/WGS_divergence/data/recip/")
d1 = read.table("mtDNAdiv_estimates.csv", sep=",", header=T)
d2 = read.table("mtDNAdiv_estimates12.csv", sep=",", header=T)

td1 = data.frame(spp=d1$spp, K80=d1$K80)
td2 = data.frame(spp=d2$NA., K80=d2$K80)

td = merge(x=td1, y=td2, by="spp")
cor(td$K80.x, td$K80.y, method="spearman")
