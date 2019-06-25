setwd("/Users/jannawilloughby/GDrive/WGS_divergence/data/")
library(scales)
library(vioplot)
library(phytools)

#read in data
data = read.table("data_june19.csv", header=T, sep=",")

####genome size####
plot(data$sp1_assembly_length, data$sp2_assembly_length)
plot(data$sp1_cvalue, data$sp2_cvalue)
plot(data$sp2_cvalue, data$sp2_assembly_length)
plot(data$sp1_cvalue, data$sp1_assembly_length)
