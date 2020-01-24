cs = unique(data$class)
m = g= NULL
for(c in cs[2:4]){
  t=data[data$class==as.character(c),]
  g = c(g, sd(t$mtW_K80, na.rm=T))
  m = c(m, sd(t$gen_K80, na.rm=T))
}

#body size rates of change
strt = 100
inc = 100
#fishes - mt
((0.322 + -0.009*(sqrt(strt)))-(0.322 + -0.009*(sqrt((strt+inc)))))^2
#fishes - nuc
((0.196 + -0.005*(sqrt(strt)))-(0.196 + -0.005*(sqrt((strt+inc)))))^2

#birds - mt
((0.237 + -0.000*(sqrt(strt)))-(0.237 + -0.000*(sqrt((strt+inc)))))^2
#birds - nuc
((0.133 + -0.001*(sqrt(strt)))-(0.133 + -0.001*(sqrt((strt+inc)))))^2

#mammals - mt
((0.317 + -0.002*(sqrt(strt)))-(0.317 + -0.002*(sqrt((strt+inc)))))^2
#birds - nuc
((1.613 + -0.001*(sqrt(strt)))-(1.613 + -0.001*(sqrt((strt+inc)))))^2


x = cbind( data$meanbd, data$gen_K80)
x = x[complete.cases(x),]
cor(x, method="spearman")
