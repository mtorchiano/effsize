## Issue #49
library(tidyverse)
# Paired samples effect size
df <- data.frame(
  id = 1:5,
  pre  = c(110, 122, 101, 120, 140),
  post = c(150, 160, 110, 140, 155)
)
df <- df %>% gather(key = "treatment", value = "value", -id)

set.seed(1234)
(df %>% rstatix::cohens_d(value ~ treatment, paired = T))$effsize
#1.75

lsr::cohensD(value ~ treatment,data = df,method = "paired")
#1.75

effsize::cohen.d(value ~ treatment, data=df, paired=T,hedges.correction=F)$estimate
#1.32

effsize::cohen.d(value ~ treatment | Subject(id), data=df, paired=T)

effsize::cohen.d(value ~ treatment | Subject(id), data=df, paired=T, within=FALSE)

#### PROBLEM SOLUTION

pairdiff = diff(df$value,lag=dim(df)[1]/2)
effsize::cohen.d(pairdiff,NA, data=df, paired=T,hedges.correction=F)$estimate
effsize::cohen.d(pairdiff~., paired=T,hedges.correction=F)$estimate

effsize::cohen.d(diff(value,lag=length(value)/2)~.,data=df, paired=T,hedges.correction=F)$estimate
effsize::cohen.d(diff(value,lag=dim(df)[1]/2)~.,data=df, paired=T,hedges.correction=F)$estimate
