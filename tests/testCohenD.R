### Cohen's d test
library(effsize)
library(tools)

assert <- function(label,condition){
  cat(label,": ")
  if(!condition){
    cat("Failed!\n")
  }else{
    cat("OK.\n")
  }
}

set.seed(52)
x = rnorm(100,mean=10)
y = rnorm(100,mean=12)
d = (c(x,y))
f = rep(c("A","B"),each=100)
eff.d = cohen.d(d,f)

assert("Two samples with large difference", eff.d$conf.int[1] < -2 & -2 < eff.d$conf.int[2]  )


eff.g = cohen.d(d,f,hedges.correction=TRUE)
assert("Two samples with Hedges G", eff.g$conf.int[1] < -2 & -2 < eff.g$conf.int[2]  )


set.seed(54)
d <- rnorm(200)
f <- rep(c(1,2),100)
assertWarning( eff.d <<- cohen.d(d ~ f) )

assert("Two samples from same population",eff.d$conf.int[1] < 0 & 0 < eff.d$conf.int[2]  )

## noncentrality t

delta = c(1.73, 1.06, 2.03, 1.40, 0.95, 1.13, 1.41, 1.73, 1.63, 1.56) - 1
set.seed(50)
a = delta
set.seed(50)
b = delta + runif(10)

eff.d = cohen.d(a,b,paired=TRUE)
assert("Paired measures",abs(eff.d$estimate)-1.42 <0.01 )


set.seed(60)
x = rnorm(10,mean=10)
y = rnorm(10,mean=12)
x[4] <- NA
eff.d = cohen.d(x,y,paired=TRUE,na.rm=TRUE)

assert("Paired measures w/NA",abs(eff.d$estimate)-1.73 <0.01 )

eff.d = cohen.d(c(0.2,0.4,0.8,0.9,NA), c(0.12,0.14,0.18,0.119,NA), paired = T, na.rm = T)

assert("Paired measures multiple w/NA",abs(eff.d$estimate)-1.35 <0.01 )

# eff.dc = cohen.d(a,b,paired=TRUE,noncentral = TRUE)
# assert("Two samples from same population",eff.d$conf.int[1] < 0 & 0 < eff.d$conf.int[2]  )
# 
# 
# set.seed(22)
# a = rnorm(35,24,sqrt(148.87))
# mean(a)
# set.seed(31)
# b = rnorm(29,16.5,sqrt(139.16))
# mean(b)
# cohen.d(a,b,noncentral = TRUE)
# cohen.d(a,b)
# 
