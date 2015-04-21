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

