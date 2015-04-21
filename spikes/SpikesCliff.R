bsearch5 <- 
  function(val, tab, L=1L, H=length(tab)) 
  { 
    n <- length(val)
    b <- cbind(L=rep(L, n), H=rep(H, n)) 
    i0 <- seq_along(val) 
    repeat { 
      M <- (b[i0,"H"] + b[i0,"L"]) %/% 2L 
      i <- tab[M] > val[i0] 
      b[i0 + i * n] <- 
        ifelse(i, M - 1L, ifelse(tab[M] < val[i0], M + 1L, M)) 
      i0 <- which(b[i0, "H"] >= b[i0, "L"]) 
      if (!length(i0)) break; 
    } 
    b[,"L"] - 1L 
  } 

a = c (1,3,4,5,6)

bsearch5(3,a)


## spike

n = 100000
b <- cbind(L=rep(L, n), H=rep(H, n))

system.time((b[i0,"H"] + b[i0,"L"]) %/% 2L) ## fastest!!!

system.time(apply(b,1,sum) %/% 2)

system.time(apply(b,1,function(r){sum(r) %/% 2}) )

## spike

x= as.integer(rnorm(10000000,10,2))
y= as.integer(rnorm(10000000,11,2))

time.asint = c()
time.trunc = c()
time.div = c()
for(i in 1:100){
  time.div = c(time.div,system.time( (x+y) %/% 2L )["user.self"])
  time.trunc = c(time.trunc,system.time( trunc((x+y)/2) )["user.self"])
  time.asint = c(time.asint,system.time( as.integer((x+y)/2) )["user.self"])
}

library(car)
qqPlot(time.div)
shapiro.test(time.div)
shapiro.test(time.asint)
shapiro.test(time.trunc)

par(mar=c(3,3,.2,.1))
boxplot(list("%/%"=time.div,"as.integer"=time.asint,"trunc"=time.trunc),
        ylim=c(0,0.12),log="y")
boxplot(list("%/%"=time.div,"as.integer"=time.asint,"trunc"=time.trunc))

wilcox.test(time.div,time.trunk)

mean.asint = mean(time.asint)
mean.div = mean(time.div)
mean.trunc = mean(time.trunc)

(mean.asint - mean.div) / mean.div * 100
(mean.div - mean.trunc) / mean.trunc * 100

## spike call

showFunction <- function(f){
  print(f)
  if(is.function(f)){
    c = sys.call() ## get the call
    cat(as.character(c[[2]]),"-> "); #find the original name
    print(f)
  }
}

cl = showFunction(partition)
cl = showFunction(showFunction)
cl = ff(x)