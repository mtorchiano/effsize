## Test for Cliff's delta
library(effsize)

context("Cohen d")


test_that("Small differences", {
  treatment = c(10,10,20,20,20,30,30,30,40,50)
  control = c(10,20,30,40,40,50)
  
  res = cliff.delta(treatment,control,use.unbiased=F,use.normal=T)
  res.dm = cliff.delta(treatment,control,return.dm=T,use.unbiased=F,use.normal=T)
#Cliff's Delta
#
#delta estimate: -0.25 (small)
#95 percent confidence interval:
#  inf        sup 
#-0.7265846  0.3890062 
  expect_equal(as.character(res$magnitude),"small")
  expect_equal(res$estimate,-0.25)
  expect_equal(res$estimate,res.dm$estimate)
  
}) 


test_that("Dominance matrix on request", {
  treatment = c(10,10,20,20,20,30,30,30,40,50)
  control = c(10,20,30,40,40,50)
  
  res = cliff.delta(treatment,control,return.dm = TRUE)
  expect_true(any(grepl("dm",names(res))))
})


# 
#   d = c(control,treatment)
#   f = rep(c("Control","_Treat"),c(length(control),length(treatment)))
#   cliff.delta(d,f,use.unbiased=F,use.normal=T)
#   
# 

test_that("Negligible differences", {
  x1 = c(10, 20, 20, 20, 30, 30, 30, 40, 50, 100)
  x2 = c(10, 20, 30, 40, 40, 50)

  res<-cliff.delta(x1, x2)
  res.dm <- cliff.delta(x1, x2,return.dm = TRUE)
  expect_equal(as.character(res$magnitude),"negligible")
  expect_equal(res$estimate,-0.06666,tolerance=0.00001)
  expect_equal(res$estimate,res.dm$estimate)
})

test_that("Non overlapping", {
  x1 = c(10, 20, 20, 20, 30, 30, 30, 40, 50, 100)
  x2 = c(10, 20, 30, 40, 40, 50) + 110
  
  res<-cliff.delta(x1, x2)
  expect_equal(as.character(res$magnitude),"large")
  
})

test_that("Consistence between naive and partitioning", {
  set.seed(1097)
  
  for(i in 1:100){
    n = round(runif(2,10,20))
    
    x1 = round(runif(n[1],1,20))
    x2 = round(runif(n[2],15,30))
    
    res<-cliff.delta(x1, x2)
    res.dm<-cliff.delta(x1, x2,return.dm = TRUE)
    
    expect_equal(res$estimate,res.dm$estimate)  
  }
})


test_that("Resonable CI for extrem cases", {
  set.seed(1097)
  
    n = round(runif(2,10,20))
    
    x1 = round(runif(n[1],1,10))
    x2 = round(runif(n[2],15,30))
    
    expect_warning( res<<-cliff.delta(x1, x2) )

    expect_false(any(is.na(res$conf.int)))
})

test_that("double factor", {
  d <- data.frame(v = c("A","B","A","C","B","C","B","B","C","B"),
                  f = rep(c("G1","G2"),each=5))
  cliff.delta(v ~ f, data=d)
  cliff.delta(d$v , d$f)
})