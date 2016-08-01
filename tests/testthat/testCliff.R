## Test for Cliff's delta
library(effsize)

context("Cohen d")


test_that("Small differences", {
  treatment = c(10,10,20,20,20,30,30,30,40,50)
  control = c(10,20,30,40,40,50)
  
  res = cliff.delta(treatment,control,use.unbiased=F,use.normal=T)
  
#Cliff's Delta
#
#delta estimate: -0.25 (small)
#95 percent confidence interval:
#  inf        sup 
#-0.7265846  0.3890062 
  expect_equal(as.character(res$magnitude),"small")
  expect_equal(res$estimate,-0.25)
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
  expect_equal(as.character(res$magnitude),"negligible")

})

