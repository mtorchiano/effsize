### Cohen's d test
library(effsize)

context("Cohen d")

test_that("Two samples with large difference", {
  set.seed(52)
  x = rnorm(100,mean=10)
  y = rnorm(100,mean=12)
  d = (c(x,y))
  f = rep(c("A","B"),each=100)
  eff.d = cohen.d(d,f)
  expect_lt(eff.d$conf.int[1], -2)
  expect_gt(eff.d$conf.int[2], -2)
})


test_that("Two samples with large difference with Hedges G", {
  set.seed(52)
  x = rnorm(100,mean=10)
  y = rnorm(100,mean=12)
  d = (c(x,y))
  f = rep(c("A","B"),each=100)
  eff.g = cohen.d(d,f,hedges.correction=TRUE)
  expect_lt(eff.g$conf.int[1], -2)
  expect_gt(eff.g$conf.int[2], -2)
})


test_that("Two samples from same population", {
  set.seed(54)
  d <- rnorm(200)
  f <- rep(c(1,2),100)
  expect_warning( eff.d <- cohen.d(d ~ f) )
  
  expect_lt(eff.d$conf.int[1], 0)
  expect_gt(eff.d$conf.int[2], 0)
})


test_that("Paired measures", {
  delta = c(1.73, 1.06, 2.03, 1.40, 0.95, 1.13, 1.41, 1.73, 1.63, 1.56) - 1
  set.seed(50)
  a = delta
  set.seed(50)
  b = delta + runif(10)

  eff.d = cohen.d(a,b,paired=TRUE)
  expect_equal(abs(eff.d$estimate),1.42, tolerance=.01)
})


test_that("Paired measures w/NA", {
  set.seed(60)
  x = rnorm(10,mean=10)
  y = rnorm(10,mean=12)
  x[4] <- NA
  eff.d = cohen.d(x,y,paired=TRUE,na.rm=TRUE)
  
  expect_equal(abs(eff.d$estimate),1.73, tolerance=.01)
})


test_that("Paired measures with multiple NA", {
  eff.d = cohen.d(c(0.2,0.4,0.8,0.9,NA), 
                  c(0.12,0.14,0.18,0.119,NA), 
                  paired = T, na.rm = T)
  
  expect_equal(abs(eff.d$estimate),1.35, tolerance=.01)
})

test_that("Non centrality parameter",{
  set.seed(22)
  a = rnorm(35,24,sqrt(148.87))
  set.seed(31)
  b = rnorm(29,16.5,sqrt(139.16))
  eff.d = cohen.d(a,b,noncentral = TRUE)
  expect_equal(as.numeric(eff.d$conf.int[1]),0.137,tolerance=0.01)
  expect_equal(as.numeric(eff.d$conf.int[2]),1.147,tolerance=0.01)
})


test_that("Two samples with large negative difference and noncentral", {
  set.seed(52)
  x = rnorm(100,mean=10)
  y = rnorm(100,mean=12)
  d = (c(x,y))
  f = rep(c("A","B"),each=100)
  eff.d = cohen.d(d,f,noncentral = TRUE)
  expect_lt(eff.d$conf.int[1], -2)
  expect_gt(eff.d$conf.int[2], -2)
})
