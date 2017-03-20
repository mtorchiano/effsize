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


try_with_time_limit <- function(expr, cpu = Inf, elapsed = Inf)
{
  y <- try({setTimeLimit(cpu, elapsed); expr}, silent = TRUE) 
  if(inherits(y, "try-error")) NULL else y 
}

# issue #23 noncentral error (infinite loop)
test_that("Two samples with large negative difference and noncentral", {
  a = c(9.81605624621576, 8.93891560898168, 9.05436620537713, 6.01771305071382, 
        8.98575772172126, 7.2595761413452, 13.7537304479165, 9.89861975363729, 
        7.37246838418731, 8.33888742066376)
  b = c(8.62798051757027, 5.25390514295654, 12.0869496333763, 7.20553319908951, 
        6.99017188584104, 8.22363490828756, 8.05485287987422, 6.50534569228974, 
        8.68871033652186, 4.62897740315003)
  res = try_with_time_limit(cohen.d(a,b,hedges.correction = FALSE, noncentral =TRUE),1) 
  expect_equal(as.numeric(res$estimate),0.63236,tolerance = .0001)
})
