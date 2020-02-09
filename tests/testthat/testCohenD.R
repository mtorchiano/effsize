### Cohen's d test
library(effsize)

try_with_time_limit <- function(expr, cpu = Inf, elapsed = Inf){
  y <- try({setTimeLimit(cpu, elapsed); expr}, silent = TRUE) 
  if(inherits(y, "try-error")) stop("Operation timed out") else y 
}


generate_data <- function(n,m,stdev){
  x <- rnorm(n,m,stdev)
  sd.adj = stdev/sd(x)
  x <- x * sd.adj
  m.adj = m - mean(x)
  x <- x + m.adj
  return(x)
}

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


test_that("Two independent samples with Hedges G, example from Borenstein et al, Introduction to Meta-Analysis", {
  set.seed(321)
  x1 = generate_data(50,103,5.5)
  x2 = generate_data(50,100,4.5)
  eff.g = cohen.d(x1,x2,hedges.correction=TRUE)
  expect_equal(eff.g$estimate, 0.5924, tolerance=.0001)
})


test_that("Two samples from same population", {
  set.seed(54)
  d <- rnorm(200)
  f <- rep(c(1,2),100)
  expect_warning( eff.d <- cohen.d(d ~ f) )
  
  expect_lt(eff.d$conf.int[1], 0)
  expect_gt(eff.d$conf.int[2], 0)
})


# test_that("Paired measures", {
#   delta = c(1.73, 1.06, 2.03, 1.40, 0.95, 1.13, 1.41, 1.73, 1.63, 1.56) - 1
#   set.seed(50)
#   a = delta
#   set.seed(50)
#   b = delta + runif(10)
# 
#   eff.d = cohen.d(a, b, paired=TRUE)
#   expect_equal(abs(eff.d$estimate),1.42, tolerance=.01)
# })


test_that("Paired measures w/NA", {
  set.seed(60)
  x = rnorm(10,mean=10)
  y = rnorm(10,mean=12)
  x[4] <- NA
  eff.d = cohen.d(x,y,paired=TRUE,na.rm=TRUE)
  
  expect_equal(abs(eff.d$estimate),2.8, tolerance=.01)
})


test_that("Paired measures with multiple NA", {
  eff.d = cohen.d(c(0.2,0.4,0.8,0.9,NA), 
                  c(0.12,0.14,0.18,0.119,NA), 
                  paired = T, na.rm = T)
  
  expect_equal(abs(eff.d$estimate),1.55, tolerance=.01)
})

test_that("Non centrality parameter",{
  # From Confidence Intervals on Effect Size, David C. Howell
  # Data from: Adams, Wright, and Lohr (1996)
  
  set.seed(22)
  Homophobic = generate_data(35,24,sqrt(148.87))
  Nonhomophobic = generate_data(29,16.5,sqrt(139.16))
  eff.d = cohen.d(Homophobic, Nonhomophobic, noncentral = TRUE)
  expect_equal(as.numeric(eff.d$estimate[1]),0.62,tolerance=0.01)
  expect_equal(as.numeric(eff.d$conf.int[1]),0.117,tolerance=0.001)
  expect_equal(as.numeric(eff.d$conf.int[2]),1.125,tolerance=0.001)
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


#issue #27 Cohen.d gives wrong value when data is not arranged by f
test_that("Order of factor values does not affect result",{
  set.seed(25)
  d.data <- data.frame(group = factor(sample(c(1,2), 20, replace = T, prob = c(.33,.67))),
                       value = rnorm(60,100,15))
  group1 <- d.data$value[d.data$group==1]
  group2 <- d.data$value[d.data$group==2]
  d.data.arranged <- d.data[order(d.data$group),]
  r1 <- cohen.d(d.data$value ~ d.data$group)
  r2 <- cohen.d(d.data.arranged$value ~ d.data.arranged$group)
  
  expect_equal(as.numeric(r1$estimate),as.numeric(r2$estimate))
})

#spin-off of issue #27 Cohen.d gives wrong value when data is not arranged by f
test_that("When inverting control and treatment effsize just change sign",{
  set.seed(7)
  group1 <- rnorm(20,100,15)
  group2 <- rnorm(40,100,15)
  r1 <- cohen.d(group1,group2)
  r2 <- cohen.d(group2,group1)
  
  expect_equal(r1$estimate,-r2$estimate)
})


# issue #28 confidence interval with non-central distribution for paired data
test_that("confidence interval with non-central distribution for paired data",{
  # data from https://www.uvm.edu/%7Edhowell/methods7/Supplements/Confidence%20Intervals%20on%20Effect%20Size.pdf
  moon.data = c(1.73, 1.06, 2.03, 1.40, 0.95, 1.13, 1.41, 1.73, 1.63, 1.56)
  g1 = rep(1,length(moon.data))
#  x2 = seq(1,1.9,by=0.1)
#  x1 = moon.data+x2
  res = cohen.d(moon.data,g1,
                   paired = TRUE,
                   noncentral = TRUE,
                   conf.level = 0.95
  )  
  expect_equal(as.numeric(res$conf.int[1]),0.4907785,tolerance = .000001)
  expect_equal(as.numeric(res$conf.int[2]),2.3358769,tolerance = .000001)
})


# issue #28 confidence interval with non-central distribution for two samples
test_that("confidence interval with non-central distribution for two samples",{
  # data from https://www.uvm.edu/%7Edhowell/methods7/Supplements/Confidence%20Intervals%20on%20Effect%20Size.pdf
  set.seed(537)
  hf = generate_data(35,24,sqrt(148.87))
  nhf = generate_data(29,16.5,sqrt(139.16))
  
  res= cohen.d(hf,nhf,noncentral=TRUE)

  expect_equal(as.numeric(res$conf.int[1]),0.117,tolerance = .001)
  expect_equal(as.numeric(res$conf.int[2]),1.126,tolerance = .001)
})


# issue #?? unpaired noncentral error (infinite loop)
test_that("Two samples paired formula paired and  noncentral", {
  data(sleep)
  res = try_with_time_limit(cohen.d(extra ~ group,
                                    data = sleep, 
                                    noncentral =TRUE),1)
  expect_equal(as.numeric(res$estimate),-0.8321,tolerance = .0001)
  
})

# issue #32 paired Hedges correction
# Using the example from 
#    Michael Borenstein, L. V. Hedges, J. P. T. Higgins and H. R. Rothstein
#    Introduction to Meta-Analysis.
test_that("Two samples paired pooled with Hedges correction", {
    set.seed(7828)
    x1 = generate_data(50,103,5.5)
    d = generate_data(50,3,5.5)
    x2 = x1 - d

    res.d = cohen.d(x1, x2, paired=TRUE)
    res.g = cohen.d(x1, x2, paired=TRUE, hedges.correction = TRUE)
  
    expect_equal(as.numeric(res.d$estimate),0.4225,tolerance = .0001)
    expect_equal(as.numeric(res.d$var),0.0131,tolerance = .0001)
    expect_equal(as.numeric(res.g$estimate),0.4160,tolerance = .0001)
    expect_equal(as.numeric(res.g$var),0.0127,tolerance = .0001)

})

# from issue #34
test_that("Two samples paired normal cohen d", {
 x1 = c(131,124,130,105,102,120,101,120)
 x2 = c(133,134,129,112,106,125,106,129)
  
  res.d = cohen.d(x2, x1, paired=TRUE)

  expect_equal(as.numeric(res.d$estimate),0.422,tolerance = .001)
})


test_that("Single sample with null mu",{
  moon.data = c(1.73, 1.06, 2.03, 1.40, 0.95, 1.13, 1.41, 1.73, 1.63, 1.56)
  
  res.d = cohen.d( ~ moon.data,data=data.frame(moon.data=moon.data))
  
  expect_equal(as.numeric(res.d$estimate),4.294,tolerance = .001)
})

test_that("Single sample with non-null mu",{
  moon.data <- c(1.73, 1.06, 2.03, 1.40, 0.95, 1.13, 1.41, 1.73, 1.63, 1.56)
  
  res.d = cohen.d( ~ moon.data, mu=1,data=data.frame(moon.data=moon.data))

  expect_equal(as.numeric(res.d$estimate),1.359,tolerance = .001)
})

test_that("Cohen pooled false",{
  Before = c(45L, 52L, 63L, 68L, 57L, 55L, 60L, 59L)
  After = c(49L, 50L, 70L, 71L, 53L, 61L, 62L, 67L)
  
  res.d = cohen.d(After, Before, pooled=FALSE)
  
  expect_equal(as.numeric(res.d$estimate),0.429,tolerance = .001)
})

test_that("Cohen numeric first arg with own class",{
  d = c (2,3,4,5,6,7)
  f = factor(rep(c("a","b"),each=3))
  
  class(d) <- "Cheese"
  res.d = cohen.d(d, f)
  
  expect_equal(as.numeric(res.d$estimate),-3,tolerance = .1)
})
