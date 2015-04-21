
##########  TEST ###################

library(effsize)

assert <- function(label,condition){
  cat(label,": ")
  if(!condition){
    cat("Failed!\n")
  }else{
    cat("OK.\n")
  }
}

test_cliff.delta <- function(){
  treatment = c(10,10,20,20,20,30,30,30,40,50)
  control = c(10,20,30,40,40,50)
  
  res = cliff.delta(treatment,control,use.unbiased=F,use.normal=T)
  
  print(res)
  
#Cliff's Delta
#
#delta estimate: -0.25 (small)
#95 percent confidence interval:
#  inf        sup 
#-0.7265846  0.3890062 
  assert("small difference", res$magnitude=="small"  )


  x1 = c(10, 20, 20, 20, 30, 30, 30, 40, 50, 100)
  x2 = c(10, 20, 30, 40, 40, 50)


  print(res<-cliff.delta(x1, x2))
  assert("Negligible difference", res$magnitude=="negligible"  )

  print(res <- cliff.delta(x1, x2, return.dm = TRUE))
  assert("Dominance matrix ", any(grepl("dm",names(res)))  )
}

test_cliff.delta();
