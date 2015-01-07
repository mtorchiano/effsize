## Data from:
##     M.Hess, J.Kromrey. 
##     Robust Confidence Intervals for Effect Sizes: 
##     A Comparative Study of Cohen’s d and Cliff’s Delta Under Non-normality and Heterogeneous Variances
##     American Educational Research Association, San Diego, April 12 – 16, 2004
##
## Available at: http://www.coedu.usf.edu/main/departments/me/documents/cohen.pdf
##
## Kristine Y. Hogarty and Jeffrey D. Kromrey
## Using SAS to Calculate Tests of Cliff’s Delta
## Proceedings of the Twenty-Fourth Annual SAS® Users Group ￼International Conferenc, Miami Beach, Florida, 1999
## Availabel at: http://www2.sas.com/proceedings/sugi24/Posters/p238-24.pdf

.bsearch.partition <- function(x, a, b=1, e=length(a) ){
  n = length(x);
  low = rep(NA,n)
  L = rep(b,n)
  H = rep(e,n)
  i0 = seq(n)
  
  repeat{
    #M = (L+H) %/% 2
    M = as.integer((L+H) / 2)
    left = x <= a[M]
    H[left] = M[left]
    L[!left] = M[!left]+1
    if( all(H<=L) ) {
      break;
    }
  }
  
  H=L
  repeat{
    below = a[H] == x
    below[is.na(below)] = FALSE
    if(!any(below)) break;
    H[below] = H[below] + 1
  }
  
  repeat{
    L.clean = L
    L.clean[L.clean<1] = NA
    above = a[L.clean] >= x
    above[is.na(above)] = FALSE
    if(!any(above)) break
    L[above] = L[above] - 1
  }
  cbind(below=L,above=H)
}


cliff.delta <- function(treatment, ... ) UseMethod("cliff.delta")

cliff.delta.default <- function( treatment, control, conf.level=.95, 
                         use.unbiased=TRUE, use.normal=FALSE, return.dm=FALSE, ...){

  if(conf.level>.9999 | conf.level<.5) stop("conf.level must be within 50% and 99.99%")
  
  n1 = length(treatment)
  n2 = length(control)

  if(!return.dm & is.factor(treatment) & is.factor(control)){ ## use factor algorithm
    algorithm="Factor"
    ft = table(treatment) / n1
    fc = table(control) / n2
    ## Note: in principle they should share the same levels
    lt = as.numeric(levels(treatment))
    lc = as.numeric(levels(control))
    
    dm_L = sign(outer(lt,lc,FUN="-"))
    
    d_i.C = dm_L %*% fc
    d_.jC = t(dm_L) %*% ft
    d =  as.numeric( t(d_i.C) %*% ft )

    d_i. = rep(d_i.C,ft*n1) 
    d_.j = rep(d_.jC,fc*n2) 
    
    SSR = t((dm_L-d)^2 %*% (fc*n2)) %*% ft *n1 
  }else{
  treatment = sort(treatment)
  control = sort(control)
  
  if(return.dm){ ## explicitly compute dominance matrix
    algorithm="Naive"
    dominance = sign(outer(treatment, control, FUN="-")) 
    row.names(dominance) = treatment
    colnames(dominance) = control
    
    d = mean(dominance)
    d_i. = apply(dominance,1,mean)
    d_.j = apply(dominance,2,mean)
    SSR = sum( (dominance-d)^2 )
  }else{ ## uses row partitioning algorithm
    algorithm="Row partitioning"
    partitions = .bsearch.partition(treatment,control)
    partitions[,2] = n2 - partitions[,2] + 1L
    partitions[partitions[,1]>n2,1] = n2
    d_i. = partitions %*% c(1L,-1L) / n2
    
    d = mean(d_i.)
    pb = sum(partitions[,1])
    pa = sum(partitions[,2])
    partitions = .bsearch.partition(control,treatment)
    partitions[,2] = n1 - partitions[,2] + 1L
    partitions[partitions[,1]>n1,2] = n1
    d_.j = partitions %*% c(-1L,1L) / n1
    
    #     d_.j = rep(NA,n2)
    #     for(i in 1:n2){
    #       d_.j = sum(partitions[,1]>=i) -sum(partitions[,1]<n2-i)
    #     }
    
    SSR = pb * (1-d)^2 + (as.double(n1)*n2-pa-pb)*d^2 + pa*(1+d)^2 
  }
  }
  ## Compute variance
  if(use.unbiased){
    # method 1: unbiased estimate:
    S_d = ( n2^2 * sum( (d_i. - d)^2) + n1^2*sum( (d_.j-d)^2) - SSR ) / (as.numeric(n1)*n2*(n1-1)*(n2-1))
  }else{
    # method 2: consistent estimate
    S_i. = sum( (d_i.-d)^2 ) / (n1-1);
    S_.j = sum( (d_.j-d)^2 ) / (n2-1);
    S_ij = SSR / ( ( n1-1)*(n2-1 ) ) ### Cliff 1996
    # The three variance terms should be:
    # 0.2669753,  0.447, and  1.0055556, respectively.
    #S_ij = SSR / ( n1*n2-1 )  ## Long et al 2003
    S_d = ( (n2-1)*S_i. + (n1-1)*S_.j + S_ij) / ( as.double(n1) * n2)
  }
  
if(use.normal){
  ## assume a normal distribution
Z = -qnorm((1-conf.level)/2)
}else{
  ## assume a Student t distribution See (Feng & Cliff, 2004) 
Z = -qt((1-conf.level)/2,n1+n2-2)
}
conf.int = c(
  ( d - d^3 - Z * sqrt(S_d) * sqrt((1-d^2)^2+Z^2*S_d )) / ( 1 - d^2+Z^2*S_d),
  ( d - d^3 + Z * sqrt(S_d) * sqrt((1-d^2)^2+Z^2*S_d )) / ( 1 - d^2+Z^2*S_d)
)
names(conf.int) = c("inf","sup")
  
  levels = c(0.147,0.33,0.474) ## effect sizes from (Hess and Kromrey, 2004)
  magnitude = c("negligible","small","medium","large")
res= 
 list(
  estimate = d,
  conf.int = conf.int,
  var = S_d,
  conf.level = conf.level,
  magnitude = magnitude[findInterval(abs(d),levels)+1],
  method = "Cliff's Delta",
  algorithm = algorithm,
  variance.estimation = if(use.unbiased){ "Unbiased"}else{"Consistent"},
  CI.distribution = if(use.normal){ "Normal"}else{"Student-t"}
  )
  if(return.dm){
    res$dm = dominance;
  }
  res$name = "delta"
  
  class(res) <- "effsize"
  return(res)
}


cliff.delta.formula <-function(formula, data=list(),conf.level=.95, 
                                use.unbiased=TRUE, use.normal=FALSE, 
                               return.dm=FALSE, ...){
#cliff.delta.formula <-function(f, data=list(),...){
  mf <- model.frame(formula=formula, data=data)
  if(dim(mf)[2]!=2){
     stop("Fomula must be a variable vs a factor")
  }
  x <- mf[[2]]
  y <- mf[[1]]
  vals = split(y,x)
  x = vals[[1]]
  y = vals[[2]]
  res = cliff.delta.default(x,y,conf.level, use.unbiased, use.normal, return.dm)
  return(res)
}

# test_cliff.delta <- function(){
#   treatment = c(10,10,20,20,20,30,30,30,40,50)
#   control = c(10,20,30,40,40,50)
#   
#   res = cliff.delta(treatment,control,use.unbiased=F,use.normal=T)
#   
#   print(res)
##Cliff's Delta
##
##delta estimate: -0.25 (small)
##95 percent confidence interval:
##  inf        sup 
##-0.7265846  0.3890062 
# }


