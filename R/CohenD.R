## Computes Cohen's d
## Cohen, J. (1988). 
## Statistical power analysis for the behavioral sciences (2nd ed.). 
## New York:Academic Press.

cohen.d <- function(d,...) UseMethod("cohen.d")

cohen.d.default = function(d,f,pooled=TRUE,paired=FALSE,na.rm=FALSE,
                   hedges.correction=FALSE,conf.level=0.95, ...){
  if( ! any(c("numeric","integer") %in% class(d))){
    stop("First parameter must be a numeric type")
  }
  if( any(c("character","factor") %in% class(f)) ){
    ## it is data and factor
    if(length(f)!=length(d)){
      stop("Data and factor must have the same length")
    }
    if( "character" %in% class(f)){
        f = factor(f)
    }  
    if(length(levels(f))!=2){
      stop("Factor should have only two levels");
      return;
    }
  }else{
    ## it is treatment and control
    treatment = d
    control = f
    d = c(treatment,control)
    f = factor(rep(c("Treatment","Control"),c(length(treatment),length(control))),
               levels=c("Treatment","Control"),ordered=T)
  }
  if(na.rm){
    nas = is.na(d) | is.na(f);
    d = d[!nas];
    f = f[!nas];
  }
  ns = table(f)
  n1 = ns[1]
  n2 = ns[2]
  m = c();
  sd = c();
  for( l in levels(f)){
    m = c(m,mean(d[f==l]));
    sd = c(sd,sd(d[f==l]));
  }
  
  delta.m = m[1] - m[2];
  
  if(pooled){
    pool_sd = sqrt(((n1-1)*sd[1]^2+(n2-1)*sd[2]^2)/(n1+n2-2))
    d = (delta.m) / pool_sd;
  }else{
    d = (delta.m) / sd(d);
  }
  df = n1+n2-2
  
  res = list()
  if(hedges.correction){
    # Hedges, L. V. & Olkin, I. (1985). Statistical methods for meta-analysis. Orlando, FL: Academic Press.
    d = d * (1 - 3 / ( 4 * (n1+n2) - 9))
    res$method = "Hedges's g"
    res$name = "g"
  }else{
    res$method = "Cohen's d"
    res$name = "d"
  }
  
  # The Handbook of Research Synthesis and Meta-Analysis (Cooper, Hedges, & Valentine, 2009)
  ## p 238
  S_d = sqrt(((n1+n2)/(n1*n2) + .5*d^2/df) * ((n1+n2)/df))
  
  Z = -qt((1-conf.level)/2,df)
  
  conf.int=c(
    d - Z*S_d,
    d + Z*S_d
  );
  names(conf.int)=c("inf","sup")
  
  levels = c(0.2,0.5,0.8)
  magnitude = c("negligible","small","medium","large")
  ## Cohen, J. (1992). A power primer. Psychological Bulletin, 112, 155-159. Crow, E. L. (1991).
  
  res$estimate = d
  res$conf.int = conf.int
  res$var = S_d
  res$conf.level = conf.level
  res$magnitude = magnitude[findInterval(abs(d),levels)+1]
#      variance.estimation = if(use.unbiased){ "Unbiased"}else{"Consistent"},
#      CI.distribution = if(use.normal){ "Normal"}else{"Student-t"}

  class(res) <- "effsize"
  return(res)
}

cohen.d.formula= function(formula, data=list(), ...){
  mf <- model.frame(formula=formula, data=data)
  if(dim(mf)[2]!=2){
    stop("Formula must be a variable vs a factor")
  }
  d <- mf[[1]]
  f <- mf[[2]]
  if( ! any(c("character","factor") %in% class(f)) ){
    warning("Cohercing rhs of formula to factor")
    f = factor(f)
  }  
  res = cohen.d.default(d,f,...)
  return(res)
}

# set.seed(52)
# x = rnorm(100,mean=10)
# y = rnorm(100,mean=12)
# d = (c(x,y))
# f = rep(c("A","B"),each=100)
# eff.d = cohen.d(d,f)
# print(eff.d)
# eff.g = cohen.d(d,f,hedges.correction=TRUE)
# print(eff.g)
# set.seed(12345)
# d <- rnorm(200)
# f <- rep(c(1,2),400)
# cohen.d(d ~ factor(f))




