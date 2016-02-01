## Computes Cohen's d
## Cohen, J. (1988). 
## Statistical power analysis for the behavioral sciences (2nd ed.). 
## New York:Academic Press.

cohen.d <- function(d,...) UseMethod("cohen.d")

cohen.d.default = function(d,f,pooled=TRUE,paired=FALSE,na.rm=FALSE,
                   hedges.correction=FALSE,conf.level=0.95,noncentral=FALSE, ...){
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
  
  if(paired & (n1!=n2)){
    stop("Paired computation requires equal number of measures.");
    return;
  }
  
  m = c();
  sd = c();
  for( l in levels(f)){
    m = c(m,mean(d[f==l]));
    sd = c(sd,sd(d[f==l]));
  }
  
  delta.m = m[1] - m[2];
  if(paired){
    dd = (delta.m) / sd(diff(d,lag=n1));
  }else
  if(pooled){
    pool_sd = sqrt(((n1-1)*sd[1]^2+(n2-1)*sd[2]^2)/(n1+n2-2))
    dd = (delta.m) / pool_sd;
  }else{
    dd = (delta.m) / sd(d);
  }
  df = n1+n2-2
  
  res = list()
  if(hedges.correction){
    # Hedges, L. V. & Olkin, I. (1985). Statistical methods for meta-analysis. Orlando, FL: Academic Press.
    dd = dd * (1 - 3 / ( 4 * (n1+n2) - 9))
    res$method = "Hedges's g"
    res$name = "g"
  }else{
    res$method = "Cohen's d"
    res$name = "d"
  }
  
  if(noncentral){
    # Based on the document:
    # David C. Howell (2010)
    # Confidence Intervals on Effect Size
    # https://www.uvm.edu/%7Edhowell/methods7/Supplements/Confidence%20Intervals%20on%20Effect%20Size.pdf
    #
    # Additional reference:
    # Cumming, G.; Finch, S. (2001) 
    # A primer on the understanding, use, and calculation of confidence intervals 
    # that are based on central and noncentral distributions. 
    # Educational and Psychological Measurement, 61, 633-649.
    #
    if(paired){
      t = mean(diff(d,lag=n1))/(sd(diff(d,lag=n1))/sqrt(n1))
      df=n1-1
    }else{
      if(pooled) s = pool_sd
      else s = sd(d)
      
      t = delta.m / sqrt(s^2*(1/n1+1/n2))
    }
    end1 = t
    while( pt(q=t,df=df,ncp=end1) > (1-conf.level)/2 ){
      end1 = end1 * 2
    }
    ncp1 = uniroot(function(x) (1-conf.level)/2-pt(q=t,df=df,ncp=x),c(-5,end1))$root
    
    end2 = t
    while( pt(q=t,df=df,ncp=end2) > (1+conf.level)/2 ){
      end2 = end2 * 2
    }
    cat("t: ",t,"  df:",df,"\n")
    #       cat("-5 -> ",pt(q=t,df=df,ncp=-5),"\n")
    #       cat(end2," -> ",pt(q=t,df=df,ncp=end2),"\n")
    ncp2 = uniroot(function(x) (1+conf.level)/2-pt(q=t,df=df,ncp=x),c(-5,end2))$root
    cat("ncp1:",ncp1,"\n")
    cat("ncp2:",ncp2,"\n")
    
    if(paired){
      conf.int=sort(c(
        ncp1/sqrt(df),
        ncp2/sqrt(df)
      ));
    }else{
      conf.int=sort(c(
        ncp1*sqrt(1/n1+1/n2),
        ncp2*sqrt(1/n1+1/n2)
      ));
    }
  }else{
  # The Handbook of Research Synthesis and Meta-Analysis 
  # (Cooper, Hedges, & Valentine, 2009)
  # p 238
  S_d = sqrt(((n1+n2)/(n1*n2) + .5*dd^2/df) * ((n1+n2)/df))
  
  Z = -qt((1-conf.level)/2,df)
  
  conf.int=c(
    dd - Z*S_d,
    dd + Z*S_d
  );
  }
  names(conf.int)=c("inf","sup")
  
  levels = c(0.2,0.5,0.8)
  magnitude = c("negligible","small","medium","large")
  ## Cohen, J. (1992). A power primer. Psychological Bulletin, 112, 155-159. Crow, E. L. (1991).
  
  res$estimate = dd
  res$conf.int = conf.int
#  res$var = S_d
  res$conf.level = conf.level
  res$magnitude = factor(magnitude[findInterval(abs(dd),levels)+1],levels = magnitude,ordered=T)
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


