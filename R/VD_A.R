# Computes Vargha and Delaney A index
# A. Vargha and H. D. Delaney. 
# A critique and improvement of the CL common language 
# effect size statistics of McGraw and Wong. 
# Journal of Educational and Behavioral Statistics, 25(2):101-132, 2000
#
# The formula to compute A has been transformed to minimize accuracy errors
# See: http://mtorchiano.wordpress.com/2014/05/19/effect-size-of-r-precision/
#

VD.A <- function(d,...) UseMethod("VD.A")

VD.A.default <- function(d,f,...){
  if( "character"%in%class(f)|"factor"%in%class(f) ){
    ## it is data and factor
    if(length(f)!=length(d)){
      stop("Data d and factor f must have the same length")
    }
    if( "character" %in% class(f)){
      f = factor(f)
    }  
    if(length(levels(f))!=2){
      stop("Factor f should have only two levels");
      return;
    }
    tc = split(d,f)
    treatment = tc[[1]]
    control = tc[[2]]
  }else{
    ## it is treatment and control
    treatment = d
    control = f
  }

  m = length(treatment)
  n = length(control)
  r = rank(c(treatment,control))
  r1 = sum(r[seq_len(m)])
  
  # Compute the measure
  # A = (r1/m - (m+1)/2)/n # formula (14) in Vargha and Delaney, 2000
  A = (2*r1 - m*(m+1)) / (2*n*m) # equivalent formula to avoid accuracy errors
  # 

  
  levels = c(0.147,0.33,0.474) ## effect sizes from Hess and Kromrey, 2004
  magnitude = c("negligible","small","medium","large")
  scaled.A=(A-0.5)*2
    
  res = list(
    method = "Vargha and Delaney A",
    name = "A",
    magnitude = factor(magnitude[findInterval(abs(scaled.A),levels)+1],levels = magnitude,ordered=T),
    estimate = A
  )
  class(res) <- "effsize"
  return(res)
}

VD.A.formula <-function(formula, data=list(), ...){
  #cliff.delta.formula <-function(f, data=list(),...){
  mf <- model.frame(formula=formula, data=data)
  if(dim(mf)[2]!=2){
    stop("Fomula must be in the form: variable ~ factor")
  }
  fac <- mf[[2]]
  resp <- mf[[1]]
  res = VD.A.default(resp,fac)
  return(res)
}
