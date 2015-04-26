##
### Print effect size information
##

print.effsize <- function(x, ...){
  cat("\n")
  cat(x$method)
  cat("\n\n")
  cat(x$name)
  cat(" estimate: ")
  cat(x$estimate)
  cat(" (")
  cat(as.character(x$magnitude))
  cat(")\n")
  if("conf.level" %in% names(x)){
    conf = x$conf.level*100
    cat(conf)
    cat(" percent confidence interval:\n")
    print(x$conf.int)
  }
}

#S3method(print, effsize)