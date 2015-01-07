\name{print.effsize}
\alias{print.effsize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prints effect size
}
\description{
Prints the results of an effect size computation
}
\usage{
\method{print}{effsize}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
the effect size result
}
  \item{...}{
further parameters are currently ignored
}
}
\details{
Shows the estimate value and, when available, the confidence interval.
}
\references{
See the main function \code{\link{cliff.delta}}.
}
\author{
Marco Torchiano \url{http://softeng.polito.it/torchiano/}
}
\note{
This is still work in progress..
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cliff.delta}} \code{\link{cohen.d}}
}