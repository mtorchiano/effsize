\name{VD.A}
\alias{VD.A}
\alias{VD.A.default}
\alias{VD.A.formula}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Vargha and Delaney A measure
}
\description{
Computes the Vargha and Delaney A effect size measure.
}
\usage{
VD.A(d, ...)

\method{VD.A}{formula}(formula,data=list(), ...)

\method{VD.A}{default}(d,f, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
a numeric vector giving either the data values (if \code{f} is a factor) or the treatment group values (if \code{f} is a numeric vector)
}
  \item{f}{
either a factor with two levels or a numeric vector of values
}
  \item{formula}{
    a formula of the form \code{y ~ f}, where \code{y} is a numeric variable giving the data values and \code{f} a factor with two levels giving the corresponding group
  }
  \item{data}{
  an optional matrix or data frame containing the variables in the formula \code{formula}. By default the variables are taken from \code{environment(formula)}.
  }
  \item{...}{
  further arguments to be passed to or from methods.
  }
}
\details{
The function computes the Vargha and Delaney A effect size measure (Vargha and Delaney, 2000).
}
\value{
A list of class \code{effsize} containing the following components:
\item{estimate}{the A statistics estimate}
\item{magnitude}{a qualitative assessment of the magnitude of effect size}
\item{method}{the method used, i.e. \code{"Vargha and Delaney A"}}
}
\references{
A. Vargha and H. D. Delaney. "A critique and improvement of the CL common language effect size statistics of McGraw and Wong." Journal of Educational and Behavioral Statistics, 25(2):101-132, 2000
}
\author{
Marco Torchiano \url{http://softeng.polito.it/torchiano/}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cliff.delta}}, \code{\link{cohen.d}}, \code{\link{print.effsize}}
}
\examples{
treatment = rnorm(100,mean=10)
control = rnorm(100,mean=12)
d = (c(treatment,control))
f = rep(c("Treatment","Control"),each=100)
## compute Vargha and Delaney A
## treatment and control
VD.A(treatment,control)
## data and factor
VD.A(d,f)
## formula interface
VD.A(d ~ f)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ effect size }
\keyword{ Vargha and Delaney }% __ONLY ONE__ keyword per line
