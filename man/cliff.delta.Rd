\name{cliff.delta}
\alias{cliff.delta}
\alias{cliff.delta.default}
\alias{cliff.delta.formula}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cliff's Delta effect size for ordinal variables
}
\description{
Computes the Cliff's Delta effect size for ordinal variables
with the related confidence interval using efficient algorithms.
}
\usage{
cliff.delta(treatment, ... )

\method{cliff.delta}{formula}(formula, data=list() ,conf.level=.95, 
                                use.unbiased=TRUE, use.normal=FALSE, 
                                return.dm=FALSE, ...)

\method{cliff.delta}{default}(d, f, treatment, control, conf.level=.95, 
                         use.unbiased=TRUE, use.normal=FALSE, 
                         return.dm=FALSE, ...)


}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{treatment}{
numeric vector or ordered factor of data values for the treatment group (see Details)
}
  \item{control}{
numeric vector or ordered factor of data values for the control group (see Details)
}
 \item{d}{
a numeric vector giving either the data values (if \code{f} is a factor) or the treatment group values (if \code{f} is a numeric vector)
}
  \item{f}{
either a factor with two levels or a numeric vector of values
}
  \item{conf.level}{
confidence level of the confidence interval
}
  \item{use.unbiased}{
a logical indicating whether to compute the delta's variance using the
"unbiased" estimate formula or the "consistent" estimate
}
  \item{use.normal}{
logical indicating whether to use the normal or Student-t distribution
for the confidence interval estimation
}
  \item{return.dm}{
logical indicating whether to return the dominance matrix. \strong{Warning}: the explicit computation of the dominance uses a sub-optimal algorithm both in terms of memory and time
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
Uses the original formula reported in (Cliff 1996).

If the dominance matrix is  required i.e. \code{return.dm=TRUE}) the full matrix is computed thus using the naive algorithm.
Otherwise, if  \code{treatment} and \code{control} are \code{factor}s then the optimized linear complexity algorithm is used, otherwise the RLE algorithm (with complexity n log n) is used.

}
\value{
A list of class \code{effsize} containing the following components:
\item{estimate}{the Cliff's delta estimate}
\item{conf.int}{the confidence interval of the delta}
\item{var}{the estimated variance of the delta}
\item{conf.level}{the confidence level used to compute the confidence interval}
\item{dm}{the dominance matrix used for computation, only if \code{return.dm} is TRUE}
\item{magnitude}{a qualitative assessment of the magnitude of effect size}
\item{method}{the method used for computing the effect size, always \code{"Cliff's Delta"}}
\item{variance.estimation}{the method used to compute the delta variance estimation, either \code{"unbiased"} or \code{"consistent"}}
\item{CI.distribution}{the distribution used to compute the confidence interval, either \code{"Normal"} or \code{"Student-t"}}

The magnitude is assessed using the thresholds provided in (Romano 2006), i.e. |d|<0.147 \code{"negligible"}, |d|<0.33 \code{"small"}, |d|<0.474 \code{"medium"}, otherwise \code{"large"}
}
\references{
Norman Cliff (1996). Ordinal methods for behavioral data analysis. Routledge.

J. Romano, J. D. Kromrey, J. Coraggio, J. Skowronek, Appropriate statistics for ordinal level data: Should we really be using t-test and cohen's d for evaluating group differences on the NSSE and other surveys?, in: Annual meeting of the Florida Association of Institutional Research, 2006.

K.Y. Hogarty and J.D.Kromrey (1999). Using SAS to Calculate Tests of Cliff's Delta. Proceedings of the Twenty-Foursth Annual SAS User Group International Conference, Miami Beach, Florida, p 238.
Available at: \url{http://www2.sas.com/proceedings/sugi24/Posters/p238-24.pdf}
}
\author{
Marco Torchiano \url{http://softeng.polito.it/torchiano/}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cohen.d}}, \code{\link{print.effsize}}
}
\examples{
## Example data from Hogarty and Kromrey (1999)
treatment <- c(10,10,20,20,20,30,30,30,40,50)
control <- c(10,20,30,40,40,50)
res = cliff.delta(treatment,control,return.dm=TRUE)
print(res)
print(res$dm)
}
\keyword{ effect size }
\keyword{Cliff}% __ONLY ONE__ keyword per line
