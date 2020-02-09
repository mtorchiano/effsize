\name{cohen.d}
\alias{cohen.d}
\alias{cohen.d.default}
\alias{cohen.d.formula}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cohen's d and Hedges g effect size
}
\description{
Computes the Cohen's d and Hedges'g effect size statistics.
}
\usage{

cohen.d(d, ...)

\method{cohen.d}{formula}(formula,data=list(),...)

\method{cohen.d}{default}(d,f,pooled=TRUE,paired=FALSE,
                   na.rm=FALSE, mu=0, hedges.correction=FALSE,
                   conf.level=0.95,noncentral=FALSE, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
a numeric vector giving either the data values (if \code{f} is a factor) or the treatment group values (if \code{f} is a numeric vector)
}
  \item{f}{
either a factor with two levels or a numeric vector of values
}
  \item{pooled}{
a logical indicating whether compute pooled standard deviation or the whole sample standard deviation. If \code{pooled=TRUE} (default) pooled sd is used, if \code{pooled=FALSE} the standard deviation of the the control group (the second argument or the one corresponding the the second level of the factor) is used instead.
}
  \item{paired}{
a logical indicating whether to consider the values as paired 
}
  \item{na.rm}{
logical indicating whether \code{NA}s should be removed before computation;
if \code{paired==TRUE} then all incomplete pairs are removed.
}
  \item{mu}{
numeric indicating there reference mean for single sample effect size.
}
  \item{hedges.correction}{
logical indicating whether apply the Hedges correction
}
  \item{conf.level}{
confidence level of the confidence interval
}
  \item{formula}{
    a formula of the form \code{y ~ f}, where \code{y} is a numeric variable giving the data values and \code{f} a factor with two levels giving the corresponding groups
  }
  \item{data}{
  an optional matrix or data frame containing the variables in the formula \code{formula}. By default the variables are taken from \code{environment(formula)}.
  }
  \item{noncentral}{
  logical indicating whether to use non-central t distributions for
  computing the confidence interval.}
  \item{...}{
  further arguments to be passed to or from methods.
  }
}
\details{
When \code{f} in the default version is a factor or a character, it must have two values and it identifies the two groups to be compared. Otherwise (e.g. \code{f} is numeric), it is considered as a sample to be compare to \code{d}.

In the formula version, if \code{f} is expected to be a factor, if that is not the case it is coherced to a factor and a warning is issued.

The function computes the value of Cohen's d statistics (Cohen 1988).
If required (\code{hedges.correction==TRUE}) the Hedges g statistics is computed instead (Hedges and Holkin, 1985).

When \code{paired} is set, the effect size is computed using the approach suggested in (Gibbons et al. 1993).

It is possible to perform a single sample effect size estimation either using a formula \code{~x} or passing \code{f=NA}. 


The computation of the CI requires the use of non-central Student-t distributions that are used when \code{noncentral==TRUE}; otherwise a central distribution is used.

Also a quantification of the effect size magnitude is performed using the thresholds define in Cohen (1992).
The magnitude is assessed using the thresholds provided in (Cohen 1992), i.e. |d|<0.2 \code{"negligible"}, |d|<0.5 \code{"small"}, |d|<0.8 \code{"medium"}, otherwise \code{"large"}

The variance of the \code{d} is computed using the conversion formula reported at page 238 of Cooper et al. (2009):

\deqn{ S^2_d = \left( \frac{n_1+n_2}{n_1 n_2} + \frac{d^2}{2 df}\right) \left( \frac{n_1+n_2}{df} \right)}{((n1+n2)/(n1*n2) + .5*d^2/df) * ((n1+n2)/df)}

}
\value{
A list of class \code{effsize} containing the following components:
\item{estimate}{the statistic estimate}
\item{conf.int}{the confidence interval of the statistic}
\item{sd}{the within-groups standard deviation}
\item{conf.level}{the confidence level used to compute the confidence interval}
\item{magnitude}{a qualitative assessment of the magnitude of effect size}
\item{method}{the method used for computing the effect size, either \code{"Cohen's d"} or \code{"Hedges' g"}}
}
\references{
Cohen, J. (1988). 
Statistical power analysis for the behavioral sciences (2nd ed.). 
New York:Academic Press.

Hedges, L. V. & Olkin, I. (1985). Statistical methods for meta-analysis. Orlando, FL: Academic Press.

Cohen, J. (1992). A power primer. Psychological Bulletin, 112, 155-159.

The Handbook of Research Synthesis and Meta-Analysis (Cooper, Hedges, & Valentine, 2009)

David C. Howell (2010). Confidence Intervals on Effect Size. Available at: \url{https://www.uvm.edu/\%7Edhowell/methods7/Supplements/Confidence\%20Intervals\%20on\%20Effect\%20Size.pdf}

Cumming, G.; Finch, S. (2001). A primer on the understanding, use, and calculation of confidence intervals that are based on central and noncentral distributions. Educational and Psychological Measurement, 61, 633-649.

Gibbons, R. D., Hedeker, D. R., & Davis, J. M. (1993). 
Estimation of effect size from a series of experiments involving paired comparisons. 
Journal of Educational Statistics, 18, 271-279.

}
\author{
Marco Torchiano \url{http://softeng.polito.it/torchiano/}
}
%\note{
%%
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cliff.delta}}, \code{\link{VD.A}}, \code{\link{print.effsize}}
}
\examples{
treatment = rnorm(100,mean=10)
control = rnorm(100,mean=12)
d = (c(treatment,control))
f = rep(c("Treatment","Control"),each=100)
## compute Cohen's d
## treatment and control
cohen.d(treatment,control)
## data and factor
cohen.d(d,f)
## formula interface
cohen.d(d ~ f)
## compute Hedges' g
cohen.d(d,f,hedges.correction=TRUE)
}
\keyword{ effect size }
\keyword{Hedges}
\keyword{Cohen}% __ONLY ONE__ keyword per line
