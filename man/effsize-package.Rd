\name{effsize-package}
\alias{effsize-package}
\alias{effsize}
\docType{package}
\title{
Efficient Effect Size Computation
}
\description{
This packages contains functions to compute effect sizes both based on means difference (Cohen's d and Hedges g), dominance matrices (Cliff's Delta) and stochastic superiority (Vargha-Delaney A).

The computation (especially for Cliff's Delta) is carried on with higly efficient algorithms.
}
\details{

The main functions are:

\code{\link{cliff.delta}}.

\code{\link{cohen.d}}.

\code{\link{VD.A}}.


\bold{Change history}

\describe{
\item{0.3.1}{Fixed a bug in \code{cohen.d} when \code{PAIRED=TRUE}, now the \code{PAIRED} parameter has no effect, it is left just for compatibility. In a future code clean-up it may be removed}
\item{0.4}{Implemented a new algorithm with improved memory and time complexity. In particular new time complexity is T = O(n1*log(n2)) vs. the previous T = O(n1*n2), and new memory complexity M = O( n1 + n2 ) vs. the previous M = O( n1 * n2). In practice now the computation becomes feasible in a "\emph{reasonable}" time.}
\item{0.4.1}{Code clean-up and optimization using vectorized binary partioning.}
\item{0.5}{Added Vargha and Delaney A and fixed minor bugs with Cohen.d.}
\item{0.5.1}{Modified the Vargha and Delaney A computation to minimize accuracy errors.}
\item{0.5.2}{Fixed bug in \code{cliff.delta}.}
\item{0.5.3}{Fixed bug in  \code{cohen.d.formula}.}
\item{0.5.4}{Fixed minor issue detected by check.}
\item{0.5.5}{Changed the effsize field magnitude to a factor value.}
\item{0.6.0}{Implemented paired computation and CI computation with non-central t-distributions for cohen.d.}
\item{0.6.1}{Added ability to specify factor vector and data vector for `cliff.delta` function (thanks to Joses W. Ho).}
\item{0.6.2}{\code{na.rm} in \code{cohen.d} removes all incomplete pairs when paired.}
\item{0.6.3}{fixed bug in \code{cohen.d} when \code{na.rm=TRUE}, minor changes in the documentation (thanks to P.Thomas)}
\item{0.6.4}{Fixed a bug related to paired\code{cohen.d} with NAs.
Minor documentation changes}
\item{0.7.0}{Refactored tests using \code{testthat} package.
Fixed a bug in \code{cliff.delta} returning inconsistent
results when the dominance matrix is returned. Fixed issue concerning CI.
Fixed bug in \code{cohen.d} when using noncentral parameter for negative effect sizes.}
\item{0.7.1}{Fixed minor bugs in  \code{cliff.delta} and \code{cohen.d}}
\item{0.7.2}{Fixed bugs in \code{cohen.d}, order of factors is now observed and CI are computed correctly}
\item{0.7.3}{Fixed bugs in \code{cohen.d}, possible endless loop, cleaned code}
\item{0.7.4}{Fixed bugs in \code{cliff.delta} when values are factors}
\item{0.7.5}{Fixed bugs in \code{cohen.d} for paired data}
\item{0.7.6}{Fixed bugs in \code{cohen.d} for CI of paired data}
\item{0.7.7}{Fixed bugs in \code{cohen.d} for non-pooled SD, plus a few pull requests on documentation}
\item{0.7.8}{Fixed bug in \code{cohen.d} wrong correct type check }
\item{0.7.9}{Fixed tests to be compatible with upcoming R 4.0, that sets stringsAsFactors to FALSE by default }
\item{0.8.0}{Added non-central CI estimation for single sample \code{cohen.d}, fixed a bug related to order of data and added a \code{subject} parameter for paired \code{cohen.d}}
}}
\author{
Marco Torchiano \url{http://softeng.polito.it/torchiano/}
}
\keyword{ package }
\keyword{ effect size }
\keyword{ Cohen }
\keyword{ Hedges }
\keyword{ Cliff }
