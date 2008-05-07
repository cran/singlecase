\name{singlecase-package}
\alias{singlecase-package}
\alias{singlecase}
\docType{package}
\title{Single-case inferential methods}
\description{The package implement frequentist and Bayesian inferential methods for use in single-case studies in which a patient is compared to a small control sample. 
}
\details{
\tabular{ll}{
Package: \tab singlecase\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2008-05-06\cr
License: \tab GPL 2 or newer\cr
LazyLoad: \tab yes\cr
}
At this stage of  development, the package contains implementation of functions for: 
(1) comparison of an individual's score on a single test with the score of a small normative or control sample: \code{crawford.t.test} and \code{SingleBayes} ; 
(2) test whether the difference between scores on two tests observed for a patient is significantly greater than the differences observed for a small control or normative sample: \code{crawford.diff.test} and \code{DiffBayes}.

\code{crawford.t.test}, and \code{crawford.diff.test} rely upon the $t$-distribution and treat summary statistics from the control group as estimates rather than parameters.  \code{SingleBayes} and \code{DiffBayes} are Bayesian inferential alternatives, based on a Monte Carlo process.  All methods may take either the raw data or summary statistics from the control population as inputs.
}
\author{Matthieu Dubois.

Maintainer: Matthieu Dubois \email{matthdub@gmail.com}, \url{http://www.code.psp.ucl.ac.be/MatthieuDubois/r_code.html}
}
\references{Crawford, J. and Garthwaite, P. (2002) Investigation of the single case in neuropsychology: Confidence limits on the abnormality of test scores and test score differences. \emph{Neuropsychologia}, \bold{40(8)}, 1196--1208.

Crawford, J. and Garthwaite, P. (2005) Testing for suspected impairments and dissociations in single-case studies in neuropsychology: Evaluation of alternatives using Monte Carlo simulations and revised tests for dissociations. \emph{Neuropsychology}, \bold{19(3)}, 318--331.

Crawford, J. and Garthwaite, P. (2007) Comparison of a single case to a control or normative sample in neuropsychology: Development of a bayesian approach. \emph{Cognitive Neuropsychology}, \bold{24(4)}, 343--372.

Crawford, J. and Howell, D. (1998) Comparing an individual's test score against norms drived from small samples. \emph{Clinical Neuropsychologist}, \bold{12(4)}, 482--486.

Crawford, J., Howell, D., and Garthwaite, P. (1998) Payne and Jones revisited: Estimating the abnormality of test score differences using a modified paired samples $t$-test. \emph{Journal of clinical and experimental neuropsychology}, \bold{20(6)}, 898--905.

}
\keyword{ package }
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
%\examples{
%~~ simple examples of the most important functions ~~
%}
