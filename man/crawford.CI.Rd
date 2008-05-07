\name{.crawford.CI}
\alias{.crawford.CI}
\title{Confidence Intervals based on a non-central t distribution}
\description{\code{.crawford.CI} is an internal function that computes 95\% confidence intervals around the rarity of a score (Crawford & Garthwaite, 2002). These confidence intervals are based on non-central $t$-distributions. }
\usage{
.crawford.CI(c,n)
}
\arguments{
  \item{c}{test statistics}
  \item{n}{number of control subjects}
}
\details{\code{c} is an observation from a non-central $t$-distribution on \code{n} - 1 degrees of freedom. Confidence intervals are based on two values of the non-centrality parameter of non-central $t$-distribution, such that the resulting non-central $t$-distribution has $c \sqrt{n}$ as its $100 \alpha / 2 $ percentile (value 1) or as its $100 (1-\alpha/2)$ (value 2) (for details, see Crawford & Garthwaite, 2002). The purpose of the \code{.crawford.CI} function is to find those values. 
}
\value{
  \item{2.5\%}{Lower bound of the confidence interval.}
  \item{97.5\%}{Upper bound of the confidence interval.}
}
\references{
Crawford, J. R., & Garthwaite, P. H. (2002). Investigation of the single case in neuropsychology: Confidence limits on the abnormality of test scores and test score differences. \emph{Neuropsychologia}, \bold{40(8)}, 1196--208. 
}
\author{Matthieu Dubois. \email{matthdub@gmail.com}, \url{http://www.code.ucl.ac.be/MatthieuDubois/r_code.html}}
\seealso{\code{\link{crawford.t.test}}, \code{\link{crawford.diff.test}}}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{univar}
