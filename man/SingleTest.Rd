\name{SingleTest}
\alias{SingleTest}
\alias{crawford.t.test}
\alias{SingleBayes}
\title{Single score comparison}
\description{\code{crawford.t.test} and \code{SingleBayes} are for comparison of an individual's score on a single test with the score of a small normative or control sample.  They provide a significance test, point estimate of the abnormality of the individual's score, and confidence limits on the abnormality.  \code{crawford t-test} is based on the $t$ distribution, \code{SingleBayes} is a Bayesian alternative.
}
\usage{
crawford.t.test(patient, controls, mean.c = 0, sd.c = 0, n = 0, na.rm = FALSE)
SingleBayes(patient, controls, mean.c = 0, sd.c = 0, n = 0, 
	na.rm = FALSE, n.simul = 1e+05)
}
\arguments{
  \item{patient}{a vector with the patient's score}
  \item{controls}{a vector with the control subjects raw scores}
  \item{mean.c}{mean of the control group}
  \item{sd.c}{standard deviation of the control group}
  \item{n}{size of the control group}
  \item{na.rm}{a logical value indicating whether \code{NA} values should be stripped before the computation proceeds}
  \item{n.simul}{a numerical value indicating the number of observations generated for the Monte Carlo estimation. Set at 100.000 by default}
}
\details{
\code{crawford.t.test} is a modified $t$-test that compare a patient's score to that of a control sample of small size.  The test treats the mean and standard deviation form the control group as estimates rather than parameters.  It provides a significance test, point estimate of the abnormality of the individual's score, and 95\% confidence limits on the abnormality. Abnormality (or rarity) is defined as the percentage of the population that would obtain a score \emph{lower} than the patient's. 
 
\code{SingleBayes} uses Bayesian Monte Carlo methods to test if a patient's score is sufficiently below the scores of controls such that the null hypothesis that it is an observation from the control population can be rejected (i.e., it is used to test for a deficit).  It also provides a point estimate of the percentage of the control population that would obtain a lower score (i.e., a point estimate of the abnormality of the score) and a 95\% credible interval for this quantity. 

Bayesian and frequentist methods provide converging results. Both methods may take either the raw data (using the argument \code{controls}) or summary statistics from the control population as inputs.  In this later case, the control mean (\code{mean.c}), standard deviation (\code{sd.c}) and sample size (\code{n}) are required. 
}
\value{
\code{crawford.t.test} returns a list with class \code{"htest"} containing all the following components. \code{SingleBayes} only returns a list with \code{p-value} and \code{rarity}:
  \item{statistic}{the one-tailed value of the $t$-statistic.}
  \item{df}{the degrees of freedom for the $t$-statistic.}
  \item{p.value}{the $p$-value of the one-tailed hypothesis that patient's score is \emph{lower} than that of the control participants. }
  \item{rarity}{a vector containing the point estimate of the rarity of the score and a 95\% interval around the rarity estimate (additional 5\% and 95\% bounds are also provided in the case of \code{SingleBayes}).}
  \item{method}{a character string indicating the test used.}
}
\references{
Crawford, J. and Howell, D. (1998) Comparing an individual's test score against norms drived from small samples. \emph{Clinical Neuropsychologist}, \bold{12(4)},482--486.

Crawford, J. and Garthwaite, P. (2002) Investigation of the single case in neuropsychology: Confidence limits on the abnormality of test scores and test score differences. \emph{Neuropsychologia}, \bold{40(8)}, 1196--208.

Crawford, J. and Garthwaite, P. (2007) Comparison of a single case to a control or normative sample in neuropsychology: Development of a bayesian approach. \emph{Cognitive Neuropsychology}, \bold{24(4)}, 343--372.

John Crawford's website: \url{http://www.abdn.ac.uk/\~psy086/dept/SingleCaseMethodology.htm}
}
\author{Matthieu Dubois. \email{matthdub@gmail.com}, \url{http://www.code.ucl.ac.be/MatthieuDubois/r_code.html}}
\seealso{\code{\link{crawford.diff.test}}, \code{\link{DiffBayes}}}
\examples{
#function provided with either raw data or summary statistics
controls <- rnorm(15,50,10)
crawford.t.test(patient=30, controls=controls)
crawford.t.test(patient=30, mean.c = mean(controls), 
	sd.c = sd(controls), n = length(controls))

#results provided by the two methods are converging
X <- crawford.t.test(patient=30, mean.c=50, sd.c=10, n=15)
X
X$rarity
SingleBayes(patient=30, mean.c=50, sd.c=10, n=15)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{univar}
\keyword{htest}
