\name{DiffTest}
\alias{DiffTest}
\alias{crawford.diff.test}
\alias{DiffBayes}
\title{Difference between patient's scores on two tasks}
\description{These functions test whether the difference between scores on two tests observed for a patient is significantly greater than the differences observed for a small control or normative sample.  \code{crawford.diff.test} implements Crawford's Unstandardized Difference Test and Crawford's Revised Standardized Difference Test (Crawford \& Garthwaite, 2005). \code{DiffBayes} corresponds to a Bayesian alternative (Crawford \& Garthwaite, 2007). }
\usage{
crawford.diff.test(patient, controls, mean.c = 0, sd.c = 0, r = 0, 
	n = 0, na.rm = FALSE, standardized=FALSE)
DiffBayes(patient, controls, mean.c = 0, sd.c = 0, r = 0, n = 0, 
	na.rm = FALSE, n.simul = 1e+05, standardized = TRUE)
}
\arguments{
  \item{patient}{a vector containing the patient's two scores}
  \item{controls}{a \emph{\code{n}*2} matrix/data-frame containing the control subjects raw scores, one column for each task}
  \item{mean.c}{a vector containing the control group's two means, one for each task}
  \item{sd.c}{a vector containing the control group's two standard deviations, one for each task}
  \item{r}{the correlation coefficient between the two tasks in the control group}
  \item{n}{size of the control group}
  \item{na.rm}{a logical value indicating whether \code{NA} values should be stripped before the computation proceeds}
  \item{n.simul}{a numerical value indicating the number of observations generated for the Monte Carlo estimation. Set at 100.000 by default}
   \item{standardized}{a logical value indicating whether the data should be standardized (\code{TRUE}, the default in function \code{DiffBayes}) or not (\code{FALSE}, the default in function \code{crawford.diff.test}).}
}
\details{These functions examine the difference between a patient's score on two tasks, relative to the difference observed between the same tasks within the control group.  \code{crawford.diff.test} implements a modified $t$-test that treat summary statistics from the control group as estimates rather than parameters.  \code{DiffBayes} uses Bayesian Monte Carlo methods as an alternative to the $t$-test.  Both methods may take either the raw data (using the argument \code{controls}) or summary statistics from the control population as inputs.  In this later case, the controls' means (\code{mean.c}), standard deviations (\code{sd.c}), correlation between tasks (\code{r}) and sample size (\code{n}) are required. 

When it is sensible to examine the raw differences between a patient's scores on two tasks against the raw differences in controls, then the argument \code{standardized=FALSE} can be applied.  In this case, the frequentist and Bayesian alternatives are converging (Crawford and Garthwaite, 2007). However, it is usually necessary to standardize the patient's scores on each task (using the data from the controls).  The \code{DiffBayes} function (with argument \code{standardized=TRUE}, the default) should be applied in these circumstances.  This test does not exhibit convergence with its frequentist alternative, the \code{crawford.diff.test} (\code{standardized=TRUE}), which is deprecated. Indeed, the Bayesian alternative has a number of advantages (see Crawford \& Garthwaite, 2007, for details) including the fact that it factors in the uncertainty over the standard deviations of the two tasks used to standardize the patient's scores.  

Those tests further include a point estimate of the rarity of the patient's score, as well as an interval estimate around this quantity (with the exception of \code{crawford.diff.test}, with  argument \code{standardized} set to \code{TRUE}). Rarity (or abnormality) of the score corresponds to the percentage of the population that would exhibit a difference between scores \emph{lower} than the patient's. 
}
\value{
\code{crawford.diff.test} returns a list with class \code{"htest"} containing all the following components. \code{DiffBayes} only returns a list with \code{p-value} and \code{rarity}:
  \item{statistic}{the one-tailed value of the test statistic.}
  \item{df}{the degrees of freedom for the $t$-statistic.}
  \item{p.value}{the two-tailed $p$-value for the test.}
  \item{rarity}{a vector containing the point estimate of the rarity of the score and a 95\% confidence interval around  the rarity estimate (additional 5\% and 95\% bounds are also provided in the case of \code{DiffBayes}).It captures the percentage of the population that will obtain a difference more \emph{extreme} than the patient's and in the same direction. }
  \item{method}{a character string indicating what type of $t$-test was performed.}
}
\references{
Crawford, J. and Garthwaite, P. (2005) Testing for suspected impairments and dissociations in single-case studies in neuropsychology: Evaluation of alternatives using Monte Carlo simulations and revised tests for dissociations. \emph{Neuropsychology}, \bold{19(3)}, 318--31.

Crawford, J. and Garthwaite, P. (2007) Comparison of a single case to a control or normative sample in neuropsychology: Development of a bayesian approach. \emph{Cognitive Neuropsychology}, \bold{24(4)} 343--372.

John Crawford's website: \url{http://www.abdn.ac.uk/\~psy086/dept/SingleCaseMethodology.htm}
}
\author{Matthieu Dubois. \email{matthdub@gmail.com}, \url{http://www.code.ucl.ac.be/MatthieuDubois/r_code.html}}
\seealso{\code{\link{crawford.t.test}}, \code{\link{SingleBayes}}}
\examples{
#Both methods can take either raw data or summary measures as arguments
controls <- as.data.frame(matrix(rnorm(50,100,10),ncol=2))
crawford.diff.test(patient=c(95,105), controls)
crawford.diff.test(patient=c(95,105), mean.c=mean(controls), 
	sd.c=sd(controls), r = cor(controls[,1],controls[,2]), 
	n = dim(controls)[1])

#In the case of unstandardized tests, the methods are converging
# Note that the default is unstandardized for crawford.t.test, 
#	and standardized for DiffBayes
controls <- as.data.frame(matrix(rnorm(50,100,10),ncol=2))
X <- crawford.diff.test(patient=c(95,105), controls)
X
X$rarity
DiffBayes(patient=c(95,105), controls=controls,standardized=FALSE)


#In the case of standardized tests, the methods are not converging
#DiffBayes has to be preferred
crawford.diff.test(patient=c(90,110), mean.c=c(100,100), sd.c=c(10,10),
	n=5, r=.6, standardized=TRUE)
DiffBayes(patient=c(90,110), mean.c=c(100,100), sd.c=c(10,10), n=5, r=.6)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{univar}
\keyword{htest}