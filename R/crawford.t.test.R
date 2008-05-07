`crawford.t.test` <-
function(patient, controls, mean.c=0, sd.c=0, n=0, na.rm=FALSE)
{
	#if no summaries are entered, they are computed
	if(missing(n)) 
	{
		n <- length(controls)
		mean.c <- mean(controls, na.rm=na.rm)
		sd.c <- sd(controls, na.rm=na.rm)
	}
	if(!sd.c > 0) stop('standard deviation is not strictly positive')
	dl <- n-1	#degrees of freedom of the  test
	
	#t.test computation
	t.obs <- (patient-mean.c) / (sd.c*sqrt((n+1)/n))
	proba.onetailed <- pt(t.obs, df=dl)
	
	#confidence intervals computation on the rarity (Crawford & Garthwaite, 2002, Neuropsychologia)
	c <- (patient-mean.c)/sd.c
	CI <- .crawford.CI(c,n)
	
	#output
	output <- list(statistic=t.obs, p.value=proba.onetailed, rarity=c(rarity=100*proba.onetailed, CI), df=dl, method=paste("Crawford modified t test with", dl, "degrees of freedom", "\n\np-value is one-tailed", sep=" "))
	class(output)<-"htest"
	return(output)
}

