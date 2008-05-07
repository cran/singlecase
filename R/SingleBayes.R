`SingleBayes` <-
function(patient,controls, mean.c=0, sd.c=0, n=0, na.rm=FALSE, n.simul=100000)
{
	#if no summaries are entered, they are computed
	if(missing(n)) 
	{
		n <- length(controls)
		mean.c <- mean(controls, na.rm=na.rm)
		sd.c <- sd(controls, na.rm=na.rm)
	}
	if(!sd.c > 0) stop('standard deviation is not strictly positive')
	df <- n-1
	var.c <- sd.c^2
	 
	#estimation of the variance (theta)
	psy <- rchisq(n.simul,df)
	theta <- df*var.c/psy
	
	#estimation of the mean (mu)
	z <- rnorm(n.simul)
	mu <- mean.c + z*sqrt(theta/n)
	
	#conditional p-value
	z.star <- (patient - mu) / sqrt(theta)
	p <- pnorm(z.star)
	
	#outputs
	pval <- mean(p)
	CI <- 100*quantile(p,c(0.025,0.05,0.95,0.975))
	output <- list(p.value=pval, rarity=c(rarity=100*pval,CI))
	output	
}

