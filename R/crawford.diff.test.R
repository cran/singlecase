`crawford.diff.test` <-
function(patient, controls, mean.c=0, sd.c=0, r=0, n=0, na.rm=FALSE, standardized=FALSE)
{
	#if no summaries are entered, they are computed
	if(missing(n)) 
	{
		if(!is.data.frame(controls)) controls <- as.data.frame(controls)
		n <- dim(controls)[1]
		mean.c <- mean(controls, na.rm=na.rm)
		sd.c <- sd(controls, na.rm=na.rm)
		
		na.method <- ifelse(na.rm,"complete.obs","all.obs")
		r <- cor(controls[,1], controls[,2], na.method)
	}
	if(!sd.c[1] > 0 || !sd.c[2] > 0) stop('standard deviations are not strictly positive')
	dl <- n-1		#degrees of freedom 
		
	#choice of the method (standardized vs. not)
	# in the case of the standardized test, the Revised version is used (see Crawford & Garthwaite, 2005, Neuropsychology)
	if(standardized)
	{
		method <- "Crawford's Revised Standardised Difference Test"
		
		#computing the test statistic psy
		alpha <- 0.05
		critical <- qt(1-alpha/2, dl)
		num <- ((patient[1]-mean.c[1])/sd.c[1]) - ((patient[2]-mean.c[2])/sd.c[2])
		deno1 <- (n+1)/n
		deno2 <- (2-2*r)
		deno3 <- (2*(1-r^2))/(n-1)
		deno4 <- ((5+critical^2)*(1-r^2)) / (2*((n-1)^2))
		deno5 <- (r*(1+critical^2)*(1-r^2)) / (2*((n-1)^2))

		psy <- num / ((deno1 * (deno2+deno3+deno4+deno5))^0.5)

		#computing the p-value
		a <- (1+r) * (1-r^2)
		b <- (1-r) * ((4*(n-1)^2) + (4*(1+r)*(n-1)) + ((1+r)*(5+r)))
		c <- -2 * num^2 * ((n*(n-1)^2)/(n+1))

		t.obs <- ((-b + ((b^2)-(4*a*c))^0.5) / (2*a))^0.5
		
		#confidence intervals around the rarity (Crawford & Garthwaite, 2002, Neuropsychologia)
		# no interval is computed in standardized mode
		CI <- c('2.5%' = NA, '97.5%' = NA)
	}
	else
	{
		method <- "Crawford's Unstandardised Difference Test"
		
		num <- (patient[1]-mean.c[1]) - (patient[2]-mean.c[2])
		var.diff <- sd.c[1]^2 + sd.c[2]^2 - 2*sd.c[1]*sd.c[2]*r 
		t.obs <- num / ((var.diff * ((n+1)/n))^0.5)	
		
		#confidence intervals around the rarity (Crawford & Garthwaite, 2002, Neuropsychologia)
		c <- num / sqrt(var.diff)
		CI <- .crawford.CI(c,n)
	}
	#computation of the p-value
	p.onetailed <- 1-pt(abs(t.obs), df=dl)
	rar <- pt(t.obs, df=dl)
	p.twotailed <- 2*p.onetailed
		
	#output
	output <- list(statistic=ifelse(standardized, psy, t.obs), p.value=p.twotailed, rarity=c(rarity=100*rar, CI), df=dl, method=paste(method, "with", dl, "degrees of freedom", "\n\ntwo-tailed p-value", sep=" "))
	class(output)<-"htest"
	return(output)
}

