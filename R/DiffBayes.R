`DiffBayes` <-
function(patient, controls, mean.c=0, sd.c=0 , r=0, n=0, na.rm=FALSE, n.simul=100000, standardized=TRUE) 
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
	
	#sum-of-squares and cross-products matrix
	s.xx <- (sd.c[1]^2) * (n-1)
	s.yy <- (sd.c[2]^2) * (n-1)
	s.xy <- sd.c[1] * sd.c[2] * r * (n-1)
	A <- matrix(c(s.xx, s.xy, s.xy, s.yy), ncol=2)
	
	#C call to compute the stats
	#dyn.load(paste("C", .Platform$dynlib.ext,sep=""))
	result <- .C("DiffBayes_C", as.integer(n), as.integer(n.simul), as.integer(standardized), as.double(patient), as.double(mean.c), as.double(as.vector(A)), out=as.double(rep(0, n.simul)), PACKAGE="singlecase")
	#dyn.unload(paste("C", .Platform$dynlib.ext,sep=""))
	p <- result$out
			
	#outputs
	pval <- mean(p)	
	pval2 <- ifelse(pval>0.5,1-pval,pval) 
	CI <- 100*quantile(p,c(0.025, 0.05, 0.95, 0.975))
	return(list(p.value=2*pval2, rarity=c(rarity=100*pval,CI)))
}

