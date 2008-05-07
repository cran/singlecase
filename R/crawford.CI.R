`.crawford.CI` <-
function(c,n)
{
	#finding the non central parameter
	f <- function(delta, pr, x, df) pt(x, df = df, ncp = delta) - pr
	deltaL <- suppressWarnings(try(uniroot(f, lower=-500, upper=500, pr = 0.025, x = c*(n^0.5), df = n-1)))
	deltaU <- suppressWarnings(try(uniroot(f, lower=-500, upper=500, pr = 0.975, x = c*(n^0.5), df = n-1)))
	CI.U <- pnorm(deltaL$root/(n^0.5)) * 100
	CI.L <- pnorm(deltaU$root/(n^0.5)) * 100
	
	output <- c('2.5%' = CI.L, '97.5%' = CI.U)
	return(output)
}