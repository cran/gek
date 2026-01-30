###                                                             ###
###                	   SUMMARY-METHOD FOR gekm    		        ###
###                                                             ###

## summary.gekm - simulates Gaussian process path given an object
##		of class gekm
## 
## @param object: gekm[1]
##		object of class gekm
## @param nsim: num[1]
##		number of simulated paths to be generated
## @param seed: num[1]
##		NULL
## @param newdata: data.frame[1]
##		
## @param ...:
##		further arguments, not used		
##
## @output:
##		invisible(x)


summary.gekm <- function(object, scale = FALSE, ...){

	
	res <- object[c("call", "terms")]
	p <- length(coef(object))

	vcov.est <- vcov(object, scale = scale)
	se <- sqrt(diag(vcov.est))
    tval <- coef(object) / se
	df <- nobs(object) - p
	
	res$coefficients <- cbind("Estimate" = coef(object), "Std. Error" = se, 
        "t value" = tval, "Pr(>|t|)" = 2 * pt(abs(tval), df, lower.tail = FALSE))
	
	res$sigma <- sigma(object, scale = scale)
	res$df <- df
	res$cov.scaled <- vcov.est
	res$covtype <- object$covtype
	res$theta <- object$theta
	
	class(res) <- "summary.gekm"
	res

}




