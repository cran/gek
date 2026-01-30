###                                                             ###
###                	 LOG-LIKELIHOOD FOR gekm	  	   	        ###
###                                                             ###

## logLik.gekm - logLik method for gekm object
## 
## @param object: gekm[1]
##		object of class gekm
## @param ...:
##		further arguments, not used	
##
## @output:
##		value of the log-likelihood


logLik.gekm <- function(object, ...){	
	nobs <- nobs(object)

	val <- -object$logLik - 0.5 * nobs * (log(2 * pi) + 1)

	attr(val, "nobs") <- nobs
	if(!is.null(object$optimizer))
		attr(val, "df") <- length(coef(object)) + length(object$theta) + 1
	else
		attr(val, "df") <- length(coef(object)) + 1	
	class(val) <- "logLik"
	return(val)
}
