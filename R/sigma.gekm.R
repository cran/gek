###                                                             ###
###                		 SIGMA-METHOD FOR gekm	     	        ###
###                                                             ###

## sigma.gekm - sigma methode for an object of class gekm
## 
## @param object: gekm[1]
##		object of class gekm
## @param scale: logi[1]
##		scale process variance
## @param ...:
##		further arguments, not used
##
## @output:
##		(scaled) process standard deviation


sigma.gekm <- function(object, scale = FALSE, ...){
	
	nobs <- nobs(object)
	p <- length(coef(object))
	
	if(scale) sqrt(nobs / (nobs - p - 2L)) * object$sigma else object$sigma
	
}