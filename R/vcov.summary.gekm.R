###                                                             ###
###				  VCOV-METHOD FOR summary.gekm	     	        ###
###                                                             ###

## vcov.summary.gekm - vcov methode for an object of class 
##		summary.gekm
## 
## @param object: gekm[1]
##		object of class gekm
## @param ...:
##		further arguments, not used
##
## @output:
##		(scaled) covariance matrix of the regression coefficients


vcov.summary.gekm <- function(object, ...){
	
	object$cov.scaled
		
}


