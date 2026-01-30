###                                                             ###
###                		 VCOV-METHOD FOR gekm	     	        ###
###                                                             ###

## vcov.gekm - vcov methode for an object of class gekm
## 
## @param object: gekm[1]
##		object of class gekm
## @param scale: logi[1]
##		scale process variance
## @param ...:
##		further arguments, not used
##
## @output:
##		(scaled) covariance matrix of the regression coefficients


vcov.gekm <- function(object, scale = FALSE, ...){
	
	x <- model.matrix(object)
	A <- backsolve(object$chol$L, x, transpose = TRUE)
	ans <- crossprod(A)

	if(object$derivatives){
		
		dx <- derivModelMatrix(object)		
		dA <- backsolve(object$chol$M, dx - crossprod(object$chol$Q, A), transpose = TRUE)

		ans <- ans + crossprod(dA)
				
	}
	

	ans <- sigma(object, scale = scale)^2 * chol2inv(chol(ans))
	colnames(ans) <- rownames(ans) <- names(coef(object))
	ans
	
}


