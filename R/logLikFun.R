
###                                                             ###
###          	 	  LOG LIKELIHOOD FUCNTION    		        ###
###                                                             ###

## logLikFun - calculates the negative condensed/concentrated log-likelihood
##
## @param theta: num[d]
##		hyperparameters of correlation function
## @param x: num[n, p]
##		model matrix
## @param y: num[n]
##		response vector
## @param dx: num[nd, p]
##		derivative of model matrix
## @param dy: num[nd]
##		derivative of response
## @param X: num[n, d]
##		matrix of inputs
## @param covtype: chr[1]
##		type of correlation function
## @param tolerance: num[1]
##		tolerance for regularization
## @param envir: env[1]
##		Environemt
##
## @output:
##		value of negative 'condensed' log-likelihood



logLikFun <- function(theta, x, y, dx, dy, X, covtype, tolerance = NULL, envir = NULL){

	derivatives <- if(is.null(dx) | is.null(dy)) FALSE else TRUE
	
	#if(any(theta < 0)) return(NA)
	corList <- blockCor(x = X, theta = theta, covtype = covtype, derivatives = derivatives)
	corChol <- blockChol(corList$K, corList$R, corList$S, tol = tolerance)
	
	logDetK <- sum(log(diag(corChol$L)))
	
	w <- backsolve(corChol$L, y, transpose = TRUE)
	A <- backsolve(corChol$L, x, transpose = TRUE)
	
	## with derivatives
	if(derivatives){
		
		dw <- backsolve(corChol$M, c(dy - t(corChol$Q) %*% w), transpose = TRUE)
		w <- c(w, dw)
		dA <- backsolve(corChol$M, dx - t(corChol$Q) %*% A, transpose = TRUE)
		A <- rbind(A, dA)

		const <- nrow(X) * (ncol(X) + 1L)
		
	## without derivatives
	}else{
		const <- nrow(X)
	}
	
	Q <- qr.Q(qr(A))
	U <- Q %*% t(Q)
	z <- w - U %*% w
	sigmaSq <- drop(crossprod(z)) / const
	
	
	## condensed/concentrated log-likelihood with derivatives
	if(derivatives){
	
		logDetN <- sum(log(diag(corChol$M)))
		res <- 0.5 * const * log(sigmaSq) + logDetK + logDetN
			
	## condensed/concentrated log-likelihood without derivatives	
	}else{
	
		res <- 0.5 * const * log(sigmaSq) + logDetK

	}
	
	if(!is.null(envir)){
		envir$L <- corChol$L
		envir$Q <- corChol$Q
		envir$M <- corChol$M
		
		envir$z <- z
		envir$sigmaSq <- sigmaSq
	}

	res
}
