###                                                             ###
###            DERIVATIVE OF  LOG LIKELIHOOD FUCNTION           ###
###                                                             ###


## logLikGrad - calculates derivative of log-likelihood at theta
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
##		value of the derivative of the negative 'condensed' log-likelihood


logLikGrad <- function(theta, x, y, dx, dy, X, covtype, tolerance = NULL, envir){

	derivatives <- if(is.null(dx) | is.null(dy)) FALSE else TRUE

	covDevList <- derivBlockCor(x = X, theta = theta, covtype = covtype, derivatives = derivatives, envir = envir)

	n <- nrow(X)
	d <- ncol(X)

	if(derivatives){

		dv <- backsolve(envir$M, envir$z[(n + 1):(n * (d + 1)), ])
		v <- backsolve(envir$L, envir$z[1:n, ] - envir$Q %*% dv)

		res <- sapply(1:d, function(h){
			
			Y <- backsolve(envir$L, cbind(covDevList$DK[[h]], t(covDevList$DR[[h]])), transpose = TRUE)
			dY <- backsolve(envir$M, cbind(covDevList$DR[[h]], covDevList$DS[[h]]) - t(envir$Q) %*% Y , transpose = TRUE)
			
			dZ <- backsolve(envir$M, dY)
			Z <- backsolve(envir$L, Y - envir$Q %*% dZ)
			
			u <- t(v) %*% (covDevList$DK[[h]] %*% v + t(covDevList$DR[[h]]) %*% dv) +
					t(dv) %*% (covDevList$DR[[h]] %*% v + covDevList$DS[[h]] %*% dv)
						
			-0.5 * (u / envir$sigmaSq - sum(diag(rbind(Z, dZ))))
		})
		
	## derivative of likelihood without derivatives
	}else{

		v <- backsolve(envir$L, envir$z)

		res <- sapply(1:d, function(h){

			Y <- backsolve(envir$L, covDevList$DK[[h]], transpose = TRUE)
			-0.5 * ((t(v) %*% covDevList$DK[[h]] %*% v) / envir$sigmaSq - 
						sum(diag(backsolve(envir$L, Y))))
		})
	}

	res
}


