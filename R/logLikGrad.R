###                                                             ###
###            DERIVATIVE OF  LOG LIKELIHOOD FUCNTION           ###
###                                                             ###


logLikGrad <- function(param, object, ...) UseMethod("logLikGrad", object)

## logLikGrad.default - calculates derivative of log-likelihood at param
##
## @param param: num[d]
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


logLikGrad.default <- function(param, object, x, y, dx = NULL, dy = NULL, 
	covtype = c("matern5_2", "matern3_2", "gaussian"), tolerance = NULL, 
	envir, ...){

	derivatives <- if(is.null(dx) | is.null(dy)) FALSE else TRUE
	covtype <- match.arg(covtype)

	covDevList <- derivBlockCor(x = object, theta = param, covtype = covtype, derivatives = derivatives, envir = envir)

	n <- nrow(object)
	d <- ncol(object)

	if(derivatives){

		dv <- backsolve(envir$M, envir$z[(n + 1):(n * (d + 1)), ])
		v <- backsolve(envir$L, envir$z[1:n, ] - envir$Q %*% dv)

		res <- sapply(1:d, function(h){
			
			Y <- backsolve(envir$L, cbind(covDevList$DK[[h]], t(covDevList$DR[[h]])), transpose = TRUE)
			dY <- backsolve(envir$M, cbind(covDevList$DR[[h]], covDevList$DS[[h]]) - crossprod(envir$Q, Y) , transpose = TRUE)
			
			dZ <- backsolve(envir$M, dY)
			Z <- backsolve(envir$L, Y - envir$Q %*% dZ)
			
			u <- crossprod(v, covDevList$DK[[h]] %*% v + crossprod(covDevList$DR[[h]], dv)) +
					crossprod(dv, covDevList$DR[[h]] %*% v + covDevList$DS[[h]] %*% dv)
						
			-0.5 * (u / envir$sigmaSq - sum(diag(rbind(Z, dZ))))
		})
		
	## derivative of likelihood without derivatives
	}else{

		v <- backsolve(envir$L, envir$z)

		res <- sapply(1:d, function(h){

			Y <- backsolve(envir$L, covDevList$DK[[h]], transpose = TRUE)
			-0.5 * ((crossprod(v, covDevList$DK[[h]] %*% v)) / envir$sigmaSq - 
						sum(diag(backsolve(envir$L, Y))))
		})
	}

	res
}




###                                                             ###
###    		 		 LOGLIKGRAD-METHOD FOR gekm 	     	    ###
###                                                             ###


## logLikGrad.gekm - determine the derivatives of a model matrix
##		for an object of class gekm
## 
## @param object: gekm[1]
##		object of class gekm
## @param ...:
##		further arguments
##
## @output:
##		derivatives of the model matrix




logLikGrad.gekm <- function(param, object, ...){

	mf <- model.frame(object)
	yname <- all.vars(formula(object), max.names = 1L)

	y <- model.response(mf)

	if(object$derivatives){
	
		dx <- derivModelMatrix(object)
		dy <- c(t(object$deriv)) 
	
	}else{
		
		dx <- dy <- NULL
		
	}
	
	X <- as.matrix(object$data[ , setdiff(names(object$data), yname), drop = FALSE])

	env <- new.env()		

	do.call("logLikFun.default", list(param = param, object = X, x = model.matrix(object), y = y,
		dx = dx, dy = dy, covtype = object$covtype, envir = env, ...))
	
    do.call("logLikGrad.default", list(param = param, object = X, x = model.matrix(object), y = y,
		dx = dx, dy = dy, covtype = object$covtype, envir = env, ...))

}

