
###                                                             ###
###          	 	  LOG LIKELIHOOD FUCNTION    		        ###
###                                                             ###


logLikFun <- function(param, object, ...) UseMethod("logLikFun", object)

## logLikFun.default - calculates the negative condensed/concentrated log-likelihood
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
##		value of negative 'condensed' log-likelihood

logLikFun.default <- function(param, object, x, y, dx = NULL, dy = NULL,
	covtype = c("matern5_2", "matern3_2", "gaussian"), tolerance = NULL, 
	envir = NULL, ...){

	derivatives <- if(is.null(dx) | is.null(dy)) FALSE else TRUE
	covtype <- match.arg(covtype)
	
	#if(any(theta < 0)) return(NA)
	corList <- blockCor(x = object, theta = param, covtype = covtype, derivatives = derivatives)
	corChol <- blockChol(corList$K, corList$R, corList$S, tol = tolerance)
	
	logDetK <- sum(log(diag(corChol$L)))
	
	w <- backsolve(corChol$L, y, transpose = TRUE)
	A <- backsolve(corChol$L, x, transpose = TRUE)
	
	## with derivatives
	if(derivatives){
		
		dw <- backsolve(corChol$M, c(dy - crossprod(corChol$Q, w)), transpose = TRUE)
		w <- c(w, dw)
		dA <- backsolve(corChol$M, dx - crossprod(corChol$Q, A), transpose = TRUE)
		A <- rbind(A, dA)

		const <- nrow(object) * (ncol(object) + 1L)
		
	## without derivatives
	}else{
		const <- nrow(object)
	}
	
	Q <- qr.Q(qr(A))
	U <- tcrossprod(Q)
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
		
		envir$K <- corList$K
		envir$R <- corList$R
		envir$S <- corList$S
		
		envir$z <- z
		envir$sigmaSq <- sigmaSq
	}

	res
}



###                                                             ###
###    		 		 LOGLIKFUN-METHOD FOR gekm 		     	    ###
###                                                             ###


## logLikFun.gekm - determine the derivatives of a model matrix
##		for an object of class gekm
## 
## @param object: gekm[1]
##		object of class gekm
## @param ...:
##		further arguments
##
## @output:
##		derivatives of the model matrix



logLikFun.gekm <- function(param, object, ...){

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
	
    do.call("logLikFun.default", list(param = param, object = X, x = model.matrix(object), y = y,
		dx = dx, dy = dy, covtype = object$covtype, ...))

}

