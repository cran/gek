
###                                                             ###
###                	  	    FIT-FUNKTION		   	  	        ###
###                                                             ###


## gekm.fit - fitting function for (Gradient-Enhanced) Kriging Model
## 
## @param x: num[n, p]
##		design matrix
## @param y: num[n]
##		vector of observations
## @param dx: num[n*d, p]
##		derivatives of design matrix
## @param dy: num[n*d]
##		vector of observed derivatives
## @param data: num[n, d]
##		data matrix
## @param covtype: chr[1]
##		covariance type
## @param nugget: logi[1]
##		nugget
## @param method: chr[1]
##		optimization methode
## @param optimizer: chr[1]
##		optimization algorithm
## @param start: num[d]
##		initial values for optimization
## @param control: list[]
##		further argument for optimizer
##
## @output:
##		list


gekm.fit <- function(x, y, dx, dy, X, covtype = c("matern5_2", "matern3_2", "gaussian"), 
	theta = NULL, tol = NULL, optimizer = c("L-BFGS-B", "NMKB"), 
	lower = NULL, upper = NULL, start = NULL, ncalls = 20, control = list(), ...){
	
	derivatives <- if(is.null(dx) | is.null(dy)) FALSE else TRUE
	d <- ncol(X)
	
	if(!is.null(theta)){

		optimizer <- NULL
		convergence <- NULL
		message <- NULL
		logLik <- logLikFun(param = theta, object = X, x = x, y = y, dx = dx, dy = dy,
			covtype = covtype, tolerance = tol)

	}else{
	
		if(is.null(lower)) lower <- rep(10^-10, d)
		if(is.null(upper)) upper <- apply(X, 2, function(x) diff(range(x)) * 2)

		if(d == 1L){
			
			opts <- optimise(logLikFun,	lower = lower, upper = upper,
						x = x, y = y, dx = dx, dy = dy, 
						object = X, covtype = covtype, tolerance = tol)
			theta <- opts$minimum
			optimizer <- "Brent"
			convergence <- 0L
			message <- NULL
			logLik <- opts$objective
			
		}else{
		
			if(is.null(start)){
				parinit <- matrix(runif(ncalls * d), nrow = d)
				parinit <- lower + parinit * (upper - lower)
			}else{
				parinit <- matrix(start, nrow = d)
				ncalls <- ncol(parinit)
			}
			
			optMethod <- switch(optimizer,
				"L-BFGS-B" = optim,
				"NMKB" = nmkb)
			
			argList <- list(fn = logLikFun, lower = lower, upper = upper,
				x = x, y = y, dx = dx, dy = dy, object = X, covtype = covtype,
				tolerance = tol, control = control)
						
			if(optimizer == "L-BFGS-B") argList <- c(argList, method = "L-BFGS-B", gr = logLikGrad, envir = new.env())
		
			res <- lapply(1:ncalls, function(i){
				try(do.call(optMethod, c(list(par = parinit[ , i]), argList)), silent = TRUE)
			})
			
			failed <- sapply(res, inherits, "try-error")
			if(all(failed)) stop("maximum likelihood estimation failed, try to specify argument 'tol'")
			else res <- res[!failed]
			
			converged <- sapply(res, function(x) x$convergence) == 0
			res <- if(any(converged)) res[converged] else res
			opts <- res[[which.min(sapply(res, function(x) x$value))]]

			theta <- opts$par
			convergence <- opts$convergence
			message <- opts$message
			logLik <- opts$value
					
		}
	}
	
	names(theta) <- ifelse(grepl("[[:blank:]]", colnames(X)), sprintf("`%s`", colnames(X)), colnames(X))
	
	corList <- blockCor(x = X, theta = theta, covtype = covtype, derivatives = derivatives)	
	corChol <- blockChol(corList$K, corList$R, corList$S, tol = tol)

	# estimate coefficients and process variance
	w <- backsolve(corChol$L, y, transpose = TRUE)
	A <- backsolve(corChol$L, x, transpose = TRUE)
	
	## with derivatives
	if(derivatives){
		dw <- backsolve(corChol$M, c(dy - crossprod(corChol$Q, w)), transpose = TRUE)
		w <- c(w, dw)
		dA <- backsolve(corChol$M, dx - crossprod(corChol$Q, A), transpose = TRUE)
		A <- rbind(A, dA)		
		const <- nrow(X) * (ncol(X) + 1L)
	# without derivatives
	}else{
		const <- nrow(X)
	}	
	
	coef <- .lm.fit(A, w)$coefficients
	names(coef) <- colnames(x)

	Q <- qr.Q(qr(A))
	U <- tcrossprod(Q)
	z <- w - U %*% w
	sigma <- sqrt(drop(crossprod(z)) / const)	
	
	res <- list("coefficients" = coef,
			"sigma" = sigma,
			"theta" = theta,
			"covtype" = covtype,
			"chol" = corChol,
			"optimizer" = optimizer,
			"convergence" = convergence,
			"message" = message,
			"logLik" = logLik)
	
	res
}

