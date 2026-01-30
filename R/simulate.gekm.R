###                                                             ###
###                	 SIMULATE-METHOD FOR gekm	     	        ###
###                                                             ###

## simulate.gekm - simulates Gaussian process path given an object
##		of class gekm
## 
## @param object: gekm[1]
##		object of class gekm
## @param nsim: num[1]
##		number of simulated paths to be generated
## @param seed: num[1]
##		NULL
## @param newdata: data.frame[1]
##		points at which a simulated value should be generated
## @param scale: logi[1]
##		scale process variance
## @param df: num[1]
##		degreees of freedom
## @param tol: num[1]
##		tolerance for condition number of correlation matrix
## @param ...:
##		further arguments, not used	
##
## @output:
##		matrix with simulated paths


simulate.gekm <- function(object, nsim = 1, seed = NULL, newdata = NULL, 
	scale = FALSE, df = NULL, tol = NULL, ...){
    
	if(!is.null(seed)) stop("argument 'seed' is not supported")
	if(is.null(newdata)) stop("'newdata' must be provided in order to simulate process paths")

	x <- model.matrix(object)
	y <- model.response(model.frame(object))

	n <- length(y)
	d <- length(object$theta)
	m <- nrow(newdata)

	tt <- terms(object, data = newdata)
	tt <- delete.response(tt)
	newx <- model.matrix(tt, newdata)

	data <- object$data
	yname <- all.vars(formula(object), max.names = 1L)
	X <- as.matrix(data[ , setdiff(names(data), yname), drop = FALSE])
	newX <- as.matrix(newdata[ , colnames(X), drop = FALSE])

	w <- y - x %*% coef(object)
	z <- backsolve(object$chol$L, w, transpose = TRUE)		

	out <- .C("corMat2", as.double(X), as.integer(n),
					as.double(newX), as.integer(m),
					as.integer(d), as.double(object$theta), as.character(object$covtype),
					ans = double(n * m))
	newK <- matrix(out$ans, n, m)
		
	k <- backsolve(object$chol$L, newK, transpose = TRUE) 
	
	fit <- drop(newx %*% coef(object) + crossprod(k, z))

	corList <- blockCor(newX, theta = object$theta, covtype = object$covtype)
	condK <- corList$K - crossprod(k)

	## with derivatives
	if(object$derivatives){
		
		dx <- derivModelMatrix(object)
		dy <- c(t(object$deriv))

		dw <- dy - dx %*% coef(object)
		dz <- backsolve(object$chol$M, dw - crossprod(object$chol$Q, z), transpose = TRUE)
		out <- .C("corMat2_dx_all",
					as.double(X), as.integer(n),
					as.double(newX), as.integer(m), as.integer(d),
					as.double(object$theta), as.character(object$covtype),
					as.double(newK), double(n * m),
					ans = double(n * d * m))
		newdK <- matrix(out$ans, n * d, m)
		
		dk <- backsolve(object$chol$M, newdK - crossprod(object$chol$Q, k), transpose = TRUE) 
		
		fit <- fit + drop(crossprod(dk, dz))
		
		condK <- condK - crossprod(dk)
	
	}
	
	corChol <- blockChol(condK, tol = tol)

	# mat <- matrix(rnorm(m * nsim), m, nsim)
	df <- if(is.null(df)) nobs(object) - length(coef(object)) else df
	mat <- matrix(rt(m * nsim, df = df), m, nsim)
	rand <- sigma(object, scale = scale) * crossprod(corChol$L, mat) * if(is.finite(df)) sqrt((df - 2L) / df) else 1L
	val <- rand + fit
    colnames(val) <- paste0("sim_", seq_len(nsim))
    val
	
}
