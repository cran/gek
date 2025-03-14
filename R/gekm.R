###                                                             ###
###               GRADIENT-ENHANCED KRIGING MODEL   	        ###
###                                                             ###

## gekm - fits a (gradient-enhanced) Kriging model
## 
## @param formula: formula[1]
##		formula woth regression functions and response
## @param data: data.frame
##		data frame with n observations and d variables + response
## @param deriv: data.frame
##		data frame with derivatives with n rows and d columns
## @param covtype: chr{[1]
##		correlation function
## @param theta: num[d]
##		vector of correlation parameters
## @param tol: num[1]
##		tolerance for the condition number of the correlation matrix
## @param optimizer: chr[1]
##		optimization algorithm for maximum likelihood estimation
## @param lower: num[d]
##		lower bound for the optimization method
## @param upper: num[d]
##		lower upper for the optimization method
## @param start: num[d]
##		initial values for optimization method
## @param ncalls: num[1]
##		number of randomly drawn initial values for optimization method
## @param control: list[1]
##		list with control arguments for the optimization method
## @param model: logi[1]
##		should model frame be returned
## @param x: logi[1]
##		should model matrix be returned
## @param y: logi[1]
##		should response be returned
## @param dx: logi[1]
##		should derivatives of model matrix be returned
## @param dy: logi[1]
##		should derivatives of response be returned
## @param ...:
##		further arguments, not used		
##
## @output:
##		object of class gekm


	
gekm <- function(formula, data, deriv,
					covtype = c("matern5_2", "matern3_2", "gaussian"),
					theta = NULL, tol = NULL,
					optimizer = c("NMKB", "L-BFGS-B"), 
					lower = NULL, upper = NULL,
					start = NULL, ncalls = 20, control = NULL,
					model = TRUE, x = FALSE, y = FALSE, dx = FALSE, dy = FALSE, ...){

	cl <- match.call()

	covtype <- match.arg(covtype)
	optimizer <- match.arg(optimizer)
	
	ret.x <- x
	ret.y <- y
	ret.dx <- dx
	ret.dy <- dy
	
	derivatives <- !missing(deriv)
	
	if(!is.data.frame(data)) stop("'data' must be a data.frame")

	mf <- model.frame(formula, data = data)

	yname <- all.vars(formula, max.names = 1L)
    y <- model.response(mf, "numeric")

	if(is.null(y)) stop("'formula' must contain a response (left-hand side)")

	# mt <- terms(formula, data = data) 
	mt <- attr(mf, "terms")
	# model matrix
	x <- model.matrix(mf, data = data)			
	
	X <- as.matrix(data[ , setdiff(names(data), yname), drop = FALSE])
	
	## with derivatives
	if(derivatives){

		if(!is.data.frame(deriv)) stop("'deriv' must be a data.frame")
		
		deriv <- get_all_vars(reformulate(sprintf("`%s`", colnames(X))), deriv)
	
		if(!all(dim(X) == dim(deriv))) stop("dimensions of 'data' and 'deriv' do not match")
		
		# derivatives of model matrix
		dx <- derivModelMatrix(formula, data)
		dy <- c(t(deriv))
			
	}else{
		dx <- dy <- NULL
	}
	

	res <- gekm.fit(x = x, y = y, dx = dx, dy = dy, X = X, 
		covtype = covtype, theta = theta, tol = tol, 
		optimizer = optimizer, lower = lower, upper = upper, start = start,
		ncalls = ncalls, control = control, ...)
	
	
	res$derivatives <- derivatives	
	res$data <- data
	if(derivatives) res$deriv <- deriv
	res$call <- cl
	res$terms <- mt
	if(model) res$model <- mf
	if(ret.x) res$x <- x
	if(ret.y) res$y <- y
	if(ret.dx & derivatives) res$dx <- dx
	if(ret.dy & derivatives) res$dy <- dy

	class(res) <- "gekm"
	res
}




###                                                             ###
###                	  PRINT-METHOD FOR gekm 	    	        ###
###                                                             ###

## print.gekm - print method for gekm
## 
## @param x: gekm[1]
##		an object of class gekm
##
## @output:
##		invisible(x)

print.gekm <- function(x, digits = 4L, ...){
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    cat("Coefficients:\n")
        print(format(x$coefficients, digits = digits), print.gap = 2L, 
            quote = FALSE)
    cat("\nStandard Deviation:", format(x$sigma, digits = digits))
	cat("\n\nCorrelation Stucture:", x$covtype)	
	cat("\nCorrelation Parameters:\n")
		print(format(x$theta, digits = digits), print.gap = 2L,
			quote = FALSE)
	cat("\n")
    invisible(x)
}



###                                                             ###
###                	 PREDICT-METHOD FOR gekm	     	        ###
###                                                             ###

## predict.gekm - predict methode for an object of class gekm
## 
## @param object: gekm[1]
##		object of class gekm
## @param newdata: num[m, d]
##		new data.frame for prediction
## @param ...:
##		further arguments, not used
##
## @output:
##		predicted values and standard deviations

predict.gekm <- function(object, newdata, ...){

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
	newX <- as.matrix(newdata[ , colnames(X)])

	w <- y - x %*% coef(object)
	z <- backsolve(object$chol$L, w, transpose = TRUE)		

	out <- .C("corMat2", as.double(X), as.integer(n),
					as.double(newX), as.integer(m),
					as.integer(d), as.double(object$theta), as.character(object$covtype),
					ans = double(n * m))
	newK <- matrix(out$ans, n, m)
		
	k <- backsolve(object$chol$L, newK, transpose = TRUE) 
	
	fit <- drop(newx %*% coef(object) + t(k) %*% z)

	A <- backsolve(object$chol$L, x, transpose = TRUE)

	## with derivatives
	if(object$derivatives){
		
		dx <- derivModelMatrix(object)
		dy <- c(t(object$deriv))

		dw <- dy - dx %*% coef(object)
		dz <- backsolve(object$chol$M, dw - t(object$chol$Q) %*% z, transpose = TRUE)
		out <- .C("corMat2_dx_all",
					as.double(X), as.integer(n),
					as.double(newX), as.integer(m), as.integer(d),
					as.double(object$theta), as.character(object$covtype),
					as.double(newK), double(n * m),
					ans = double(n * d * m))
		newdK <- matrix(out$ans, n * d, m)
		
		dk <- backsolve(object$chol$M, newdK - t(object$chol$Q) %*% k, transpose = TRUE) 
		
		fit <- fit + drop(t(dk) %*% dz)
		
		dA <- backsolve(object$chol$M, dx - t(object$chol$Q) %*% A, transpose = TRUE)

		B <- chol(t(A) %*% A + t(dA) %*% dA)
		
		b <- backsolve(B, t(newx - t(k) %*% A - t(dk) %*% dA), transpose = TRUE)
		
		sd.fit <- object$sigma * sqrt(pmax(1L - colSums(k * k) - colSums(dk * dk) + colSums(b * b), 0L))
	
	}else{
	
		B <- chol(t(A) %*% A)
		b <- backsolve(B, t(newx - t(k) %*% A), transpose = TRUE)

		sd.fit <- object$sigma * sqrt(pmax(1L - colSums(k * k) + colSums(b * b), 0L))
	}
	
	list("fit" = fit, "sd.fit" = sd.fit)

}




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
##		
##
## @output:
##		matrix with simulated paths


simulate.gekm <- function(object, nsim = 1, seed = NULL, newdata = NULL, tol = NULL, ...){
    
	if(!is.null(seed)) stop("argument 'seed' is not supported")
	
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
	newX <- as.matrix(newdata[ , colnames(X)])

	w <- y - x %*% coef(object)
	z <- backsolve(object$chol$L, w, transpose = TRUE)		

	out <- .C("corMat2", as.double(X), as.integer(n),
					as.double(newX), as.integer(m),
					as.integer(d), as.double(object$theta), as.character(object$covtype),
					ans = double(n * m))
	newK <- matrix(out$ans, n, m)
		
	k <- backsolve(object$chol$L, newK, transpose = TRUE) 
	
	fit <- drop(newx %*% coef(object) + t(k) %*% z)

	corList <- blockCor(newX, theta = object$theta, covtype = object$covtype)
	condK <- corList$K - t(k) %*% k

	## with derivatives
	if(object$derivatives){
		
		dx <- derivModelMatrix(object)
		dy <- c(t(object$deriv))

		dw <- dy - dx %*% coef(object)
		dz <- backsolve(object$chol$M, dw - t(object$chol$Q) %*% z, transpose = TRUE)
		out <- .C("corMat2_dx_all",
					as.double(X), as.integer(n),
					as.double(newX), as.integer(m), as.integer(d),
					as.double(object$theta), as.character(object$covtype),
					as.double(newK), double(n * m),
					ans = double(n * d * m))
		newdK <- matrix(out$ans, n * d, m)
		
		dk <- backsolve(object$chol$M, newdK - t(object$chol$Q) %*% k, transpose = TRUE) 
		
		fit <- fit + drop(t(dk) %*% dz)
		
		condK <- condK - t(dk) %*% dk
	
	}
	
	corChol <- blockChol(condK, tol = tol)

	mat <- matrix(rnorm(m * nsim), m, nsim)
	rand <- object$sigma * t(corChol$L) %*% mat

	val <- rand + fit
    colnames(val) <- paste0("sim_", seq_len(nsim))
    val
	
}





###                                                             ###
###                	 FORMULA-METHOD FOR gekm 	    	        ###
###                                                             ###

## formula.gekm - formula method for gekm object
## 
## @param x: gekm[1]
##		object of class gekm
##
## @output:
##		extracted formula


formula.gekm <- function(x, ...){
	formula(x$terms)
}

