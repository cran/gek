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
	
	
	res$nobs <- if(derivatives) length(y) * (length(res$theta) + 1L) else length(y)
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







