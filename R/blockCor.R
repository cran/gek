
###                                                             ###
###           CORRELATIONMATIX WITH DERIVATIVES		 	        ###
###                                                             ###

blockCor <- function(x, ...) UseMethod("blockCor")


## blockCor - calculates block correlation matrix with or 
##		without derivatives
##
## @param x: num[n, d]
##		matrix with observations 
## @param theta: num[d]
##		vector of hyperparameters
## @param covtype: chr[1]
##		correlation function to be used
## @param derivatives: logi[1]
##		should derivatives be included
##
## @output:
##		list of block matrices
## 		Ktilde <- cbind(rbind(K, R), rbind(t(R), S))

blockCor.default <- function(x, theta, covtype = c("matern5_2", "matern3_2", "gaussian"), derivatives = FALSE, ...){

	if(is.data.frame(x)) x <- as.matrix(x)
	if(!is.numeric(x)) stop("is.numeric(x) is not TRUE")

	covtype <- match.arg(covtype)

	n <- nrow(x)
	d <- ncol(x)
	
	if(length(theta) != d) stop("incompatible dimensions, length of 'theta' must be equal to ncol(x)")
	
	if(any(theta < .Machine$double.neg.eps)) stop("at least one entry of 'theta' is zero or negative")
	
	out <- .C("corMat", as.double(x), as.integer(n),
					as.integer(d), as.double(theta), as.character(covtype),
					ans = double(n * n))
	K <- matrix(out$ans, n, n)
	
	## with derivatives
	if(derivatives){

		out <- .C("corMat_dx_all", as.double(x), as.integer(n), 
					as.integer(d), as.double(theta), as.character(covtype),
					as.double(K), double(n * n),
					ans = double(n * n * d))
		R <- matrix(out$ans, n * d, n)

		out <- .C("corMat_dxdy_all", as.double(x), as.integer(n), 
					as.integer(d), as.double(theta), as.character(covtype),
					as.double(K), ans = double(n * n * d * d))
		S <- matrix(out$ans, n * d, n * d)
		
	## without derivatives
	}else{
		R <- NULL
		S <- NULL
	}
	
	## fullK <- cbind(rbind(K, R), rbind(t(R), S))
	return(list("K" = K, "R" = R, "S" = S))
}


blockCor.gekm <- function(x, ...){
	
	dots <- list(...)
	args <- c("theta", "covtype", "derivatives")
	
	if(length(dots)){
		m <- pmatch(names(dots), args, 0L)
		dots <- c(dots, x[args[-m]])
	}else{
		dots <- x[args]
	}
	
	yname <- all.vars(formula(x), max.names = 1L)	
	X <- as.matrix(x$data[ , setdiff(names(x$data), yname), drop = FALSE])

	do.call("blockCor.default", c(list(x = X), dots))

}
