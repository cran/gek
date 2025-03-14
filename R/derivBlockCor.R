###                                                             ###
###        		 DERIVATIVES OF CORRELATION MATRIX		        ###
###                                                             ###


## derivBlocCor - calculates the derivatives of the correlation
##		matrix with derivatives with respect to the hyperparameters
##		of the correlation function 
##
## @param x: num[n, d]
##		matrix with inputs
## @param theta: num[d]
##		hyperparameters of the correlation function
## @param covtype: chr[1]
##		type of correlation function
## @param derivatives: logi[1]
##		with derivatives?
##
## @output:
##		list with derivatives of correlation matrices



derivBlockCor <- function(x, theta, covtype = c("matern5_2", "matern3_2", "gaussian"), derivatives = FALSE, envir = NULL){

	if(is.data.frame(x)) x <- as.matrix(x)
	if(!is.numeric(x)) stop("is.numeric(x) is not TRUE")

	covtype <- match.arg(covtype)

	n <- nrow(x)
	d <- ncol(x)
	
	if(length(theta) != d) stop("incompatible dimensions, length of 'theta' must be equal to ncol(x)")
	
	if(any(theta < .Machine$double.neg.eps)) stop("at least one entry of 'theta' is zero or negative")
	
	DKs <- lapply(1:d, function(h){
			out <- .C("corMat_dp", as.double(x), as.integer(n),
					as.double(theta), as.character(covtype), 
					as.integer(h), as.double(envir$K),
					ans = double(n * n))
			matrix(out$ans, n, n)
	})
	
	
	## with derivatives
	if(derivatives){
		
		DRs <- lapply(1:d, function(h){
			out <- .C("corMat_dxdp", as.double(x), as.integer(n),
					as.integer(d), as.double(theta), as.character(covtype),
					as.integer(h), as.double(envir$R),
					ans = double(n * n * d))
			matrix(out$ans, n * d, n)
		})

		DSs <- lapply(1:d, function(h){
			out <- .C("corMat_dxdydp", as.double(x), as.integer(n),
						as.integer(d), as.double(theta), as.character(covtype),
						as.integer(h), as.double(envir$R), as.double(envir$S),
						ans = double(n * n * d * d))
			matrix(out$ans, n * d, n * d)
		})
	
	## without derivatives
	}else{
		DRs <- NULL
		DSs <- NULL
	}

	## Output: Einzelteile der Block-Matrix
	return(list("DK" = DKs, "DR" = DRs, "DS" = DSs))
}


