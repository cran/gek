###                                                             ###
###     	         DERIVATIVES OF MODEL MATRIX 		        ###
###                                                             ###


derivModelMatrix <- function(object, ...) UseMethod("derivModelMatrix")

## derivModelMatrix.default - determine the derivatives of a model matrix
## 
## @param formula: object[1]
##		formula
## @param data: data.frame[n, d]
##		data.frame
## @param ...:
##		further arguments
##
## @output:
##		derivatives of the model matrix

derivModelMatrix.default <- function(object, data, ...){
	
	# mf <- model.frame(formula, data)
	if(!is.data.frame(data)) stop("'data' must be a data.frame")
	if(any(sapply(data, is.factor))) stop("'data' contains at least one 'factor'")
	
	tt <- if(missing(data))
		terms(object)
	else terms(object, data = data)
	# tt <- delete.response(tt)
	
	xs <- if(attr(tt, "response")){
		setdiff(names(data), all.vars(object, max.names = 1L))
	}else{
		names(data)
	}
	# xs <- names(data) #all.vars(reformulate(labels(terms(mf))))
	# xs <- setdiff(names(data), all.vars(formula, max.names = 1L))
	
	t1 <- labels(tt)

	# if(!all(sapply(data, is.numeric))) stop("")
		
	t2 <- gsub("I\\(([^\\)]*)\\)", "\\1", t1)
	t2 <- gsub(":", "*", t2)
	
	if(attr(tt, "intercept")) t2 <- c("1", t2)
	
	## regression functions
	hs <- sapply(t2, function(x) as.formula(paste("~", x)))
	
	## gradients
	grads <- lapply(hs, stats::deriv, namevec = xs)
	
	## Jacobian matrices
	J <- function(xk){
		sapply(grads, function(gr) attr(eval(gr, envir = xk), "gradient"))
	}
	ans <- do.call("rbind", lapply(1:nrow(data), function(i) J(data[i, , drop = FALSE])))
			
	## attributes
	attr(ans, "assign") <- if(1 %in% colnames(ans)) 0:(ncol(ans) - 1) else 1:ncol(ans)
	colnames(ans) <- if(1 %in% colnames(ans)) c("(Intercept)", t1) else t1
	rownames(ans) <- 1:nrow(ans)
	
	ans
} 


###                                                             ###
###    		     DERIVATIVES OF MODEL MATRIX FOR gekm	        ###
###                                                             ###


## derivModelMatrix.gekm - determine the derivatives of a model matrix
##		for an object of class gekm
## 
## @param object: gekm[1]
##		object of class gekm
## @param ...:
##		further arguments
##
## @output:
##		derivatives of the model matrix

derivModelMatrix.gekm <- function(object, ...){

    if (!object$derivatives) 
        warning("derivatives = FALSE in 'object'")
    do.call("derivModelMatrix.default", list(object = formula(object), 
        data = object$data))

}

###                                                             ###
###     	     	   MODEL MATRIX FOR gekm			        ###
###                                                             ###



## model.matrix.gekm - determines the model matrix for an object
##		of class gekm
## 
## @param object: gekm[1]
##		object of class gekm
## @param ...: 
##		arguments passed to model.matrix.default
##
## @output:
##		model matrix

model.matrix.gekm <- function(object, ...){
	
	mf <- model.frame(object)
	model.matrix.default(object, mf, ...)
		
}
