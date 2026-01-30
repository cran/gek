###                                                             ###
###                	 PREDICT-METHOD FOR gekm	     	        ###
###                                                             ###

## predict.gekm - predict methode for an object of class gekm
## 
## @param object: gekm[1]
##		object of class gekm
## @param newdata: num[m, d]
##		new data.frame for prediction
## @param sd.fit: logi[1]
##		estimate root mean squared error of prediction?
## @param scale: logi[1]
##		scale process standard deviation?
## @param df: num[1]
##		degreees of freedom
## @param interval: chr[1]
##		calculate confidence intervals?
## @param level: num[1]
##		confidence level
## @param ...:
##		further arguments, not used
##
## @output:
##		predicted values, standard deviations and confidence intervals

predict.gekm <- function(object, newdata, sd.fit = TRUE, scale = TRUE, 
	df = Inf, interval = c("none", "confidence"), level = 0.95, ...){

	interval <- match.arg(interval)

	x <- model.matrix(object)
	y <- model.response(model.frame(object))

	n <- length(y)
	d <- length(object$theta)
	if(missing(newdata)) newdata <- object$data
	m <- nrow(newdata)

	tt <- terms(object, data = newdata)
	# tt <- delete.response(tt)
	newx <- model.matrix(tt, newdata)
	newdx <- derivModelMatrix(tt, newdata)
	
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
	
	## with derivatives
	if(object$derivatives){

		out <- .C("corMat2_dx_all",
				as.double(X), as.integer(n),
				as.double(newX), as.integer(m), as.integer(d),
				as.double(object$theta), as.character(object$covtype),
				as.double(newK), double(n * m),
				ans = double(n * d * m))
		newdK <- matrix(out$ans, n * d, m)

		dx <- derivModelMatrix(object)
		dy <- c(t(object$deriv))

		dw <- dy - dx %*% coef(object)
		dz <- backsolve(object$chol$M, dw - crossprod(object$chol$Q, z), transpose = TRUE)		
		dk <- backsolve(object$chol$M, newdK - crossprod(object$chol$Q, k), transpose = TRUE) 
		
		fit <- fit + drop(crossprod(dk, dz))

	}	
	
	
	# compute standard deviation of prediction
	if(sd.fit || interval == "confidence"){

		A <- backsolve(object$chol$L, x, transpose = TRUE)
		
		if(object$derivatives){
			
			dA <- backsolve(object$chol$M, dx - crossprod(object$chol$Q, A), transpose = TRUE)
			B <- chol(crossprod(A) + crossprod(dA))
			b <- backsolve(B, t(newx - crossprod(k, A) - crossprod(dk, dA)), transpose = TRUE)

			rmse <- pmax(sqrt(1L - colSums(k * k) - colSums(dk * dk) + colSums(b * b)), 0L)
			
		}else{
		
			B <- chol(crossprod(A))
			b <- backsolve(B, t(newx - crossprod(k, A)), transpose = TRUE)

			rmse <- pmax(sqrt(1L - colSums(k * k) + colSums(b * b)), 0L)
			
		}
		
		
		if(scale){
		
			p <- ncol(x)
			df <- if(object$derivatives) n + n * d - p else n - p
			
			sigma <- object$sigma * sqrt((df + p) / df)
		
		}else{
			
			sigma <- object$sigma
		
		}
		
		if(interval == "confidence"){
		
			alpha <- qt((1 - level) / 2, df)
			fit <- cbind(fit, fit + sigma * rmse %o% c(alpha, -alpha))
            colnames(fit) <- c("fit", "lower", "upper")
			
		}
		
	}
	

	gradient <- sapply(1:d, function(i){
		ind <- seq(i, n * d, d)
		dc <- backsolve(object$chol$L, newdK[ind, , drop = FALSE], transpose = TRUE)
		matrix(newdx %*% coef(object), ncol = d, byrow = TRUE)[ , i, drop = FALSE] - t(dc) %*% z
	})
	gradient <- matrix(gradient, ncol = d)



	newdx %*% coef(object)
	
	z
	
	## with derivatives
	if(object$derivatives){
		
		dx <- derivModelMatrix(object)
		dy <- c(t(object$deriv))

		dw <- dy - dx %*% coef(object)
		dz <- backsolve(object$chol$M, dw - crossprod(object$chol$Q, z), transpose = TRUE)		
		dk <- backsolve(object$chol$M, newdK - crossprod(object$chol$Q, k), transpose = TRUE) 
		
		fit <- fit + drop(crossprod(dk, dz))
		
		out <- .C("corMat2_dxdy_all",
					as.double(X), as.integer(n),
					as.double(newX), as.integer(m), as.integer(d),
					as.double(object$theta), as.character(object$covtype),
					as.double(newK),
					ans = double(d * n * d * m))
		newddK <- matrix(-out$ans, n * d, d * m)

		dddk <- backsolve(object$chol$M, newddK - t(object$chol$Q) %*% ddk, transpose = TRUE) 
		grad <- gradient + drop(t(dddk) %*% dz)

		grad <- gradient - sapply(1:d, function(i){
			ind <- seq(i, n * d, d)
			ind2 <- seq(i, m * d, d) # (1:m) + (i - 1) * m 
			dc <- backsolve(object$chol$L, newdK[ind, , drop = FALSE], transpose = TRUE)
			dddk <- backsolve(object$chol$M[ind, ind, drop = FALSE], newddK[ind , (1:m) + (i - 1) * m, drop = FALSE] - t(object$chol$Q)[ind, , drop = FALSE] %*% -dc, transpose = TRUE)
			t(dddk) %*% dz[ind, , drop = FALSE]			
		})
		
		# corList <- blockCor(object)
		
		
		# ind <- seq(1, length(newdata$x), length = 10)
		# ind <- 1:length(newdata$x)
		# plot(newdata$x, fit)
		# plot(newdata$x, fit1, col = "green")
		# tangents(newdata$x, fit1, grad[ , 1], col = "red")
		# points(newdata$x1, fit, col = "red")
		# tangents(newdata$x, fit1, gradient[ , 1], col = "red")
		
		dA <- backsolve(object$chol$M, dx - crossprod(object$chol$Q, A), transpose = TRUE)
		B <- chol(crossprod(A) + crossprod(dA))
		b <- backsolve(B, t(newx - crossprod(k, A) - crossprod(dk, dA)), transpose = TRUE)

		rmse <- object$sigma * sqrt(pmax(1L - colSums(k * k) - colSums(dk * dk) + colSums(b * b), 0L))
	
	}else{
	
		B <- chol(crossprod(A))
		b <- backsolve(B, t(newx - crossprod(k, A)), transpose = TRUE)

		rmse <- object$sigma * sqrt(pmax(1L - colSums(k * k) + colSums(b * b), 0L))
	}
	
	if(sd.fit) list("fit" = fit, "sd.fit" = sigma * rmse) else fit

}




predict.gekm <- function(object, newdata, sd.fit = TRUE, scale = FALSE, 
	df = NULL, interval = c("none", "confidence"), level = 0.95, ...){

	interval <- match.arg(interval)

	x <- model.matrix(object)
	y <- model.response(model.frame(object))

	n <- length(y)
	d <- length(object$theta)
	if(missing(newdata)) newdata <- object$data
	m <- nrow(newdata)

	tt <- terms(object, data = newdata)
	tt <- delete.response(tt)
	newx <- model.matrix(tt, newdata)
	# newdx <- derivModelMatrix(tt, newdata)
	
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
	
	## with derivatives
	if(object$derivatives){

		out <- .C("corMat2_dx_all",
				as.double(X), as.integer(n),
				as.double(newX), as.integer(m), as.integer(d),
				as.double(object$theta), as.character(object$covtype),
				as.double(newK), double(n * m),
				ans = double(n * d * m))
		newdK <- matrix(out$ans, n * d, m)

		dx <- derivModelMatrix(object)
		dy <- c(t(object$deriv))

		dw <- dy - dx %*% coef(object)
		dz <- backsolve(object$chol$M, dw - crossprod(object$chol$Q, z), transpose = TRUE)		
		dk <- backsolve(object$chol$M, newdK - crossprod(object$chol$Q, k), transpose = TRUE) 
		
		fit <- fit + drop(crossprod(dk, dz))

	}	
	
	# compute standard deviation of prediction
	if(sd.fit || interval == "confidence"){

		A <- backsolve(object$chol$L, x, transpose = TRUE)
		
		if(object$derivatives){
			
			dA <- backsolve(object$chol$M, dx - crossprod(object$chol$Q, A), transpose = TRUE)
			B <- chol(crossprod(A) + crossprod(dA))
			b <- backsolve(B, t(newx - crossprod(k, A) - crossprod(dk, dA)), transpose = TRUE)

			rmse <- sqrt(pmax(1L - colSums(k * k) - colSums(dk * dk) + colSums(b * b), 0L))
			
		}else{
		
			B <- chol(crossprod(A))
			b <- backsolve(B, t(newx - crossprod(k, A)), transpose = TRUE)

			rmse <- sqrt(pmax(1L - colSums(k * k) + colSums(b * b), 0L))
			
		}
		

		rmse <- sigma(object, scale = scale) * rmse
		
		if(interval == "confidence"){

			df <- if(is.null(df)) nobs(object) - length(coef(object)) else df
			alpha <- qt((1 - level) / 2, df)
			fit <- cbind(fit, fit + rmse %o% c(alpha, -alpha))
            colnames(fit) <- c("fit", "lower", "upper")
			
		}
		
	}

	if(sd.fit) list("fit" = fit, "sd.fit" = rmse) else fit

}
