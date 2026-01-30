###                                                             ###
###            	 LEAVE-ONE-OUT-METHOD FOR gekm	     		    ###
###                                                             ###

loo <- function(object, ...) UseMethod("loo")


###                                                             ###
###            	 LEAVE-ONE-OUT-METHOD FOR gekm	     		    ###
###                                                             ###


## loo.gekm - leave-one-out cross validation for an object
##		of class gekm
## 
## @param object: gekm[1]
##		object of class gekm
##		
## @param ...:
##		further arguments, not used		
##
## @output:
##		loo prediction and standard deviation

loo.gekm <- function(object, reestim = TRUE, sd.fit = TRUE, scale = FALSE, 
	df = NULL, interval = c("none", "confidence"), level = 0.95, ...){

	interval <- match.arg(interval)

	x <- model.matrix(object)
	y <- model.response(model.frame(object))

	if(object$derivatives){
		dx <- derivModelMatrix(object)
		dy <- c(t(object$deriv))
	}
	
	n <- length(y)
	d <- length(object$theta)
	p <- ncol(x)

	sigma <- sigma(object, scale = scale)

	if(!reestim){

		fit <- variance <- numeric(n)
	
		corList <- blockCor(object)	
		eps <- attr(object$chol, "eps")
		K <- corList$K + diag(eps, n) # crossprod(object$chol$L)
		R <- corList$R
		S <- corList$S + diag(eps, n * d)

		for(i in 1:n){
		
			Li <- blockChol(K[-i, -i])$L
			
			wi <- backsolve(Li, y[-i], transpose = TRUE)
			Ai <- backsolve(Li, x[-i, , drop = FALSE], transpose = TRUE)
			
			zi <- wi - Ai %*% coef(object)
			ki <- backsolve(Li, K[-i, i], transpose = TRUE)
			
			fit[i] <- x[i, ] %*% coef(object) + crossprod(ki, zi)

			
			if(object$derivatives){
			
				index <- ((i - 1) * d + 1):(i * d)
				
				corChol <- blockChol(K[-i, -i], R[-index, -i], S[-index, -index])
				Qi <- corChol$Q
				Mi <- corChol$M
						
				dwi <- backsolve(Mi, c(dy[-index] - crossprod(Qi, wi)), transpose = TRUE)
				dAi <- backsolve(Mi, dx[-index, , drop = FALSE] - crossprod(Qi, Ai), transpose = TRUE)
	
				zi <- wi - Ai %*%coef(object)
				dzi <- dwi - dAi %*% coef(object)
				dki <- backsolve(Mi, R[-index, i] - crossprod(Qi, ki), transpose = TRUE)
		
				fit[i] <- x[i, ] %*% coef(object) + crossprod(ki, zi) + crossprod(dki, dzi)	
				
				Bi <- chol(crossprod(Ai) + crossprod(dAi))
				ci <- backsolve(Bi, c(x[i, ] - crossprod(ki, Ai) - crossprod(dki, dAi)), transpose = TRUE)
					
				variance[i] <- 1L - crossprod(ki) - crossprod(dki) + crossprod(ci)
				
			}else{
			
				Bi <- chol(crossprod(Ai))
				ci <- backsolve(Bi, c(x[i, ] - crossprod(ki, Ai)), transpose = TRUE)
					
				variance[i] <- 1L - crossprod(ki) + crossprod(ci)

			}
		}
	
		sd.loo <- sigma * sqrt(variance)
		sd.loo <- pmax(sd.loo, 0L)

	}else{	
	
		A <- backsolve(object$chol$L, x, transpose = TRUE)

		if(object$derivatives){
					
			dA <- backsolve(object$chol$M, dx - crossprod(object$chol$Q, A), transpose = TRUE)
			B <- chol(crossprod(A) + crossprod(dA))

			Linv <- chol2inv(object$chol$L)
			Minv <- chol2inv(object$chol$M) # M^{-1} t(M)^{-1}
			Qinv <- backsolve(object$chol$L, -object$chol$Q %*% Minv)
			Kinv <- Linv + tcrossprod(tcrossprod(Qinv, object$chol$M))
			
			z <- Kinv %*% x + Qinv %*% dx
			dz <- crossprod(Qinv, x) + Minv %*% dx
	
			o <- backsolve(B, t(z), transpose = TRUE)
			do <- backsolve(B, t(dz), transpose = TRUE)
			
			Q <- cbind(rbind(Kinv, t(Qinv)), rbind(Qinv, Minv)) - crossprod(cbind(o, do))
			Qy <- (Q %*% c(y, dy)) / sigma^2

			fit <- sd.loo <- numeric(n)

			for(i in 1:n){
			
				index <- n + ((i - 1) * d + 1):(i * d)
				Qi <- chol(Q[c(i, index), c(i, index)])
				Vi <- chol2inv(Qi) * sigma^2
				
				residual <- Vi %*% Qy[c(i, index)]
				fit[i] <- y[i] - residual[1] # (c(y, dy)[c(i, index)] - residual)[1]
				sd.loo[i] <- sqrt(Vi[1, 1])
				
			}

		}else{
			
			B <- chol(crossprod(A))

			Kinv <- chol2inv(object$chol$L)

			o <- backsolve(B, t(Kinv %*% x), transpose = TRUE)
			Q <- Kinv - crossprod(o)
			# Q <- Kinv - Kinv %*% H %*% solve(t(H) %*% Kinv %*% H) %*% t(H) %*% Kinv

			sd.loo <- sigma / sqrt(diag(Q))
			residual <- (Q %*% y) / diag(Q) 
			fit <- c(y - residual)
		
		}

	}
	
	if(interval == "confidence"){
		
		nobs <- if(object$derivatives) nobs(object) - d - 1 else nobs(object) - 1
		df <- if(is.null(df)) nobs - p else df
		alpha <- qt((1 - level) / 2, df)
		fit <- cbind(fit, fit + sd.loo %o% c(alpha, -alpha))
        colnames(fit) <- c("fit", "lower", "upper")
	
	}

	if(sd.fit) list(fit.loo = fit, sd.loo = sd.loo) else fit
	
}

