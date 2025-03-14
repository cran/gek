###                                                             ###
###                BLOCK CHOLESKY DECOMPOSITION 	  	        ###
###                                                             ###

## blockChol - calculates block Cholesky decomposition
##
## @param K: num[n, n]
##		covariance matrix K
## @param R: num[nd, n]
##		block matrix R (optional)
## @param S: num[nd, nd]
##		block matrix S (optional)
## @param tol: num[1]
##		tolerance for log(kappa(A + eps I)) <= tol
##
## @output:
##		list of block matrices

blockChol <- function(K, R = NULL, S = NULL, tol = NULL){
	
	# regularization (Ranjan et al., 2011)
	if(!is.null(tol)){
	
		if(is.null(R) | is.null(S)) A <- K
		else A <- cbind(rbind(K, R), rbind(t(R), S))
		s <- svd(A, nu = 0L, nv = 0L)$d
		cond <- max(s) / min(s[s > 0])
		eps <- max(c(max(s) * (cond - exp(tol)) / cond / (exp(tol) - 1), 0))
	
	}else{
		eps <- 0L
	}
	
	## Cholesky decomposition of K
	L <- chol(K + diag(rep(eps, nrow(K))))
	# K == t(L) %*% L
			
	## with derivatives
	if(!is.null(R) & !is.null(S)){
	
		Q <- backsolve(L, t(R), transpose = TRUE)
		## Schur complement
		N <- S - t(Q) %*% Q
		
		## Cholesky decomposition of the Schur complement
		M <- chol(N + diag(rep(eps, nrow(S))))
		# N == t(M) %*% M
		
	## without derivatives
	}else{
		Q <- NULL
		M <- NULL
	}
	
	## Output: 
	## fullL <- cbind(rbind(L, matrix(0, d * n, n)), rbind(Q, M))
	## A <- t(fullL) %*% fullL
	res <- list("L" = L, "Q" = Q, "M" = M)
	attr(res, "eps") <- eps
	
	res
}
