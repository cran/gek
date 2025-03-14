###                                                             ###
###					2-DIMENSIONAL TESTFUNCTIONS	     	        ###
###                                                             ###



###				 	 	BRANIN-HOO FUNCTION		  		        ###

## input domain: [-5, 10] x [0, 15]
## global minimum: f(x) = 0.3997887
## at x = (-pi, 12.2756), x = (pi, 2.475) and x = (9.42478, 2.475)

branin <- function(x){
	if(is.vector(x)) x <- matrix(x, nrow = 1)
	d <- ncol(x)
	if(d != 2) stop("Branin-Hoo function is only defined for 2 dimensions")
	
	x1 <- x[ , 1]
	x2 <- x[ , 2]	

	(x2 - 5.1 * x1^2 / (4 * pi^2) + 5 * x1 / pi - 6)^2 + 10 * (1 - 1/(8 * pi)) * cos(x1) + 10
}

braninGrad <- function(x){
	if(is.vector(x)) x <- matrix(x, nrow = 1)
	d <- ncol(x)
	if(d != 2) stop("Branin-Hoo function is only defined for 2 dimensions")
	
	x1 <- x[ , 1]
	x2 <- x[ , 2]	
	
	res <- cbind((-10.2 * x1 / (4 * pi^2) + 5 / pi) * 2 * (x2 - 5.1 * x1^2 / (4 * pi^2) + 5 * x1 / pi - 6) - 10 * (1 - 1/(8 * pi)) * sin(x1),
		2 * (x2 - 5.1 * x1^2 / (4 * pi^2) + 5 * x1 / pi - 6))
	colnames(res) <- colnames(x)
	res
}



###				  THREE-HUMP CAMELBACK FUNCTION		  	        ###


## input domain: [-5, 5] x [-5, 5]
## input domain: [-2, 2] x [-2, 2]
## global minimum: f(x) = 0
## at x = (0, 0)

camel3 <- function(x){
	if(is.vector(x)) x <- matrix(x, nrow = 1)
	d <- ncol(x)
	if(d != 2) stop("Three-hump camel function is only defined for 2 dimensions")
	
	x1 <- x[ , 1]
	x2 <- x[ , 2]
	2 * x1^2 - 1.05 * x1^4 + x1^6 / 6 + x1 * x2 + x2^2
}

camel3Grad <- function(x){
	if(is.vector(x)) x <- matrix(x, nrow = 1)
	d <- ncol(x)
	if(d != 2) stop("Three-hump camel function is only defined for 2 dimensions")
	
	x1 <- x[ , 1]
	x2 <- x[ , 2]
	res <- cbind(4 * x1 - 4.2 * x1^3 + x1^5 + x2, x1 + 2 * x2)
	colnames(res) <- colnames(x)
	res	
}


###				  SIX-HUMP CAMELBACK FUNCTION		  	        ###

## input domain: [-3, 3] x [-2, 2]
## input domain: [-2, 2] x [-1, 1]
## global minimum: f(x) = -1.0316
## at x = (0.0898, -0.7126) and x = (-0.0898, 0.7126)

camel6 <- function(x){
	if(is.vector(x)) x <- matrix(x, nrow = 1)
	d <- ncol(x)
	if(d != 2) stop("Six-hump camel function is only defined for 2 dimensions")
	
	x1 <- x[ , 1]
	x2 <- x[ , 2]	
	(4 - 2.1 * x1^2 + x1^4 / 3) * x1^2 + x1 * x2 + (-4 + 4 * x2^2) * x2^2
}

camel6Grad <- function(x){
	if(is.vector(x)) x <- matrix(x, nrow = 1)
	d <- ncol(x)
	if(d != 2) stop("Six-hump camel function is only defined for 2 dimensions")
	
	x1 <- x[ , 1]
	x2 <- x[ , 2]	
	res <- cbind(8 * x1 - 8.4 * x1^3 + 2 * x1^5 + x2, x1 - 8 * x2 + 16 * x2^3)
	colnames(res) <- colnames(x)
	res
}



###				 	 	HIMMELBLAU'S FUNCTION	  		        ###

## input domain: [-5, 5]^2
## local maximum: f(x) = -181.617
## at x = (-0.270845, -0.923039)

himmelblau <- function(x){
	if(is.vector(x)) x <- matrix(x, nrow = 1)
	d <- ncol(x)
	if(d != 2) stop("Himmelblau's function is only defined for 2 dimensions")
	
	x1 <- x[ , 1]
	x2 <- x[ , 2]	
	
	(x1^2 + x2 - 11)^2 + (x1 + x2^2 - 7)^2
}

himmelblauGrad <- function(x){
	if(is.vector(x)) x <- matrix(x, nrow = 1)
	d <- ncol(x)
	if(d != 2) stop("Himmelblau's function is only defined for 2 dimensions")
	
	x1 <- x[ , 1]
	x2 <- x[ , 2]	
	
	res <- cbind(4 * x1 * (x1^2 + x2 - 11) + 2 * (x1 + x2^2 - 7),
		2 * (x1^2 + x2 - 11) + 4 * x2 * (x1 + x2^2 - 7))
	colnames(res) <- colnames(x)
	res
}




###                                                             ###
###				MULTIDIMENSIONAL TESTFUNCTIONS	     	        ###
###                                                             ###


###					ROSENBROCK (BANANA) FUNCTION	  	        ###

## input domain: [-5, 10]^d or [-2.048, 2.048]^d
## global minimum: f(x) = 0
## at x = (1, ..., 1) 

banana <- function(x){
	if(is.vector(x)) stop("Rosenbrock-Banana function is only defined for dimensions > 1")	

	d <- ncol(x)	
	
	rowSums(100 * (x[ , 2:d, drop = FALSE] - x[ , 1:(d - 1), drop = FALSE]^2)^2 + (1 - x[ , 1:(d - 1), drop = FALSE])^2)
}

bananaGrad <- function(x){
	if(is.vector(x)) stop("Rosenbrock-Banana function is only defined for dimensions > 1")

	d <- ncol(x)	
	
	if(d == 2){
		res <- cbind(-400 * (x[ , 2] - x[ , 1]^2) * x[ , 1] + 2 * (x[ , 1] - 1), 
						200 * (x[ , 2] - x[ , 1]^2))
	}else{
		res <- cbind(-400 * (x[ , 2] - x[ , 1]^2) * x[ , 1] + 2 * (x[ , 1] - 1),
			200 * (x[ , 2:(d - 1)] - x[ , 1:(d - 2)]^2) - 
			400 * x[ , 2:(d - 1)] * (x[ , 3:d] - x[ , 2:(d - 1)]^2) + 
			2 * (x[ , 2:(d - 1)] - 1),
			200 * (x[ , d] - x[ , d - 1]^2))	
	}
	
	colnames(res) <- colnames(x)
	res
}




###						SPHERE FUNCTION	   			  	        ###

## input domain: [-5.12, 5.12]^d
## global minimum: f(x) = 0
## at x = (0, ..., 0) 

sphere <- function(x){
	if(is.vector(x)) x <- matrix(x, ncol = 1)	
	
	rowSums(x^2)
}

sphereGrad <- function(x){
	if(is.vector(x)) x <- matrix(x, ncol = 1)
	
	2 * x
}



###						BENT CIGAR FUNCTION	   			  	    ###

## input domain: [-100, 100]^d
## global minimum: f(x) = 0
## at x = (0, ..., 0) 

cigar <- function(x){
	if(is.vector(x)) x <- matrix(x, ncol = 1)	
	d <- ncol(x)

	if(d == 1) return(x^2)
	x[ , 1]^2 + 10^6 * rowSums(x[ , 2:d, drop = FALSE]^2)
}

cigarGrad <- function(x){
	if(is.vector(x)) x <- matrix(x, ncol = 1)	
	d <- ncol(x)

	if(d == 1) res <- 2 * x
	else res <- cbind(2 * x[ , 1], 2 * 10^6 * x[ , 2:d])
	
	colnames(res) <- colnames(x)
	res
}


###						RASTRIGIN FUNCTION	  	  	    	    ###

## input domain: [-5.12, 5.12]^d
## global minimum: f(x) = 0
## at x = (0, ..., 0) 


rastrigin <- function(x){
	if(is.vector(x)) x <- matrix(x, ncol = 1)	
	d <- ncol(x)
	
	10 * d + rowSums(x^2 - 10 * cos(2 * pi * x))
}

rastriginGrad <- function(x){
	if(is.vector(x)) x <- matrix(x, ncol = 1)	
	
	2 * x + 20 * pi * sin(2 * pi * x)
}


###						SCHWEFEL-FUNCTION	  		  	        ###

## input domain: [-500, 500]^d
## global minimum: f(x) = 0
## at x = (420.9687, ..., 420.9687) 

schwefel <- function(x){
	if(is.vector(x)) x <- matrix(x, ncol = 1)	
	d <- ncol(x)
	
	418.982887272433799807913601398 * d - rowSums(x * sin(sqrt(abs(x))))
}

schwefelGrad <- function(x){
	if(is.vector(x)) x <- matrix(x, ncol = 1)	
	
	res <- -sin(sqrt(abs(x))) - x^2 * cos(sqrt(abs(x))) / (2 * abs(x)^(3/2))
	res[x == 0] <- 0
	res
}



###					STYBLINSKI-TANG FUNCTION	  	  	        ###

## input domain: [-5, 5]^d
## global minimum: f(x) = -29.16599d
## at x = (-2.903534, ..., -2.903534) 


styblinski <- function(x){
	if(is.vector(x)) x <- matrix(x, ncol = 1)	
	
	0.5 * rowSums(x^4 - 16 * x^2 + 5 * x)
}

styblinskiGrad <- function(x){
	if(is.vector(x)) x <- matrix(x, ncol = 1)	
	
	2 * x^3 - 16 * x + 2.5
}


###						 QING FUNCTION	  			  	        ###

## input domain: [-500, 500]^d
## global minimum: f(x) = 0
## at x = (-sqrt(1), ..., -sqrt(d))
## and x =  (sqrt(1), ..., sqrt(d))


qing <- function(x){
	if(is.vector(x)) x <- matrix(x, ncol = 1)	
	d <- ncol(x)
	n <- nrow(x)
	
	rowSums((x^2 - tcrossprod(rep(1, n), 1:d))^2)
	
}

qingGrad <- function(x){
	if(is.vector(x)) x <- matrix(x, ncol = 1)	
	d <- ncol(x)
	n <- nrow(x)	
	
	4 * x * (x^2 - tcrossprod(rep(1, n), 1:d)) 
}



###						 GRIEWANK FUNCTION	  		  	        ###

## input domain: [-600, 600]^d
## global minimum: f(x) = 0
## at x = (0,...,0)


griewank <- function(x){
	if(is.vector(x)) x <- matrix(x, ncol = 1)	
	d <- ncol(x)
	n <- nrow(x)
	
	mat <- sqrt(tcrossprod(rep(1, n), 1:d))
	rowSums(x^2 / 4000) - apply(cos(x / mat), 1, prod) + 1
}

griewankGrad <- function(x){
	if(is.vector(x)) x <- matrix(x, ncol = 1)	
	d <- ncol(x)
	n <- nrow(x)	
	
	mat <- sqrt(tcrossprod(rep(1, n), 1:d))
	x / 2000 + sin(x / mat) / mat * tcrossprod(apply(cos(x / mat), 1, prod), rep(1, d)) / cos(x / mat)
}


###                                                             ###
###					UNCERTAINTY TESTFUNCTIONS	     	        ###
###                                                             ###


###						BOREHOLE-FUNCTION	  		  	        ###

## input domain: (x[1], ..., x[8])
## rw %in% [0.05, 0.15]		radius of borehole (m)
## r %in% [100, 50000]		radius of influence (m)
## Tu %in% [63070, 115600]	transmissivity of upper aquifer (m2/yr)
## Hu %in% [990, 1110]		potentiometric head of upper aquifer (m)
## Tl %in% [63.1, 116]		transmissivity of lower aquifer (m2/yr)
## Hl %in% [700, 820]		potentiometric head of lower aquifer (m)
## L %in% [1120, 1680]		length of borehole (m)
## Kw %in% [9855, 12045] 	hydraulic conductivity of borehole (m/yr)
##
## input distributions for uncertainty analysis
## rw ~ N(0.10, 0.0161812)
## r ~ Lognormal(7.71, 1.0056)   
## Tu ~ Uniform[63070, 115600]
## Hu ~ Uniform[990, 1110]
## Tl ~ Uniform[63.1, 116]
## Hl ~ Uniform[700, 820]
## L ~ Uniform[1120, 1680]
## Kw ~ Uniform[9855, 12045]
## remark: N(mu, sigma) not sigma^2
##
## output/response: 
## water flow rate in m^2/yr

borehole <- function(x){
	if(is.vector(x)) x <- matrix(x, nrow = 1)	
	d <- ncol(x)
	if(d != 8) stop("Borehole function is only defined for 8 dimensions")
	
	r_w <- x[ , 1]
	r  <- x[ , 2]
	T_u <- x[ , 3]
	H_u <- x[ , 4]
	T_l <- x[ , 5]
	H_l <- x[ , 6]
	L  <- x[ , 7]
	K_w <- x[ , 8]

	numerator <- 2 * pi * T_u * (H_u - H_l)
	denominator <- log(r / r_w) * (1 + 2 * L * T_u / (log(r / r_w) * r_w^2 * K_w) + T_u / T_l)

	numerator / denominator
}


boreholeGrad <- function(x){
	if(is.vector(x)) x <- matrix(x, nrow = 1)	
	d <- ncol(x)
	if(d != 8) stop("Borehole function is only defined for 8 dimensions")
	
	r_w <- x[ , 1]
	r  <- x[ , 2]
	T_u <- x[ , 3]
	H_u <- x[ , 4]
	T_l <- x[ , 5]
	H_l <- x[ , 6]
	L  <- x[ , 7]
	K_w <- x[ , 8]
	
	denominator <- log(r / r_w) * r_w^2 * K_w * (T_u + T_l) + 2 * L * T_l * T_u
	
	
	cbind("r_w" = 2 * pi * r_w * K_w * T_l * T_u * (H_u - H_l) * (4 * L * T_l * T_u + r_w^2 * K_w * T_u + r_w^2 * K_w * T_l) / denominator^2,
		"r" = -2 * pi * r_w^4 * K_w^2 * T_l * T_u * (T_l + T_u) * (H_u - H_l) / (r * denominator^2),
		"T_u" = 2 * pi * r_w^4 * K_w^2 * T_l^2 * (H_u - H_l) * log(r / r_w) / denominator^2,
		"H_u" = 2 * pi * r_w^2 * K_w * T_l * T_u / denominator,
		"T_l" = 2 * pi * r_w^4 * K_w^2 * T_u^2 * (H_u - H_l) * log(r / r_w) / denominator^2,
		"H_l" = -2 * pi * r_w^2 * K_w * T_l * T_u / denominator,
		"L" = -4 * pi * r_w^2 * K_w * T_l^2 * T_u^2 * (H_u - H_l) / denominator^2,
		"K_w" = 4 * pi * r_w^2 * L * T_l^2 * T_u^2 * (H_u - H_l) / denominator^2)
}



###						 SULFUR-FUNCTION	  		  	        ###

## input distributions for uncertainty analysis
## rw ~ N(0.10, 0.0161812)
## r ~ Lognormal(7.71, 1.0056)   
## Tu ~ Uniform[63070, 115600]
## Hu ~ Uniform[990, 1110]
## Tl ~ Uniform[63.1, 116]
## Hl ~ Uniform[700, 820]
## L ~ Uniform[1120, 1680]
## Kw ~ Uniform[9855, 12045]
## remark: N(mu, sigma) not sigma^2
##
## output/response: 
## water flow rate in m^2/yr

sulfur <- function(x, S_0 = 1366, A = 5.1e+14){
	if(is.vector(x)) x <- matrix(x, nrow = 1)	
	d <- ncol(x)
	if(d != 9) stop("Sulfur function is only defined for 9 dimensions")
	
	Q <- x[ , 1]
	Y  <- x[ , 2]
	L <- x[ , 3]
	Psi <- x[ , 4]
	beta <- x[ , 5]
	f_Psi <- x[ , 6]
	T <- x[ , 7]
	A_c <- x[ , 8]
	R_s <- x[ , 9]

	-1/2 * S_0^2 * A_c * T^2 * R_s^2 * beta * Psi * f_Psi * 3 * Q * Y * L / A
}




sulfurGrad <- function(x, S_0 = 1366, A = 5.1e+14){
	if(is.vector(x)) x <- matrix(x, nrow = 1)	
	d <- ncol(x)
	if(d != 9) stop("Sulfur function is only defined for 9 dimensions")
	
	Q <- x[ , 1]
	Y  <- x[ , 2]
	L <- x[ , 3]
	Psi <- x[ , 4]
	beta <- x[ , 5]
	f_Psi <- x[ , 6]
	T <- x[ , 7]
	A_c <- x[ , 8]
	R_s <- x[ , 9]
		
	cbind("Q" = -1/2 * S_0^2 * A_c * T^2 * R_s^2 * beta * Psi * f_Psi * 3 * Y * L / A,
		"Y" = -1/2 * S_0^2 * A_c * T^2 * R_s^2 * beta * Psi * f_Psi * 3 * Q * L / A,
		"L" = -1/2 * S_0^2 * A_c * T^2 * R_s^2 * beta * Psi * f_Psi * 3 * Q * Y / A,
		"Psi" = -1/2 * S_0^2 * A_c * T^2 * R_s^2 * beta * f_Psi * 3 * Q * Y * L / A,
		"beta" = -1/2 * S_0^2 * A_c * T^2 * R_s^2 * Psi * f_Psi * 3 * Q * Y * L / A,
		"f_Psi" = -1/2 * S_0^2 * A_c * T^2 * R_s^2 * beta * Psi * 3 * Q * Y * L / A,
		"T" = -S_0^2 * A_c * T * R_s^2 * beta * Psi * f_Psi * 3 * Q * Y * L / A,
		"1 - A_c" = -1/2 * S_0^2 * T^2 * R_s^2 * beta * Psi * f_Psi * 3 * Q * Y * L / A,
		"1 - R_s" = -S_0^2 * A_c * T^2 * R_s * beta * Psi * f_Psi * 3 * Q * Y * L / A)
}



###						SHORT COLUMN-FUNCTION		  	        ###

## input distributions for uncertainty analysis 
## Y ~ Lognormal(5, 0.5)	yield stress
## M ~ N(2000, 400)			bending moment
## P ~ N(500, 100)			axial force
## remark: 	N(mu, sigma) not sigma^2
## 			The input variables are uncorrelated,
##			except for a correlation coefficient of 0.5 between M and P. 
##
## output/response: 
## 		limit state function

short <- function(x, b = 5, h = 15){
	if(is.vector(x)) x <- matrix(x, nrow = 1)	
	d <- ncol(x)
	if(d != 3) stop("Short column function is only defined for 3 dimensions")

	Y <- x[ , 1]
	M <- x[ , 2]
	P <- x[ , 3]
	
	1 - 4 * M / (b * h^2 * Y) - P^2 / (b^2 * h^2 * Y^2)
}


shortGrad <- function(x, b = 5, h = 15){
	if(is.vector(x)) x <- matrix(x, nrow = 1)	
	d <- ncol(x)
	if(d != 3) stop("Short column function is only defined for 3 dimensions")

	Y <- x[ , 1]
	M <- x[ , 2]
	P <- x[ , 3]
		
	cbind("Y" = 2 * (2 * b * M * Y + P^2) / (b^2 * h^2 * Y^3), 
		"M" = -4 / (b * h^2 * Y),
		"P" = -2 * P / (b^2 * h^2 * Y^2))
}



###					STEEL COLUMN-FUNCTION			  	        ###

## input distributions for uncertainty analysis 
## Fs ~ Lognormal(400, 35) 		yield stress (MPa)
## P1 ~ N(500000, 50000) 		dead weight load (N)
## P2 ~ Gumbel(600000, 90000)  	variable load (N)
## P3 ~ Gumbel(600000, 90000)	variable load (N)
## B ~ Lognormal(b, 3) 			flange breadth (mm)
## D ~ Lognormal(t, 2) 			flange thickness (mm)
## H ~ Lognormal(h, 5) 			profile height (mm)
## F0 ~ N(30, 10) 				initial deflection (mm)
## E ~ Weibull(210000, 4200)    Young's modulus (MPa)
## remark: 	N(mu, sigma) not sigma^2
## 			The input variables are uncorrelated,
##			except for a correlation coefficient of 0.5 between M and P. 
##
## output/response: 
## trade-off between cost and reliability for a steel column

steel <- function(x, L = 7500){
	if(is.vector(x)) x <- matrix(x, nrow = 1)	
	d <- ncol(x)
	if(d != 9) stop("Steel column function is only defined for 9 dimensions")

	Fs <- x[ , 1]
	P1  <- x[ , 2]
	P2 <- x[ , 3]
	P3 <- x[ , 4]
	B <- x[ , 5]
	D <- x[ , 6]
	H  <- x[ , 7]
	F0 <- x[ , 8]
	E <- x[ , 9]
	
	P <- P1 + P2 + P3
	Eb <- pi^2 * E * B * D * H^2 / (2 * L^2)
	
	Fs - P * (1 / (2 * B * D) + F0 * Eb / (B * D * H * (Eb - P)))
}


steelGrad <- function(x, L = 7500){
	if(is.vector(x)) x <- matrix(x, nrow = 1)	
	d <- ncol(x)
	if(d != 9) stop("Steel column function is only defined for 9 dimensions")

	Fs <- x[ , 1]
	P1  <- x[ , 2]
	P2 <- x[ , 3]
	P3 <- x[ , 4]
	B <- x[ , 5]
	D <- x[ , 6]
	H  <- x[ , 7]
	F0 <- x[ , 8]
	E <- x[ , 9]

	P <- P1 + P2 + P3
	Eb <- pi^2 * E * B * D * H^2 / (2 * L^2)
	
	cbind("Fs" = 1, 
		"P1" = -(1 / (2 * B * D) + F0 * Eb / (B * D * H * (Eb - P))) - P * F0 * Eb / (B * D * H * (Eb - P)^2),
		"P2" = -(1 / (2 * B * D) + F0 * Eb / (B * D * H * (Eb - P))) - P * F0 * Eb / (B * D * H * (Eb - P)^2),
		"P3" = -(1 / (2 * B * D) + F0 * Eb / (B * D * H * (Eb - P))) - P * F0 * Eb / (B * D * H * (Eb - P)^2),
		"B" = -P / (B^2 * D) * (-0.5 + F0 * Eb / (H * (Eb - P)) - F0 * Eb * (2 * Eb - P) / (H * (Eb - P)^2)),
		"D" = -P / (B * D^2) * (-0.5 + F0 * Eb / (H * (Eb - P)) - F0 * Eb * (2 * Eb - P) / (H * (Eb - P)^2)),
		"H" = P * F0 * Eb * (Eb + P) / (B * D * H^2 * (Eb - P)^2),
		"F0" = -P * Eb / (B *D * H * (Eb - P)),
		"E" = P^2 * F0 * Eb / (B * D * H * E * (Eb - P)^2))
}


