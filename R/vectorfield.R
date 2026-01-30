###                                                             ###
###					  VECTORFIELD-FUNCTION	     	 	        ###
###                                                             ###

## vectorfield - drawing a vector filed to an existing plot
## 
## @param x: num[n]
##		numeric vector with x-coordinates
## @param gradient: num[n] 
##		numeric vector with y-coordinates
## @param scale: num[1]
##		a scaling factor for the length of the arrows to be drawn
## @param max.len: num[1]
##		scaling factor for the length of the tangents to be drawn
## @param min.len: num[1]
##		
## @param ...:
##		further arguments to be passed to arrows()
##
## @output:
##		invisible(NULL)


vectorfield <- function(x, gradient, scale = 1, max.len = 0.1, min.len = 0.001, ...){
	
	if(ncol(x) != 2L | ncol(gradient) != 2L)
		stop("'x' and 'gradient' must have 2 columns")
	if(nrow(x) != nrow(gradient))
		stop("'x' and 'gradient' must have the same number of rows")
	if(min.len > max.len)
		stop("'min.len' must be less than 'max.len'")
	n <- nrow(x)	
	
	args.arrows <- list(...)
	
	min.dist <- min(dist(x))
	grad.len <- sqrt(rowSums(gradient^2))
	off <- min.dist / max(grad.len)

	x0 <- x[ , 1]
	y0 <- x[ , 2]
	x1 <- x0 + scale * off * gradient[ , 1]
	y1 <- y0 + scale * off * gradient[ , 2]


	if("length" %in% names(args.arrows)){
		if(!length(args.arrows[["length"]]) %in% c(1, n))
			stop("invalid length of 'length'")
		if(length(args.arrows[["length"]]) == 1) 
			len <- rep(args.arrows[["length"]], n)
		else
			len <- args.arrows[["length"]]
		args.arrows[["length"]] <- NULL
	}else{
		a <- min(grad.len)
		b <- max(grad.len)
		c <- min.len
		d <- max.len
		len <- (d - c) / (b - a) * grad.len + 
			(b * c - a * d) / (b - a)
	}

		
	arg.arrows <- list(angle = 15, code = 2)
	arg.arrows[names(args.arrows)] <- args.arrows	

	for(i in seq(along = x0)){
		do.call("arrows", c(arg.arrows, x0 = x0[i], y0 = y0[i], 
			x1 = x1[i], y1 = y1[i], length = scale * len[i]))
	}
	
	invisible(NULL)
	
}
