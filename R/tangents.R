###                                                             ###
###						TANGENTS-FUNCTION	     	 	        ###
###                                                             ###

### Function for drawing tangent lines to an existing plot 

tangents <- function(x, y, slope, length = 1, ...){

	w <- par("pin")[1] / diff(par("usr")[1:2])
	h <- par("pin")[2] / diff(par("usr")[3:4])
	asp <- w / h
	
	intercept <- y - x * slope
	eps <- 0.5 * length / sqrt(1 + (slope / asp)^2) 
		
	segments(x - eps, slope * (x - eps) + intercept, 
		x + eps, slope * (x + eps) + intercept, ...)
		
	invisible(NULL)
	
}
