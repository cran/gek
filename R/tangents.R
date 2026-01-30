###                                                             ###
###						TANGENTS-FUNCTION	     	 	        ###
###                                                             ###

## tangents - drawing tangent lines to an existing plot
## 
## @param x: num[n]
##		numeric vector with x-coordinates
## @param y: num[n] 
##		numeric vector with y-coordinates
## @param slope: num[n]
##		numeric vector with slopes at x
## @param length: num[1]
##		scaling factor for the length of the tangents to be drawn
## @param ...:
##		further arguments to be passed to segments()
##
## @output:
##		invisible(NULL)

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
