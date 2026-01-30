
###                                                             ###
###                	   PLOT-METHOD FOR gekm	     		        ###
###                                                             ###

## plot.gekm - simulates Gaussian process path given an object
##		of class gekm
## 
## @param x: gekm[1]
##		object of class gekm
## @param y: 
##		not used
## @param main: chr[1]
##		main title
## @param ylim: num[2]
##		range of y-axis
## @param panel.first: expr[1]
##		expressions to add to plot
## @param add: logi[1]
##		add to existing plot
## @param reestim: logi[1]
##		re-estimate regression coefficients
## @param sd.fit: logi[1]
##		estimated standard deviation, i.e. root mean squared error of prediction
## @param scale: logi[1]
##		scale estimated process variance
## @param df: num[1]
##		degrees of freedom
## @param add.interval: logi[1]
##		add confidence interval
## @param level: num[1]
##		confidence level
## @param args.arrows: list[1]
##		arguments to be passed to arrows
## @param ...:
##		further arguments to be passed to plot.default() or points()
##		
## @output:
##		invisible(res)


plot.gekm <- function(x, y = NULL, main = "Leave-One-Out",
	ylim = NULL, panel.first = abline(0, 1),
	add = FALSE, reestim = TRUE, scale = FALSE, df = NULL,
	add.interval = FALSE, level = 0.95, args.arrows = NULL, ...){

	response <- model.response(model.frame(x))

	res <- loo(x, reestim = reestim, sd.fit = FALSE, scale = scale, df = df,
		interval = if(add.interval) "confidence" else "none", level = level)
	prediction <- if(add.interval) res[ , "fit"] else res
	
	if(add.interval){
		
		lower <- res[ , "lower"]
		upper <- res[ , "upper"]
	
		arg.arrows <- list(x0 = response, y0 = lower, y1 = upper,
			length = 0.05, angle = 90, code = 3)
		arg.arrows[names(args.arrows)] <- args.arrows	
	
	}
		
	if(add){
	
		if(add.interval) do.call("arrows", arg.arrows)		
		points(response, prediction, ...)
	
	}else{
	
		if(add.interval) ylim <- if(is.null(ylim)) range(lower, upper) else ylim
		
		plot.default(x = response, y = prediction, main = main, ylim = ylim,
			panel.first = eval(c(panel.first, {
				if(add.interval)
					do.call("arrows", arg.arrows)
			})), ...)
			
	}
	
	invisible(res)
}
