###                                                             ###
###                		CONFINT-METHOD FOR gekm	     	        ###
###                                                             ###

## confint.gekm - confint methode for an object of class gekm
## 
## @param object: gekm[1]
##		object of class gekm
## @param parm: chr[m] or int[m]
##		a character or interger vector specifing the relevant 
##		regression coefficients
## @param level: num[1]
##		confidence level
## @param ...:
##		further arguments, not used
##
## @output:
##		confidence intervals for the regression coefficients


confint.gekm <- function(object, parm, level = 0.95, scale = FALSE, ...){
	cf <- coef(object)
	pnames <- names(cf)
	if(missing(parm)) parm <- pnames
	else if(is.numeric(parm)) 
		parm <- pnames[parm]
	p <- (1 - level) / 2
	a <- c(p, 1 - p)

	df <- nobs(object) - length(cf)
	alpha <- qt(a, df)

	ses <- sqrt(diag(vcov(object, scale = scale)))[parm]

	ci <- cf[parm] + ses %o% alpha
	colnames(ci) <- paste(format(100 * a, trim = TRUE, scientific = FALSE, digits = 3L), "%")

	ci
}

