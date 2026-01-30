###                                                             ###
###                 PRINT-METHOD FOR summary.gekm    	        ###
###                                                             ###

## print.summary.gekm - print method for summary.gekm
## 
## @param x: summary.gekm[1]
##		object of class summary.gekm
## @param digits: num[1]
##		number of digits to be printed
## @param ...:
##		further arguments passed to printCoefmat
##
## @output:
##		invisible(x)



print.summary.gekm <- function(x, digits = 4L, ...){

    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
	cat("Coefficients:\n")
	printCoefmat(coef(x), digits = digits, ...)
	cat("\nStandard Deviation:", format(x$sigma, digits = digits))
	cat("\n\nCorrelation Stucture:", x$covtype)	
	cat("\nCorrelation Parameters:\n")
		print(format(x$theta, digits = digits), print.gap = 2L,
			quote = FALSE)
	cat("\n")
	
	invisible(x)

}
