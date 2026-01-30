###                                                             ###
###                	  PRINT-METHOD FOR gekm 	    	        ###
###                                                             ###

## print.gekm - print method for gekm
## 
## @param x: gekm[1]
##		an object of class gekm
## @param digits: num[1]
##		number of digits to be printed
## @param scale: logi[1]
##		scale estimated process standard deviation
## @param ...:
##		further arguments, not used		
##
## @output:
##		invisible(x)

print.gekm <- function(x, digits = 4L, scale = FALSE, ...){

    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    cat("Coefficients:\n")
        print(format(coef(x), digits = digits), print.gap = 2L, 
            quote = FALSE)
    cat("\nStandard Deviation:", format(sigma(x, scale = scale), digits = digits))
	cat("\n\nCorrelation Stucture:", x$covtype)	
	cat("\nCorrelation Parameters:\n")
		print(format(x$theta, digits = digits), print.gap = 2L,
			quote = FALSE)
	cat("\n")
    invisible(x)

}
