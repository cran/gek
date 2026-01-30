###                                                             ###
###                	 FORMULA-METHOD FOR gekm 	    	        ###
###                                                             ###

## formula.gekm - formula method for gekm object
## 
## @param x: gekm[1]
##		object of class gekm
## @param ...:
##		further arguments, not used	
##
## @output:
##		extracted formula


formula.gekm <- function(x, ...){
	formula(terms(x))
}
