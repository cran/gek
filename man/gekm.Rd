\encoding{UTF-8}

\name{gekm}
\alias{gekm}
\alias{print.gekm}

\title{
Fitting (Gradient-Enhanced) Kriging Models
}

\description{
Estimation of a Kriging model with or without derivatives.
}

\usage{
gekm(formula, data, deriv, covtype = c("matern5_2", "matern3_2", "gaussian"),
	theta = NULL, tol = NULL, optimizer = c("NMKB", "L-BFGS-B"), 
	lower = NULL, upper = NULL, start = NULL, ncalls = 20, control = NULL,
	model = TRUE, x = FALSE, y = FALSE, dx = FALSE, dy = FALSE, \dots)
	
\method{print}{gekm}(x, digits = 4L, \dots)
}

\arguments{
  \item{formula}{
a \code{formula} that defines the regression functions.
Note that only formula expressions for which the derivations are contained in the derivatives table of \code{\link{deriv}} are supported in the gradient-enhanced Kriging model.
In addition, formulas containing \code{\link{I}} also work for the gradient-enhanced Kriging model, although this function is not included in the derivatives table of \code{\link{deriv}}.
See \code{\link{derivModelMatrix}} for some examples of the trend specification.
}
  \item{data}{
a \code{\link{data.frame}} with named columns of \code{n} training points of dimension \code{d}. Note, all variables contained in \code{data} are used for the construction of the correlation matrix without and with derivatives.
}
  \item{deriv}{
an optional \code{\link{data.frame}} with the derivatives, whose columns should be named like those of \code{data}. If not specified, a Kriging model without derivatives is estimated.
}
  \item{covtype}{
a \code{\link{character}} to specify the covariance structure to be used. One of \code{matern5_2}, \code{matern3_2} or \code{gaussian}. Default is \code{matern5_2}.
}
  \item{theta}{
a \code{\link{numeric}} vector of length \code{d} for the hyperparameters (optional). If not given, hyperparameters will be estimated via maximum likelihood.
}
  \item{tol}{
a tolerance for the conditional number of the correlation matrix, see \code{\link{blockChol}} for details. Default is \code{NULL}, i.e. no regularization is applied.
}
  \item{optimizer}{
an optional \code{\link{character}} that characterizes the optimization algorithms to be used for maximum likelihood estimation. See \sQuote{Details}.
}
  \item{lower}{
an optional lower bound for the optimization of the correlation parameters.
}
  \item{upper}{
an optional upper bound for the optimization of the correlation parameters.
}
  \item{start}{
an optional \code{\link{vector}} of inital values for the optimization of the correlation parameters.
}
  \item{ncalls}{
an optional \code{\link{integer}} that defines the number of randomly selected initial values for the optimization.
}
  \item{control}{
a \code{\link{list}} of control parameters for the optimization routine. See \code{\link{optim}} or \code{\link[dfoptim]{nmkb}}.
}
  \item{model}{
\code{\link{logical}}. Should the model frame be returned? Default is \code{TRUE}.
}
  \item{x}{
\code{\link{logical}}. Should the model matrix be returned? Default is \code{FALSE}.
}
  \item{y}{
\code{\link{logical}}. Should the response vector be returned? Default is \code{FALSE}.
}
  \item{dx}{
\code{\link{logical}}. Should the derivative of the model matrix be returned? Default is \code{FALSE}.
}
  \item{dy}{
\code{\link{logical}}. Should the derivatives of the response be returned? Default is \code{FALSE}.
}
  \item{digits}{
number of digits to be used for the \code{print} method.
}
	\item{\dots}{
further arguments, currently not used.
}}

\details{

Parameter estimation is performed via maximum likelihood.
The \code{optimizer} argument can be used to select one of the optimization algorithms \code{"NMBK"} or \code{"L-BFGS-B"}.
In the case of the \code{"L-BFGS-B"}, analytical gradients of the \dQuote{concentrated} log likelihood are used.
For one-dimensional problems, \code{\link{optimize}} is called and the algorithm selected via \code{optimizer} is ignored.
}

\value{
\code{gekm} returns an object of \code{\link{class}} \code{"gekm"} whose underlying structure is a list containing at least the following components:

\item{coefficients}{
the estimated regression coefficients.
}
  \item{sigma}{
the estimated process standard deviation.
}
  \item{theta}{
the (estimated) correlation parameters.
}
  \item{covtype}{
the name of the correlation function.
}
  \item{chol}{
(the components of) the upper triangular matrix of the Cholesky decomposition of the correlation matrix.
}
  \item{optimizer}{
the optimization algorithm.
}
  \item{convergence}{
the convergence code.
}
  \item{message}{
information from the optimizer.
}
  \item{logLik}{
the value of the negative \dQuote{concentrated} log likelihood at the estimated parameters.
}
  \item{derivatives}{
\code{TRUE} if a gradient-enhanced Kriging model was adapted, otherwise \code{FALSE}.
}
  \item{data}{
the \code{\link{data.frame}} that was specified via the \code{data} argument.
}
  \item{deriv}{
if \code{derivatives = TRUE}, the \code{\link{data.frame}} with the derivatives that was specified via the \code{deriv} argument.
}
  \item{call}{
the matched call.
}
  \item{terms}{
the \code{\link{terms}} object used.
}
  \item{model}{
if requested (the default), the model frame used.
}
  \item{x}{
if requested, the model matrix.
}
\item{y}{
if requested, the response vector.
}
  \item{dx}{
if requested, the derivatives of the model matrix.
}
  \item{dx}{
if requested, the vector of derivatives of the response.
}
}

\references{
Cressie, N. A. C. (1993). \emph{Statistics for Spartial Data}. John Wiley & Sons. \doi{10.1002/9781119115151}.

Koehler, J. and Owen, A. (1996). Computer Experiments. In Ghosh, S. and Rao, C. (eds.), \emph{Design and Analysis of Experiments}, volume 13 of \emph{Handbook of Statistics}, pp. 261–308. Elsevier Science. \doi{10.1016/S0169-7161(96)13011-X}.

Krige, D. G. (1951). A Statistical Approach to Some Basic Mine Valuation Problems on the Witwatersrand. \emph{Journal of the Southern African Institute of Mining and Metallurgy}, \bold{52}(6):199--139. 

Laurent, L., Le Riche, R., Soulier, B., and Boucard, PA. (2019). An Overview of Gradient-Enhanced Metamodels with Applications. \emph{Archives of Computational Methods in Engineering}, \bold{26}(1):61--106. \doi{10.1007/s11831-017-9226-3}.

Martin, J. D. and Simpson, T. W. (2005). Use of Kriging Models to Approximate Deterministic Computer Models. \emph{AIAA Journal}, \bold{43}(4):853--863. \doi{10.2514/1.8650}.

Morris, M., Mitchell, T., and Ylvisaker, D. (1993). Bayesian Design and Analysis of Computer Experiments: Use of Derivatives in Surface Prediction. \emph{Technometrics}, \bold{35}(3):243--255. \doi{10.1080/00401706.1993.10485320}.

Oakley, J. and O'Hagan, A. (2002). Bayesian Inference for the Uncertainty Distribution of Computer Model Outputs. \emph{Biometrika}, \bold{89}(4):769--784. \doi{10.1093/biomet/89.4.769}.

O'Hagan, A., Kennedy, M. C., and Oakley, J. E. (1999). Uncertainty Analysis and Other Inference Tools for Complex Computer Codes. In \emph{Bayesian Statistics 6}, Ed. J. M. Bernardo, J. O. Berger, A. P. Dawid and A .F. M. Smith, 503--524. Oxford University Press. 

O'Hagan, A. (2006). Bayesian Analysis of Computer Code Outputs: A Tutorial. \emph{Reliability Engineering & System Safet}, \bold{91}(10):1290--1300. \doi{10.1016/j.ress.2005.11.025}.

Park, J.-S. and Beak, J. (2001). Efficient Computation of Maximum Likelihood Estimators in a Spatial Linear Model with Power Exponential Covariogram. \emph{Computers & Geosciences}, \bold{27}(1):1--7. \doi{10.1016/S0098-3004(00)00016-9}.

Ranjan, P., Haynes, R. and Karsten, R. (2011). A Computationally Stable Approach to Gaussian Process Interpolation of Deterministic Computer Simulation Data. \emph{Technometrics}, \bold{53}:366--378. \doi{10.1198/TECH.2011.09141}.

Rasmussen, C. E. and Williams, C. K. I. (2006). \emph{Gaussian Processes for Machine Learning}. The MIT Press. \url{https://gaussianprocess.org/gpml/}.

Ripley, B. D. (1981). \emph{Spatial Statistics}. John Wiley & Sons. \doi{10.1002/0471725218}. 

Sacks, J., Welch, W. J., Mitchell, T. J., and Wynn, H. P. (1989). Design and Analysis of Computer Experiments. \emph{Statistical Science}, \bold{4}(4):409--423. \doi{10.1214/ss/1177012413}.

Santner, T. J., Williams, B. J., and Notz, W. I. (2018). \emph{The Design and Analysis of Computer Experiments}. 2nd edition. Springer-Verlag.

Stein, M. L. (1999). \emph{Interpolation of Spatial Data: Some Theory for Kriging}. Springer Series in Statistics. Springer-Verlag. \doi{10.1007/978-1-4612-1494-6}.

}


\author{
Carmen van Meegen
}

\seealso{
\code{\link{predict.gekm}} for prediction at new data points based on a model of class \code{"gekm"}.

\code{\link{simulate.gekm}} for simulation of process paths conditional on a model of class \code{"gekm"}.
}

\examples{
## 1-dimensional example: Oakley and O’Hagan (2002)

# Define test function and its gradient
f <- function(x) 5 + x + cos(x)
fGrad <- function(x) 1 - sin(x)

# Generate coordinates and calculate slopes
x <- seq(-5, 5, length = 5)
y <- f(x)
dy <- fGrad(x)
dat <- data.frame(x, y)
deri <- data.frame(x = dy)

# Fit Kriging model
km.1d <- gekm(y ~ x, data = dat, covtype = "gaussian", theta = 1)
km.1d

# Fit Gradient-Enhanced Kriging model
gekm.1d <- gekm(y ~ x, data = dat, deriv = deri, covtype = "gaussian", theta = 1)
gekm.1d


## 2-dimensional example: Morris et al. (1993)

# List of inputs with their distributions and their respective ranges
inputs <- list("r_w" = list(dist = "norm", mean =  0.1, sd = 0.0161812, min = 0.05, max = 0.15),
	"r" = list(dist = "lnorm", meanlog = 7.71, sdlog = 1.0056, min = 100, max = 50000),
	"T_u" = list(dist = "unif", min = 63070, max = 115600),
	"H_u" = list(dist = "unif", min = 990, max = 1110),
	"T_l" = list(dist = "unif", min = 63.1, max = 116),
	"H_l" = list(dist = "unif", min = 700, max = 820),
	"L" = list(dist = "unif", min = 1120, max = 1680),
	# for a more nonlinear, nonadditive function, see Morris et al. (1993)
	"K_w" = list(dist = "unif", min = 1500, max = 15000))

# Generate design
design <- data.frame("r_w" = c(0, 0.268, 1),
	"r" = rep(0, 3),
	"T_u" = rep(0, 3),
	"H_u" = rep(0, 3),
	"T_l" = rep(0, 3),
	"H_l" = rep(0, 3),
	"L" = rep(0, 3),
	"K_w" = c(0, 1, 0.268))

# Function to transform design onto input range
transform <- function(x, data){
	for(p in names(data)){
		data[ , p] <- (x[[p]]$max - x[[p]]$min) * data[ , p] + x[[p]]$min
	}
	data
}

# Function to transform derivatives
deriv.transform <- function(x, data){
	for(p in colnames(data)){
		data[ , p] <- data[ , p] * (x[[p]]$max - x[[p]]$min)  
	}
	data
}

# Generate outcome and derivatives
design.trans <- transform(inputs, design)
design$y <- borehole(design.trans)
deri.trans <- boreholeGrad(design.trans)
deri <- data.frame(deriv.transform(inputs, deri.trans))

# Design and data
cbind(design[ , c("r_w", "K_w", "y")], deri[ , c("r_w", "K_w")])

# Fit gradient-enhanced Kriging model with Gaussian correlation function
mod <- gekm(y ~ 1, data = design[ , c("r_w", "K_w", "y")], 
	deriv = deri[ , c("r_w", "K_w")], covtype = "gaussian")
mod

## Compare results with Morris et al. (1993):

# Estimated hyperparameters
# in Morris et al. (1993): 0.429 and 0.467
1 / (2 * mod$theta^2)
# Estimated intercept
# in Morris et al. (1993): 69.15
coef(mod)
# Estimated standard deviation
# in Morris et al. (1993): 135.47
mod$sigma
# Posterior mean and standard deviation at (0.5, 0.5)
# in Morris et al. (1993): 69.4 and 2.7
predict(mod, data.frame("r_w" = 0.5, "K_w" = 0.5))
# Posterior mean and standard deviation at (1, 1)
# in Morris et al. (1993): 230.0 and 19.2
predict(mod, data.frame("r_w" = 1, "K_w" = 1))

## Graphical comparison: 

# Generate a 21 x 21 grid for prediction
n_grid <- 21
x <- seq(0, 1, length.out = n_grid)
grid <- expand.grid("r_w" = x, "K_w" = x)
pred <- predict(mod, grid)

# Compute ground truth on (transformed) grid
newdata <- data.frame("r_w" = grid[ , "r_w"], 
	"r" = 0, "T_u" = 0, "H_u" = 0, 
	"T_l" = 0, "H_l" = 0, "L" = 0,
	"K_w" = grid[ , "K_w"])
newdata <- transform(inputs, newdata)
truth <- borehole(newdata)

# Contour plots of predicted and actual output
par(mfrow = c(1, 2), oma = c(3.5, 3.5, 0, 0.2), mar = c(0, 0, 1.5, 0))
contour(x, x, matrix(pred$fit, nrow = n_grid, ncol = n_grid, byrow = TRUE),
	levels = c(seq(10, 50, 10), seq(100, 250, 50)),
	main = "Predicted output")
points(design[ , c("K_w", "r_w")], pch = 16)
contour(x, x, matrix(truth, nrow = n_grid, ncol = n_grid, byrow = TRUE), 
	levels = c(seq(10, 50, 10), seq(100, 250, 50)),
	yaxt = "n", main = "Ground truth")
points(design[ , c("K_w", "r_w")], pch = 16)
mtext(side = 1, outer = TRUE, line = 2.5, "Normalized hydraulic conductivity of borehole")
mtext(side = 2, outer = TRUE, line = 2.5, "Normalized radius of borehole")
}

\keyword{models}
