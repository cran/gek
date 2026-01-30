\encoding{UTF-8}

\name{logLikGrad}
\alias{logLikGrad}
\alias{logLikGrad.default}
\alias{logLikGrad.gekm}

\title{
Log-Likelihood Function
}

\description{
Calculates the gradient of the negative \dQuote{concentrated} log-likelihood.
}

\usage{
logLikGrad(param, object, \dots)

\method{logLikGrad}{default}(param, object, x, y, dx = NULL, dy = NULL, 
	covtype = c("matern5_2", "matern3_2", "gaussian"),
	tolerance = NULL, envir, \dots)
\method{logLikGrad}{gekm}(param, object, \dots)
}

\arguments{
 \item{param}{
a numeric \code{\link{vector}} corresponding to the values of the correlation parameters at which the negative \dQuote{concentrated} log-likelihood is to be evaluated.
}
   \item{object}{
a numeric \code{\link{matrix}} containing the  or an object of class \code{gekm}.
}
  \item{x}{
a \code{\link{model.matrix}}.
}  
   \item{y}{
a \code{\link{vector}} of response values.
}
   \item{dx}{
the derivatives of the model matrix \code{x}, see \code{\link{derivModelMatrix}} for details. Default is \code{NULL}.
}
   \item{dy}{
the \code{\link{vector}} of derivatives of the response \code{y}. Default is \code{NULL}.
}
   \item{covtype}{
a \code{\link{character}} to specify the correlation structure to be used. One of \code{matern5_2}, \code{matern3_2} or \code{gaussian}. Default is \code{matern5_2}, see \code{\link{blockCor}} for details.
}
   \item{tolerance}{
a tolerance for the conditional number of the correlation matrix, see \code{\link{blockChol}} for details. Default is \code{NULL}, i.e. no regularization is applied.
}
   \item{envir}{
an \code{\link{environment}} with auxiliary variables obtained from a corresponding \code{\link{call}} to \code{\link{logLikFun}}.
}
  \item{\dots}{
arguments to be passed to the default method.
}
}

\value{
The value of the negative \dQuote{concentrated} log-likelihood at \code{param}. 
}

\references{
Park, J.-S. and Beak, J. (2001). Efficient Computation of Maximum Likelihood Estimators in a Spatial Linear Model with Power Exponential Covariogram. \emph{Computers & Geosciences}, \bold{27}(1):1--7. \doi{10.1016/S0098-3004(00)00016-9}.

Rasmussen, C. E. and Williams, C. K. I. (2006). \emph{Gaussian Processes for Machine Learning}. The MIT Press. \url{https://gaussianprocess.org/gpml/}.

Santner, T. J., Williams, B. J., and Notz, W. I. (2018). \emph{The Design and Analysis of Computer Experiments}. 2nd edition. Springer-Verlag.

Zimmermann, R. (2015). On the Condition Number Anomaly of Gaussian Correlation Matrices. \emph{Linear Algebra and its Applications}, \bold{466}:512--526. \doi{10.1016/j.laa.2014.10.038}.
}

\author{
Carmen van Meegen
}

\seealso{
\code{\link{logLikFun}} for computing the value of the negative \dQuote{concentrated} log-likelihood.

\code{\link{gekm}} for fitting a (gradient-enhanced) Kriging model.
}

\examples{
## 2-dimensional example

# Generate coordinates and calculate slopes
x1 <- seq(-1.75, 1.75, length = 3)
x2 <- seq(-0.75, 0.75, length = 3)
X <- expand.grid(x1 = x1, x2 = x2)
y <- camel6(X)
dy <- camel6Grad(X)
dat <- data.frame(X, y)
deri <- data.frame(dy)

# Fit (gradient-enhanced) Kriging model
km.2d <- gekm(y ~ 1, data = dat, covtype = "gaussian", optimizer = "L-BFGS-B")
gekm.2d <- gekm(y ~ 1, data = dat, deriv = deri, covtype = "gaussian", optimizer = "L-BFGS-B")

# Compute negative 'concentrated' log-likelihood values
n.grid <- 20
theta1.grid <- seq(0.5, 4, length = n.grid)
theta2.grid <- seq(0.5, 2, length = n.grid)
params <- expand.grid(theta1 = theta1.grid, theta2 = theta2.grid)

logLik.km.2d <- apply(params, 1, logLikFun, km.2d)
logLik.gekm.2d <- apply(params, 1, logLikFun, gekm.2d)

logLikGrad.km.2d <- t(apply(params, 1, logLikGrad, km.2d))
logLikGrad.gekm.2d <- t(apply(params, 1, logLikGrad, gekm.2d))

# Plot negative 'concentrated' log-likelihood
par(mfrow = c(1, 2), oma = c(3.6, 3.5, 1.5, 0.2), mar = c(0, 0, 1.5, 0))
contour(theta1.grid, theta2.grid, matrix(logLik.km.2d, nrow = n.grid, ncol = n.grid), 
		nlevels = 50, main = "Kriging")
vectorfield(params, logLikGrad.km.2d, col = 4, lwd = 2, length = 0.1)
points(km.2d$theta[1], km.2d$theta[2], col = "red", pch = 16)

contour(theta1.grid, theta2.grid, matrix(logLik.gekm.2d, nrow = n.grid, ncol = n.grid), 
		nlevels = 50, main = "GEK", yaxt = "n")
points(gekm.2d$theta[1], gekm.2d$theta[2], col = "red", pch = 16)
vectorfield(params, logLikGrad.gekm.2d, col = 4, lwd = 2, length = 0.1)

title(main = "Negative 'concentrated' log-likelihood", outer = TRUE)
mtext(side = 1, outer = TRUE, line = 2.5, expression(theta[1]))
mtext(side = 2, outer = TRUE, line = 2.5, expression(theta[2]))
}
