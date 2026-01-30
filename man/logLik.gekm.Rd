\encoding{UTF-8}

\name{logLik.gekm}
\alias{logLik.gekm}

\title{
Log-Likelihood of a gekm Object
}
\description{
Returns the log-likelihood of a \code{\link{gekm}} object.
}
\usage{
\method{logLik}{gekm}(object, \dots)
}
\arguments{
  \item{object}{
an object of class \code{gekm}.
}
  \item{\dots}{
not used.
}
}
\value{
The log-likelihood value of the model evaluated at the estimated coefficients. 
}
\references{
Oakley, J. and O'Hagan, A. (2002). Bayesian Inference for the Uncertainty Distribution of Computer Model Outputs. \emph{Biometrika}, \bold{89}(4):769--784. \doi{10.1093/biomet/89.4.769}.

Park, J.-S. and Beak, J. (2001). Efficient Computation of Maximum Likelihood Estimators in a Spatial Linear Model with Power Exponential Covariogram. \emph{Computers & Geosciences}, \bold{27}(1):1--7. \doi{10.1016/S0098-3004(00)00016-9}.

Rasmussen, C. E. and Williams, C. K. I. (2006). \emph{Gaussian Processes for Machine Learning}. The MIT Press. \url{https://gaussianprocess.org/gpml/}.

Santner, T. J., Williams, B. J., and Notz, W. I. (2018). \emph{The Design and Analysis of Computer Experiments}. 2nd edition. Springer-Verlag.

Zimmermann, R. (2015). On the Condition Number Anomaly of Gaussian Correlation Matrices. \emph{Linear Algebra and its Applications}, \bold{466}:512--526. \doi{10.1016/j.laa.2014.10.038}.
}

\author{
Carmen van Meegen
}

\seealso{
\code{\link{gekm}} for fitting a (gradient-enhanced) Kriging model.
}

\examples{
## 1-dimensional example

# Define test function and its gradient from Oakley and O’Hagan (2002)
f <- function(x) 5 + x + cos(x)
fGrad <- function(x) 1 - sin(x)

# Generate coordinates and calculate slopes
x <- seq(-5, 5, length = 5)
y <- f(x)
dy <- fGrad(x)
dat <- data.frame(x, y)
deri <- data.frame(x = dy)

# Fit (gradient-enhanced) Kriging model
km.1d <- gekm(y ~ x, data = dat, covtype = "gaussian", theta = 1)
gekm.1d <- gekm(y ~ x, data = dat, deriv = deri, covtype = "gaussian", theta = 1)

# Extract log-likelihood value 
logLik(km.1d)
logLik(gekm.1d)
}
