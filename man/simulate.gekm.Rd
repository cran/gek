\encoding{UTF-8}

\name{simulate.gekm}
\alias{simulate.gekm}

\title{
Simulates Conditional Process Paths
}
\description{
Simulates process paths conditional on a fitted \code{\link{gekm}} model.
}
\usage{
\method{simulate}{gekm}(object, nsim = 1, seed = NULL, newdata = NULL, tol = NULL, \dots)
}
\arguments{
  \item{object}{
an object of class \code{"gekm"}.
}
  \item{nsim}{
number of simulated process paths. Default is \code{1}.
}
  \item{seed}{
argument is not supported.
}  
  \item{newdata}{
a \code{\link{data.frame}} containing the points where the process in \code{object} should be evaluated.
}
  \item{tol}{
a tolerance for the conditional number of the conditional correlation matrix of \code{newdata}, see \code{\link{blockChol}} for details. Default is \code{NULL}, i.e. no regularization is applied.
}
  \item{\dots}{
further arguments, not used.
}
}
\value{
  \item{val}{a \code{\link{matrix}} with \code{nrow(newdata)} rows and \code{nsim} columns of simulated response values at the points of \code{newdata}.
  Each column represents one conditional simulated process path.}
}

\references{
Cressie, N. A. C. (1993). \emph{Statistics for Spartial Data}. John Wiley & Sons. \doi{10.1002/9781119115151}.

Ripley, B. D. (1981). \emph{Spatial Statistics}. John Wiley & Sons. \doi{10.1002/0471725218}. 
}

\author{
Carmen van Meegen
}


\seealso{
\code{\link{gekm}} for fitting a (gradient-enhanced) Kriging model.

\code{\link{predict.gekm}} for prediction at new data points based on a model of class \code{"gekm"}.
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

# Fit Kriging model
km.1d <- gekm(y ~ x, data = dat, covtype = "gaussian", theta = 1)

# Fit Gradient-Enhanced Kriging model
gekm.1d <- gekm(y ~ x, data = dat, deriv = deri, covtype = "gaussian", theta = 1)

# Generate new data for prediction and simulation
newdat <- data.frame(x = seq(-6, 6, length = 600))

# Prediction for both models
pred.km.1d <- predict(km.1d, newdat)
pred.gekm.1d  <- predict(gekm.1d, newdat)

# Simulate process path conditional on fitted models
set.seed(1)
n <- 50
sim.km.1d <- simulate(km.1d, nsim = n, newdata = newdat, tol = 35)
sim.gekm.1d <- simulate(gekm.1d, nsim = n, newdata = newdat, tol = 35)

par(mfrow = c(1, 2), oma = c(3.5, 3.5, 0, 0.2), mar = c(0, 0, 1.5, 0))
matplot(newdat$x, sim.km.1d, type = "l", lty = 1, col = 2:8, lwd = 1, 
	ylim = c(-1, 12), main = "Kriging")
lines(newdat$x, pred.km.1d$fit + qnorm(0.975) * pred.km.1d$sd, lwd = 1.5)
lines(newdat$x, pred.km.1d$fit - qnorm(0.975) * pred.km.1d$sd, lwd = 1.5)
points(x, y, pch = 16, cex = 1)

matplot(newdat$x, sim.gekm.1d, type = "l", lty = 1, col = 2:8, 
	lwd = 1, ylim = c(-1, 12), main = "GEK", yaxt = "n")
lines(newdat$x, pred.gekm.1d$fit + qnorm(0.975) * pred.gekm.1d$sd, lwd = 1.5)
lines(newdat$x, pred.gekm.1d$fit - qnorm(0.975) * pred.gekm.1d$sd, lwd = 1.5)
points(x, y, pch = 16, cex = 1)

mtext(side = 1, outer = TRUE, line = 2.5, "x")
mtext(side = 2, outer = TRUE, line = 2.5, "f(x)")
}
