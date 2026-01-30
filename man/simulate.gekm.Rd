\encoding{UTF-8}

\name{simulate.gekm}
\alias{simulate.gekm}

\title{
Simulation of Conditional Process Paths
}
\description{
Simulates process paths conditional on a fitted \code{\link{gekm}} object.
}
\usage{
\method{simulate}{gekm}(object, nsim = 1, seed = NULL, newdata = NULL, 
	scale = FALSE, df = NULL, tol = NULL, \dots)
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
a \code{\link{data.frame}} containing the points at which the process path should be realized. The column names must be identical to those in the \code{data} used to construct the \code{"gekm"} object.
}
  \item{scale}{
\code{\link{logical}}. Should the estimated process standard deviation be scaled? Default is \code{FALSE}, see \code{\link{sigma.gekm}} for details.
}
  \item{df}{
degrees of freedom of the \eqn{t} distribution. Default is \code{NULL}, see \code{\link{predict.gekm}} for details.
}
  \item{tol}{
a tolerance for the conditional number of the conditional correlation matrix of \code{newdata}, see \code{\link{blockChol}} for details. Default is \code{NULL}, i.e. no regularization is applied.
}
  \item{\dots}{
further arguments, not used.
}
}

\details{
By setting \code{df = Inf}, paths of a Gaussian process are simulated.
}

\value{
  A \code{\link{matrix}} with \code{nrow(newdata)} rows and \code{nsim} columns of simulated response values at the points of \code{newdata}.
  Each column represents one conditional simulated process path.
}

\references{
Cressie, N. A. C. (1993). \emph{Statistics for Spartial Data}. John Wiley & Sons. \doi{10.1002/9781119115151}.

Oakley, J. and O'Hagan, A. (2002). Bayesian Inference for the Uncertainty Distribution of Computer Model Outputs. \emph{Biometrika}, \bold{89}(4):769--784. \doi{10.1093/biomet/89.4.769}.

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
df <- NULL
scale <- FALSE
pred.km.1d <- predict(km.1d, newdat, sd.fit = FALSE, interval = "confidence", 
 df = df, scale = scale)
pred.gekm.1d  <- predict(gekm.1d, newdat, sd.fit = FALSE, interval = "confidence", 
 df = df, scale = scale)

# Simulate process paths conditional on fitted models
set.seed(1)
n <- 500
sim.km.1d <- simulate(km.1d, nsim = n, newdata = newdat, tol = 35, df = df, scale = scale)
sim.gekm.1d <- simulate(gekm.1d, nsim = n, newdata = newdat, tol = 35, df = df, scale = scale)

par(mfrow = c(1, 2), oma = c(3.5, 3.5, 0, 0.2), mar = c(0, 0, 1.5, 0))
matplot(newdat$x, sim.km.1d, type = "l", lty = 1, col = 2:8, lwd = 1, 
	ylim = c(-1, 12), main = "Kriging")
matplot(newdat$x, pred.km.1d, type = "l", lwd = 2, add = TRUE,
	col = "black", lty = 1)
points(x, y, pch = 16, cex = 1, col = "red")

matplot(newdat$x, sim.gekm.1d, type = "l", lty = 1, col = 2:8, 
	lwd = 1, ylim = c(-1, 12), main = "GEK", yaxt = "n")
matplot(newdat$x, pred.gekm.1d, type = "l", lwd = 2, add = TRUE,
	col = "black", lty = 1)
points(x, y, pch = 16, cex = 1, col = "red")

mtext(side = 1, outer = TRUE, line = 2.5, "x")
mtext(side = 2, outer = TRUE, line = 2.5, "f(x)")

# Compare predicted means and standard deviations from predict() and simulate()
pred.km.1d <- predict(km.1d, newdat, sd.fit = TRUE, df = df, scale = scale)
pred.gekm.1d  <- predict(gekm.1d, newdat, sd.fit = TRUE, df = df, scale = scale)

# Predicted means
plot(newdat$x, pred.km.1d$fit, type = "l", lty = 1, lwd = 1, 
	ylim = c(-1, 12), main = "Kriging")
lines(newdat$x, rowMeans(sim.km.1d), col = 4)
points(x, y, pch = 16, cex = 1, col = "red")

plot(newdat$x, pred.gekm.1d$fit, type = "l", lty = 1, lwd = 1, 
	ylim = c(-1, 12), main = "GEK", yaxt = "n")
lines(newdat$x, rowMeans(sim.gekm.1d), col = 4)
points(x, y, pch = 16, cex = 1, col = "red")

mtext(side = 1, outer = TRUE, line = 2.5, "x")
mtext(side = 2, outer = TRUE, line = 2.5, "f(x)")

# Standard deviation
plot(newdat$x, pred.km.1d$sd.fit, type = "l", lty = 1, lwd = 1, 
	ylim = c(0, 0.8), main = "Kriging")
lines(newdat$x, apply(sim.km.1d, 1, sd), col = 4)
points(x, rep(0, 5), pch = 16, cex = 1, col = "red")

plot(newdat$x, pred.gekm.1d$sd.fit, type = "l", lty = 1, lwd = 1, 
	ylim = c(0, 0.8), main = "GEK", yaxt = "n")
lines(newdat$x, apply(sim.gekm.1d, 1, sd), col = 4)
points(x, rep(0, 5), pch = 16, cex = 1, col = "red")

mtext(side = 1, outer = TRUE, line = 2.5, "x")
mtext(side = 2, outer = TRUE, line = 2.5, "standard deviation")


}
