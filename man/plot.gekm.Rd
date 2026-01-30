\encoding{UTF-8}

\name{plot.gekm}
\alias{plot.gekm}

\title{
Plot Method for a gekm Object
}
\description{
Visualization of the leave-one-out cross-validation results of a \code{gekm} model.
}

\usage{
\method{plot}{gekm}(x, y, main = "Leave-One-Out", ylim = NULL, panel.first = abline(0, 1),
	add = FALSE, reestim = TRUE, scale = FALSE, df = NULL, add.interval = FALSE, 
	level = 0.95, args.arrows = NULL, \dots)
}

\arguments{
  \item{x}{
an object of class \code{"gekm"}.
}
  \item{y}{
not used.
}
  \item{main}{
main title for the plot.
}
  \item{ylim}{
limits for the y-axis.
}
  \item{panel.first}{
an \code{\link{expression}} to be evaluated before the actual plot. Default is \code{abline(0, 1)}.
}
  \item{add}{
\code{\link{logical}}. Should results be added to an already existing plot? Default is \code{FALSE}.
}
  \item{reestim}{
\code{\link{logical}}. Should the regression coefficients be re-estimated? Default is \code{TRUE}, see \code{\link{loo.gekm}} for details.
}
  \item{scale}{
\code{\link{logical}}. Should the estimated process variance be scaled? Default is \code{FALSE}, see \code{\link{sigma.gekm}} for details.
}
  \item{df}{
degrees of freedom for the \eqn{t} distribution. Default is \code{NULL}, see \code{\link{predict.gekm}} for details.
}
  \item{add.interval}{
\code{\link{logical}}. Should confidence intervals be added? Default is \code{FALSE}.
}
  \item{level}{
confidence level for calculating confidence intervals. Default is \code{0.95}.
}
  \item{args.arrows}{
a \code{\link{list}} with further arguments to be passed to \code{\link{arrows}}. Only used if \code{add.interval = TRUE}.
}
  \item{\dots}{
further arguments to be passed to \code{\link{plot.default}} or \code{\link{points}}.
}
}


\details{
For further details on the arguments \code{scale} see
}

\value{
Returns the predicted values of the leave-one-out cross-validation invisibly.
If \code{add.interval = TRUE}, the lower and upper bounds of the confidence intervals are also returned.
}


\references{
Bachoc, F. (2013). Cross Validation and Maximum Likelihood Estimations of Hyper-parameters of Gaussian Processes with Model Misspecification. \emph{Computational Statistics and Data Analysis}, \bold{66}:55--69. \doi{10.1016/j.csda.2013.03.016}.

Dubrule, O. (1983). Cross Validation of Kriging in a Unique Neighborhood. \emph{Mathematical Geology}, \bold{15}:687--699. \doi{10.1007/BF01033232}.

Martin, J. D. and Simpson, T. W. (2005). Use of Kriging Models to Approximate Deterministic Computer Models. \emph{AIAA Journal}, \bold{43}(4):853--863. \doi{10.2514/1.8650}.
}

\author{
Carmen van Meegen
}


\seealso{
\code{\link{gekm}} for fitting a (gradient-enhanced) Kriging model.

\code{\link{loo}} for leave-one-out cross-validation.

\code{\link{branin}} for the Branin-Hoo function.

\code{\link{arrows}} for drawing arrows.
}


\examples{
## 2-dimensional example: Branin-Hoo function

# Generate a grid for training
n <- 4
x1 <- seq(-5, 10, length = n)
x2 <- seq(0, 15, length = n)
x <- expand.grid(x1 = x1, x2 = x2)
y <- branin(x)
dy <- braninGrad(x)
dat <- data.frame(x, y)
deri <- data.frame(dy)

# Fit (gradient-enhanced) Kriging model
km.2d <- gekm(y ~ .^2, data = dat)
gekm.2d <- gekm(y ~ .^2, data = dat, deriv = deri)

# Plot leave-one-out cross-validation results
plot(km.2d)
plot(km.2d, panel.first = grid())
plot(km.2d, panel.first = {grid(); abline(0, 1, col = 8)})
plot(km.2d, add.interval = TRUE)
plot(km.2d, add.interval = TRUE, pch = 16, col = 4)
plot(km.2d, add.interval = TRUE, pch = 16, col = 4,
	panel.first = {grid(); abline(0, 1)}, 
	args.arrows = list(col = 4, length = 0))

plot(km.2d, pch = 1, col = 4, cex = 1.2, lwd = 2)
plot(gekm.2d, pch = 4, col = 2, cex = 1.2, lwd = 2, add = TRUE)
legend("topleft", legend = c("Kriging", "GEK"), col = c(4, 2), pch = c(1, 4), pt.lwd = 2)

par(mfrow = c(1, 2), oma = c(3.6, 3.5, 1.5, 0.2), mar = c(0, 0, 1.5, 0))
res <- plot(km.2d, col = 7, pch = 16, add.interval = TRUE, main = "Kriging",
	scale = TRUE, panel.first = {grid(); abline(0, 1, col = 8)})
res
plot(gekm.2d, col = 3, pch = 16, add.interval = TRUE, scale = TRUE, main = "GEK",
	ylim = range(res), yaxt = "n", panel.first = {grid(); abline(0, 1, col = 8)})
title(main = "Leave-One-Out", outer = TRUE)
mtext(side = 1, outer = TRUE, line = 2.5, "response")
mtext(side = 2, outer = TRUE, line = 2.5, "prediction")
}


\keyword{aplot}
