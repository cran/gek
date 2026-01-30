\name{tangents}
\alias{tangents}

\title{Add Tangent Lines to a Plot}

\description{
Draw tangent lines to an existing plot.
}

\usage{
tangents(x, y, slope, length = 1, \dots)
}

\arguments{
\item{x, y}{
	coordinate vectors of points \code{x} and function values \code{y}.
}
\item{slope}{
	vector of slopes at the points \code{x}.
}
\item{length}{
	desired length of tangent lines, see \sQuote{Details}.
}
\item{\dots}{
	further \link[=par]{graphical parameters} to be passed to \code{\link{segments}}.
}
}

\details{
The length of the tangent lines is scaled according to the current aspect ratio of the existing plot.
}

\references{
Oakley, J. and O'Hagan, A. (2002). Bayesian Inference for the Uncertainty Distribution of Computer Model Outputs. \emph{Biometrika}, \bold{89}(4):769--784. \doi{10.1093/biomet/89.4.769}.
}

\author{
Carmen van Meegen
}

\seealso{
\code{\link{segments}} for drawing line segments between pairs of points.
}

\examples{
# Define test function and its gradient from Oakley and O'Hagan (2002)
f <- function(x) 5 + x + cos(x)
fGrad <- function(x) 1 - sin(x)

# Generate coordinates and calculate slopes
x <- seq(-5, 5, length = 5)
y <- f(x)
dy <- fGrad(x)

# Draw curve and tangent lines
curve(f(x), from = -6, to = 6)
tangents(x, y, dy, length = 2, lwd = 2, col = 2:6)
points(x, y, pch = 16)
}

\keyword{aplot}
