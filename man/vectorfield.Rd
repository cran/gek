\name{vectorfield}
\alias{vectorfield}

\title{Add a Vector Field to a Plot}

\description{
Draw a vector field to an existing plot.
}

\usage{
vectorfield(x, gradient, scale = 1, max.len = 0.1, min.len = 0.001, \dots)
}

\arguments{
\item{x}{
	a \code{\link{matrix}} with \code{2} columns containing the coordinates.
}
\item{gradient}{
	a \code{\link{matrix}} with \code{2} columns columns containing the corresponding gradients at the locations in \code{x}.
}
\item{scale}{
	a scaling factor for the arrows to be drawn.
}
\item{max.len}{
	the maximum length of the edges of the arrow head.
}
\item{min.len}{
	the minimum length of the edges of the arrow head.
}
\item{\dots}{
	further \link[=par]{graphical parameters} to be passed to \code{\link{arrows}}.
}
}



\author{
Carmen van Meegen
}

\seealso{
\code{\link{arrows}} for drawing arrows between pairs of points.
}

\examples{

x <- seq(-2, 2, 0.2)
n <- length(x)
X <- expand.grid(X = x, Y= x)

f <- function(x,y) x * exp(-x^2-y^2)
df <- deriv(~ x * exp(-x^2-y^2), c("x", "y"), function(x, y){})
res <- df(X[, 1], X[,2])
grad <- attr(res, "gradient")

contour(x, x, matrix(res, n, n), asp = 1)
vectorfield(X, grad)
contour(x, x, matrix(res, n, n), asp = 1)
vectorfield(X, grad, col = 4, code = 1, scale = 2)
contour(x, x, matrix(res, n, n), asp = 1)
vectorfield(X, grad, col = 4, lwd = 2, scale = 2)
contour(x, x, matrix(res, n, n), asp = 1)
vectorfield(X, grad, col = 1, angle = 20, scale = 2)
contour(x, x, matrix(res, n, n), asp = 1)
vectorfield(X, grad, col = 4, angle = 20, length = 0.1)
contour(x, x, matrix(res, n, n), asp = 1)
vectorfield(X, grad, col = 4, angle = 20, length = 0.1, scale = 2)
contour(x, x, matrix(res, n, n), asp = 1)
vectorfield(X, grad, col = 3, min.len = 0.1, max.len = 0.15, scale = 2)

}

\keyword{aplot}
