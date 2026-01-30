\name{sphere}
\alias{sphere}
\alias{sphereGrad}

\title{Sphere Function}

\description{
The sphere function is defined by
\deqn{f_{\rm sphere}(x_1, ..., x_d) = \sum_{k = 1}^{d} x_k^2}
with \eqn{x_k \in [-5.12, 5.12]} for \eqn{k = 1, ..., d}.
}

\usage{
sphere(x)
sphereGrad(x)
}

\arguments{
\item{x}{
a numeric \code{\link{vector}} or a numeric \code{\link{matrix}} with \code{n} rows and \code{d} columns.
If a \code{\link{vector}} is passed, the 1-dimensional version of the sphere function is calculated.}
}

\details{
The gradient of the sphere function is
\deqn{\nabla f_{\rm sphere}(x_1, \dots, x_d) = \begin{pmatrix} 2 x_1 \\ \vdots \\ 2 x_d \end{pmatrix}.}

The sphere function has one global minimum \eqn{f_{\rm sphere}(x^{\star}) = 0} at \eqn{x^{\star} = (0, \dots, 0)}. 
}

\value{
\code{sphere} returns the function value of the sphere function at \code{x}.

\code{sphereGrad} returns the gradient of the sphere function at \code{x}.
}

\references{
Plevris, V. and Solorzano, G. (2022). A Collection of 30 Multidimensional Functions for Global Optimization Benchmarking. \emph{Data}, \bold{7}(4):46. \doi{10.3390/data7040046}. 

Surjanovic, S. and Bingham, D. (2013). Virtual Library of Simulation Experiments: Test Functions and Datasets. \url{https://www.sfu.ca/~ssurjano/} (retrieved January 19, 2024). 
}

\author{
Carmen van Meegen
}

\seealso{
\code{\link{testfunctions}} for further test functions.

\code{\link{tangents}} for drawing tangent lines.
}

\examples{
# 1-dimensional sphere function with tangents
curve(sphere(x), from = -5, to = 5)
x <- seq(-4.5, 4.5, length = 5)
y <- sphere(x)
dy <- sphereGrad(x)
tangents(x, y, dy, length = 2, lwd = 2, col = "red")
points(x, y, pch = 16)

# Contour plot of sphere function 
n.grid <- 15
x1 <- x2 <- seq(-5.12, 5.12, length.out = n.grid)
y <- outer(x1, x2, function(x1, x2) sphere(cbind(x1, x2)))
contour(x1, x2, y, xaxs = "i", yaxs = "i", nlevels = 25, xlab = "x1", ylab = "x2")

contour(x1, x2, y, #xaxs = "i", yaxs = "i",
nlevels = 25, xlab = "x1", ylab = "x2")

x <- expand.grid(x1, x2)
gradient <- sphereGrad(x)

vectorfield(x, gradient, col = 4, scale = 1.1)

# Perspective plot of sphere function
col.pal <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow",
	"#FF7F00", "red", "#7F0000"))
colors <- col.pal(100)
y.facet.center <- (y[-1, -1] + y[-1, -n.grid] + y[-n.grid, -1] + y[-n.grid, -n.grid])/4
y.facet.range <- cut(y.facet.center, 100)
persp(x1, x2, y, phi = 30, theta = -315, expand = 0.75, ticktype = "detailed", 
	col = colors[y.facet.range])
}

