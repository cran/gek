\name{griewank}
\alias{griewank}
\alias{griewankGrad}

\title{Griewank Function}

\description{
Griewank function is defined by
\deqn{f_{\rm griewank}(x_1, ..., x_d) = \sum_{k = 1}^{d} \frac{x_k^2}{4000} - \prod_{k = 1}^d \cos\left(\frac{x_k}{\sqrt{k}}\right) + 1}
with \eqn{x_k \in [-600, 600]} for \eqn{k = 1, ..., d}.
}

\usage{
griewank(x)
griewankGrad(x)
}

\arguments{
\item{x}{
a numeric \code{\link{vector}} or a numeric \code{\link{matrix}} with \code{n} rows and \code{d} columns.
If a \code{\link{vector}} is passed, the 1-dimensional version of the Griewank function is calculated.
}
}

\details{
The gradient of Griewank function is
\deqn{\nabla f_{\rm griewank}(x_1, ..., x_d) = \begin{pmatrix} \frac{x_1}{2000} + \frac{1}{\sqrt{1}} \sin\left( \frac{x_1}{\sqrt{1}}\right) \prod_{k = 2}^d \cos\left(\frac{x_k}{\sqrt{k}}\right) \\ \vdots \\ \frac{x_d}{2000} + \frac{1}{\sqrt{d}} \sin\left( \frac{x_d}{\sqrt{d}}\right) \prod_{k = 1}^{d-1} \cos\left(\frac{x_k}{\sqrt{k}}\right) \end{pmatrix}.}

Griewank function has one global minimum \eqn{f_{\rm griewank}(x^{\star}) = 0} at \eqn{x^{\star} = (0,\dots, 0)}.
}

\value{
\code{griewank} returns the function value of Griewank function at \code{x}.

\code{griewankGrad} returns the gradient of Griewank function at \code{x}.
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
# 1-dimensional Griewank function with tangents
curve(griewank(x), from = -5, to = 5, n = 200)
x <- seq(-5, 5, length = 5)
y <- griewank(x)
dy <- griewankGrad(x)
tangents(x, y, dy, length = 2, lwd = 2, col = "red")
points(x, y, pch = 16)

# Contour plot of Griewank function 
n.grid <- 50
x1 <- x2 <- seq(-5, 5, length.out = n.grid)
y <- outer(x1, x2, function(x1, x2) griewank(cbind(x1, x2)))
contour(x1, x2, y, xaxs = "i", yaxs = "i", nlevels = 25, xlab = "x1", ylab = "x2")

# Perspective plot of Griewank function
col.pal <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow",
	"#FF7F00", "red", "#7F0000"))
colors <- col.pal(100)
y.facet.center <- (y[-1, -1] + y[-1, -n.grid] + y[-n.grid, -1] + y[-n.grid, -n.grid])/4
y.facet.range <- cut(y.facet.center, 100)
persp(x1, x2, y, phi = 30, theta = -315, expand = 0.75, ticktype = "detailed", 
	col = colors[y.facet.range])
}

