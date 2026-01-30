\name{schwefel}
\alias{schwefel}
\alias{schwefelGrad}

\title{Schwefel Function}

\description{
The Schwefel function is defined by
\deqn{f_{\rm schwefel}(x_1, ..., x_d) = 418.9829 d - \sum_{k = 1}^{d} x_k \sin\left(\sqrt{|x_k|}\right)}
with \eqn{x_k \in [-500, 500]} for \eqn{k = 1, ..., d}.
}

\usage{
schwefel(x)
schwefelGrad(x)
}

\arguments{
\item{x}{
a numeric vector of length \code{d} or a numeric matrix with \code{n} rows and \code{d} columns.
}
}

\details{
The gradient of the Schwefel function is
\deqn{\nabla f_{\rm schwefel}(x_1, \dots, x_d) = \begin{pmatrix} -\sin\left( \sqrt{|x_1|} \right) - \frac{x_1^2 \cos\left(\sqrt{|x_1|}\right)}{2 |x_1|^{\frac{3}{2}}} \\ \vdots \\ -\sin\left( \sqrt{|x_d|} \right) - \frac{x_d^2 \cos\left(\sqrt{|x_d|}\right)}{2  |x_d|^{\frac{3}{2}}} \end{pmatrix}.}

The Schwefel function has one global minimum \eqn{f_{\rm schwefel}(x^{\star}) = 0} at \eqn{x^{\star} = (420.968746, \dots, 420.968746)}.
}

\value{
\code{schwefel} returns the function value of the Schwefel function at \code{x}.

\code{schwefelGrad} returns the gradient of the Schwefel function at \code{x}.
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
# 1-dimensional Schwefel function with tangents
curve(schwefel(x), from = -500, to = 500, n = 500)
x <- seq(-450, 450, length = 5)
y <- schwefel(x)
dy <- schwefelGrad(x)
tangents(x, y, dy, length = 200, lwd = 2, col = "red")
points(x, y, pch = 16)

# Contour plot of Schwefel function 
n.grid <- 75
x1 <- x2 <- seq(-500, 500, length.out = n.grid)
y <- outer(x1, x2, function(x1, x2) schwefel(cbind(x1, x2)))
contour(x1, x2, y, xaxs = "i", yaxs = "i", nlevels = 25, xlab = "x1", ylab = "x2")

# Perspective plot of Schwefel function
col.pal <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow",
	"#FF7F00", "red", "#7F0000"))
colors <- col.pal(100)
y.facet.center <- (y[-1, -1] + y[-1, -n.grid] + y[-n.grid, -1] + y[-n.grid, -n.grid])/4
y.facet.range <- cut(y.facet.center, 100)
persp(x1, x2, y, phi = 30, theta = -315, expand = 0.75, ticktype = "detailed", 
	col = colors[y.facet.range])
}

