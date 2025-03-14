\name{rastrigin}
\alias{rastrigin}
\alias{rastriginGrad}

\title{Rastrigin Function}

\description{
Rastrigin function is defined by
\deqn{f_{\rm rastrigin}(x_1, ..., x_d) = 10 d + \sum_{k = 1}^{d} (x_k^2 - 10\cos(2\pi x_k))}
with \eqn{x_k \in [-5.12, 5.12]} for \eqn{k = 1, ..., d}.
}

\usage{
rastrigin(x)
rastriginGrad(x)
}

\arguments{
\item{x}{
a numeric \code{\link{vector}} or a numeric \code{\link{matrix}} with \code{n} rows and \code{d} columns.
If a \code{\link{vector}} is passed, the 1-dimensional version of the Rastrigin function is calculated.
}
}

\details{
The gradient of Rastrigin function is
\deqn{\nabla f_{\rm rastrigin}(x_1, ..., x_d) = \begin{pmatrix} 2x_1 + 20 \sin(2\pi x_1) \\ \vdots \\ 2x_d + 20 \sin(2\pi x_d))\end{pmatrix}.}

Rastrigin function has one global minimum \eqn{f_{\rm rastrigin}(x^{\star}) = 0} at \eqn{x^{\star} = (0,\dots, 0)}.
}

\value{
\code{rastrigin} returns the function value of Rastrigin function at \code{x}.

\code{rastriginGrad} returns the gradient of Rastrigin function at \code{x}.
}

\references{
Plevris, V. and Solorzano, G. (2022). A Collection of 30 Multidimensional Functions for Global Optimization Benchmarking. \emph{Data}, \bold{7}(4):46. \doi{10.3390/data7040046}. 

Surjanovic, S. and Bingham, D. (2013). Virtual Library of Simulation Experiments: Test Functions and Datasets. \url{https://www.sfu.ca/~ssurjano/} (retrieved January 19, 2024). 
}

\author{
Carmen van Meegen
}

\seealso{
\code{\link{tangents}} for drawing tangent lines.
}

\examples{
# 1-dimensional Rastrigin function with tangents
curve(rastrigin(x), from = -5, to = 5, n = 200)
x <- seq(-4.5, 4.5, length = 5)
y <- rastrigin(x)
dy <- rastriginGrad(x)
tangents(x, y, dy, length = 2, lwd = 2, col = "red")
points(x, y, pch = 16)

# Contour plot of Rastrigin function 
n.grid <- 100
x1 <- x2 <- seq(-5.12, 5.12, length.out = n.grid)
y <- outer(x1, x2, function(x1, x2) rastrigin(cbind(x1, x2)))
contour(x1, x2, y, xaxs = "i", yaxs = "i", nlevels = 25, xlab = "x1", ylab = "x2")

# Perspective plot of Rastrigin function
col.pal <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow",
	"#FF7F00", "red", "#7F0000"))
colors <- col.pal(100)
y.facet.center <- (y[-1, -1] + y[-1, -n.grid] + y[-n.grid, -1] + y[-n.grid, -n.grid])/4
y.facet.range <- cut(y.facet.center, 100)
persp(x1, x2, y, phi = 30, theta = -315, expand = 0.75, ticktype = "detailed", 
	col = colors[y.facet.range])
}

