\name{cigar}
\alias{cigar}
\alias{cigarGrad}

\title{Bent Cigar Function}

\description{
The Bent cigar function is defined by
\deqn{f_{\rm cigar}(x_1, ..., x_d) = x_1^2 + 10^6 \sum_{k = 2}^{d} x_k^2}
with \eqn{x_k \in [-100, 100]} for \eqn{k = 1, ..., d}.
}

\usage{
cigar(x)
cigarGrad(x)
}

\arguments{
\item{x}{
a numeric vector of length \code{2} or a numeric matrix with \code{n} rows and \code{2} columns.
}
}

\details{
The gradient of the bent cigar function is
\deqn{\nabla f_{\rm cigar}(x_1, ..., x_d) = \begin{pmatrix} 2x_1 \\ 20^6x_2\\ \vdots \\  20^6x_d\end{pmatrix}.}

The bent cigar function has one global minimum \eqn{f_{\rm cigar}(x^{\star}) = 0} at \eqn{x^{\star} = (1,\dots, 1)}.
}

\value{
\code{cigar} returns the function value of the bent cigar function at \code{x}.

\code{cigarGrad} returns the gradient of the bent cigar function at \code{x}.
}

\references{
Plevris, V. and Solorzano, G. (2022). A Collection of 30 Multidimensional Functions for Global Optimization Benchmarking. \emph{Data}, \bold{7}(4):46. \doi{10.3390/data7040046}. 
}

\author{
Carmen van Meegen
}

\seealso{
\code{\link{tangents}} for drawing tangent lines.
}

\examples{
# 1-dimensional Cigar function with tangents
curve(cigar(x), from = -5, to = 5, n = 200)
x <- seq(-4.5, 4.5, length = 5)
y <- cigar(x)
dy <- cigarGrad(x)
tangents(x, y, dy, length = 2, lwd = 2, col = "red")
points(x, y, pch = 16)

# Contour plot of Cigar function 
n.grid <- 50
x1 <- x2 <- seq(-100, 100, length.out = n.grid)
y <- outer(x1, x2, function(x1, x2) cigar(cbind(x1, x2)))
contour(x1, x2, y, xaxs = "i", yaxs = "i", nlevels = 25, xlab = "x1", ylab = "x2")

# Perspective plot of Cigar function
col.pal <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow",
	"#FF7F00", "red", "#7F0000"))
colors <- col.pal(100)
y.facet.center <- (y[-1, -1] + y[-1, -n.grid] + y[-n.grid, -1] + y[-n.grid, -n.grid])/4
y.facet.range <- cut(y.facet.center, 100)
persp(x1, x2, y, phi = 30, theta = -315, expand = 0.75, ticktype = "detailed", 
	col = colors[y.facet.range])
}

