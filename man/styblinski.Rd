\encoding{UTF-8}

\name{styblinski}
\alias{styblinski}
\alias{styblinskiGrad}

\title{Styblinski-Tang Function}

\description{
Styblinski-Tang function is defined by
\deqn{f_{\rm styblinski}(x_1, ..., x_d) = \frac{1}{2} \sum_{k = 1}^{d} (x_k^4 - 16x_k^2 + 5x_k)}
with \eqn{x_k \in [-5, 5]} for \eqn{k = 1, ..., d}.
}

\usage{
styblinski(x)
styblinskiGrad(x)
}

\arguments{
\item{x}{
a numeric \code{\link{vector}} or a numeric \code{\link{matrix}} with \code{n} rows and \code{d} columns.
If a \code{\link{vector}} is passed, the 1-dimensional version of the Rastrigin function is calculated.
}
}

\details{
The gradient of Styblinski-Tang function is
\deqn{\nabla f_{\rm styblinski}(x_1, ..., x_d) = \begin{pmatrix} 2x_1^3 - 16 x_1 + 2.5\\ \vdots \\ 2 x_d^3 - 16 x_d + 2.5 \end{pmatrix}.}

Styblinski-Tang function has one global minimum \eqn{f_{\rm styblinski}(x^{\star}) = -39.16599d} at
\eqn{x^{\star} = (-2.903534,\dots, -2.903534)}.
}

\value{
\code{styblinski} returns the function value of Styblinski-Tang function at \code{x}.

\code{styblinskiGrad} returns the gradient of Styblinski-Tang function at \code{x}.
}

\references{
Jamil, M. and Yang, X.-S. (2013). A Literature Survey of Benchmark Functions for Global Optimization Problems. \emph{International Journal of Mathematical Modelling and Numerical Optimisation}, \bold{4}(2):150-–194. \doi{10.1504/IJMMNO.2013.055204}.

Surjanovic, S. and Bingham, D. (2013). Virtual Library of Simulation Experiments: Test Functions and Datasets. \url{https://www.sfu.ca/~ssurjano/} (retrieved January 19, 2024). 
}

\author{
Carmen van Meegen
}

\seealso{
\code{\link{tangents}} for drawing tangent lines.
}

\examples{
# 1-dimensional Styblinski-Tang function with tangents
curve(styblinski(x), from = -5, to = 5)
x <- seq(-4.5, 4.5, length = 5)
y <- styblinski(x)
dy <- styblinskiGrad(x)
tangents(x, y, dy, length = 2, lwd = 2, col = "red")
points(x, y, pch = 16)

# Contour plot of Styblinski-Tang function 
n.grid <- 50
x1 <- seq(-5, 5, length.out = n.grid)
x2 <- seq(-5, 5, length.out = n.grid)
y <- outer(x1, x2, function(x1, x2) styblinski(cbind(x1, x2)))
contour(x1, x2, y, xaxs = "i", yaxs = "i", nlevels = 25, xlab = "x1", ylab = "x2")

# Perspective plot of Styblinski-Tang function
col.pal <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow",
	"#FF7F00", "red", "#7F0000"))
colors <- col.pal(100)
y.facet.center <- (y[-1, -1] + y[-1, -n.grid] + y[-n.grid, -1] + y[-n.grid, -n.grid])/4
y.facet.range <- cut(y.facet.center, 100)
persp(x1, x2, y, phi = 30, theta = -315, expand = 0.75, ticktype = "detailed", 
	col = colors[y.facet.range])
}

