\encoding{UTF-8}

\name{camel3}
\alias{camel3}
\alias{camel3Grad}

\title{Three-Hump Camel Function}

\description{
The three-hump camel function is defined by
\deqn{f_{\rm camel3}(x_1, x_2) = 2 x_1^2 - 1.05 x_1^4 + \frac{x_1^6}{6} + x_1 x_2 + x_2^2}
with \eqn{x_1, x_2 \in [-5, 5]}.
}

\usage{
camel3(x)
camel3Grad(x)
}

\arguments{
\item{x}{
a numeric \code{\link{vector}} of length \code{2} or a numeric \code{\link{matrix}} with \code{n} rows and \code{2} columns.
}
}

\details{
The gradient of the three-hump camel function is
\deqn{\nabla f_{\rm camel3}(x_1, x_2) = \begin{pmatrix} 4 x_1 - 4.2 x_1^3 + x_1^5 + x_2 \\ x_1 + 2 x_2 \end{pmatrix}.}

The three-hump camel function has one global minimum \eqn{f_{\rm camel3}(x^{\star}) = 0} at \eqn{x^{\star} = (0, 0)}.
}

\value{
\code{camel3} returns the function value of the three-hump camel function at \code{x}.

\code{camel3Grad} returns the gradient of the three-hump camel function at \code{x}.
}

\references{
Jamil, M. and Yang, X.-S. (2013). A Literature Survey of Benchmark Functions for Global Optimization Problems. \emph{International Journal of Mathematical Modelling and Numerical Optimisation}, \bold{4}(2):150-–194. \doi{10.1504/IJMMNO.2013.055204}.

Surjanovic, S. and Bingham, D. (2013). Virtual Library of Simulation Experiments: Test Functions and Datasets. \url{https://www.sfu.ca/~ssurjano/} (retrieved January 19, 2024). 
}

\author{
Carmen van Meegen
}

\examples{
# Contour plot of three-hump camel function 
n.grid <- 50
x1 <- x2 <- seq(-2, 2, length.out = n.grid)
y <- outer(x1, x2, function(x1, x2) camel3(cbind(x1, x2)))
contour(x1, x2, y, xaxs = "i", yaxs = "i", nlevels = 25, xlab = "x1", ylab = "x2")

# Perspective plot of three-hump camel function
col.pal <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow",
	"#FF7F00", "red", "#7F0000"))
colors <- col.pal(100)
y.facet.center <- (y[-1, -1] + y[-1, -n.grid] + y[-n.grid, -1] + y[-n.grid, -n.grid])/4
y.facet.range <- cut(y.facet.center, 100)
persp(x1, x2, y, phi = 30, theta = -315, expand = 0.75, ticktype = "detailed", 
	col = colors[y.facet.range])
}
