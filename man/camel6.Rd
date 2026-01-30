\encoding{UTF-8}

\name{camel6}
\alias{camel6}
\alias{camel6Grad}

\title{Six-Hump Camel Function}

\description{
The six-hump camel function is defined by
\deqn{f_{\rm camel6}(x_1, x_2) = \left(4 - 2.1 x_1^2 + \frac{x_1^4}{3}\right) x_1^2 + x_1 x_2 + (-4 + 4 x_2^2) x_2^2}
with \eqn{x_1 \in [-3, 3]} and \eqn{x_2 \in [-2, 2]}.
}

\usage{
camel6(x)
camel6Grad(x)
}

\arguments{
\item{x}{
a numeric \code{\link{vector}} of length \code{2} or a numeric \code{\link{matrix}} with \code{n} rows and \code{2} columns.
}
}

\details{
The gradient of the six-hump camel function is
\deqn{\nabla f_{\rm camel6}(x_1, x_2) = \begin{pmatrix} 8 x_1 - 8.4 x_1^3 + 2 x_1^5 + x_2 \\ x_1 - 8 x_2 + 16 x_2^3 \end{pmatrix}.}

The six-hump camel function has two global minima \eqn{f_{\rm camel6}(x^{\star}) = -1.031628} at \eqn{x^{\star} = (0.0898, -0.7126)} and \eqn{x^{\star} = (-0.0898, 0.7126)}.
}

\value{
\code{camel6} returns the function value of the six-hump camel function function at \code{x}.

\code{camel6Grad} returns the gradient of the six-hump camel function function at \code{x}. 
}

\references{
Jamil, M. and Yang, X.-S. (2013). A Literature Survey of Benchmark Functions for Global Optimization Problems. \emph{International Journal of Mathematical Modelling and Numerical Optimisation}, \bold{4}(2):150-–194. \doi{10.1504/IJMMNO.2013.055204}.

Surjanovic, S. and Bingham, D. (2013). Virtual Library of Simulation Experiments: Test Functions and Datasets. \url{https://www.sfu.ca/~ssurjano/} (retrieved January 19, 2024). 
}

\author{
Carmen van Meegen
}

\seealso{
\code{\link{testfunctions}} for further test functions.
}

\examples{
# Contour plot of six-hump camel function 
n.grid <- 50
x1 <- seq(-2, 2, length.out = n.grid)
x2 <- seq(-1, 1, length.out = n.grid)
y <- outer(x1, x2, function(x1, x2) camel6(cbind(x1, x2)))
contour(x1, x2, y, xaxs = "i", yaxs = "i", nlevels = 25, xlab = "x1", ylab = "x2")

# Perspective plot of six-hump camel function
col.pal <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow",
	"#FF7F00", "red", "#7F0000"))
colors <- col.pal(100)
y.facet.center <- (y[-1, -1] + y[-1, -n.grid] + y[-n.grid, -1] + y[-n.grid, -n.grid])/4
y.facet.range <- cut(y.facet.center, 100)
persp(x1, x2, y, phi = 30, theta = -315, expand = 0.75, ticktype = "detailed", 
	col = colors[y.facet.range])
}

