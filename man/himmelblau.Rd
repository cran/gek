\encoding{UTF-8}

\name{himmelblau}
\alias{himmelblau}
\alias{himmelblauGrad}

\title{Himmelblau's Function}

\description{
Himmelblau's function is defined by
\deqn{f_{\rm himmelblau}(x_1, x_2) = (x_1^2 + x_2 - 11)^2 + (x_1 + x_2^2 -7)^2}
with \eqn{x_1, x_2 \in [-5, 5]}.
}

\usage{
himmelblau(x)
himmelblauGrad(x)
}

\arguments{
\item{x}{
a numeric \code{\link{vector}} of length \code{2} or a numeric \code{\link{matrix}} with \code{n} rows and \code{2} columns.
}
}

\details{
The gradient of Himmelblau's function is
\deqn{\nabla f_{\rm himmelblau}(x_1, x_2) = \begin{pmatrix} 4 x_1 (x_1^2 + x_2 - 11) + 2 (x_1 + x_2^2 - 7) \\ 2 (x_1^2 + x_2 - 11) + 4 x_2 (x_1 + x_2^2 - 7) \end{pmatrix}.}

Himmelblau's function has four global minima \eqn{f_{\rm himmelblau}(x^{\star}) = 0} at \eqn{x^{\star} = (3, 2)}, 
\eqn{x^{\star} = (-2.805118, 3.131312)}, \eqn{x^{\star} = (-3.779310, -3.283186)} and \eqn{x^{\star} = (3.584428, -1.848126)}.
}

\value{
\code{himmelblau} returns the function value of Himmelblau's function at \code{x}.

\code{himmelblauGrad} returns the gradient of Himmelblau's function at \code{x}.
}

\references{
Himmelblau, D. (1972). Applied Nonlinear Programming. McGraw-Hill. ISBN 0-07-028921-2.

Jamil, M. and Yang, X.-S. (2013). A Literature Survey of Benchmark Functions for Global Optimization Problems. \emph{International Journal of Mathematical Modelling and Numerical Optimisation}, \bold{4}(2):150-–194. \doi{10.1504/IJMMNO.2013.055204}.

Surjanovic, S. and Bingham, D. (2013). Virtual Library of Simulation Experiments: Test Functions and Datasets. \url{https://www.sfu.ca/~ssurjano/} (retrieved January 19, 2024). 
}

\author{
Carmen van Meegen
}

\examples{
# Contour plot of Himmelblau's function 
n.grid <- 50
x1 <- x2 <- seq(-5, 5, length.out = n.grid)
y <- outer(x1, x2, function(x1, x2) himmelblau(cbind(x1, x2)))
contour(x1, x2, y, xaxs = "i", yaxs = "i", nlevels = 25, xlab = "x1", ylab = "x2")

# Perspective plot of Himmelblau's function
col.pal <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow",
	"#FF7F00", "red", "#7F0000"))
colors <- col.pal(100)
y.facet.center <- (y[-1, -1] + y[-1, -n.grid] + y[-n.grid, -1] + y[-n.grid, -n.grid])/4
y.facet.range <- cut(y.facet.center, 100)
persp(x1, x2, y, phi = 30, theta = -315, expand = 0.75, ticktype = "detailed", 
	col = colors[y.facet.range])
}
