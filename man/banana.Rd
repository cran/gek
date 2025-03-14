\encoding{UTF-8}

\name{banana}
\alias{banana}
\alias{bananaGrad}

\title{Rosenbrock's Banana Function}

\description{
Rosenbrock's banana function is defined by
\deqn{f_{\rm banana}(x_1, ..., x_d) = \sum_{k = 1}^{d - 1} (100 (x_{k+1} - x_k^2)^2 + (x_k - 1)^2)}
with \eqn{x_k \in [-5, 10]} for \eqn{k = 1, ..., d} and \eqn{d \geq 2}.
}

\usage{
banana(x)
bananaGrad(x)
}

\arguments{
\item{x}{
a numeric \code{\link{vector}} of length \code{d} or a numeric \code{\link{matrix}} with \code{n} rows and \code{d} columns, where \code{d} must be greater than \code{1}.
}
}

\details{
The gradient of Rosenbrock's banana function is
\deqn{\nabla f_{\rm banana}(x_1, ..., x_d) = \begin{pmatrix} -400 (x_2 - x_1)^2 x_1 + 2 (x_1 - 1) \\ 200 (x_2 - x_1)^2 - 400 x_2 (x_3 - x_2^2) + 2 (x_2 - 1) \\ \vdots \\ 200 (x_{d-1} - x_{d-2})^2 - 400 x_{d-1} (x_d - x_{d-1}^2) + 2 (x_{d-1} - 1) \\ 200 (x_d - x_{d -1}^2)\end{pmatrix}.}

Rosenbrock's banana function has one global minimum \eqn{f_{\rm banana}(x^{\star}) = 0} at \eqn{x^{\star} = (1,\dots, 1)}.
}

\value{
\code{banana} returns the function value of Rosenbrock's banana function at \code{x}.

\code{bananaGrad} returns the gradient of Rosenbrock's banana function at \code{x}.
}

\references{
Jamil, M. and Yang, X.-S. (2013). A Literature Survey of Benchmark Functions for Global Optimization Problems. \emph{International Journal of Mathematical Modelling and Numerical Optimisation}, \bold{4}(2):150-–194. \doi{10.1504/IJMMNO.2013.055204}.

Rosenbrock, H. H. (1960). An Automatic Method for Finding the Greatest or least Value of a Function. \emph{The Computer Journal}, \bold{3}(3):175--184. \doi{10.1093/comjnl/3.3.175}.

Surjanovic, S. and Bingham, D. (2013). Virtual Library of Simulation Experiments: Test Functions and Datasets. \url{https://www.sfu.ca/~ssurjano/} (retrieved January 19, 2024). 
}

\author{
Carmen van Meegen
}

\examples{
# Contour plot of Rosenbrock's banana function 
n.grid <- 50
x1 <- seq(-2, 2, length.out = n.grid)
x2 <- seq(-1, 3, length.out = n.grid)
y <- outer(x1, x2, function(x1, x2) banana(cbind(x1, x2)))
contour(x1, x2, y, xaxs = "i", yaxs = "i", nlevels = 25, xlab = "x1", ylab = "x2")

# Perspective plot of Rosenbrock's banana function
col.pal <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow",
	"#FF7F00", "red", "#7F0000"))
colors <- col.pal(100)
y.facet.center <- (y[-1, -1] + y[-1, -n.grid] + y[-n.grid, -1] + y[-n.grid, -n.grid])/4
y.facet.range <- cut(y.facet.center, 100)
persp(x1, x2, y, phi = 30, theta = -315, expand = 0.75, ticktype = "detailed", 
	col = colors[y.facet.range])
}

