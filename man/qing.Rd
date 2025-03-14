\encoding{UTF-8}

\name{qing}
\alias{qing}
\alias{qingGrad}

\title{Qing Function}

\description{
Qing function is defined by
\deqn{f_{\rm qing}(x_1, ..., x_d) = \sum_{k = 1}^{d} (x_i^2 - i)^2}
with \eqn{x_k \in [-500, 500]} for \eqn{k = 1, ..., d}.
}

\usage{
qing(x)
qingGrad(x)
}

\arguments{
\item{x}{
a numeric \code{\link{vector}} or a numeric \code{\link{matrix}} with \code{n} rows and \code{d} columns.
If a \code{\link{vector}} is passed, the 1-dimensional version of the Rastrigin function is calculated.
}
}

\details{
The gradient of Qing function is
\deqn{\nabla f_{\rm qing}(x_1, ..., x_d) = \begin{pmatrix} 4 x_1 (x_1^2 - 1)  \\  \vdots \\ 4 x_d (x_d^2 - d) \end{pmatrix}.}

Qing function has \eqn{2^d} global minimum \eqn{f_{\rm qing}(x^{\star}) = 0} at \eqn{x^{\star} = (\pm \sqrt{1},\dots, \pm \sqrt{d})}.
}

\value{
\code{qing} returns the function value of Qing function at \code{x}.

\code{qingGrad} returns the gradient of Qing function at \code{x}.
}

\references{
Qing, A. (2006). Dynamic Differential Evolution Strategy and Applications in Electromagnetic Inverse Scattering Problems. \emph{IEEE Transactions on Geoscience and Remote Sensing}, \bold{44}(1):116-–125. \doi{10.1109/TGRS.2005.859347}.

Jamil, M. and Yang, X.-S. (2013). A Literature Survey of Benchmark Functions for Global Optimization Problems. \emph{International Journal of Mathematical Modelling and Numerical Optimisation}, \bold{4}(2):150-–194. \doi{10.1504/IJMMNO.2013.055204}.
}

\author{
Carmen van Meegen
}

\examples{
# 1-dimensional Qing function with tangents
curve(qing(x), from = -1.7, to = 1.7)
x <- seq(-1.5, 1.5, length = 5)
y <- qing(x)
dy <- qingGrad(x)
tangents(x, y, dy, length = 1, lwd = 2, col = "red")
points(x, y, pch = 16)

# Contour plot of Qing function 
n.grid <- 50
x1 <- seq(-2, 2, length.out = n.grid)
x2 <- seq(-2, 2, length.out = n.grid)
y <- outer(x1, x2, function(x1, x2) qing(cbind(x1, x2)))
contour(x1, x2, y, xaxs = "i", yaxs = "i", nlevels = 25, xlab = "x1", ylab = "x2")

# Perspective plot of Qing function
col.pal <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow",
	"#FF7F00", "red", "#7F0000"))
colors <- col.pal(100)
y.facet.center <- (y[-1, -1] + y[-1, -n.grid] + y[-n.grid, -1] + y[-n.grid, -n.grid])/4
y.facet.range <- cut(y.facet.center, 100)
persp(x1, x2, y, phi = 30, theta = -315, expand = 0.75, ticktype = "detailed", 
	col = colors[y.facet.range])

}

