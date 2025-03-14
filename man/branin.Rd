\name{branin}
\alias{branin}
\alias{braninGrad}

\title{Branin-Hoo Function}

\description{
The Branin-Hoo function is defined by
\deqn{f_{\rm branin}(x_1, x_2) = \left(x_2 - \frac{5.1}{4 \pi^2}x_1^2 + \frac{5}{\pi}x_1 - 6\right)^2 + 10 \left(1-\frac{1}{8\pi}\right)\cos(x_1) + 10}
with \eqn{x_1 \in [-5, 10]} and \eqn{x_2 \in [0, 15]}.
}

\usage{
branin(x)
braninGrad(x)
}

\arguments{
\item{x}{
a numeric \code{\link{vector}} of length \code{2} or a numeric \code{\link{matrix}} with \code{n} rows and \code{2} columns.
}
}

\details{
The gradient of the Branin-Hoo function is
\deqn{\nabla f_{\rm branin}(x_1, x_2) = \begin{pmatrix} 2 \left(x_2 - \frac{5.1 x_1^2}{4 \pi^2} + \frac{5 x_1}{\pi} - 6\right) \left(-10.2 \frac{x_1}{4\pi^2} + \frac{5}{\pi}\right) - 10 \left(1 - \frac{1}{8\pi}\right) \sin(x_1) \\ 2 \left( x_2 - \frac{5.1 x_1^2}{4 \pi^2} + \frac{5 x_1}{\pi} - 6\right)\end{pmatrix}.}

The Branin-Hoo function has three global minima \eqn{f_{\rm branin}(x^{\star}) = 0.397887} at \eqn{x^{\star} = (-\pi, 12.275)}, \eqn{x^{\star} = (\pi, 2.275)} and \eqn{x^{\star} = (9.42478, 2.475)}. 
}

\value{
\code{branin} returns the function value of the Branin-Hoo function at \code{x}.

\code{braninGrad} returns the gradient of the Branin-Hoo function at \code{x}.
}

\references{
Branin, Jr., F. H. (1972). Widely Convergent Method of Finding Multiple Solutions of Simultaneous Nonlinear Equations. \emph{IBM Journal of Research and Development}, \bold{16}(5):504--522.

Surjanovic, S. and Bingham, D. (2013). Virtual Library of Simulation Experiments: Test Functions and Datasets. \url{https://www.sfu.ca/~ssurjano/} (retrieved January 19, 2024). 
}

\author{
Carmen van Meegen
}



\examples{
# Contour plot of Branin-Hoo function 
n.grid <- 50
x1 <- seq(-5, 10, length.out = n.grid)
x2 <- seq(0, 15, length.out = n.grid)
y <- outer(x1, x2, function(x1, x2) branin(cbind(x1, x2)))
contour(x1, x2, y, xaxs = "i", yaxs = "i", nlevels = 25, xlab = "x1", ylab = "x2")

# Perspective plot of Branin-Hoo function
col.pal <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow",
	"#FF7F00", "red", "#7F0000"))
colors <- col.pal(100)
y.facet.center <- (y[-1, -1] + y[-1, -n.grid] + y[-n.grid, -1] + y[-n.grid, -n.grid])/4
y.facet.range <- cut(y.facet.center, 100)
persp(x1, x2, y, phi = 30, theta = -315, expand = 0.75, ticktype = "detailed", 
	col = colors[y.facet.range])
}

