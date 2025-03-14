\encoding{UTF-8}

\name{blockCor}
\alias{blockCor}
\alias{blockCor.default}
\alias{blockCor.gekm}

\title{
Correlation Matrix without or with Derivatives
}
\description{
Calculation of a correlation matrix with or without derivatives according to the specification of a correlation structure. 
}
\usage{
blockCor(x, \dots)

\method{blockCor}{default}(x, theta, covtype = c("matern5_2", "matern3_2", "gaussian"), 
	derivatives = FALSE, \dots)
\method{blockCor}{gekm}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a numeric \code{\link{matrix}} or an object of class \code{\link{gekm}}.
}
  \item{theta}{
\code{\link{numeric}} vector of length \code{d} for the hyperparameters.
}
  \item{covtype}{
\code{\link{character}} specifying the correlation function to be used. Must be one of \code{"matern5_2"}, \code{"matern3_2"} or \code{"gaussian"}.
See \sQuote{Details}.}
  \item{derivatives}{
\code{\link{logical}}, if \code{TRUE} the first and second partial derivatives of the correlation matrix are calculated, otherwise not.
}
  \item{\dots}{
  further arguments passed to or from other methods.
}
}
\details{
The correlation matrix with derivatives is defined as a block matrix of the form
\deqn{\begin{pmatrix} K & R^{\rm T} \\ R & S \end{pmatrix}.}

Three correlation functions are currently implemented in \code{blockCor}:
\itemize{
	\item{Matérn 5/2 kernel with \code{covtype = "matern5_2"}: \deqn{\phi(x, x'; \theta) = \prod_{k = 1}^d \left(1 + \frac{\sqrt{5}|x_k - x_k'|}{\theta_k} + \frac{5 |x_k - x_k'|^2}{3 \theta_k^2}\right) \exp\left( -\frac{\sqrt{5}|x_k - x_k'|}{\theta_k}\right)}}
	\item{Matérn 3/2 kernel with \code{covtype = "matern3_2"}: \deqn{\phi(x, x'; \theta) = \prod_{k = 1}^d \left(1 + \frac{\sqrt{3} |x_k - x_k'|}{\theta_k} \right) \exp\left( -\frac{\sqrt{3} |x_k - x_k'|}{\theta_k}\right)}}
	\item{Gaussian kernel with \code{covtype = "gaussian"}: \deqn{\phi(x, x'; \theta) = \prod_{k = 1}^d \exp\left( -\frac{(x_k - x_k')^2}{2\theta_k^2}\right)}}
}

}
\value{
\code{blockCor} returns a list with the following components:
\item{K}{The correlation matrix without derivatives.}
\item{R}{If \code{derivatives = TRUE}, the correlation matrix with first partial derivatives, otherwise \code{NULL}.}
\item{S}{If \code{derivatives = TRUE}, the correlation matrix with second partial derivatives, otherwise \code{NULL}.}

The components of the list can be combined in the following form to constructed the complete correlation matrix with derivatives:
\code{cbind(rbind(K, R), rbind(t(R), S))}.
}
\references{
Koehler, J. and Owen, A. (1996). Computer Experiments. In Ghosh, S. and Rao, C. (eds.), \emph{Design and Analysis of Experiments}, volume 13 of \emph{Handbook of Statistics}, pp. 261–308. Elsevier Science. \doi{10.1016/S0169-7161(96)13011-X}.

Rasmussen, C. E. and Williams, C. K. I. (2006). \emph{Gaussian Processes for Machine Learning}. The MIT Press. \url{https://gaussianprocess.org/gpml/}.

Stein, M. L. (1999). \emph{Interpolation of Spatial Data: Some Theory for Kriging}. Springer Series in Statistics. Springer-Verlag. \doi{10.1007/978-1-4612-1494-6}.
}
\author{
Carmen van Meegen
}

\seealso{
\code{\link{blockChol}} for block Cholesky decomposition.

\code{\link{tangents}} for drawing tangent lines.
}


\examples{
# Some examples of correlation matrices:

x <- matrix(1:4, ncol = 1)
blockCor(x, theta = 1)
blockCor(x, theta = 1, covtype = "gaussian")
blockCor(x, theta = 1, covtype = "gaussian", derivatives = TRUE)
blockCor(x, theta = 1, covtype = "matern3_2", derivatives = TRUE)

res <- blockCor(x, theta = 2, covtype = "gaussian", derivatives = TRUE)
cbind(rbind(res$K, res$R), rbind(t(res$R), res$S))


# Illustration of correlation functions and their derivatives:

x <- seq(0, 1, length = 100)
X <- matrix(x, ncol = 1)
gaussian <- blockCor(X, theta = 0.25, covtype = "gaussian", derivatives = TRUE)
matern5_2 <- blockCor(X, theta = 0.25, covtype = "matern5_2", derivatives = TRUE)
matern3_2 <- blockCor(X, theta = 0.25, covtype = "matern3_2", derivatives = TRUE)

# Correlation functions and first partial derivatives:
index <- c(10, 20, 40, 80)
par(mar = c(5.1, 5.1, 4.1, 2.1))

# Matern 3/2
plot(x, matern3_2$K[1, ], type = "l", xlab = expression(group("|", x - x*minute, "|")),
	ylab = expression(phi(x, x*minute, theta == 0.25)), lwd = 2)
tangents(x[index], matern3_2$K[index, 1], matern3_2$R[index, 1], 
	length = 0.15, lwd = 2, col = 2)
points(x[index], matern3_2$K[index, 1], pch = 16)

# Matern 5/2
lines(x, matern5_2$K[1, ], lwd = 2, col = 3)
tangents(x[index], matern5_2$K[index, 1], matern5_2$R[index, 1], 
	length = 0.15, lwd = 2, col = 2)
points(x[index], matern5_2$K[index, 1], pch = 16)

# Gaussian
lines(x, gaussian$K[1, ], lwd = 2, col = 4)
tangents(x[index], gaussian$K[index, 1], gaussian$R[index, 1], 
	length = 0.15, lwd = 2, col = 2)
points(x[index], gaussian$K[index, 1], pch = 16)

legend("topright", lty = 1, lwd = 2, col = c(1, 3, 4), bty = "n",
	legend = c("Matern 3/2", "Matern 5/2", "Gaussian"))


# First and second partial derivatives of correlation functions:
index <- c(5, 10, 20, 40, 80)
par(mar = c(5.1, 6.1, 4.1, 2.1))

# Gaussian
plot(x, matern3_2$R[1, ], type = "l", xlab = expression(group("|", x - x*minute, "|")),
	ylab = expression(frac(partialdiff * phi(x, x*minute, theta == 0.25), 
		partialdiff * x * minute)), lwd = 2)
tangents(x[index], matern3_2$R[1, index], matern3_2$S[index, 1], 
	length = 0.15, lwd = 2, col = 2)
points(x[index], matern3_2$R[1, index], pch = 16)

# Matern 5/2
lines(x, matern5_2$R[1, ], lwd = 2, col = 3)
tangents(x[index], matern5_2$R[1, index], matern5_2$S[index, 1], 
	length = 0.15, lwd = 2, col = 2)
points(x[index], matern5_2$R[1, index], pch = 16)

# Matern 3/2
lines(x, gaussian$R[1, ], lwd = 2, col = 4)
tangents(x[index], gaussian$R[1, index], gaussian$S[index, 1], 
	length = 0.15, lwd = 2, col = 2)
points(x[index], gaussian$R[1, index], pch = 16)

legend("topright", lty = 1, lwd = 2, col = c(1, 3, 4), bty = "n",
	legend = c("Matern 3/2", "Matern 5/2", "Gaussian"))
}

