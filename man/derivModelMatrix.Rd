\name{derivModelMatrix}
\alias{derivModelMatrix}
\alias{derivModelMatrix.default}
\alias{derivModelMatrix.gekm}
\alias{model.matrix.gekm}


\title{
Derivatives of Model Matrix
}
\description{
Determine the derivatives of a model matrix.
}
\usage{
derivModelMatrix(object, \dots)

\method{derivModelMatrix}{default}(object, data, \dots)
\method{derivModelMatrix}{gekm}(object, \dots)
\method{model.matrix}{gekm}(object, \dots)
}

\arguments{
  \item{object}{
an object of an appropriate \code{\link{class}}. For the default method, a \code{\link{formula}} defining the regression functions.
}
  \item{data}{
a \code{\link{data.frame}} with named columns.
}
  \item{\dots}{
further arguments, yet not used.
}
}
\details{
\code{derivModelMatrix} makes use of the function \code{\link{deriv}}. 
Accordingly, the calculation of derivatives is only possible for functions that are contained in the derivatives table of \code{\link{deriv}}.

Note, in contrast to \code{\link{model.matrix}}, \code{\link{factor}s} are not supported. 
}
\value{
The derivatives of the model (or design) matrix.

As in \code{\link{model.matrix}} there is an attribute \code{"assign"}.
}

\author{
Carmen van Meegen
}

\seealso{
\code{\link{deriv}} for more details on supported arithmetic operators and functions.

\code{\link{model.matrix}} for construction of a design (or model) matrix.
}

\examples{
## Several examples for the derivatives of a model matrix

dat <- data.frame(x1 = seq(-2, 2, length.out = 5))

model.matrix(~ 1, dat)
derivModelMatrix(~ 1, dat)

model.matrix(~ ., dat)
derivModelMatrix(~ ., dat)

model.matrix(~ . - 1, dat)
derivModelMatrix(~ . - 1, dat)

model.matrix(~ sin(x1) + I(x1^2), dat)
derivModelMatrix(~ sin(x1) + I(x1^2), dat)

dat <- cbind(dat, x2 = seq(1, 5, length.out = 5))

model.matrix(~ 1, dat)
derivModelMatrix(~ 1, dat)

model.matrix(~ .^2, dat)
derivModelMatrix(~ .^2, dat)

model.matrix(~ log(x2), dat)
derivModelMatrix(~ log(x2), dat)

model.matrix(~ x1:x2, dat)
derivModelMatrix(~ x1:x2, dat)

model.matrix(~ I(x1^2) * I(x2^3), dat)
derivModelMatrix(~ I(x1^2) * I(x2^3), dat)

model.matrix(~ sin(x1) + cos(x2) + atan(x1 * x2), dat)
derivModelMatrix(~ sin(x1) + cos(x2) + atan(x1 * x2), dat)
}
