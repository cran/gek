\encoding{UTF-8}

\name{confint.gekm}
\alias{confint.gekm}


\title{Confint Method for a gekm Object}

\description{
Determines confidence intervals for the estimated regression coefficients.
}

\usage{
\method{confint}{gekm}(object, parm, level = 0.95, scale = FALSE, \dots)
}

\arguments{
  \item{object}{
an object of class \code{"gekm"}.
}
  \item{parm}{
a \code{\link{vector}} of numbers or names specifying the parameters for which the confidence intervals are to be calculated. By default, all parameters are considered.
}
  \item{level}{
confidence level for calculating confidence intervals. Default is \code{0.95}.
}
  \item{scale}{
\code{\link{logical}}. Should the estimated process variance be scaled? Default is \code{FALSE}, see \code{\link{sigma.gekm}} for details.
}
  \item{\dots}{
further arguments, currently not used.
}
}



\value{
A \code{\link{matrix}} with the lower and upper bounds of the confidence intervals for each parameter.
}


\references{
Koehler, J. and Owen, A. (1996). Computer Experiments. In Ghosh, S. and Rao, C. (eds.), \emph{Design and Analysis of Experiments}, volume 13 of \emph{Handbook of Statistics}, pp. 261–308. Elsevier Science. \doi{10.1016/S0169-7161(96)13011-X}.

Santner, T. J., Williams, B. J., and Notz, W. I. (2018). \emph{The Design and Analysis of Computer Experiments}. 2nd edition. Springer-Verlag.
}

\author{
Carmen van Meegen
}


\seealso{
\code{\link{gekm}} for fitting a (gradient-enhanced) Kriging model.

\code{\link{coef}} for extracting the (matrix of) coefficients.

\code{\link{vcov}} for calculating the covaraince matrix of the regression coefficients.
}


\examples{
## 1-dimensional example: Oakley and O’Hagan (2002)

# Define test function and its gradient 
f <- function(x) 5 + x + cos(x)
fGrad <- function(x) 1 - sin(x)

# Generate coordinates and calculate slopes
x <- seq(-5, 5, length = 5)
y <- f(x)
dy <- fGrad(x)
dat <- data.frame(x, y)
deri <- data.frame(x = dy)

# Fit gradient-enhanced Kriging model
gekm.1d <- gekm(y ~ ., data = dat, deriv = deri, covtype = "gaussian", theta = 1)

# Determine confidence intervals
confint(gekm.1d)
confint(gekm.1d, scale = TRUE)
confint(gekm.1d, parm = "x", scale = TRUE)
confint(gekm.1d, parm = 1, scale = TRUE)
}



