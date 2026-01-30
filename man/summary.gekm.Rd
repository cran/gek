\encoding{UTF-8}

\name{summary.gekm}
\alias{summary.gekm}
\alias{print.summary.gekm}

\title{Summary Method for a gekm Object}

\description{
Summarizing (Gradient-Enhanced) Kriging Models.
}

\usage{
\method{summary}{gekm}(object, scale = FALSE, \dots)

\method{print}{summary.gekm}(x, digits = 4L, \dots)
}

\arguments{
  \item{object}{
an object of class \code{"gekm"}.
}
  \item{x}{
an object of class \code{"summary.gekm"}.
}
  \item{scale}{
\code{\link{logical}}. Should the estimated process standard deviation be scaled? Default is \code{FALSE}, see \code{\link{sigma.gekm}} for details.
}
  \item{digits}{
number of digits to be used for the \code{print} method.
}
  \item{\dots}{
further arguments passed to \code{\link{printCoefmat}} in the \code{print} method.
}}


\value{The \code{summary} method for an object of class \code{"gekm"} returns a list with the following components:

\item{call}{
the matched call of object.
}
\item{terms}{
the \code{\link{terms}} object used.
}
\item{coefficients}{
a \code{\link{matrix}} with the estimated regression coefficients.
}
\item{sigma}{
the estimated (scaled) process standard deviation.
}
\item{df}{
degrees of freedom, i.e. the number of observations used to fit the model minus the number of regression coefficients.
}
\item{cov.scaled}{
the (scaled) covariance matrix of the estimated regression coefficients.
}
\item{covtype}{
the name of the correlation function.
}
\item{theta}{
the (estimated) correlation parameteres.
}
}

\references{
Morris, M., Mitchell, T., and Ylvisaker, D. (1993). Bayesian Design and Analysis of Computer Experiments: Use of Derivatives in Surface Prediction. \emph{Technometrics}, \bold{35}(3):243--255. \doi{10.1080/00401706.1993.10485320}.

Oakley, J. and O'Hagan, A. (2002). Bayesian Inference for the Uncertainty Distribution of Computer Model Outputs. \emph{Biometrika}, \bold{89}(4):769--784. \doi{10.1093/biomet/89.4.769}.

Park, J.-S. and Beak, J. (2001). Efficient Computation of Maximum Likelihood Estimators in a Spatial Linear Model with Power Exponential Covariogram. \emph{Computers & Geosciences}, \bold{27}(1):1--7. \doi{10.1016/S0098-3004(00)00016-9}.

Rasmussen, C. E. and Williams, C. K. I. (2006). \emph{Gaussian Processes for Machine Learning}. The MIT Press. \url{https://gaussianprocess.org/gpml/}.

Ripley, B. D. (1981). \emph{Spatial Statistics}. John Wiley & Sons. \doi{10.1002/0471725218}. 

Sacks, J., Welch, W. J., Mitchell, T. J., and Wynn, H. P. (1989). Design and Analysis of Computer Experiments. \emph{Statistical Science}, \bold{4}(4):409--423. \doi{10.1214/ss/1177012413}.

Santner, T. J., Williams, B. J., and Notz, W. I. (2018). \emph{The Design and Analysis of Computer Experiments}. 2nd edition. Springer-Verlag.

Stein, M. L. (1999). \emph{Interpolation of Spatial Data: Some Theory for Kriging}. Springer Series in Statistics. Springer-Verlag. \doi{10.1007/978-1-4612-1494-6}.

Zimmermann, R. (2015). On the Condition Number Anomaly of Gaussian Correlation Matrices. \emph{Linear Algebra and its Applications}, \bold{466}:512-–526. \doi{10.1016/j.laa.2014.10.038}.
}


\author{
Carmen van Meegen
}

\seealso{
\code{\link{gekm}} for fitting a (gradient-enhanced) Kriging model.

\code{\link{coef}} for extracting the (matrix of) coefficients.

\code{\link{vcov}} for calculating the covaraince matrix of the regression coefficients.

\code{\link{confint}} for computing confidence intervals for the regression coefficients.
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

# Fit (gradient-enhanced) Kriging model
km.1d <- gekm(y ~ . + I(x^2), data = dat, covtype = "gaussian", theta = 1)
gekm.1d <- gekm(y ~ . + I(x^2), data = dat, deriv = deri, covtype = "gaussian", theta = 1)

# Model summaries
summary(km.1d)
summary(gekm.1d)
summary(gekm.1d, scale = TRUE)
}

\keyword{models}
