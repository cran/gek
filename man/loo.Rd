\encoding{UTF-8}

\name{loo}
\alias{loo}
\alias{loo.gekm}


\title{
Leave-One-Out Cross-Validation
}

\description{
Calculation of the leave-one-out prediction, standard deviation and confidence intervals of a \code{gekm} object.
}

\usage{
\method{loo}{gekm}(object, reestim = TRUE, sd.fit = TRUE, scale = FALSE, 
	df = NULL, interval = c("none", "confidence"), level = 0.95, \dots)
}

\arguments{
  \item{object}{
an object of class \code{"gekm"}.
}
  \item{reestim}{
\code{\link{logical}}. Should the regression coefficients be re-estimated? Default is \code{TRUE}.
}
  \item{sd.fit}{
\code{\link{logical}}. Should the standard deviation of the prediction, i.e., the root mean squared error, be computed? Default is \code{TRUE}.
}
  \item{scale}{
\code{\link{logical}}. Should the estimated process variance be scaled? Default is \code{FALSE}, see \code{\link{sigma.gekm}} for details.
}
  \item{df}{
degrees of freedom for the \eqn{t} distribution. Default is \code{NULL}, see \sQuote{Details}.
}
  \item{interval}{
a \code{\link{character}} that specifies the type of interval calculation.
}
  \item{level}{
confidence level for calculating confidence intervals. Default is \code{0.95}.
}
  \item{\dots}{
further arguments, currently not used.
}
}

\details{
For \code{reestim = TRUE} (default), the formulas form Dubrule (1983) are used.
These enable a faster calculation of the leave-one-out prediction and the associated standard deviation, especially for a large number of observations.
However, with few observations, the re-estimated regression coefficients may differ considerably from those based on the entire data set.
Note that the process variance and correlation parameters are not re-estimated.
}


\value{
The \code{loo} method of class \code{"gekm"} returns a \code{\link{vector}} of leave-one-out predictions, if \code{sd.fit = FALSE} and \code{interval = "none"}.
As with \code{\link{predict.gekm}}, setting \code{sd.fit = FALSE} and \code{interval = "confidence"} generates a \code{\link{matrix}} with the leave-one-out predicted values and the lower and upper limits of the confidence intervals.
For \code{sd.fit = TRUE}, a \code{\link{list}} with the following components is returned:
  \item{fit.loo}{either a \code{\link{vector}} or a \code{\link{matrix}}, as described above.}
  \item{sd.loo}{leave-one-out predicted standard deviation.}
}


\references{
Bachoc, F. (2013). Cross Validation and Maximum Likelihood Estimations of Hyper-parameters of Gaussian Processes with Model Misspecification. \emph{Computational Statistics and Data Analysis}, \bold{66}:55--69. \doi{10.1016/j.csda.2013.03.016}.

Dubrule, O. (1983). Cross Validation of Kriging in a Unique Neighborhood. \emph{Mathematical Geology}, \bold{15}:687--699. \doi{10.1007/BF01033232}.

Martin, J. D. and Simpson, T. W. (2005). Use of Kriging Models to Approximate Deterministic Computer Models. \emph{AIAA Journal}, \bold{43}(4):853--863. \doi{10.2514/1.8650}.
}

\author{
Carmen van Meegen
}


\seealso{
\code{\link{gekm}} for fitting a (gradient-enhanced) Kriging model.

\code{\link{predict.gekm}} for prediction at new data points based on a model of class \code{"gekm"}. 

\code{\link{plot.gekm}} for plotting the results of a leave-one-out cross-validation.
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
gekm.1d <- gekm(y ~ x, data = dat, deriv = deri, covtype = "gaussian", theta = 1)

# Perform leave-one-out cross-validation
loo(gekm.1d)
loo(gekm.1d, sd.fit = FALSE)
loo(gekm.1d, sd.fit = FALSE, reestim = FALSE)
loo(gekm.1d, sd.fit = TRUE, scale = TRUE)
loo(gekm.1d, sd.fit = TRUE, reestim = FALSE, scale = TRUE)
loo(gekm.1d, sd.fit = FALSE, interval = "confidence")
loo(gekm.1d, sd.fit = TRUE, interval = "confidence")
}



