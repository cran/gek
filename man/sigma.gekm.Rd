\encoding{UTF-8}

\name{sigma.gekm}
\alias{sigma.gekm}

\title{Extract Process Standard Deviation}

\description{
Extract the estimated process standard deviation of a Kriging model with or without derivatives.
}

\usage{
\method{sigma}{gekm}(object, scale = FALSE, \dots)
}

\arguments{
  \item{object}{
an object of class \code{"gekm"}.
}
  \item{scale}{
\code{\link{logical}}. Should the estimated process standard deviation be scaled? Default is \code{FALSE}, see \sQuote{Details}.
}
  \item{\dots}{
further arguments, currently not used.
}}

\details{
By default, the process variance \eqn{\sigma^2} is estimated using the maximum likelihood estimator, which uses \code{nobs} in the denominator, where \code{nobs} is the total number of observations used to fit the model.
Note for gradient-enhanced Kriging: \code{nobs = n + n * d} with \code{n} and \code{d} being the number of response values and inputs, respectively. 

Setting \code{scale = TRUE} replaces the denominator \code{nobs} with \code{nobs - p - 2}, where \code{p} is the number of regression coefficients.
If the correlation parameters are known and weak priors are assumed for the hyperparameters (the regression coefficients and the process variance), i.e., \eqn{f(\beta, \sigma^2) \propto \sigma^{-2}}, 
this leads to the Bayesian estimator of the process variance.
}

\value{The (scaled) estimated process standard deviation.}

\references{
Morris, M., Mitchell, T., and Ylvisaker, D. (1993). Bayesian Design and Analysis of Computer Experiments: Use of Derivatives in Surface Prediction. \emph{Technometrics}, \bold{35}(3):243--255. \doi{10.1080/00401706.1993.10485320}.

Oakley, J. and O'Hagan, A. (2002). Bayesian Inference for the Uncertainty Distribution of Computer Model Outputs. \emph{Biometrika}, \bold{89}(4):769--784. \doi{10.1093/biomet/89.4.769}.

O'Hagan, A. (1991). Bayes-Hermite Quadrature. \emph{Journal of Statistical Planning an Inference}, \bold{29}(3):245--260. \doi{10.1016/0378-3758(91)90002-V}.

O'Hagan, A., Kennedy, M. C., and Oakley, J. E. (1999). Uncertainty Analysis and Other Inference Tools for Complex Computer Codes. In \emph{Bayesian Statistics 6}, Ed. J. M. Bernardo, J. O. Berger, A. P. Dawid and A .F. M. Smith, 503--524. Oxford University Press. 

Park, J.-S. and Beak, J. (2001). Efficient Computation of Maximum Likelihood Estimators in a Spatial Linear Model with Power Exponential Covariogram. \emph{Computers & Geosciences}, \bold{27}(1):1--7. \doi{10.1016/S0098-3004(00)00016-9}.

Santner, T. J., Williams, B. J., and Notz, W. I. (2018). \emph{The Design and Analysis of Computer Experiments}. 2nd edition. Springer-Verlag.

Zimmermann, R. (2015). On the Condition Number Anomaly of Gaussian Correlation Matrices. \emph{Linear Algebra and its Applications}, \bold{466}:512-–526. \doi{10.1016/j.laa.2014.10.038}.
}


\author{
Carmen van Meegen
}

\seealso{
\code{\link{gekm}} for fitting a (gradient-enhanced) Kriging model.
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
km.1d <- gekm(y ~ x, data = dat, covtype = "gaussian", theta = 1)
gekm.1d <- gekm(y ~ x, data = dat, deriv = deri, covtype = "gaussian", theta = 1)

# Extact estimated process standard deviation
sigma(km.1d)
sigma(gekm.1d)					
sigma(gekm.1d, scale = TRUE)
}

\keyword{models}
