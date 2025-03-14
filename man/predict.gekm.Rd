\encoding{UTF-8}

\name{predict.gekm}
\alias{predict.gekm}

\title{
Predict Method for (Gradient-Enhanced) Kriging Model Fits
}
\description{
Predicted values and standard deviations based on a \code{gekm} model.
}

\usage{
\method{predict}{gekm}(object, newdata, \dots)
}

\arguments{
  \item{object}{
an object of class \code{"gekm"}.
}
  \item{newdata}{
a \code{\link{data.frame}} containing the points where to perform predictions. 
}
  \item{\dots}{
further arguments, currently not used.
}
}

\value{
  \item{fit}{predicted mean computed at newdata.}
  \item{sd.fit}{predicted standard deviation of predicted mean.}
}


\references{
Cressie, N. A. C. (1993). \emph{Statistics for Spartial Data}. John Wiley & Sons. \doi{10.1002/9781119115151}.

Koehler, J. and Owen, A. (1996). Computer Experiments. In Ghosh, S. and Rao, C. (eds.), \emph{Design and Analysis of Experiments}, volume 13 of \emph{Handbook of Statistics}, pp. 261–308. Elsevier Science. \doi{10.1016/S0169-7161(96)13011-X}.

Krige, D. G. (1951). A Statistical Approach to Some Basic Mine Valuation Problems on the Witwatersrand. \emph{Journal of the Southern African Institute of Mining and Metallurgy}, \bold{52}(6):199--139. 

Laurent, L., Le Riche, R., Soulier, B., and Boucard, PA. (2019). An Overview of Gradient-Enhanced Metamodels with Applications. \emph{Archives of Computational Methods in Engineering}, \bold{26}(1):61--106. \doi{10.1007/s11831-017-9226-3}.

Martin, J. D. and Simpson, T. W. (2005). Use of Kriging Models to Approximate Deterministic Computer Models. \emph{AIAA Journal}, \bold{43}(4):853--863. \doi{10.2514/1.8650}.

Morris, M., Mitchell, T., and Ylvisaker, D. (1993). Bayesian Design and Analysis of Computer Experiments: Use of Derivatives in Surface Prediction. \emph{Technometrics}, \bold{35}(3):243--255. \doi{10.1080/00401706.1993.10485320}.

Oakley, J. and O'Hagan, A. (2002). Bayesian Inference for the Uncertainty Distribution of Computer Model Outputs. \emph{Biometrika}, \bold{89}(4):769--784. \doi{10.1093/biomet/89.4.769}.

O'Hagan, A., Kennedy, M. C., and Oakley, J. E. (1999). Uncertainty Analysis and Other Inference Tools for Complex Computer Codes. In \emph{Bayesian Statistics 6}, Ed. J. M. Bernardo, J. O. Berger, A. P. Dawid and A .F. M. Smith, 503--524. Oxford University Press. 

O'Hagan, A. (2006). Bayesian Analysis of Computer Code Outputs: A Tutorial. \emph{Reliability Engineering & System Safet}, \bold{91}(10):1290--1300. \doi{10.1016/j.ress.2005.11.025}.

Park, J.-S. and Beak, J. (2001). Efficient Computation of Maximum Likelihood Estimators in a Spatial Linear Model with Power Exponential Covariogram. \emph{Computers & Geosciences}, \bold{27}(1):1--7. \doi{10.1016/S0098-3004(00)00016-9}.

Ranjan, P., Haynes, R. and Karsten, R. (2011). A Computationally Stable Approach to Gaussian Process Interpolation of Deterministic Computer Simulation Data. \emph{Technometrics}, \bold{53}:366--378. \doi{10.1198/TECH.2011.09141}.

Rasmussen, C. E. and Williams, C. K. I. (2006). \emph{Gaussian Processes for Machine Learning}. The MIT Press. \url{https://gaussianprocess.org/gpml/}.

Ripley, B. D. (1981). \emph{Spatial Statistics}. John Wiley & Sons. \doi{10.1002/0471725218}. 

Sacks, J., Welch, W. J., Mitchell, T. J., and Wynn, H. P. (1989). Design and Analysis of Computer Experiments. \emph{Statistical Science}, \bold{4}(4):409--423. \doi{10.1214/ss/1177012413}.

Santner, T. J., Williams, B. J., and Notz, W. I. (2018). \emph{The Design and Analysis of Computer Experiments}. 2nd edition. Springer-Verlag.

Stein, M. L. (1999). \emph{Interpolation of Spatial Data: Some Theory for Kriging}. Springer Series in Statistics. Springer-Verlag. \doi{10.1007/978-1-4612-1494-6}.
}

\author{
Carmen van Meegen
}


\seealso{
\code{\link{gekm}} for fitting a (gradient-enhanced) Kriging model.

\code{\link{tangents}} for drawing tangent lines.
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

# Fit Kriging model
km.1d <- gekm(y ~ x, data = dat, covtype = "gaussian", theta = 1)

# Fit gradient-enhanced Kriging model
gekm.1d <- gekm(y ~ x, data = dat, deriv = deri, covtype = "gaussian", theta = 1)

# Generate new data for prediction
newdat <- data.frame(x = seq(-6, 6, length = 100))

# Prediction for both models
pred.km.1d <- predict(km.1d, newdat)
pred.gekm.1d  <- predict(gekm.1d, newdat)

# Draw curves
curve(f(x), from = -6, to = 6, lwd = 2)
lines(newdat$x, pred.km.1d$fit, type = "l", col = rgb(0.5, 0.5, 0.5), lwd = 2, lty = 2)
lines(newdat$x, pred.gekm.1d$fit, type = "l", col = rgb(0, 0.5, 1), lwd = 2, lty = 2)

# Add pointswise confidence intervals
lower.km.1d <- pred.km.1d$fit - qnorm(0.975) * pred.km.1d$sd
upper.km.1d <- pred.km.1d$fit + qnorm(0.975) * pred.km.1d$sd
lower.gekm.1d <- pred.gekm.1d$fit - qnorm(0.975) * pred.gekm.1d$sd
upper.gekm.1d <- pred.gekm.1d$fit + qnorm(0.975) * pred.gekm.1d$sd
 
polygon(c(newdat$x, rev(newdat$x)),	c(lower.km.1d, rev(upper.km.1d)),
	col = rgb(0.5, 0.5, 0.5, alpha = 0.1), border = NA)
polygon(c(newdat$x, rev(newdat$x)),	c(lower.gekm.1d, rev(upper.gekm.1d)),
	col = rgb(0, 0.5, 1, alpha = 0.1), border = NA)

# Add tangent lines and observations
tangents(x, y, dy, lwd = 2, col = "red")
points(x, y, pch = 16)

# Add legend
legend("topleft", lty = c(1, 2, 2), col = c("black", rgb(0.5, 0.5, 0.5), rgb(0, 0.5, 1)),
	legend = c("Test function", "Kriging", "GEK"), lwd = 2, bty = "n")


## 2-dimensional example: Branin-Hoo function

# Generate a grid for training
n <- 5
x1 <- seq(-5, 10, length = n)
x2 <- seq(0, 15, length = n)
x <- expand.grid(x1 = x1, x2 = x2)
y <- branin(x)
dy <- braninGrad(x)
dat <- data.frame(x, y)
deri <- data.frame(dy)

# Fit Kriging model
km.2d <- gekm(y ~ ., data = dat)
km.2d

# Fit gradient-enhanced Kriging model
gekm.2d <- gekm(y ~ ., data = dat, deriv = deri)
gekm.2d

# Generate new data for prediction
n.grid <- 50
x1.grid <- seq(-5, 10, length = n.grid)
x2.grid <- seq(0, 15, length = n.grid)
newdat <- expand.grid(x1 = x1.grid, x2 = x2.grid)

# Prediction for both models and actual outcome
pred.km.2d <- predict(km.2d, newdat)
pred.gekm.2d  <- predict(gekm.2d, newdat)
truth <- outer(x1.grid, x2.grid, function(x1, x2) branin(cbind(x1, x2)))

# Contour plots of predicted and actual output
par(mfrow = c(1, 3), oma = c(3.5, 3.5, 0, 0.2), mar = c(0, 0, 1.5, 0))
contour(x1.grid, x2.grid, truth, nlevels = 30, 
	levels = seq(0, 300, 10), main = "Branin-Hoo")
points(x, pch = 16)
contour(x1.grid, x2.grid, matrix(pred.km.2d$fit, nrow = n.grid, ncol = n.grid),
	levels = seq(0, 300, 10), main = "Kriging", yaxt = "n")
points(x, pch = 16)
contour(x1.grid, x2.grid, matrix(pred.gekm.2d$fit, nrow = n.grid, ncol = n.grid),
	levels = seq(0, 300, 10), main = "GEK", yaxt = "n")
points(x, pch = 16)
mtext(side = 1, outer = TRUE, line = 2.5, expression(x[1]))
mtext(side = 2, outer = TRUE, line = 2, expression(x[2]))

# Contour plots of predicted variance
par(mfrow = c(1, 2), oma = c(3.5, 3.5, 0, 0.2), mar = c(0, 0, 1.5, 0))
contour(x1.grid, x2.grid, matrix(pred.km.2d$sd.fit^2, nrow = n.grid, ncol = n.grid),
	main = "Kriging variance")
points(x, pch = 16)
contour(x1.grid, x2.grid, matrix(pred.gekm.2d$sd.fit^2, nrow = n.grid, ncol = n.grid),
	main = "GEK variance", yaxt = "n")
points(x, pch = 16)
mtext(side = 1, outer = TRUE, line = 2.5, expression(x[1]))
mtext(side = 2, outer = TRUE, line = 2, expression(x[2]))
}



