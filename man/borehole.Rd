\name{borehole}
\alias{borehole}
\alias{boreholeGrad}

\title{Borehole Function}

\description{
The borehole function is defined by
\deqn{f_{\rm borehole}(x) = \frac{2 \pi T_u (H_u - H_l)}{\log(r/r_w)\left(1 + \frac{2 L T_u}{\log(r / r_w) r_w^2 K_w} + \frac{T_u}{T_l}\right)}}
with \eqn{x = (r_w, r, T_u, H_u, T_l, H_l, L, K_w)}.
}

\usage{
borehole(x)
boreholeGrad(x)
}

\arguments{
\item{x}{
a numeric \code{\link{vector}} of length \code{8} or a numeric \code{\link{matrix}} with \code{n} rows and \code{8} columns.
}
}

\details{
The borehole function calculates the water flow rate \eqn{\rm [m^3/yr]} through a borehole.

\tabular{lrrl}{
	Input \tab Domain \tab Distribution \tab Description \cr
	\eqn{r_w} \tab \eqn{[0.05, 0.15]} \tab \eqn{\mathcal{N}(0.1, 0.0161812)} \tab radius of borehole in \eqn{\rm m} \cr
	\eqn{r} \tab \eqn{[100, 50\,000]} \tab \eqn{\mathcal{LN}(7.71, 1.0056)} \tab radius of influence in \eqn{\rm m} \cr
	\eqn{T_u} \tab \eqn{[63070, 115600]} \tab \eqn{\mathcal{U}(63070, 115600)} \tab transmissivity of upper aquifer in \eqn{\rm m^2/yr} \cr
	\eqn{H_u} \tab \eqn{[990, 1100]} \tab \eqn{\mathcal{U}(990, 1110)} \tab potentiometric head of upper aquifer in \eqn{\rm m} \cr
	\eqn{T_l} \tab \eqn{[63.1, 116]} \tab \eqn{\mathcal{U}(63.1, 116)} \tab transmissivity of lower aquifer in \eqn{\rm m^2/yr} \cr
	\eqn{H_l} \tab \eqn{[700, 820]} \tab \eqn{\mathcal{U}(700, 820)} \tab potentiometric head of lower aquifer in \eqn{\rm m} \cr
	\eqn{L} \tab \eqn{[1120, 1680]} \tab \eqn{\mathcal{U}(1120, 1680)} \tab length of borehole in \eqn{\rm m} \cr
	\eqn{K_w} \tab \eqn{[9855, 12045]} \tab \eqn{\mathcal{U}(9855, 12045)} \tab hydraulic conductivity of borehole in \eqn{\rm m/yr} \cr
}
Note, \eqn{\mathcal{N}(\mu, \sigma)} represents the normal distribution with expected value \eqn{\mu} and standard deviation \eqn{\sigma}
and \eqn{\mathcal{LN}(\mu, \sigma)} is the log-normal distribution with mean \eqn{\mu} and standard deviation \eqn{\sigma} of the logarithm.
Further, \eqn{\mathcal{U}(a,b)} denotes the continuous uniform distribution over the interval \eqn{[a,b]}.
}

\value{
\code{borehole} returns the function value of borehole function at \code{x}.

\code{boreholeGrad} returns the gradient of borehole function at \code{x}.
}

\references{
Harper, W. V. and Gupta, S. K. (1983). Sensitivity/Uncertainty Analysis of a Borehole Scenario Comparing Latin Hypercube Sampling and Deterministic Sensitivity Approaches.
BMI/ONWI-516, Office of Nuclear Waste Isolation, Battelle Memorial Institute, Columbus, OH.

Morris, M., Mitchell, T., and Ylvisaker, D. (1993). Bayesian Design and Analysis of Computer Experiments: Use of Derivatives in Surface Prediction. \emph{Technometrics}, \bold{35}(3):243--255. \doi{10.1080/00401706.1993.10485320}.
}

\author{
Carmen van Meegen
}


\seealso{
\code{\link{gekm}} for another example.
}


\examples{
# List of inputs with their distributions and their respective ranges
inputs <- list("r_w" = list(dist = "norm", mean =  0.1, sd = 0.0161812, min = 0.05, max = 0.15),
	"r" = list(dist = "lnorm", meanlog = 7.71, sdlog = 1.0056, min = 100, max = 50000),
	"T_u" = list(dist = "unif", min = 63070, max = 115600),
	"H_u" = list(dist = "unif", min = 990, max = 1110),
	"T_l" = list(dist = "unif", min = 63.1, max = 116),
	"H_l" = list(dist = "unif", min = 700, max = 820),
	"L" = list(dist = "unif", min = 1120, max = 1680),
	# for a more nonlinear, nonadditive function, see Morris et al. (1993)
	"K_w" = list(dist = "unif", min = 1500, max = 15000))

# Function for Monte Carlo simulation
samples <- function(x, N = 10^5){
	switch(x$dist,
		"norm" = rnorm(N, x$mean, x$sd),
		"lnorm" = rlnorm(N, x$meanlog, x$sdlog),
		"unif" = runif(N, x$min, x$max))
}

# Uncertainty distribution of the water flow rate
set.seed(1)
X <- sapply(inputs, samples)
y <- borehole(X)
hist(y, breaks = 50, xlab = expression(paste("Water flow rate ", group("[", m^3/yr, "]"))), 
	main = "", freq = FALSE)
}