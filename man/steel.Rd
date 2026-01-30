\name{steel}
\alias{steel}
\alias{steelGrad}

\title{Steel Column Function}

\description{
The steel column function is defined by
\deqn{f_{\rm steel}(x) = F_S - P \left(\frac{1}{2BD} + \frac{F_0 E_b}{BDH(E_b - P)} \right),}
with \eqn{P = P_1 + P_2 + P_3}, \eqn{E_b = \frac{\pi^2 EBDH^2}{2L^2}} and \eqn{x = (F_S, P_1, P_2, P_3, B, D, H, F_0, E)}.
}

\usage{
steel(x, L = 7500)
steelGrad(x, L = 7500)
}

\arguments{
\item{x}{
a numeric \code{\link{vector}} of length \code{9} or a numeric \code{\link{matrix}} with \code{n} rows and \code{9} columns.
}
\item{L}{
length in \eqn{\rm mm} of the steel column. Default is \code{7500}.
}
}

\details{
The steel column function describes the limite state function of a steel column with uncertain parameters.

\tabular{lrrrrl}{
	Input \tab Distribution \tab Mean \tab Standard Deviation \tab Description \cr
	\eqn{F_S} \tab \eqn{\mathcal{LN}} \tab \eqn{400} \tab \eqn{35} \tab yield stress in \eqn{\rm MPa} \cr
	\eqn{P_1} \tab \eqn{\mathcal{N}} \tab \eqn{500000} \tab \eqn{50000} \tab dead weight load in \eqn{\rm N} \cr
	\eqn{P_2} \tab \eqn{\mathcal{G}} \tab \eqn{600000} \tab \eqn{90000} \tab variable load in \eqn{\rm N} \cr
	\eqn{P_3} \tab \eqn{\mathcal{G}} \tab \eqn{600000} \tab \eqn{90000} \tab variable load in \eqn{\rm N} \cr
	\eqn{B} \tab \eqn{\mathcal{LN}} \tab \eqn{b} \tab \eqn{3} \tab flange breadth in \eqn{\rm mm} \cr
	\eqn{D} \tab \eqn{\mathcal{LN}} \tab \eqn{t} \tab \eqn{2} \tab flange thickness in \eqn{\rm mm} \cr
	\eqn{H} \tab \eqn{\mathcal{LN}} \tab \eqn{h} \tab \eqn{5} \tab profile height in \eqn{\rm mm} \cr
	\eqn{F_0} \tab \eqn{\mathcal{N}} \tab \eqn{30} \tab \eqn{10} \tab initial deflection in \eqn{\rm mm} \cr
	\eqn{E} \tab \eqn{\mathcal{W}} \tab \eqn{210000} \tab \eqn{4200} \tab Young's modulus in \eqn{\rm MPa} \cr
}
Here, \eqn{\mathcal{N}} is the normal distribution and \eqn{\mathcal{LN}} is the log-normal distribution.
Further, \eqn{\mathcal{G}} represents the Gumbel distribution and \eqn{\mathcal{W}} denotes the Weibull distribution.
}

\value{
\code{steel} returns the function value of steel column function at \code{x}.

\code{steelGrad} returns the gradient of steel column function at \code{x}.
}

\references{
Kuschel, N. and Rackwitz, R. (1997). Two Basic Problems in Reliability-Based Structural Optimization. \emph{Mathematical Methods of Operations Research}, \bold{46}(3):309--333.

Surjanovic, S. and Bingham, D. (2013). Virtual Library of Simulation Experiments: Test Functions and Datasets. \url{https://www.sfu.ca/~ssurjano/} (retrieved January 19, 2024). 
}


\author{
Carmen van Meegen
}

\seealso{
\code{\link{testfunctions}} for further test functions.
}
