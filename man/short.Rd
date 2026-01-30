\name{short}
\alias{short}
\alias{shortGrad}

\title{Short Column Function}

\description{
The short column function is defined by
\deqn{f_{\rm short}(x) = 1 - \frac{4M}{b h^2 Y} - \frac{P^2}{b^2 h^2 Y^2}}
with \eqn{x = (Y, M, P)}.
}

\usage{
short(x, b = 5, h = 15)
shortGrad(x, b = 5, h = 15)
}

\arguments{
\item{x}{
a numeric \code{\link{vector}} of length \code{3} or a numeric \code{\link{matrix}} with \code{n} rows and \code{3} columns.
}
\item{b}{
width of the cross-section in \eqn{\rm mm} of the short column. Default is \code{5}.
}
\item{h}{
depth of the cross-section in \eqn{\rm mm} of the short column. Default is \code{15}.
}
}

\details{
The short column function describes the limite state function of a short column with uncertain material properties and loads.

\tabular{lrrrrl}{
	Input \tab Distribution \tab Mean \tab Standard deviation \tab Description \cr
	\eqn{Y} \tab \eqn{\mathcal{LN}} \tab \eqn{5} \tab \eqn{0.5} \tab yield stress in \eqn{\rm MPa} \cr
	\eqn{M} \tab \eqn{\mathcal{N}} \tab \eqn{2000} \tab \eqn{400} \tab bending moment in \eqn{\rm MNm} \cr
	\eqn{P} \tab \eqn{\mathcal{N}} \tab \eqn{500} \tab \eqn{100} \tab axial force in \eqn{\rm MPa} \cr
}
The bending moment and the axial force are correlated with \eqn{\rm{Cor}}\eqn{(M, P) = 0.5}.
Note, \eqn{\mathcal{N}} represents the normal distribution and \eqn{\mathcal{LN}} is the log-normal distribution.
}

\value{
\code{short} returns the function value of short column function at \code{x}.

\code{shortGrad} returns the gradient of short column function at \code{x}.
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
