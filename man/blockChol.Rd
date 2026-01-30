\name{blockChol}
\alias{blockChol}

\title{Block Cholesky Decomposition}

\description{
\code{blockChol} calculates the block Cholesky decomposition of a partitioned matrix of the form
\deqn{A = \begin{pmatrix} K & R^{\rm T} \\ R & S \end{pmatrix}.}
}

\usage{
blockChol(K, R = NULL, S = NULL, tol = NULL)
}

\arguments{
  \item{K}{
	a real symmetric positive-definite square submatrix.
}
  \item{R}{
	an (optional) rectangular submatrix.
}
  \item{S}{
	an (optional) real symmetric positive-definite square submatrix.
}
  \item{tol}{
	an (optional) numeric tolerance, see \sQuote{Details}.
}
}

\details{
To obtain the block Cholesky decomposition
\deqn{\begin{pmatrix} K & R^{\rm T} \\ R & S \end{pmatrix} = \begin{pmatrix} L^{\rm T} & 0 \\ Q^{\rm T} & M^{\rm T} \end{pmatrix} \begin{pmatrix} L & Q \\ 0 & M \end{pmatrix}}
the following steps are performed:
\enumerate{
	\item{Calculate \eqn{K = L^{\rm T}L} with upper triangular matrix \eqn{L}.}
	\item{Solve \eqn{L^{\rm T}Q = R^{\rm T}} via forward substitution.}
	\item{Compute \eqn{N = S - Q^{\rm T}Q} the Schur complement of the block \eqn{K} of the matrix \eqn{A}.}
	\item{Determine \eqn{N = M^{\rm T}M} with upper triangular matrix \eqn{M}.}
}

The upper triangular matrices \eqn{L} and \eqn{M} in step 1 and 4 are obtained by \code{\link{chol}}.
Forward substitution in step 2 is carried out with \code{\link{backsolve}} and the option \code{transpose = TRUE}.

If \code{tol} is specified a regularization of the form \eqn{A_{\epsilon} = A + \epsilon I} is conducted.
Here, \code{tol} is the upper bound for the logarithmic condition number of \eqn{A_{\epsilon}}.
Then
\deqn{\epsilon = \max\left\{ \frac{\lambda_{\max}(\kappa(A) - e^{\code{tol}})}{\kappa(A)(e^{\code{tol}} - 1)}, 0 \right\}}
is chosen as the minimal "nugget" that is added to the diagonal of \eqn{A} to ensure \eqn{\log(\kappa(A_{\epsilon})) \leq} \code{tol}.

Within \pkg{gek} this function is used to calculate the block Cholesky decomposition of the correlation matrix with derivatives.
Here \code{K} is the Kriging correlation matrix. 
\code{R} is the matrix containing the first partial derivatives and \code{S} consists of the second partial derivatives of the correlation matrix \code{K}.}

\value{
\code{blockChol} returns a list with the following components:
	\item{L}{the upper triangular factor of the Cholesky decomposition of \code{K}.}
	\item{Q}{the solution of the triangular system \code{t(L) \%*\% Q == t(R)}.}
	\item{M}{the upper triangular factor of the Cholesky decomposition of the Schur complement \code{N}.}
If \code{R} or \code{S} are not specified, \code{Q} and \code{M} are set to \code{NULL},
i.e. only the Cholesky decomposition of \code{K} is calculated. 

The attribute \code{"eps"} gives the minimum \dQuote{nugget} that is added to the diagonal. 
}

\references{
Chen, J., Jin, Z., Shi, Q., Qiu, J., and Liu, W. (2013). Block Algorithm and Its Implementation for Cholesky Factorization.

Gustavson, F. G. and Jonsson, I. (2000). Minimal-storage high-performance Cholesky factorization via blocking and recursion. \emph{IBM Journal of Research and Development}, \bold{44}(6):823--850. \doi{10.1147/rd.446.0823}.

Ranjan, P., Haynes, R. and Karsten, R. (2011). A Computationally Stable Approach to Gaussian Process Interpolation of Deterministic Computer Simulation Data. \emph{Technometrics}, \bold{53}:366--378. \doi{10.1198/TECH.2011.09141}.
}

\author{
	Carmen van Meegen
}

\section{Warning}{
As in \code{\link{chol}} there is no check for symmetry.
}


\seealso{
\code{\link{chol}} for the Cholesky decomposition.

\code{\link{backsolve}} for backward substitution.

\code{\link{blockCor}} for computing a correlation matrix with derivatives.
}


\examples{
# Construct correlation matrix
x <- matrix(seq(0, 1, length = 5), ncol = 1)
res <- blockCor(x, theta = 1, covtype = "gaussian", derivatives = TRUE)
A <- cbind(rbind(res$K, res$R), rbind(t(res$R), res$S))

# Cholesky decomposition of correlation matix without derivatives
cholK <- blockChol(res$K)
cholK
cholK$L == chol(res$K)

# Cholesky decomposition of correlation matix with derivatives
cholA <- blockChol(res$K, res$R, res$S) 
cholA <- cbind(rbind(cholA$L, matrix(0, ncol(cholA$Q), nrow(cholA$Q))), 
	rbind(cholA$Q, cholA$M))
cholA
cholA == chol(A)

# Cholesky decomposition of correlation matix with derivatives with regularization
res <- blockCor(x, theta = 2, covtype = "gaussian", derivatives = TRUE)
A <- cbind(rbind(res$K, res$R), rbind(t(res$R), res$S))
try(blockChol(res$K, res$R, res$S))
blockChol(res$K, res$R, res$S, tol = 35)
}

\keyword{algebra}
\keyword{array}


