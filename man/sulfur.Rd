\name{sulfur}
\alias{sulfur}
\alias{sulfurGrad}

\title{Sulfur Model Function}

\description{
The sulfur function is defined by
\deqn{f_{\rm sulfur}(x) = -\frac{1}{2} S_0^2 (1 - A_c) T^2 (1 - R_s)^2 \bar{\beta} \Psi_e f_{\Psi_e} \frac{3 Q Y L}{A}}
with \eqn{x = (Q, Y, L, \Psi_e, \bar{\beta}, f_{\Psi_e}, T, 1-A_c, 1-R_s)}.
}

\usage{
sulfur(x, S_0 = 1366, A = 5.1e+14)
sulfurGrad(x, S_0 = 1366, A = 5.1e+14)
}

\arguments{
\item{x}{
a numeric \code{\link{vector}} of length \code{9} or a numeric \code{\link{matrix}} with \code{n} rows and \code{9} columns.
}
\item{S_0}{
solar constant in \eqn{\rm W/m^2}. Default is \code{1366}.
}
\item{A}{
surface area of the earth in \eqn{\rm m^2}. Default is \code{5.1e+14}.
}
}

\details{
The sulfur model function calculates the direct radiative forcing by sulfate aerosols \eqn{\rm [W/m^2]}.

\tabular{lrrl}{
	Input \tab Central value \tab Uncertainty factor \tab Description \cr
	\eqn{Q} \tab \eqn{71} \tab \eqn{1.15} \tab source strength of anthropogenic sulfur in \eqn{\rm Tg/yr} \cr
	\eqn{Y} \tab \eqn{0.5} \tab \eqn{1.5}  \tab fraction of \eqn{\rm SO_2} oxidized to \eqn{\rm SO_4^=} \cr
	\eqn{L} \tab \eqn{5.5} \tab \eqn{1.5}  \tab average lifetime of atmospheric \eqn{\rm SO_4^=} in \eqn{\rm days} \cr
	\eqn{\Psi_e} \tab \eqn{5} \tab \eqn{1.4}  \tab aerosol mass scattering efficiency in \eqn{\rm m^2/g} \cr
	\eqn{\bar{\beta}} \tab \eqn{0.3} \tab \eqn{1.3}  \tab fraction of light scattering into upward hemisphere \cr
	\eqn{f_{\Psi_e}} \tab \eqn{1.7} \tab \eqn{1.2}  \tab fractional increase in aerosol scattering efficiency due to hygroscopic growth \cr
	\eqn{T} \tab \eqn{0.76} \tab \eqn{1.2}  \tab atmospheric transmittance above aerosol layer \cr
	\eqn{1-A_c} \tab \eqn{0.39} \tab \eqn{1.1}  \tab fraction of earth not covered by cloud \cr
	\eqn{1-R_s} \tab \eqn{0.85} \tab \eqn{1.1}  \tab surface coalbedo \cr
}
The inputs are all log-normally distributed.
}

\value{
\code{sulfur} returns the function value of sulfur function at \code{x}.

\code{sulfurGrad} returns the gradient of sulfur function at \code{x}.
}

\references{
Charlson, R. J., Schwartz, S. E., Hales, J. M., Cess, R. D., Coakley, Jr., J. A., Hansen, J. E., and Hoffman, D. J. (1992). Climate Forcing by Anthropogenic Aerosols. \emph{Science}, \bold{255}:423--430. \doi{10.1126/science.255.5043.423}.

Penner, J. E., Charlson, R. J., Hales, J. M., Laulainen, N. S., Leifer, R., Novakov, T., Ogren, J., Radke, L. F., Schwartz, S. E., and Travis, L. (1994). Quantifying and Minimizing Uncertainty of Climate Forcing by Anthropogenic Aerosols. \emph{Bulletin of the American Meteorological Society}, \bold{75}(3):375--400. \doi{10.1175/1520-0477(1994)075<0375:QAMUOC>2.0.CO;2}.

Tatang, M. A., Pan, W., Prinn, R. G., and McRae, G. J. (1997). An Efficient Method for Parametric Uncertainty Analysis of Numerical Geophysical Model. \emph{Journal of Geophysical Research Atmospheres}, \bold{102}(18):21925--21932. \doi{10.1029/97JD01654}.
}

\author{
Carmen van Meegen
}

