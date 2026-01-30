\encoding{UTF-8}

\name{testfunctions}
\alias{testfunctions}

\title{Test Functions in \pkg{gek}}

\description{
Overview of test functions and their associated gradients available in \pkg{gek}.
}

\section{2-dimensional test functions for optimization}{
\itemize{
\item{\code{\link{branin}}: Branin-Hoo function}
\item{\code{\link{camel3}}: Three-hump camel function}
\item{\code{\link{camel6}}: Six-hump camel function}
\item{\code{\link{himmelblau}}: Himmelblaus's function}
}
}

\section{Multi-dimensional test functions for optimization}{
\itemize{
\item{\code{\link{banana}}: Rosenbrock's Banana function}
\item{\code{\link{cigar}}: Bent Cigar function}
\item{\code{\link{griewank}}: Griewank function}
\item{\code{\link{qing}}: Qing function}
\item{\code{\link{rastrigin}}: Rastrigin function}
\item{\code{\link{schwefel}}: Schwefel function}
\item{\code{\link{sphere}}: Sphere function}
\item{\code{\link{styblinski}}: Styblinski-Tang function}	
}	
}

\section{Test functions for uncertainty quantification}{
\itemize{
\item{\code{\link{borehole}}: Borehole function}
\item{\code{\link{steel}}: Steel column function}
\item{\code{\link{short}}: Short column function}
\item{\code{\link{sulfur}}: Sulfur model function}
}
}

\references{
Branin, Jr., F. H. (1972). Widely Convergent Method of Finding Multiple Solutions of Simultaneous Nonlinear Equations. \emph{IBM Journal of Research and Development}, \bold{16}(5):504--522.

Jamil, M. and Yang, X.-S. (2013). A Literature Survey of Benchmark Functions for Global Optimization Problems. \emph{International Journal of Mathematical Modelling and Numerical Optimisation}, \bold{4}(2):150-–194. \doi{10.1504/IJMMNO.2013.055204}.

Harper, W. V. and Gupta, S. K. (1983). Sensitivity/Uncertainty Analysis of a Borehole Scenario Comparing Latin Hypercube Sampling and Deterministic Sensitivity Approaches.
BMI/ONWI-516, Office of Nuclear Waste Isolation, Battelle Memorial Institute, Columbus, OH.

Himmelblau, D. (1972). Applied Nonlinear Programming. McGraw-Hill. ISBN 0-07-028921-2.

Kuschel, N. and Rackwitz, R. (1997). Two Basic Problems in Reliability-Based Structural Optimization. \emph{Mathematical Methods of Operations Research}, \bold{46}(3):309--333.

Morris, M., Mitchell, T., and Ylvisaker, D. (1993). Bayesian Design and Analysis of Computer Experiments: Use of Derivatives in Surface Prediction. \emph{Technometrics}, \bold{35}(3):243--255. \doi{10.1080/00401706.1993.10485320}.

Plevris, V. and Solorzano, G. (2022). A Collection of 30 Multidimensional Functions for Global Optimization Benchmarking. \emph{Data}, \bold{7}(4):46. \doi{10.3390/data7040046}. 

Rosenbrock, H. H. (1960). An Automatic Method for Finding the Greatest or least Value of a Function. \emph{The Computer Journal}, \bold{3}(3):175--184. \doi{10.1093/comjnl/3.3.175}.

Surjanovic, S. and Bingham, D. (2013). Virtual Library of Simulation Experiments: Test Functions and Datasets. \url{https://www.sfu.ca/~ssurjano/} (retrieved January 19, 2024). 
}

\author{
Carmen van Meegen
}

\seealso{
\code{\link{banana}},
\code{\link{borehole}},
\code{\link{branin}},
\code{\link{camel3}},
\code{\link{camel6}},
\code{\link{cigar}},
\code{\link{griewank}},
\code{\link{himmelblau}},
\code{\link{qing}},
\code{\link{rastrigin}},
\code{\link{schwefel}},
\code{\link{short}},
\code{\link{sphere}},
\code{\link{steel}},
\code{\link{styblinski}},
\code{\link{sulfur}}
}
