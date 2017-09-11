\name{SobolevKernel}
\alias{SobolevKernel}
\alias{eigenval}
\alias{eigenvect}
\alias{derivatives}
\docType{data}
\title{
An example of Sobolev kernel
}
\description{
Eigenvalues, eigenfunctions (with derivatives) of a Sobolev kernel on the
domain \eqn{D = [0,1]} with smoothness parameter \eqn{\sigma = 8}.
}
\usage{data("SobolevKernel")}
\format{
  A set of variables associated to the Sobolev kernel
  \describe{
    \item{\code{eigenval}}{vector. \code{J}-dimensional vector of the first \code{J}
    eigenvalues of the kernel.}
    \item{\code{eigenvect}}{matrix. \code{m} \eqn{\times} \code{J} matrix of the
    evaluation of the \code{J} eigenfunctions on the time domain \code{seq(0, 1, length = m)}, where \code{m = 50}.}
    \item{\code{derivatives}}{matrix. \code{(m-1)} \eqn{\times} \code{J} matrix of the
    punctual evaluation of the derivatives of the \code{J} eigenfunctions.}}
}
\details{
The parameter \code{J} is selected such that:
\deqn{
\sum_{j=1}^J \theta_j \geq \textrm{thres} \sum_{j=1}^\infty \theta_j.
}
with \code{thres} in this case equal to 0.99.
}
\examples{
data(SobolevKernel)
}
\keyword{datasets}
