\name{simulation}
\alias{simulation}
\alias{B_true}
\alias{T_domain}
\alias{X}
\alias{Y_fd}
\alias{Y_full}
\alias{Y_matrix}
\alias{eps}
\docType{data}
\title{
Function-on-Scalar linear model data, used to test FLAME method.
}
\description{
Data of a Function-on-Scalar linear model, it includes the design matrix,
the time domain, the response functions (both as pointwise evaluations on the time domain,
as projections  on the kernel basis, and as  \link[=fda]{fd} objects),  the true coefficients
of the Function-on-Scalar model and the
Gaussian errors.
}
\usage{data("simulation")}
\format{
  A set of variables associated to the Function-on-Scalar linear model.
  \itemize{
    \item{\code{X}}{ matrix. \code{N} \eqn{\times} \code{I} design matrix.}
    \item{\code{Y_matrix}}{ matrix. \code{J} \eqn{\times} \code{N} matrix containing the  \code{N}
    response functions projected  on the \code{J}-dimensional kernel basis.}
      \item{\code{T_domain}}{ vector. \code{m}-length vector of the time domain where the
    response functions are evaluated.}
     \item{\code{Y_full}}{ matrix. \code{N} \eqn{\times} \code{m} matrix containing the response functions evaluated on the \code{T_domain} grid.}
    \item{\code{Y_fd}}{ fd object. \code{N}-length \link[=fda]{fd} object of the repsonse functions. The same response functions as \code{Y_matrix} and \code{Y_full}, but projected on a 20 elements cubic Bspline basis.}
    \item{\code{B_true}}{ matrix. \code{J} \eqn{\times} \code{I0} matrix of the projection on the kernel basis of the \code{I0} significant predictor coefficients. }
    \item{\code{eps}}{ matrix. \code{N} \eqn{\times} \code{m} matrix of the random errors
    added to the model. Evaluation of the \code{N} functions on the \code{m}-dimensional
    time domain.}
  }
}
\details{
It contains data on a high-dimensional simulation setting
with \code{N = 500} and \code{I = 1000} to highlight
both the efficiency of FLAME in  estimation and in variable selection. Only
\code{I0 = 10} predictors, in fact, are meaningful for the response, the others have a null
effect on the \eqn{Y}'s.

The predictor matrix \code{X} is the standardized version of a matrix randomly sampled  from a \code{N} dimension Gaussian distribution with 0 average and identity covariance matrix. The true coefficients \eqn{\beta} are sampled from a Matern process with 0 average and parameters
\eqn{(\nu = 2.5, \textrm{range}  = 1/4, \sigma^2=1)}.

Observations \eqn{y(t)} are, then,  obtained as the sum of the contribution of all the predictors and a random noise \code{eps}, a 0-mean  Matern process with parameters \eqn{(\nu = 1.5, \textrm{range} = 1/4, \sigma^2=1)}.  Functions are sampled on a \code{m = 50} points grid.

See \link[=covMaterniso]{covMaterniso} for details on the covariance structure of
coefficients and errors.
}

\examples{
data(simulation)
}
\keyword{datasets}
