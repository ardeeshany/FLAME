#' Representation on the time domain of a function defined by the coefficients of the kernel basis
#'
#' It computes the pointwise evaluation of a function (or a set of
#' functions) on the time domain, given its (their)  projection on the kernel basis.
#'
#' @param y matrix. \code{J} \eqn{\times} \code{N} matrix containing in
#'  column \eqn{n} the coefficients of the projection of the function \eqn{y_n}
#'  on the \code{J} eigenfunctions of the kernel.
#' @param eigenfun matrix. \code{m} \eqn{\times} \code{J} matrix containing
#' in each column the
#' point-wise evaluation of the eigenfunctions on the
#' kernel
#'
#' @return \code{N} \eqn{\times} \code{m} matrix containing in the row \eqn{n} the pointwise
#' evaluation of the function \eqn{y_n} on the domain \eqn{D} of length \code{m}.
#' @export
#'
#' @examples
#' data(SobolevKernel)
#' data(simulation)
#' projection_domain(Y_matrix, eigenvect) # projection of the data
#' # on the time domain seq(0, 1, length = 50)
#'
projection_domain <- function(y, eigenfun)
{
    y_mat <- eigenfun %*% y
    return(t(y_mat))
}
