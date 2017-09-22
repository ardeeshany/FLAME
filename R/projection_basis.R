#' Projection on the kernel basis
#'
#' It computes the  projection of a set of
#'  functions on the kernel basis defined by a set of eigenvalues and
#' eigenfunctions
#'
#' @param y matrix. \code{N} \eqn{\times} \code{m} matrix containing in
#'  row \eqn{n} the point-wise evaluation of the function \eqn{y_n}.
#' @param eigenvect matrix. \code{m} \eqn{\times} \code{J} matrix
#' containing in  column \eqn{j} the point-wise evaluation of an eigenfunction
#' \eqn{v_j} of the kernel evaluated on an equispaced time domain grid.
#' @param M_integ scalar, integer. Number of points of the domain \eqn{D}
#' of the functions (\code{m}) divided by the width of the domain \eqn{|D|}.
#' Constant value used to compute the integrals.
#'
#' @return
#' \code{J} \eqn{\times} \code{N} matrix (\eqn{Y}) containing in position
#' \eqn{(j,n)} the integral of the product
#' of the eigenfunction \eqn{v_j} with the function \eqn{y_n}, i.e. the coefficient
#' correspondent to \eqn{v_j} of the projection of \eqn{y_n} on the kernel basis.
#' \deqn{
#' Y[i,j] = \int_{D} y_n(t)  v_j(t) dt
#' }
#' @export
#'
#' @examples
#' data(simulation)
#' data(SobolevKernel)
#' summary(T_domain)
#' M_integ <- length(T_domain)/diff(range(T_domain))
#' projection_basis(Y_full, eigenvect, M_integ ) # projection on the J dimensional
#' # basis of the Y functions.


## ?@A: Again in definition, we should have Number of points of the domain D of the functions (m) "minus one" divided by the width of the domain |D|.

projection_basis <- function(y, eigenvect, M_integ)
{
    N <- dim(y)[1]
    y_mat_inprods <- matrix(NA, dim(eigenvect)[2],N)
    for (i in 1:N)
    {
        for(j in 1:dim(eigenvect)[2])
        {
            y_mat_inprods[j,i] <- sum(eigenvect[,j]*y[i,])/M_integ
        }
    }
    return(y_mat_inprods)
}
