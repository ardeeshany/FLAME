#' Computation of the K-norm of a set of functions
#'
#' It computes the norm associated to the kernel K of a set of functions.
#' The functions are represented as their projection on a
#' kernel basis.
#'
#' @param M matrix. \code{J} \eqn{\times} \code{N} matrix containing in the column
#' \eqn{n} the coefficients of the projection of a function \eqn{y_n} on the
#' \code{J} eigenfunctions of the kernel.
#' @param eigenval vector. \code{J}-length vector of the eigenvalues of the
#' kernel.
#'
#' @return vector of length \code{N} containing the K-norm of the functions
#' @export
#'
#' @examples
#' data(SobolevKernel)
#' data(simulation)
#' norm_K = sum(norm_matrix_K(Y_matrix, eigenval))
#'
norm_matrix_K <- function(M, eigenval)
{
    if (is.null(dim(M)))
    {
        M <- matrix(M, length(M),1)
    }
    norm <- apply(M,2, function(x){sqrt(sum(x^2/eigenval))} )
    return(norm)
}
