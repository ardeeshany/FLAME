#' Computation of the H-norm of the derivatives of a set of functions
#'
#' It computes the H norm of the derivatives of a set of functions
#' represented as their projection on a
#' basis of the space H.
#'
#' @param M matrix. \code{J} \eqn{\times} \code{N} matrix containing in the column
#' \eqn{n} the coefficients of the projection of a function \eqn{y_n} on the
#' \code{J}-dimensional basis functions.
#' @param pairwise_derivatives matrix. \code{J} \eqn{\times} \code{J} matrix
#' containing  the integral of the pairwise product of the basis functions. It can
#' be computed with the \link[=compute_pairwise_integrals]{compute_pairwise_integrals}
#' function.
#'
#' @return vector of length \code{N} containing the H-norm of the derivatives of
#' the functions
#' @export
#'
#' @examples
#' data(SobolevKernel)
#' T_domain <- seq(0, 1, length = 50)
#' pairwise_derivartives <- compute_pairwise_integrals(derivatives, T_domain)
#' sum(norm_matrix_H_derivatives(t(Y_matrix), pairwise_derivartives))

norm_matrix_H_derivatives <- function(M, pairwise_derivatives)
{
    if (is.null(dim(M)))
    {
        M <- matrix(M, length(M),1)
    }
    norm <- sqrt(diag(M %*% pairwise_derivatives %*% t(M)))
    return(norm)
}
