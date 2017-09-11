#' Computation of the H-norm of a set of functions
#'
#' It computes the H norm of a set of functions represented as their projection on the
#' basis of the space H.
#'
#' @param M matrix. \code{J} \eqn{\times} \code{N} matrix containing in the column
#' \eqn{n} the coefficients of the projection of a function \eqn{y_n} on the
#' \code{J}-dimensional basis.
#'
#' @return vector of length \code{N} containing the H-norm of the functions
#' @export
#'
#' @examples
#' data(simulation)
#' length(norm_matrix_H(Y_matrix))
#' sum(norm_matrix_H(Y_matrix))
#'
norm_matrix_H <- function(M)
{
    if (is.null(dim(M)))
    {
        M <- matrix(M, length(M),1)
    }
    norm <- apply(M,1, function(x){sqrt(sum(x^2))} )
    return(norm)
}
