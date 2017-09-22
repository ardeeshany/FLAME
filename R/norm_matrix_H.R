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

## @A: For each column, it does not give the norm H of the column! It actually compute the norm H of
## the function whose projection on basis is the column.

norm_matrix_H <- function(M)
{
    if (is.null(dim(M)))
    {
        M <- matrix(M, length(M),1)
    }

    ## ?@A: Shouldn't it be norm <- apply(M,2, function(x){sqrt(sum(x^2))})? or I guess, the M should be a N \times J matrix, not J*N !
    norm <- apply(M,1, function(x){sqrt(sum(x^2))} )
    return(norm)
}
