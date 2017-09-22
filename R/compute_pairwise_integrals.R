#' Integral of the pairwise product of functions
#'
#' Given a set of functions, represented as their pointwise evaluation on a time grid
#' and stored in the columns of a matrix, it computes the integral of the pairwise
#' product of the functions
#'
#' @param matrix matrix. \code{m} \eqn{\times} \code{J} matrix containing in the column
#' \eqn{j} the pointwise representation of the function \eqn{f_j} on the \code{T_grid} domain.
#' @param T_grid vector. \code{m}-length vector of the time domain
#'
#' @return matrix \code{J} \eqn{\times} \code{J} containing in each element \eqn{(i,j)}
#' the integral of \eqn{f_i \cdot f_j} on the domain defined by \code{T_grid}.
#'
#' @export
#' @examples
#' data(SobolevKernel)
#' T_domain <- seq(0, 1, length = 50)
#' pairwise_derivartives <- compute_pairwise_integrals(derivatives, T_domain)
#'

## ?@A: Because of T_grid[2]-T_grid[1] as M_integ, It works just for equispace T_grid, isn't it?
compute_pairwise_integrals <- function(matrix, T_grid)
{
    matrix_derivatives <- apply(matrix, 2, function(y)
    {
        t(as.vector(apply(matrix,2,
                          function(x)
                          {
                              sum(x*y*(T_grid[2]-T_grid[1]))
                          }
        )
        ))
    }
    )

    return(matrix_derivatives)
}


