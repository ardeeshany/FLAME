#' Definition of the kernel for the Sobolev Space
#'
#' @description Definition of the eigenfunctions and eigenvalues of the
#' kernel for the Sobolev Space \eqn{H^1((a,b))}.
#'
#' @param a scalar. left end point of the domain.
#' @param b scalar. right end point of the domain.
#' @param m scalar, integer. number of points of the interval domain.
#' @param sigma scalar. weight \eqn{\sigma} associated to the derivative in
#'  the norm associated to the kernel. See Details
#' for the explicit definition of the norm.
#' @param plot.eigen bool. if \code{TRUE} the cumulative sum of the
#' eigenvalues of the kernel is plotted. Default is \code{FALSE}.
#'
#' @return list containing
#' \itemize{
#' \item \code{vectors} matrix. \code{m} \eqn{\times} \code{m} matrix containing
#' the eigenvectors of the kernel. Each column contains the evaluation of an
#' eigenfunction on the domain \code{seq(a, b, length = m)}.
#' \item \code{values} vector. \code{m} length vector containing the eigenvalues
#' of the kernel
#' }
#'
#' @details The norm associated to the Sobolev kernel, dependent on the
#' smoothing parameter \eqn{\sigma} is
#' \deqn{
#' \|f\|^2 = \|f\|^2_{L^2} + 1/\sigma \|f^{\prime}\|^2_{L^2}
#' }
#' The function \link[=sobolev_kernel_generation]{sobolev_kernel_generation}
#' is implicitly called in the \link[=generation_kernel]{generation_kernel}
#' function when \code{type} parameter is \code{'sobolev'}. See the Vignette for
#' the explicit definition of the kernel.
#' @export
#'
#' @examples
#' sobolev_kernel_generation(a = 0, b = 1, m = 100,
#'                sigma = 1, plot.eigen = FALSE)
#'


## @A: sobolev_kernel_generation generates the eigenvalues and eigenvectors of the Sigma matrix not the sobolev kernel.
## sobolev kernel has the same eigenvectors but its elgenval is sobolev_kernel_generation$values/m
## more detail: when we want to find the eigenval/eigenvec of a kernel function, say K(t,s)
## we should find a \lambda and v such that \int{ K(t,s) v(s) ds} = \lambda v(t)
## numerically we can write the integral as \sum{j=1}^m K(t,s_j)v(s_j)*1/M = \lambda v(t) for any t
## or even by matrix [K(t_i_s_j)][v(s_1),...,v(s_m)]*1/m=\lambda [v(t_1),...,v(t_m)]
## Since by linear algebra we can find the eigenval/eigenvec of the matrix [K(t_i_s_j)], called \lambda_new and v_new
## we can see v=v_new but \lambda=\lambda_new*1/m

sobolev_kernel_generation <- function( a, b, m,
                           sigma, plot.eigen = FALSE )
{

    ## ?@A: shouldn't it be return(sigma * cosh(sigma*(b-s))*cosh(sigma*(t-a))/sinh(sigma*(sigma)))?

        K1<-function(t,s){
            if(t < s)
                return(sigma * cosh(sigma*(b-s))*cosh(sigma*(t-a))/sinh(sigma*(b-a)))
            else
                return(sigma * cosh(sigma*(b-t))*cosh(sigma*(s-a))/sinh(sigma*(b-a)))
        }




       ## @A: pts is the grid on t and s
        pts<-seq(0,1,length=m)
        Sigma<-matrix(nrow=m,ncol=m)
        for(i in 1:m){
            for(j in i:m){
                Sigma[i,j] = K1(pts[i],pts[j])
                Sigma[j,i] = Sigma[i,j]
            }
        }

        E_Sig<-eigen(Sigma)
        VAL <- E_Sig$values
        if (plot.eigen)
            plot(cumsum(VAL)/sum(VAL))
        V<-E_Sig$vectors

        return(list(vectors = V, values = VAL )) # per each column an eigenvector.


}
