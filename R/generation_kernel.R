#' Definition of the eigenfunctions and eigenvalues of the kernel
#'
#' Given a particular type of kernel, to be chosen among (\code{gaussian},
#' \code{exponential} and \code{sobolev}), it returns the
#' evaluation of the eigenfunctions of the kernel on the grid \code{domain}
#' and the correspondent eigenvalues.
#'
#' @param type string. Type of kernel. Three possible choices implemented:
#' \code{gaussian},  \code{exponential} and \code{sobolev}. For the other
#' types of kernel, define manually the eigenfunctions and eigenvectors
#' of the kernel.
#' @param parameter scalar. Value of the characteristic parameter of the kernel.
#' It is the \eqn{\sigma}
#' parameter of the Gaussian and the Exponential kernel, as introduced in \code{kernlab};
#' and the \eqn{\sigma} parameter of  the Sobolev
#' kernel as in the \link[=sobolev_kernel_generation]{sobolev_kernel_generation} function.
#'
#' @param domain  vector. \code{m}-length vector of the
#' abscissa grid  for the definition of the kernel.
#' @param thres scalar. Threshold to identify the significant
#' eigenvalues of the kernel.  The number of significant eigennvalues \code{J} is
#' the minimum \eqn{J} s.t.
#' \deqn{
#' \sum_{j = 1}^J \theta_j \geq \textrm{thres} \sum_{j = 1}^{\infty} \theta_j.
#' } Default is 0.99.
#'
#' @param return.derivatives bool. If \code{TRUE} the function returns also
#' the matrix of the evaluation of
#' the derivatives of the eigenfunctions on the time domain.
#'  Default is \code{FALSE}.
#'
#' @return list containing
#'  \itemize{
#'  \item \code{eigenvect} \code{m} \eqn{\times} \code{J} matrix of
#'  the eigenfunctions of the kernel evaluated on the \code{domain}.
#'  \item \code{eigenval} \code{J}-length  vector of the
#'  eigenvalues of the kernel
#'  \item \code{derivatives}. if \code{return.derivatives = TRUE}.
#'  \code{derivatives} is the (\code{m-1}) \eqn{\times} \code{J} matrix of the derivatives of
#'  the eigenfunctions evaluated on the time domain.
#'}
#' @details Here the list of the kernel defined in this function
#' \itemize{
#' \item \code{gaussian}
#' \deqn{
#' k(x, x') = \exp(-\sigma \| x- x'\|^2)
#' }
#' \item \code{exponential}
#' \deqn{
#' k(x, x') = \exp(-\sigma \| x- x'\|)
#' }
#' \item \code{sobolev}, the kernel associated to the norm in the \eqn{H^1} space
#' \deqn{
#' \| f \|^2 = \int_{D} f(t)^2 dt + \frac{1}{\sigma} \int_{D} f'(t)^2 dt
#' }
#' where \eqn{D} is the one-dimensional \code{domain} and \eqn{f'} is the first derivative of the function.
#' }
#'
#' @export
#'
#' @import kernlab
#'
#' @examples
#' # definition of the kernel
#' type_kernel <-  'sobolev'
#' param_kernel <-  8
#' T_domain <- seq(0, 1, length = 50)
#' kernel_here <- generation_kernel ( type = type_kernel,
#'                                    parameter = param_kernel,
#'                                    domain = T_domain,
#'                                    thres = 0.99,
#'                                    return.derivatives = TRUE)
#'
#' eigenvalues <- kernel_here$eigenval
#' eigenvectors <- kernel_here$eigenvect
#' der <- kernel_here$derivatives
#'

generation_kernel <- function(type = 'sobolev', parameter = NULL, domain, thres = 0.99,
                              return.derivatives = FALSE)
{

    ## ?@A: shouldn't it be M_integ <- (length(domain)-1)/diff(range(domain))?

    M_integ <- length(domain)/diff(range(domain))

    if (!(type %in% c('sobolev', 'exponential', 'gaussian')))
    {
        stop ("error: not defined kernel, please define the set
              of eigenfunctions and eigenvectors manually")
    }
    if (length(parameter) != 1 )
    {
        stop ("please provide the parameters of the kernel. See the help page
              of the function for the parameter definition")
    }

    # here the defintion of the eigenfunctions and eigenvalues of the
    # three tipes of kernel.
    if (type =='sobolev')
    {
        kernel_def <- sobolev_kernel_generation(a = domain[1], b = domain[length(domain)],
                                     m = length(domain), sigma = parameter,
                                     plot.eigen = FALSE)

    ## @A:when we want to find the eigenval/eigenvec of a kernel function, say K(t,s)
    ## we should find a \lambda and v such that \int{ K(t,s) v(s) ds} = \lambda v(t)
    ## numerically we can write the integral as \sum{j=1}^m K(t,s_j)v(s_j)*1/M = \lambda v(t) for any t
    ## or even by matrix [K(t_i_s_j)][v(s_1),...,v(s_m)]*1/m=\lambda [v(t_1),...,v(t_m)]
    ## Since by linear algebra we can find the eigenval/eigenvec of the matrix [K(t_i_s_j)], called \lambda_new and v_new
    ## we can see v=v_new but \lambda=\lambda_new*1/m

    ## ?@A: Shouln't it be kernel_def$vectors <- kernel_def$vectors
        kernel_def$values <- kernel_def$values/M_integ
        kernel_def$vectors <- kernel_def$vectors*sqrt(M_integ)
    }
    if (type == 'exponential')
    {

    ## @A: in "kernlab" package, the "laplacedot" makes the exponential kernel function
        rbfkernel <- laplacedot(sigma = parameter)
        mat <- kernelMatrix(rbfkernel, domain)

        ## ?@A: Shouln't it be vectors = eigen(mat)$vectors?
        kernel_def <- list(vectors = eigen(mat)$vectors*sqrt(M_integ),
                           values = eigen(mat)$values/M_integ)
    }
    if (type == 'gaussian')
    {

    ## @A: in "kernlab" package, the "rbfdot" makes the gaussian kernel function

        rbfkernel <- rbfdot(sigma = parameter)
        mat <- kernelMatrix(rbfkernel, domain)

        ## ?@A: Shouln't it be vectors = eigen(mat)$vectors?
         kernel_def <- list(vectors = eigen(mat)$vectors*sqrt(M_integ),
                           values = eigen(mat)$values/M_integ)
    }

    # isolate the number of significant eigenvalues
    num_eigen <- which((cumsum(kernel_def$values)/sum(kernel_def$values))>thres)[1]
    eigen_chosen <- 1:num_eigen

    # defintion of the eigenvalues and eigenvectros of the kernel

    autoval <- kernel_def$values[eigen_chosen]
    autovett <- kernel_def$vectors[,eigen_chosen]

    if (return.derivatives == TRUE)
    {



    ## @A: since v'(x)=(v(x+h)-v(x))/h so we take each column and make v' for it.
    ## It is obvious when we have m points for each v, we will have m-1 points for each v'

    ## ?@A: Shouldn't it be diff_autovett <- apply(autovett, 2, diff)/(M_integ)?
       diff_autovett <- apply(autovett, 2, diff)/(domain[2]-domain[1])
        return(list(eigenvect = autovett, eigenval = autoval,
                    derivatives = diff_autovett))
    }else{
        return(list(eigenvect = autovett, eigenval = autoval))
    }

}
