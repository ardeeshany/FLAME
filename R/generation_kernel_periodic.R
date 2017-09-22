#' Definition of the eigenfunctions and eigenvalues of the periodic kernel
#'
#' @description Given the period \code{p} and the smoothing parameter \eqn{\sigma}, it returns the
#' evaluation of the eigenfunctions of the kernel on the grid \code{domain}
#' and the correspondent eigenvalues.
#'
#' @param period scalar. Period of the kernel.
#' @param parameter scalar. Parameter to tune the smoothness level of the kernel.
#' The \eqn{\sigma} parameter introduced in Details.
#' @param domain vector. \code{m}-length vector for the abscissa grid
#' of the kernel.
#' @param thres scalar. Threshold for the identification of the significant
#' eigenvalues of the kernel. The number of significant eigennvalues \code{J} is
#' the minimum \eqn{J} s.t.
#' \deqn{
#' \sum_{j = 1}^J \theta_j \geq \textrm{thres} \sum_{j = 1}^{\infty} \theta_j.
#' } Default is 0.99.
#' @param return.derivatives bool. If \code{TRUE} the function returns the matrix of
#' the derivatives of the selected eigenfunctions evaluated on the time domain.
#'  Default is \code{FALSE}.
#'
#' @return list containing
#'  \itemize{
#'  \item \code{eigenvect} \code{m} \eqn{\times} \code{J} matrix of
#'  the eigenfunctions of the kernel evaluated on the \code{domain}.
#'  \item \code{eigenval} \code{J}-length  vector of the
#'  eigenvalues of the kernel
#'  \item \code{derivatives}. if \code{return.derivatives = TRUE}.
#'  \code{m-1} \eqn{\times} \code{J} matrix of the derivatives of
#'  the eigenfunctions.
#'}
#'
#' @details The periodic kernel of period \eqn{p} defined in this function is
#' \deqn{
#' K(x, y) = \sigma^2 \exp \left\{  -2/\sigma \sin^2\left(\frac{\pi | x - y|}{p} \right)\right\}.
#' }
#' where \eqn{\sigma} is the smoothing parameter.
#' @export
#'
#' @examples
#' param_kernel <- 8
#' T_domain <- seq(0, 1, length = 50)
#' kernel_here <- generation_kernel_periodic(period = 1/2,
#'                                           parameter = param_kernel,
#'                                           domain = T_domain,
#'                                           thres = 1-10^{-16},
#'                                           return.derivatives = TRUE)
#' names(kernel_here)
#'
generation_kernel_periodic <- function(period = NULL, parameter = NULL, domain, thres = 0.99,
                              return.derivatives = FALSE)
{

    ## ?@A: Shouldn't it be M_integ <- (length(domain)-1)/diff(range(domain))?
    M_integ <- length(domain)/diff(range(domain))

    if (length(period) != 1 )
    {
        stop ("please provide the period of the kernel.") # period is as proportion of the domain
    }

    if (length(parameter) != 1 )
    {
        stop ("please provide the parameters of the kernel. See the help page
              of the function for the parameter definition")
    }


    dist_periodic <- function(x, y, sigma, p)
    {
        sigma^2* exp(- (2*sin(pi*abs(x-y)/(p))^2)/(sigma))
    }

    ## @A: difference between apply, lapply and sapply:
    ## apply  - When you want to apply a function to the rows or columns of a matrix
    ## lapply - When you want to apply a function to each element of a list in turn and get a list back.
    ## sapply - When you want to apply a function to each element of a list in turn, but you want a vector back, rather than a list.

    generate_matrix <- function(domain, sigma, p){
    ## @A: it means we apply the function(x) on each element of "domain" but the output is a vector.
        sapply(domain, function(x){
            sapply(domain, function(y) {dist_periodic(x,y, sigma,p)})
        })
    }


   ## @A: the following is another way to define generate_matrix but takes more time
   ## generate_matrix=function(domain,sigma,p){
   #      G2=matrix(NA,length(domain),length(domain))
   #          for(i in 1:length(domain)){
   #          for(j in 1:length(domain)){
   #            G2[i,j]=  dist_periodic(domain[i],domain[j], sigma,p)
   #          }
   #      }
   # return(G2)
   #       }


    kernel_matrix <- generate_matrix(domain, sigma = parameter, p = period)

    kernel_def <- eigen(kernel_matrix)

    ## ?@A: Again, shouldn't it be kernel_def$vectors <- kernel_def$vectors
    kernel_def$values <- kernel_def$values/M_integ
    kernel_def$vectors <- kernel_def$vectors*sqrt(M_integ)

    # isolate the number of significant eigenvalues
    num_eigen <- which((cumsum(kernel_def$values)/sum(kernel_def$values))>thres)[1]
    eigen_chosen <- 1:num_eigen

    # defintion of the eigenvalues and eigenvectros of the kernel

    autoval <- kernel_def$values[eigen_chosen]
    autovett <- kernel_def$vectors[,eigen_chosen]


    if (return.derivatives == TRUE)
    {
        diff_autovett <- apply(autovett, 2, diff)/(domain[2]-domain[1])
        return(list(eigenvect = autovett, eigenval = autoval,
                    derivatives = diff_autovett))
    }else{
        return(list(eigenvect = autovett, eigenval = autoval))
    }

    }
