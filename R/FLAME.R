#' Global FLAME method: from the definition of the kernel to the estimation
#'
#' It computes  FLAME  for the Function-on-Scalar
#' regression problem. From a set of functions -stored in a matrix or in
#' an \link[=fda]{fd} object- and a set of predictors, FLAME identifies the
#' set of meaningful predictors and their smooth representation in the
#' kernel space.
#'
#' @param Y fd object or list. \code{N} functional responses of the
#' Function-on-Scalar regression probelm. It can be an \link[=fda]{fd} object or
#' a list with 2 elements: \code{data} and \code{time_domain}. \code{time_domain}
#' is the \code{m}-length vector of the absicissa grid of the functions and \code{data}
#' is the \code{N} \eqn{\times} \code{m} matrix of the point-wise evaluation of
#' the response functions.
#' @param X matrix. \code{N} \eqn{\times} \code{I} design matrix. It has
#' standardized columns.
#' @param type_kernel string. Four possible choices are implemented.  \code{gaussian},
#' \code{exponential}, \code{sobolev} or \code{periodic}. For other
#' kernels, define manually the eigenfunctions and eigenvectors and then use
#' the \link[=estimation_beta]{estimation_beta} function. Defualt is \code{sobolev}.
#' @param param_kernel scalar. Value of the characteristic smoothing parameter of the kernel.
#' It is the \eqn{\sigma}
#' parameter of the Gaussian and the Exponential kernel, as introduced
#'  in \link[kernlab]{rbfdot} and \link[kernlab]{laplacedot} functions;
#'  the \eqn{\sigma} parameter of  the Sobolev
#' kernel as in the \link[=sobolev_kernel_generation]{sobolev_kernel_generation} function or
#' the \eqn{\sigma} paramter of the periodic kernel of the
#' \link[=generation_kernel_periodic]{generation_kernel_periodic} function. Default is
#' \code{8}.
#' @param thres_eigen scalar. Threshold to identify the significant
#' eigenvalues of the kernel. The number of significant eigennvalues \code{J} is the
#' minimum \eqn{J} s.t.
#' \deqn{
#' \sum_{j = 1}^{J} \theta_j \geq \textrm{thres\_eigen} \sum_{j = 1}^{\infty} \theta_j.
#' }
#' Default is \code{0.99}.
#' @param period_kernel scalar. Period of the kernel. In case
#'  of \code{type_kernel = "periodic"}, it is a mandatory parameter with no default.
#' If other types of kernel are chosen, it is ignored.
#' @param NoI scalar. integer, maximum number of iterations in the
#' Coordinate-Descent loop. Default is \code{10}.
#' @param thres_CD scalar. tolerance in the increment of the K-norm of the estimation
#' to stop the Coordinate-Descent loop.
#' Default is \code{0.01}
#' @param number_non_zeros scalar. integer,
#' threshold on the number of non zeros parameters to be detected. It is the
#' kill switch parameter. See the Vignette for further details. Default is \code{NULL}
#' meaning that no kill switch paramter is imposed.
#' @param ratio_lambda scalar. ratio to compute the minimum value of lambda. The
#' maximum \eqn{\lambda} (\eqn{\lambda_{\textrm{max}}}) is computed as the minimum value which makes all the coefficients
#' equal to zero. And the minimum is the product \code{ratio_lambda}\eqn{\times \lambda_{\max}}. Default is \code{0.01}.
#' @param number_lambda scalar. integer, length of the grid for the
#' \eqn{\lambda} parameter. Default is \code{100}.
#' @param proportion_training_set scalar. value in (0,1), the
#' proportion for the training set for the Cross Validation.
#' Defualt is \code{0.75}.
#' @param verbose  bool. If \code{TRUE} the progression of the algorithm in the
#' adaptive and non adaptive step is shown. If \code{FALSE} no output is shown.
#' Default is \code{FALSE}.
#'
#' @return  list containing:
#'  \itemize{
#'  \item \code{beta} fd object or list. \code{I} functional coefficients
#'  estimated by FLAME. If \code{Y} is an \link[=fda]{fd} object. then
#'  also \code{beta} is  an \link[=fda]{fd} object, and \code{beta} is
#'  the projection of the estimations on a 10 elements cubic Bspline basis.
#'  If \code{Y} is a list, then also \code{beta} is a
#'  list with 2 elements: \code{data} and \code{time_domain}. \code{time_domain}
#' is the \code{m}-length domain grid, while \code{data} is \code{I} \eqn{\times} \code{m} matrix of the point-wise evaluation of
#' the estimated coefficients.
#'  \item \code{predictors} vector of the indices of
#'  the non-zero estimated predictors.}
#'
#' @export
#'
#' @import fda
#'
#' @examples
#' \dontrun{
#' data(simulation)
#' data(SobolevKernel)
#' time <- proc.time()
#' FLAME_estimation <- FLAME()
#' duration <- proc.time()-time
#' duration
#' names(FLAME_estimation)
#' }
#'
FLAME <- function(Y, X, type_kernel = 'sobolev', param_kernel = 8,
                  thres_eigen = 0.99, period_kernel = NULL,
                  NoI = 10, thres_CD = 0.01,
                  number_non_zeros = NULL, ratio_lambda = 0.01,
                  number_lambda = 100, proportion_training_set = 0.75,
                  verbose = FALSE)
{
    # the function allows the user to pass as a Y a fd object or a list
    # containing the time domain (time_grid) and the punctual evaluation
    # of data (as a matrix)
    if (class(Y) =='fd')
    {
        T_domain = sort(c(Y$basis$rangeval, Y$basis$params))
        Y_full <- t(eval.fd(T_domain, Y, Lfdobj=0, returnMatrix=TRUE))
    }else
    {
        if (sum(names(Y) == c('time_domain')) +
            sum(names(Y) == c('data')) == 2)
        {
            T_domain <- Y$time_domain
            Y_full <- Y$data
            if (dim(Y_full)[2] != length(T_domain))
            {
                if (dim(Y_full)[1] != length(T_domain))
                {
                    warning('The data matrix should be Nobs x Tdomain, so it is transposed.')
                    Y_full <- t(Y_full)
                }else
                {
                    stop('Number of columns different from the length of the domain.' )
                }
            }
        }
    }

    M_integ <- length(T_domain)/diff(range(T_domain))

    if (type_kernel == 'periodic')
    {
        if (is.null(period_kernel))
        {
            stop('periodic kernel chosen, but no period provided')
        }

        kernel_here <- generation_kernel_periodic(period = period_kernel,
                                                  parameter = param_kernel,
                                                  domain = T_domain,
                                                  thres = thres_eigen,
                                                  return.derivatives = TRUE)
        eigenval <- kernel_here$eigenval
        eigenvect <- -kernel_here$eigenvect
        derivatives <- kernel_here$derivatives
    }else
    {
        if ((type_kernel != "exponential") & (type_kernel != "sobolev") & (type_kernel != "gaussian") )
        {
            stop('provide a valid value for the type of kernel')
        }
        kernel_here <- generation_kernel(type = type_kernel,
                                         parameter = param_kernel,
                                         domain = T_domain,
                                         thres = thres_eigen,
                                         return.derivatives = TRUE)
        eigenval <- kernel_here$eigenval
        eigenvect <- kernel_here$eigenvect
        derivatives <- kernel_here$derivatives
    }

    # project on the kernel basis

    Y_matrix <- projection_basis(Y_full, eigenvect, M_integ)
    # and run the estimation
    # check dimenstions
    if (dim(X)[1] != dim(Y_matrix)[2])
    {
        stop('Number of rows of X doesn\'t coincide with the number of individuals')
    }

    if (verbose)
        {
        cat(paste('Estimation: identification of the contribution of the ',
              dim(X)[2], ' predictors on the ', dim(Y_matrix)[2], ' functions', sep =''))}

    if (is.null(number_non_zeros))
    {
        number_non_zeros = dim(X)[2]
    }
    FLAME_est <- estimation_beta(X = X, # design matrix
                             Y = Y_matrix, # response functions projected on the kernel basis
                             eigenval = eigenval, # basis
                             NoI = NoI, # max. num. iterations coordinate descent
                             thres = thres_CD, # thres for the CV
                             number_non_zeros = number_non_zeros, # kill switch parameter
                             ratio_lambda = ratio_lambda, # ratio for the min. lambda
                             number_lambda = number_lambda, # num. elements of the grid for lambda
                             proportion_training_set = proportion_training_set, # training set
                             verbose = FALSE) # no show of all the iterations

    beta_on_time_grid <- projection_domain(FLAME_est$beta, eigenvect)

    if (class(Y) == 'fd')
    {

        basis_beta <- create.bspline.basis(norder=4,
                                      rangeval = range(T_domain),
                                      nbasis = 10)

        beta <- Data2fd(T_domain, t(beta_on_time_grid),
                        basisobj = basis_beta, nderiv = 0)
    }else
    {
        beta <- vector('list', 2)
        beta[[1]] <- T_domain
        beta[[2]] <- beta_on_time_grid
        names(beta) <- c('time_domain', 'data')
    }

    return(list(beta = beta, predictors = FLAME_est$predictors))

}
