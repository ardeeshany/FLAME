#' Compute the Covariance of a Matern process
#'
#' @param nu scalar. Value s.t. 2\code{nu} is integer. It is a smoothing
#' paramter \eqn{\nu} of the covariance: the larger \code{nu}, the smoother the process.
#'
#' @param rho scalar. Non negative range parameter \eqn{\rho} of the covariance.
#'
#' @param sigma scalar. Non negative standard deviation parameter \eqn{s} of the covariance.
#'
#' @param x vector. Vector containing the
#'  abscissa points where the Covariance is defined.
#'
#' @return \code{d} \eqn{\times}\code{d} covariance matrix.
#'
#' @export covMaterniso
#'
#'
#' @description  Definition of the Matern Covariance function, if \eqn{2 \nu} is integer. For
#' \eqn{\nu =} \code{0.5} the function is also known as the exponential covariance  or the
#' Ornstein-Uhlenbeck covariance in one dimension. In general the explicit formula
#' for the Matern Covariance for one dimensional processes, with \eqn{\nu} s.t. \eqn{\nu = p + 0.5} and \eqn{p} integer is:
#' \deqn{
#' k(d) = s^2  \frac{p!}{(2p)!} exp\left(-\frac{\sqrt{2 \nu} d}{\rho}\right) \sum_{i = 0}^p \frac{(p+i)!}{(p-i)! i!} \left( \frac{\sqrt{8 \nu} d}{\rho} \right)^{p-i}
#' }
#' with \eqn{d} the distance between two points \eqn{x} and \eqn{y} in R.
#'
#' @examples
#' # Defintion of a Gaussian process
#' # with Matern Covariance
#'
#' # Time domain of the Gaussian Process
#' M <- 50
#' T_domain <- seq(0, 1, length = M)
#'
#' # paramteters of the Matern Covariance
#' nu_alpha <- 2.5
#' range <- 1/4
#' variance <- 1
#'
#' # mean of the process
#' mu_alpha <- rep(0,M)
#'
#' # covariance structure
#' Sig_alpha <- covMaterniso(nu_alpha, rho = range, sigma = sqrt(variance), T_domain)
#'
#' # definition of the process
#'
#' # alpha <- mvrnorm(mu=mu_alpha, Sigma=Sig_alpha, n=1) # if MASS is inslalled
#'

covMaterniso <- function(nu, rho, sigma, x)
{
    if((2*nu)%%1!=0)
    {
        stop ('nu is not s.t. 2*nu is integer. You need a more general definition of the Covariance Function')
    }

    if (rho <= 0)
    {
        stop('rho must be a positive parameter')
    }

    if (sigma <= 0)
    {
        stop('sigma must be a positive parameter')
    }

    mat <- sapply(x, function(v){sapply(x, function(y){k_dist(abs(v-y), nu, rho, sigma)})})
    return(mat)
}

k_dist <- function(d, nu, rho, sigma)
{
    p = nu - 0.5
    SumToP <- 0
    for (i in 0:p)
    {
        SumToP <- SumToP + (factorial(p+i)/(factorial(p-i) * factorial(i))) * (sqrt(8*nu) * d / rho)^(p-i)
    }
    K = SumToP * sigma^2 * exp(-sqrt(2*nu)*d/rho) * factorial(p) / factorial(2*p)
    return(K)
}
