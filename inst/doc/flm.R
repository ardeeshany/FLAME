## ----global_options, include = FALSE-------------------------------------
knitr::opts_chunk$set(fig.widh = 12, fig.height = 8, fig.path='Figs/', echo = FALSE,
                      warning = FALSE, message = FALSE)

## ----load_packages, echo = FALSE-----------------------------------------
library(flm)
require("MASS")

## ----defintion_kernel,  echo=TRUE, eval = TRUE, fig.height= 6, fig.width=6, fig.pos='c', fig.cap = '\\label{eigenvalues_sobolev} Plot of the first 4 eigenfunctions (left panel) and derivatives (right panel) of the Sobolev kernel with parameter $\\sigma = 8$. The correspondent explained variance is on the top of the plot.'----
type_kernel <-  'sobolev'
param_kernel <-  8
M <- 50
T_domain <-  seq(0, 1, length = M) # time point grid.
thres <- 0.99 # thresold for the eigenvalues. 
kernel_here <- generation_kernel(type = type_kernel,
                                 parameter = param_kernel,
                                 domain = T_domain,
                                 thres = 0.99,
                                 return.derivatives = TRUE)
eigenval <- kernel_here$eigenval
eigenvect <- kernel_here$eigenvect
derivatives <- kernel_here$derivatives

layout(mat = matrix(c(1,2,1,3),  nrow = 2, ncol = 2),
       heights= c(0.6,3.5),respect = TRUE, widths = c(4,4))
# legend
par(oma=c(0, 0, 0, 0), mar =c(0,0,2,0))
plot(c(1,2,3),c(0,1,1),col=0, xlab ='', ylab='', axes = FALSE)
legend('center',  legend = round(eigenval[1:4]/sum(eigenval), 2), 
       xpd = TRUE, horiz = TRUE,
       inset = c(0,  0.98), bty = "n", pch = '-', pt.cex = 2,
       col = rainbow(4), cex = 1)
title( main = expression(paste('Sobolev kernel, ', sigma, ' = 8')),
       cex.main = 1.2, font.main = 2)

# plot of the eigenfunctions
par(mar=c(4,4, 1, 2) + 0.1)
matplot(T_domain, eigenvect[,1:4], type = 'l', lwd = 2, xlab = 'time', 
        ylab ='eigenfunctions', 
        lty = 1, cex = 1, cex.main = 1, cex.lab = 1, 
        cex.axis = 1, col =rainbow(4), ylim = c(-2,2))

# plot of the derivatives
par(mar=c(4,4,1, 2) + 0.1)
matplot(T_domain[-50], derivatives[,1:4], type = 'l', lwd = 2, xlab = 'time', 
        ylab ='derivatives', 
        lty = 1, cex = 1, cex.main = 1, cex.lab = 1, 
        cex.axis = 1, col =rainbow(4))

## ----periodic_kernel, echo=TRUE, eval = TRUE, fig.height= 6, fig.width=6, fig.pos='c', fig.cap = '\\label{eigenvalues_periodic} Plot of the first 4 eigenfunctions (left panel) and derivatives (right panel) of the periodic kernel with period $p = 1/2$. The correspondent explained variance is on the top of the plot.'----

kernel_here <- generation_kernel_periodic(period = 1/2,
                                          parameter = param_kernel,
                                          domain = T_domain,
                                          thres = 1-10^{-16},
                                          return.derivatives = TRUE)
eigenval <- kernel_here$eigenval
eigenvect <- -kernel_here$eigenvect
derivatives <- kernel_here$derivatives

layout(mat = matrix(c(1,2,1,3),  nrow = 2, ncol = 2),
       heights= c(0.6,3.5),respect = TRUE, widths = c(4,4))
# legend
par(oma=c(0, 0, 0, 0), mar =c(0,0,2,0))
plot(c(1,2,3),c(0,1,1),col=0, xlab ='', ylab='', axes = FALSE)
legend('center',  legend = round(eigenval[1:4]/sum(eigenval), 3), 
       xpd = TRUE, horiz = TRUE,
       inset = c(0,  0.98), bty = "n", pch = '-', pt.cex = 2,
       col = rainbow(4), cex = 1)
title( main = expression(paste('Periodic kernel, period = 1/2')),
       cex.main = 1.2, font.main = 2)

# plot of the eigenfunctions
par(mar=c(4,4, 1, 2) + 0.1)
matplot(T_domain, eigenvect[,1:4], type = 'l', lwd = 2, xlab = 'time', 
        ylab ='eigenfunctions', 
        lty = 1, cex = 1, cex.main = 1, cex.lab = 1, 
        cex.axis = 1, col =rainbow(4), ylim = c(-2,2))

# plot of the derivatives
par(mar=c(4,4,1, 2) + 0.1)
matplot(T_domain[-50], derivatives[,1:4], type = 'l', lwd = 2, xlab = 'time', 
        ylab ='derivatives', 
        lty = 1, cex = 1, cex.main = 1, cex.lab = 1, 
        cex.axis = 1, col =rainbow(4))

## ----set_seed, echo = FALSE----------------------------------------------
set.seed(16589)

## ----defintion_data,  echo=TRUE, eval = TRUE, fig.height= 3, fig.width=9, fig.pos='c', fig.cap = '\\label{lin_mod} Random generation of data. From the left, 10  coefficients $\\beta^{*}(t)$, 20 random errors $\\varepsilon(t)$ and the correspondent 20 response functions $Y_n(t)$.'----
N <-  500
I <- 1000
I0 <- 10

# defintion of the time domain
m <- 50 # total number of points
T_domain <-  seq(0, 1, length = m) # time points, length = m
M_integ <- length(T_domain)/diff(range(T_domain)) # coefficient for the
# computation of the integrals

# defintion of the design matrix X, in this specific case the
# covariance matrix C is the identity matrix
mu_x <- rep(0, I)
C <- diag(I)
X <-  mvrnorm(n=N, mu=mu_x, Sigma=C)
X <- scale(X) # normalization

# defintion of the coefficients 
nu_beta <- 2.5 
range <-1/4 
variance <- 1 
hyp <- c(log(range), log(variance)/2) # set of parameters for the
# Matern Covariance operator of beta
mu_beta <- rep(0,m) # mean of the beta
Sig_beta <- covMaterniso(nu_beta, rho = range, sigma = sqrt(variance) , T_domain)
beta <- mvrnorm(mu=mu_beta, Sigma=Sig_beta, n=I0) # generation of the
# I0 significant coefficients

# defintion of the random errors
nu_eps <- 1.5
mu_eps <- rep(0, m)
Sig_eps <- covMaterniso(nu_eps, rho = range, sigma = sqrt(variance), T_domain)
eps <- mvrnorm(mu=mu_eps, Sigma=Sig_eps, n=N) # generation of the N
# random errors

I_X <- sort(sample(1:I, I0)) # index of the I0 significant predictors

Y_true <- X[,I_X]%*%beta
Y_full <- X[,I_X]%*%beta + eps # Y_n observations

par(mfrow = c(1,3), mar=c(4,4,3,2) + 0.1)

# plot of the true beta coefficients
matplot(T_domain, t(beta), type = 'l', lwd = 2, xlab = 'time', 
        ylab ='beta', 
        lty = 1, cex = 2,  cex.lab = 1.5,
        cex.axis = 1.5, col =rainbow(10))
title(main = expression(paste(beta, '*',(t),' coefficients', sep ='')),
      cex.main = 1.5, font.main = 2)

# plot of the errors
matplot(T_domain, t(eps[1:20,]), type = 'l', lwd = 2, xlab = 'time', 
        ylab ='beta', 
        lty = 1, cex = 2,  cex.lab = 1.5, 
        cex.axis = 1.5, col =rainbow(20))
title(main = expression(paste(epsilon(t),' coefficients', sep ='')),
      cex.main = 1.5, font.main = 2)

# plot of the response functions
matplot(T_domain, t(Y_full[1:20,]), type = 'l', lwd = 2, xlab = 'time', 
        ylab ='beta', 
        lty = 1, cex = 2, cex.lab = 1.5, 
        cex.axis = 1.5, col =rainbow(20))
title(main = expression(paste(Y(t),' functions', sep ='')),
      cex.main = 1.5, font.main = 2)

## ----kernel_and_projection, echo = TRUE----------------------------------

# defintion of the kernel
type_kernel <-  'sobolev'
param_kernel <-  8
m <- 50
T_domain <-  seq(0, 1, length = m) # time point grid.
thres <- 0.99 # thresold for the eigenvalues. 
kernel_here <- generation_kernel(type = type_kernel,
                                 parameter = param_kernel,
                                 domain = T_domain,
                                 thres = 0.99,
                                 return.derivatives = TRUE)
eigenval <- kernel_here$eigenval
eigenvect <- kernel_here$eigenvect
derivatives <- kernel_here$derivatives

# preojection on the kernel basis of y and beta
Y_matrix <- projection_basis(Y_full, eigenvect, M_integ)
B_true <- projection_basis(beta, eigenvect, M_integ)

matrix_beta_true_full <- matrix(0, dim(B_true)[1], I)
matrix_beta_true_full[,I_X] <- B_true

## ----defintion_derivatives, echo=TRUE------------------------------------
B_true_der <- t(kernel_here$derivatives %*% B_true)

Y_true_der <- X[,I_X]%*%B_true_der

## ----estimation, echo = TRUE---------------------------------------------

time <- proc.time()
FLAME <- estimation_beta(X = X, # design matrix
                         Y = Y_matrix, # response functions projected on the kernel basis
                         eigenval = eigenval, # basis
                         NoI = 10, # max. num. iterations coordinate descent
                         thres = 0.1, # stopping threshold for the coordinate descent
                         number_non_zeros = I0*2, # kill switch parameter
                         ratio_lambda = 0.01, # ratio for the min. lambda
                         number_lambda = 100, # num. elements of the grid for lambda
                         proportion_training_set = 0.75, # training set
                         verbose = FALSE) # no show of all the iterations
duration <- proc.time()-time
duration

## ----names_FLAME, echo = TRUE--------------------------------------------
names(FLAME)

## ----projection_on_the_domain, echo=TRUE, eval = TRUE, fig.height= 3.5, fig.width=6.5, fig.pos='c', fig.cap = '\\label{fig:estimation_FLAME} In the left panel the plot of the simulated $\\beta^{*}$ coefficients, while in the right panel the FLAME coefficients $\\hat{\\beta}$ are shown.'----
beta_on_time_grid <- projection_domain(FLAME$beta, eigenvect)
y_on_grid_estimated <- X %*% beta_on_time_grid

par(mfrow = c(1,2), mar=c(4,4,3,2) + 0.1)

# plot of the true beta coefficients
matplot(T_domain, t(beta), type = 'l', lwd = 2, xlab = 'time', 
        ylab ='beta', 
        lty = 1, cex = 1,  cex.lab = 1,
        cex.axis = 1, col =rainbow(10), ylim = c(-3, 2))
title(main = expression(paste(beta, '*', ' coefficients', sep ='')),
      cex.main = 1.2, font.main = 2)

# plot of the estimated beta
matplot(T_domain, t(beta_on_time_grid[FLAME$predictors,]), type = 'l', 
        lwd = 2, xlab = 'time', 
        ylab ='beta', 
        lty = 1, cex = 1,  cex.lab = 1, 
        cex.axis = 1, col =rainbow(10), ylim = c(-3, 2))
title(main = expression(paste('estimated ',beta,' coefficients', sep ='')),
      cex.main = 1.2, font.main = 2)


## ----TP, echo = TRUE-----------------------------------------------------

I_X
FLAME$predictors

true_positives <- length(which(I_X %in% FLAME$predictors))
true_positives # number of significant predictors correctly identified
false_positives <- length(which(!(FLAME$predictors %in% I_X)))
false_positives # number of non significant predictors wrongly picked by the algorithm

## ----analysis_results, echo = TRUE---------------------------------------
beta_der_on_grid_estimated<- kernel_here$derivatives %*% FLAME$beta

prediction_error <- sum(apply(Y_true - y_on_grid_estimated,
                              1,
                              function(x)
                              {
                               sqrt((2*sum(x^2)-x[1]^2-x[length(x)]^2)/(M_integ*2))
                              }
                              )
                        )
prediction_error

estimated_y_der_grid <- X %*% t(beta_der_on_grid_estimated)
prediction_error_der <- sum(apply(Y_true_der - estimated_y_der_grid,
                                  1,
                                  function(x)
                                  {
                                   sqrt((2*sum(x^2)-x[1]^2-x[length(x)]^2)/((M_integ-1)*2))
                                  }
                                  )
                            )
prediction_error_der

norm_K_beta <- sum(norm_matrix_K(matrix_beta_true_full - FLAME$beta, eigenval)^2)
norm_K_beta


## ----estimation_from_fd, echo = TRUE-------------------------------------
class(Y_fd)
estimation_auto <-  FLAME(Y_fd, # fd object for the response
                          X, # predictors matrix
                          number_non_zeros = 20)
# default choice for the kernel is Sobolev with sigma = 8,
names(estimation_auto)
class(estimation_auto$beta)
estimation_auto$predictors

