#' FLAME for the Function-on-Scalar regression coefficients.
#'
#' This function computes  FLAME for the Function-on-Scalar
#' regression problem,
#' returning the significant parameters and
#' the estimated coefficients in the kernel basis.
#'
#' @param X matrix. \code{N} \eqn{\times} \code{I} design matrix. It has
#' standardized columns.
#' @param Y matrix. \code{J} \eqn{\times} \code{N} matrix containing the
#' coefficients of the projection of the \code{N} response functions on the kernel
#' basis. In particular column \eqn{n} contains the basis projection of \eqn{y_n}.
#' @param eigenval vector. \code{J}-length vector containing the eigenvalues
#' of the kernel.
#' @param NoI scalar. integer, maximum number of iterations in the
#' Coordinate-Descent loop.
#' @param thres scalar. tolerance on the K-norm of the increment of the estimation
#' to stop the Coordinate-Descent loop
#' @param number_non_zeros scalar. integer,
#' threshold on the number of non zeros parameters to be detected. It is the
#' kill switch parameter presented in the Vignette.
#' @param ratio_lambda scalar. ratio to compute the minimum value of lambda. The
#' maximum \eqn{\lambda_{\max}} is computed as the minimum value which makes all the coefficients
#' equal to zero. The minimum is the product \code{ratio_lambda}\eqn{\times \lambda_{\max}}.
#' @param number_lambda scalar. integer, length of the grid for the
#' \eqn{\lambda} parameter.
#' @param proportion_training_set scalar. value in (0,1), the
#' proportion for the training set for the Cross Validation.
#' @param verbose bool. If \code{TRUE} the progression of the algorithm in the
#' adaptive and non adaptive step is shown. If \code{FALSE} no output is shown.
#' Default is \code{FALSE}.
#'
#' @return list containing:
#'  \itemize{
#'  \item \code{beta} matrix. \code{J} \eqn{\times} \code{I} matrix of the final estimated
#'  coefficients in the kernel basis.
#'  \item \code{beta_no_adaptive}  matrix \code{J} \eqn{\times} \code{I} matrix of
#'  the coefficients estimated at the end of the non-adaptive step.
#'  \item \code{predictors} vector of the indices of the non-zero estimated predictors.
#'  \item \code{predictors_no_adaptive}  vector of the indices of the non-zeros
#'  predictor estimated after the non-adaptive step.}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(simulation)
#' data(SobolevKernel)
#' I0 <- dim(B_true)[2]
#' time <- proc.time()
#' FLAME_estimation <- estimation_beta(X = X,
#'                                  Y = Y_matrix,
#'                                  eigenval = eigenval,
#'                                  NoI = 10,
#'                                  thres = 0.1,
#'                                  number_non_zeros = I0*2,
#'                                  ratio_lambda = 0.01,
#'                                  number_lambda = 100,
#'                                  proportion_training_set = 0.75,
#'                                  verbose = FALSE)
#' duration <- proc.time()-time
#' duration
#' names(FLAME_estimation)
#' }
#'
#' @import Rcpp
#' @import RcppArmadillo
#' @import nloptr
#' @useDynLib flm, .registration = TRUE
#'

estimation_beta <- function(X, Y, eigenval,
                 NoI, thres, number_non_zeros, ## @A: These three variables are stopping criteria
                 ratio_lambda, # lambda,
                 number_lambda, proportion_training_set,
                 verbose = FALSE)
{
#
    # @param BIC_adaptive_step bool. if \code{TRUE} the \eqn{\lambda} paramter
    # of the adaptive step is choosen with the Bayesian Index Criteria (BIC); othewise
    # Cross-Validation is used.
    BIC_adaptive_step = FALSE
    # @param lambda scalar. value of \eqn{\lambda}.
    # set this paramter to \code{numeric()} to ignore it and run the
    # algorithm to identify the proper value of \eqn{\lambda}. The vector
    # is defined with \code{ratio_lambda} and \code{number_lambda} parameters.
    # The best value of \eqn{\lambda} is detected with Cross-Validation.
    # if \code{lambda} is different from \code{numeric()}
    # further parameters are ignored

    # estimation_norm_R is the function used to compute the norm of
    # beta. It is numerically solved with the NLOPTR optimization tool.


############################################################################
#
# Part 1: Defining the function for numerically estimate norm beta in K
#
############################################################################
### %% Need to modify for AFSSEN

    estimation_norm_R <- function(lambda, omega, B, tau)
    {

    ## @A: the function we want to optimize in FLAME is:
    ## 1=\sum_{j=1}^{\infty} \dfrac{tau_j < \beta_j , v_j >^2 }{ tau_j \norm{\beta_j}_K + \lambda omega_j}^2


    ## @A: nloptr is an R package for nonlinear optimization.

        ## ?@A: why do you optimiza |1-1/sum| instead of |1-sum|?

             optimization_function <- function(x, lambda, omega, B, tau)
        {
            numerat <- B^2 * tau
            denom <- (tau * x + lambda * omega)^2
            tot <- 1 - 1 / (sum(numerat/denom))
            return(abs(tot))
        }

        ## ?@A: Are these initial values are fixed or we need to change them in different situation?
        opts= list("algorithm"= "NLOPT_LN_COBYLA", "xtol_rel"=1.0e-16)


        ## ?@A: why it ended up with semicolon (;)?
        ott_model <- nloptr(0, optimization_function, opts=opts, lb=0, omega=omega,lambda=lambda, B=B, tau=tau);

        return(ott_model$solution)

    }

############################################################################
#
# Part 2: RUN non-adaptive version to estimte weights for adaptive step
#
############################################################################



#######################################################
#
# Part 2.1: estimation_first
#
#######################################################

    function_computation <- estimation_norm_R

    weights = rep(1, dim(X)[2])

    # First step: all the weights are set to 1
    # first run to identify the possible values of lambda
    # (from lambda_max -all the predictors set to 0- to the the value for which
    # number_non_zeros predicotors are different from 0)


    num_lambda_NONad <- number_lambda

    #num_lambda_NONad <- round(number_lambda/5) # for this first step (the most
    # computationally expensive) we can reduce the number of lambda examined.
    # Then (in the Second: Adaptive Step)
    # we will consider the whole set of lambdas to define the final
    # estimation

    # num_lambda_NONad <- number_lambda # no reduction

    lambda <- numeric()

    if (verbose) print('Non-adaptive step')

    ## ?@A: what is the matrix(, 1, 0) for Beta here?

    ## @A: we can import the vector of lambda in "lambda" or same as here put it numeric() (means nothing)
    ## and ask the definition_beta function generate the log equispaced lambda function with initializing ratio_lambda and number_lambda.
     estimation_first <- definition_beta(X, Y, matrix(, 1, 0),
                                        eigenval, weights, NoI,
                                        function_computation, thres,
                                        number_non_zeros, lambda,
                                        ratio_lambda, num_lambda_NONad,
                                        BIC = 0, verbose = FALSE)


     ## @A: estimation_first$Lambda_vect is a vector {lambda_max,....,lambda*} the lambdas were run in definition_beta.
     ## lambda* can be lambda_min or some lambda larger than it.

     lambda_first <- estimation_first$Lambda_vect

#######################################################
#
# Part 2.2: estimation_first_CV and finding lamda which implies minimum CV
#
#######################################################

    # Corss Valiation to identify the optimum value of lambda among
    # the set of the possible values previously
    # identified (the optimum lambda will be
    # the ones which minimizes the cross validation error)

    if (verbose) print('Validation on the test set: identification of lambda')

    subset <- c(rep(2, ceiling(dim(X)[1]*proportion_training_set)),
                rep(1, ceiling(dim(X)[1]*(1-proportion_training_set))))
    # proportion_training_set is the percentage of the data to use to
    # fit the model (2), the remaining part to compute the prediction
    # error:left_out (1)


    ## ?@A: there exists a set.seed for lable of training set!

    set.seed(16589)
    random_groups <- subset[sample(1:dim(X)[1])] # definition of the

    # training and test set

    i <- 1

    left_out <- which(random_groups==i)

    # fitting of the model with the proportion_training_set% of the data and
    # computation of the CV error
    estimation_first_CV <- definition_beta_CV(X[-left_out,], Y[, -left_out],
                                              X[left_out,], Y[, left_out],
                                              eigenval, weights, NoI,
                                              function_computation,
                                              thres, number_non_zeros,
                                              lambda_first, verbose = FALSE)

   # print(estimation_first_CV$error)

    ## @A: estimation_first_CV$error is the \norm{Y_pred - Y_test}_H for different lambda.
    ## @A: lambda_first_selected is not the optimum lambda. It is the index of the lambda.
    lambda_first_selected <- which.min(estimation_first_CV$error) # optimum
    # lambda. It minimizes the CV error

#######################################################
#
# Part 2.3: estimation_first_definite, weights_new and
#
#######################################################

    if (verbose) { print(paste("Non adaptive step: final estimation with the optimum lambda")) }

   # if (estimation_first_CV$reached_last==1)
    #{
        # if the CV error identifies the last value of lambda,
        # we take as estimation the result of the previous run
  #      estimation_first_definite <- estimation_first
 #   }
   # if (estimation_first_CV$reached_last==0){
        # otherwise we compute again the estimation with
        # to the optimum value of lambda

    ## @A: Since lambda_first=estimation_first$Lambda_vect, is a vector {lambda_max,....,lambda*} the lambdas were run in definition beta.
    ## So we import the lambda_start, we do not need to initialize the num_lambda and ratio_lambda, so we put them 0. It can be anything else.

    ## @A: here the weights are still 1.

    ## @A: we reduced the range of lambda from [lambda_max, lambda_min_by_CV]



    ## ?@A: In part 2.1 we calculated definition_beta for lambda={\lambda_max,....,\lambda_min}. In
    ## part 2.2 we find the lambda_first_selected  which is the label of tha lambda gives the min C.V.
    ## But again here calculate definition_beta for lambda={\lambda_max , ... , \lambda_{lambda_first_selected}} which
    ## the result should be a subset of results of part 2.1. Don't you think that? It does not
    ## affect on final result but makes the code slower.

        estimation_first_definite <- definition_beta(X, Y,
                                                     matrix(, 1, 0), #estimation_first_CV$Beta, #,
                                                     eigenval, weights, NoI,
                                                     function_computation, thres,
                                                     number_non_zeros,
                                                     lambda_first[1:lambda_first_selected],
                                                     0, 0, BIC = 0, verbose = FALSE)
#    }

    # adaptive step to improve the estimation of the
    # predicotors and to make a further selection.

    if(verbose) {print (paste("Adaptive Step: update of the estimation of the",
                 length(estimation_first_definite$Pred),
                 "non zeros predictors identified."))}

    # isolation of the significant predictors fitted by the Non-Adaptive step.


    ## @A: estimation_first_definite$Pred is just the label of nonzero predictors

    if (length(estimation_first_definite$Pred)==0)
    {
        print('No significant predictors indentified.')
        result <- list(beta=NULL,
                       beta_no_adaptive=NULL,
                       predictors=NULL,
                       predictors_no_adaptive=NULL)
        return(result)
    }else{
        if(length(estimation_first_definite$Pred)==1)
        {

            ## @A: In the following matrices, we have just one column because the length(estimation_first_definite$Pred)==1
             beta_selected <- matrix(estimation_first_definite$Beta[,estimation_first_definite$Pred],
                                    length(estimation_first_definite$Beta[,estimation_first_definite$Pred]), 1)
            X2_data <- matrix(X[,estimation_first_definite$Pred], length(X[,estimation_first_definite$Pred]),1)

        }else
        {
            beta_selected <- estimation_first_definite$Beta[,estimation_first_definite$Pred]
            X2_data <- X[,estimation_first_definite$Pred]
        }

        # defintion of the weights
        weights_new <- 1/norm_matrix_K(beta_selected, eigenval)


        ## @A: For adaptive step, we just consider the non zero predictors and delete the zero ones.
        ## you can see the length of weights_new can be less than the initial weights.

############################################################################
#
# Part 3: RUN adaptive step
#
############################################################################


#######################################################
#
# Part 3.1: estimation_second
#
#######################################################

        # the optimum value of lambda can be identified both with Cross-Validation
        # and BIC. (depending on the BIC_adaptive_step input parameter)


        ## @A: Here we changed X to X2_data. It's because we don't want to consider the zero predictors.

         if (BIC_adaptive_step == TRUE) # BIC
        {
            estimation_second <- definition_beta(X2_data, Y, matrix(, 1, 0),
                                                 eigenval, weights_new, NoI,
                                                 function_computation, thres,
                                                 number_non_zeros, numeric(),
                                                 ratio_lambda, number_lambda, BIC = 1, verbose = FALSE)
            lambda_second <- estimation_second$Lambda_vect

        }

        if (BIC_adaptive_step == FALSE) # Cross-Validation
        {

            estimation_second <- definition_beta(X2_data, Y, matrix(, 1, 0),
                                                 eigenval, weights_new, NoI,
                                                 function_computation, thres,
                                                 number_non_zeros, numeric(),
                                                 ratio_lambda, number_lambda, BIC = 0, verbose = FALSE)

            lambda_second <- estimation_second$Lambda_vect

#######################################################
#
# Part 3.2: estimation_second_CV and finding lamda which implies minimum CV
#
#######################################################

             # Cross Validation procedure
            if (verbose) print('Validation on the test set: identification of lambda')

            subset <- c(rep(2, ceiling(dim(X2_data)[1]*proportion_training_set)), rep(1, ceiling(dim(X2_data)[1]*(1-proportion_training_set)))) # 75% test set (2), 25% training:left_out (1)
            set.seed(16589)
            random_groups <- subset[sample(1:dim(X2_data)[1])]

            i <- 1 # training set is 1

            left_out <- which(random_groups==i)
            estimation_second_CV <- definition_beta_CV(X2_data[-left_out,], Y[, -left_out],
                                                       X2_data[left_out,], Y[, left_out],
                                                       eigenval, weights_new, NoI,
                                                       function_computation,
                                                       thres, number_non_zeros,
                                                       lambda_second, verbose = FALSE)

            #print(estimation_second_CV$error)

            ## @A: lambda_second_selected is the index of the lambda which implies the min C.V.
            lambda_second_selected <- which.min(estimation_second_CV$error)

            if (lambda_second_selected == lambda_second[length(lambda_second)])
            {
                warning('Last value of the gird of lambda selected by Cross Validation. Change the grid for lambda and run FLAME again!')
            }

            # print(lambda_second_selected)
            # plot(lambda_first,estimation_first_CV$error, pch=19, main='First step: selection of the value of lambda')
            # points(lambda_first[lambda_first_selected], estimation_first_CV$error[lambda_first_selected], pch=19, col=2)

            if(verbose)
            {
                print(paste("Non adaptive step: final estimation with the optimum lambda")) }

            # if (estimation_second_CV$reached_last==1)
            #    {
            # if the CV error identifies as optimum the last value of lambda,
            # we take as estimation the result of the previous run
            #     estimation_second_definite <- estimation_second
            #   }
            #  if (estimation_second_CV$reached_last==0){
            # otherwise we compute again the estimation with the
            # optimum value of lambda

#######################################################
#
# Part 3.3:  estimation_second_definite and finding lamda which implies minimum CV
#
#######################################################

            ## ?@A: Same as part 2.3 I have the problem here.
            ## In part 3.1 we calculated definition_beta for lambda={\lambda_max,....,\lambda_min}. In
            ## part 3.2 we find the lambda_second_selected  which is the label of tha lambda gives the min C.V.
            ## But again here we are gonna calculate the definition_beta for lambda={\lambda_max , ... , \lambda_{lambda_first_selected}} which
            ## the result should be a subset of results of part 3.1. Don't you think that? It does not affect on final result but does not change anything.

              estimation_second_definite <- definition_beta(X2_data, Y,
                                                          matrix(, 1, 0), #estimation_second_CV$Beta,
                                                          eigenval, weights_new,
                                                          NoI, function_computation,
                                                          thres, number_non_zeros,
                                                          lambda_second[1:lambda_second_selected],
                                                          0, 0, BIC = 0, verbose = FALSE)
            # }

            estimation_second <- estimation_second_definite


        }

#######################################################
#
# Part 3.3:  The output values
#
#######################################################

        predictors_2=estimation_second$Pred # final set of
        # predictors different from 0 (among the ones isolated in the
        # First: Non-Adaptive step)
        beta_2=estimation_second$Beta[,predictors_2] # estimated betas,
        # it is the matrix of the coefficients of the basis expansion
        # of the betas (with respect to the basis defined by the eigenvectors
        # of the kernel)


        ## @A: Pay attention: predictors_2 are the labels of nonzro predictors among the
        ## ones isolated in the First: Non-Adaptive step. So for finding the actual labels among the
        ## whole predictors we need to define "predictor_def" as follows

        predictor_def <- estimation_first_definite$Pred[predictors_2] # index of
        # the non-zero predictors computed in the estimation

        beta_def <- matrix(0, dim(estimation_first_definite$Beta)[1],
                           dim(estimation_first_definite$Beta)[2])

        ## @A: "beta_def" is a J*I matrix which the nonzero ones estimated by estimation_second and
        ## rest of the columns are zero.

        beta_def[, predictor_def] <- beta_2 # final matrix of the estimated betas
        # (still as coefficients of the kernel basis)

        if (verbose) {print(paste("Total number of non zeros predictor estimated is ", length(predictors_2)))}

          result <- list(beta=beta_def,
                       beta_no_adaptive=estimation_first_definite$Beta,
                       predictors=predictor_def,
                       predictors_no_adaptive=estimation_first_definite$Pred)
        return(result)
    }

}

## End
