#include <math.h>
#include <cstdio>
//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <functional>


// [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::plugins(cpp11)]]
// // [[Rcpp::export]]
class linear_model {

private:
    arma::mat Y, X, B, E, B_ls;
    unsigned int num_pred;
    unsigned int num_data;
    unsigned int num_zeros_pred;
    double lambda_max;

public:

    // c++ function which calls the R function to perform the optimization and
    // compute the norm of beta_j
    double computation_norm_cpp( double &lambda, double &omega, arma::vec &B_j,
                                 arma::vec &tau,  Rcpp::Function f)
    {
        return( Rcpp::as<double>(f(lambda, omega, B_j, tau)) );
    }

    // norm in K of K(x), given the components in the kernel basis (B) and the
    // eigenvalues of the kernel (tau)
    double norm_K_Kx(arma::vec &B, arma::vec &tau)
    {
        return(sqrt(sum(B % B % tau)));
    }

    // K norm of the matrix B_here (square-root of the sum of the squares of the
    // K norms of the columns of B_here). tau is the vector of the eigenvalues of
    // the kernel
    double norm_K_matrix(arma::mat &B_here, arma::vec tau)
    {
        double norm_matrix=0;
        for (unsigned int i=0; i<B_here.n_cols; i++)
        {
            norm_matrix += sum(B_here.col(i)%B_here.col(i)/tau);
        }
        return (sqrt(norm_matrix));
    }

    // square of the H norm of the matrix Y_here (sum of the squares of the
    // H norms of the columns of Y_here)
    double square_norm_H_matrix(arma::mat &Y_here)
    {
        double norm_matrix=0;
        for (unsigned int i=0; i<Y_here.n_cols; i++)
        {
            norm_matrix += sum(Y_here.col(i)%Y_here.col(i));
        }
        return (norm_matrix);
    }

    // square of the K norm of the matrix Y_here (sum of the squares of the
    // K norms of the columns of Y_here). tau is the vector of the eigenvalues of
    // the kernel
    double square_norm_K_matrix(arma::mat &Y_here, arma::vec &tau)
    {
        double norm_matrix=0;
        for (unsigned int i=0; i<Y_here.n_cols; i++)
        {
            norm_matrix += sum(Y_here.col(i)%Y_here.col(i)/tau);
        }
        return (norm_matrix);
    }

    // computation of the BIC given as input
    // - Y_here: kxN matrix; the matrix of the true values of the observations in the kernel basis
    // - X_here: Nxp matrix; predictors' matrix
    // - B_here: kxp matrix; estimated coefficients
    // - p: int; number of non-zero predictors identified
    // - N: int; number of observations
    // - tau: k vector; eigenvalues of the kernel
    double computation_BIC (arma::mat &Y_here, arma::mat &X_here, arma::mat &B_here,
                            unsigned int p, unsigned int N, arma::vec &tau)
    {
        arma::mat difference = Y_here - B_here*X_here.t();
        // we can compute the RSS both with the H and K norm.
        // double error_H = square_norm_H_matrix(difference); //H norm
        double error_H = square_norm_K_matrix(difference, tau); //K norm
        //std::cout<<error_H<<N<<p<<std::endl;
        double BIC_computed =  N*log(error_H/N) + log(N)*p;
        return BIC_computed;

    }

    // constructor of the linear_model element given the matrix of the predictors
    // X_data (Nxp) and the matrix of the observations Y_data (kxN)
    linear_model(const arma::mat &X_data, const arma::mat &Y_data){
        Y=Y_data;
        X=X_data;
        E=Y; // inital estimation of the error equal to the data, since beta is initalized with 0
        num_data = X.n_rows;
        B_ls = E * X / num_data; // first LS estimation of beta
        num_pred= X.n_cols;
        num_zeros_pred=0;
        B.zeros(B_ls.n_rows, B_ls.n_cols);
    }

    // constructor of the linear_model element given the matrix of the predictors
    // X_data (Nxp), the matrix of the observations Y_data (kxN) and the first
    // estimation of the coefficients Beta (kxp)
    linear_model(const arma::mat &X_data, const arma::mat &Y_data, const arma::mat &Beta){
        Y=Y_data;
        X=X_data;
        num_data = X.n_rows;
        B = Beta;
        E= Y - B*X.t();
        B_ls = E * X / num_data;
        num_pred= X.n_cols;
        num_zeros_pred=0;
    }

    // output function to print the number of the predictors and the predictor matrix
    void predictors()
    {
        Rcpp::Rcout<<"matrix of predictors" <<X <<"number of predictors "<< num_pred;
    }
    
    //make the betas 0
    void set_B_zeros()
    {
        B.zeros();
    }

    // definition of the vector of all the possible values of lambda,
    // given the number of lambdas (num_lambda) and the ratio between the maximum
    // and the minimum (ratio_lambda).
    // in particular the vector is defined from lambda_max
    // (computes as the minimum value for which all the predictors are
    // guaranteed to be 0) to lambda_min = lambda_max*raito_lambda
    // not equispaced vector
    arma::vec definition_lambda(double ratio_lambda, int num_lambda, arma::vec omega, arma::vec tau)
    {

        arma::vec lambda_max_vec(num_pred);
        arma::vec lambda(num_lambda);
        for (unsigned int j=0; j<num_pred; j++)
        {
            arma::vec B_temp=B_ls.col(j);
            lambda_max_vec(j) = norm_K_Kx(B_temp, tau)/omega(j);
            //std::cout<<"j="<<j<<" lambda "<<lambda_max_vec(j)<<std::endl;
        }
        lambda_max = 2*max(lambda_max_vec); //CAMPPPPP
        // std::cout<<"maximum lambda "<<lambda_max<<std::endl;
        double lambda_min= lambda_max*ratio_lambda;
        // std::cout<<"defintion of a set of lambda form "<<lambda_max<<" to "<<
        //  lambda_min<< " of length "<< num_lambda<<std::endl;

        double tau_max = log10(lambda_max);
        double tau_min = log10(lambda_min);

        arma::vec tau_vect(num_lambda);
        for (unsigned int i=0; i<num_lambda; i++)
        {
            tau_vect(i)=tau_max - i*(tau_max-tau_min)/(num_lambda -1);

          //  std::cout<<tau_vect(i)<<" "<<std::endl;
        }

        lambda = arma::exp10(tau_vect);

        // uncomment these lines to define an equispaced vector

        //     for (unsigned int i=0; i<num_lambda; i++)
        //     {
        //       lambda(i) = lambda_max - i*(lambda_max-lambda_min)/(num_lambda-1);
        //     }

        return lambda;
    }

    // MAIN FUNCTION TO COMPUTE THE ESTIMATION,
    // called from define_beta
    Rcpp::List estimate( arma::vec &tau, arma::vec &omega_vect,
                         const int &N_iter, Rcpp::Function computation_norm, double &thres,
                         unsigned int &non_zeros_pred, arma::vec &lambda,
                         unsigned int BIC,
                         bool verbose)
    {


      //  B.zeros();
        arma::mat B_subset(B.n_rows, B.n_cols-1);
        arma::mat X_subset(X.n_rows, X.n_cols-1);

        arma::vec X_j(num_data);
        arma::vec B_j(B.n_rows);

        arma::mat E_glob(E.n_rows, E.n_cols);
        arma::mat B_old(B.n_rows, B.n_cols);

        double lambda_iter;
        arma::vec vector_non_zeros_pred(num_pred);
        double num_zeros_pred;
        int number_lambda_computed=lambda.n_elem;


        arma::vec BIC_computed(lambda.n_elem); //to store the BIC of each iteration

        arma::mat Beta_def_BIC;
        arma::vec Pred_def_BIC;
        unsigned int number_non_zeros_def_BIC;
        double lambda_BIC;
        double BIC_value=1000000; //initialization of the best BIC value


        for (unsigned int l=0; l<lambda.n_elem; l++) //loop on lambda
        {
            lambda_iter=lambda(l);

            if (verbose) {Rcpp::Rcout << "Lambda = "<< lambda_iter<< std::endl;}

            for (unsigned int k=0; k<N_iter; k++)
            {
                if (verbose) {Rcpp::Rcout << "*";}

                //  if (k==0 && l==0)
                 // {
                  //         E_glob=E;
                 // }else{
                     E_glob=Y-B*X.t(); //initialization with the previous estimation
                    
                 // }
                   
                B_old=B;

                num_zeros_pred=0;
                vector_non_zeros_pred.zeros(num_pred,1);

                for (unsigned int j=0; j<num_pred; j++)
                {
                    double omega=omega_vect(j);

                    X_j=X.col(j);
                    B_j= B_old.col(j);

                    arma::vec B_tilde(B_j.n_elem);

                    //if(k != 0)
                //    {
                        E= E_glob + B_j * X_j.t();

                        B_tilde = E * X_j / num_data;
                  //  }else
                //    {
                  //      B_tilde = B_j;
                //    }


                    if (norm_K_Kx(B_tilde, tau) <= lambda_iter*omega)
                    {
                        B_j.zeros();
                        num_zeros_pred ++;

                    }else
                    {
                        vector_non_zeros_pred(j-num_zeros_pred)=j;

                        double norm_B_j;

                        norm_B_j = computation_norm_cpp(lambda_iter, omega, B_tilde,
                                                        tau, computation_norm); // numerical optimization to compute the norm
                        B_j = tau % B_tilde * norm_B_j / (tau * norm_B_j + lambda_iter*omega );

                    }
                    B.col(j) = B_j;
                } //end of the predictors loop

                arma::mat error_beta_matrix = B_old-B;
                // check on the improuvment of the estimation
                if (norm_K_matrix(error_beta_matrix, tau) < thres)
                {
                    // std::cout<<std::endl<<"Not significant change of Beta => Stop iteration at iter = "<<k+1<<std::endl;
                    k=N_iter;
                }


            } //end of the N_iter loop


            vector_non_zeros_pred.resize(num_pred-num_zeros_pred);

            if (verbose)
            {
                Rcpp::Rcout<<std::endl<<"number of non zero predictors fitted with lambda "<<lambda_iter<<
                    " is " <<num_pred-num_zeros_pred<<std::endl;
            }
            
            /*
            if (BIC==1) // if we need to compute the BIC to identify the best value of lambda
            {
                unsigned int p=num_pred-num_zeros_pred;
                //  std::cout<<"number predictors"<<p<<" number data"<<num_data<<std::endl;
                BIC_computed(l) = computation_BIC(Y,X,B,p,num_data, tau);
                //    std::cout<<"BIC "<< BIC_computed(l)<<std::endl;
                // std::cout<<"BIC"<<BIC_computed(l)<<std::endl;
                if (l==0)
                {
                    Beta_def_BIC = B;
                    Pred_def_BIC=vector_non_zeros_pred+1;
                    number_non_zeros_def_BIC=num_pred-num_zeros_pred;
                    lambda_BIC=lambda(l);
                    BIC_value=BIC_computed(l);
                }
                   if (l!=0)
                {
                if ( (BIC_computed(l) < BIC_computed(l-1)) && (BIC_computed(l-1) < BIC_value) )
                {
                Beta_def_BIC=B;
                Pred_def_BIC=vector_non_zeros_pred+1;
                number_non_zeros_def_BIC=num_pred-num_zeros_pred;
                lambda_BIC=lambda(l);
                } else
                {
                if (BIC_computed(l) >= BIC_computed(l-1))
                {
                if (BIC_computed(l-1)< BIC_value )
                {
                BIC_value = BIC_computed(l-1);
                }
                }
                }
                } 

                if(l!=0)
                {
                    if (BIC_computed(l)<BIC_value)
                    {
                        Beta_def_BIC=B;
                        Pred_def_BIC=vector_non_zeros_pred+1;
                        number_non_zeros_def_BIC=num_pred-num_zeros_pred;
                        lambda_BIC=lambda(l);
                        BIC_value = BIC_computed(l);
                    }
                }



            } */
        
         //    std::cout<<"number of zeros predictor with lambda "<<lambda_iter<<" is "<< num_zeros_pred<< std::endl;

            if( num_pred-num_zeros_pred >= non_zeros_pred )
            {
                // std::cout << "Number of non zeros predictors reached!";
                number_lambda_computed=l+1;
                // std::cout<<"number lamba analyzed "<<number_lambda_computed<<std::endl;
                l = lambda.n_elem ;
            }

        } //end of the lambda loop (l)

        //   std::cout<<"last value of lambda "<< lambda_iter<<" which identifies "<< num_zeros_pred<< " zeros predictors "<<std::endl;

        arma::vec lambda_subset=lambda.rows(0,number_lambda_computed-1);

        if (BIC==1)
        {
            // uncomment these lines to print the list of the non-zeros predictors
            /*  std::cout<<Pred_def_BIC.n_elem<<"predictors";
            for (unsigned int i = 0 ; i< Pred_def_BIC.n_elem; i++ )
            {
            std::cout<<Pred_def_BIC(i);
            }
            std::cout<<std::endl; */

            arma::vec BIC_computed_subset=BIC_computed.rows(0,number_lambda_computed-1);
            return(Rcpp::List::create(
                    Rcpp::_["Beta"] = Beta_def_BIC,
                    Rcpp::_["Pred"] = Pred_def_BIC, // R codifies c++ + 1
                    Rcpp::_["Number_non_zeros"] = number_non_zeros_def_BIC,
                    Rcpp::_["estimated_lambda"] = lambda_BIC,
                    Rcpp::_["Lambda_vect"] = lambda_subset,
                    Rcpp::_["BIC_vector"]=BIC_computed_subset,
                    Rcpp::_["optimum_BIC"] = BIC_value,
                    Rcpp::_["optimum_lambda"] = lambda_BIC)
            );
        }else
        {
            return(Rcpp::List::create(
                    Rcpp::_["Beta"] = B,
                    Rcpp::_["Pred"] = vector_non_zeros_pred + 1, // R codifies c++ + 1
                    Rcpp::_["Number_non_zeros"] = num_pred-num_zeros_pred,
                    Rcpp::_["estimated_lambda"] = lambda_iter,
                    Rcpp::_["Lambda_vect"] = lambda_subset)
            );
        }


    };


    Rcpp::List estimate_CV(const arma::mat X_test, const arma::mat Y_test, arma::vec &tau, arma::vec &omega_vect,
                           const int &N_iter, Rcpp::Function computation_norm, double &thres,
                           unsigned int &non_zeros_pred, arma::vec &lambda, bool verbose)
    {

        //B.zeros();

        arma::mat B_subset(B.n_rows, B.n_cols-1);
        arma::mat X_subset(X.n_rows, X.n_cols-1);

        arma::mat Beta_best;
        
        unsigned int updated_beta=0;

        arma::vec X_j(num_data);
        arma::vec B_j(B.n_rows);

        arma::mat E_glob(E.n_rows, E.n_cols);
        arma::mat B_old(B.n_rows, B.n_cols);

        double lambda_iter;
        arma::vec vector_non_zeros_pred(num_pred);
        double num_zeros_pred;
        int number_lambda_computed=lambda.n_elem;

        arma::vec error_lambda(number_lambda_computed);
        
        
        double error_def = 100000000;

        for (unsigned int l=0; l<lambda.n_elem; l++)
        {
            lambda_iter=lambda(l);

            if (verbose) {Rcpp::Rcout<< "Lambda = "<< lambda_iter<< std::endl;}

            for (unsigned int k=0; k<N_iter; k++)
            {
                  if (verbose) {Rcpp::Rcout << "*";}

                  //  if (k==0 && l==0)
                //    {
                  //    E_glob=E;
                   //}else{
                        E_glob = Y - B*X.t();
             // }
                B_old=B;

                num_zeros_pred=0;
                vector_non_zeros_pred.zeros(num_pred,1);

                for (unsigned int j=0; j<num_pred; j++)
                {
                    // std::cout<<"     Predictor "<<j+1<<std::endl;

                    double omega=omega_vect(j);

                    X_j=X.col(j);
                    B_j= B.col(j);

                    E = E_glob + B_j * X_j.t();

                    arma::vec B_tilde = E * X_j / num_data;


                    if (norm_K_Kx(B_tilde, tau) <= lambda_iter*omega)
                    {
                        B_j.zeros();
                        num_zeros_pred ++;

                    }else
                    {
                        vector_non_zeros_pred(j-num_zeros_pred)=j;

                        double norm_B_j;

                        norm_B_j = computation_norm_cpp(lambda_iter, omega, B_tilde, tau, computation_norm);
                        B_j = tau % B_tilde * norm_B_j / (tau * norm_B_j + lambda_iter*omega );

                    }
                    B.col(j) = B_j;
                } //end of the predictors loop

                arma::mat error_beta_matrix = B_old-B;
                if (norm_K_matrix(error_beta_matrix, tau) < thres)
                {
                    //    std::cout<<std::endl<<"Not significant change of Beta => Stop iteration at iter = "<<k+1<<std::endl;
                    k=N_iter-1;
                }

            } //end of the N_iter loop

            if (verbose)
            {
                Rcpp::Rcout<<std::endl<<"number of non zero predictors fitted with lambda "<<lambda_iter<<
                    " is " <<num_pred-num_zeros_pred<<std::endl;
            }
            
            arma::mat Y_pred_subset = B * X_test.t();
            arma::mat difference = Y_pred_subset-Y_test;
            error_lambda(l) = norm_K_matrix(difference, tau);
            //error_lambda(l) = square_norm_H_matrix(difference);
            if(verbose) {Rcpp::Rcout<<"error K norm "<< error_lambda(l) <<std::endl;}
            

           if (l!=0)
            {
                //   std::cout<<error_lambda(l)<<error_lambda(l-1);
                if (( error_lambda(l) < error_def ) && ( updated_beta == 0) )
                {
                    Beta_best=B;
                    error_def = error_lambda(l);
                }else
                {
                    //  std::cout<<"not updated beta"<<std::endl;
                    //      updated_beta=1;
                }
            }

        } //end of the lambda loop

        unsigned int reached_last;
        if (updated_beta==1) // if it is 1, you have stopped the update of Beta before since you have reached a min
        {
            reached_last = 0;
        }else{
            reached_last = 1;
        }

        return Rcpp::List::create(
            Rcpp::_["Beta"] = Beta_best,
            Rcpp::_["error"] = error_lambda,
            Rcpp::_["reached_last"] = reached_last
        );


    }


    };

// function to compute the estimation of Beta.
// called from R and explained in the file FLAME.R

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List definition_beta(const arma::mat X, const arma::mat Y, arma::mat Beta, arma::vec tau,
                           arma::vec omega_vect, const int N_iter,
                           Rcpp::Function computation_norm, double thres, unsigned int non_zeros_pred,
                           arma::vec lambda_start,
                           double ratio_lambda, int num_lambda, unsigned int BIC, bool verbose)
{
    linear_model LM(X,Y);
    
    if (Beta.n_cols!=0)
    {
        Rcpp::Rcout<<"Beta read from input"<<std::endl;
        linear_model LM_2(X,Y, Beta);
        LM = LM_2;
    }

    arma::vec lambda_estimation;

    if (lambda_start.is_empty())
    {
        lambda_estimation=LM.definition_lambda(ratio_lambda, num_lambda, omega_vect, tau);
    }else
    {
        //Rcpp::Rcout<<"lambda is not empty => ratio_lambda and num_lambda ignored!"<<std::endl;
        lambda_estimation =  lambda_start;
    }


  /*  if (Beta.n_cols==0)
    {
         LM.set_B_zeros();
    }
*/

    Rcpp::List estimation;
    // call the main function of the class linear_model
    estimation = LM.estimate(tau,  omega_vect, N_iter, computation_norm, thres, non_zeros_pred, lambda_estimation, BIC, verbose);
    return(estimation);
}


// function to identify the best value of lambda with the cross validation method.
// called from R and described in FLAME.R
// N.B. the value of lambda_start MUST be a vector: the vector of all the
// values of lambda!!! Not empty!

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List definition_beta_CV(const arma::mat X_train, const arma::mat Y_train, const arma::mat X_test,const arma::mat Y_test,
                              arma::vec tau, arma::vec omega_vect, const int N_iter,
                              Rcpp::Function computation_norm, double thres, unsigned int non_zeros_pred,
                              arma::vec lambda_start, bool verbose)
{
    linear_model LM(X_train, Y_train);
    
   // LM.set_B_zeros();
    
    arma::vec lambda_estimation  = lambda_start;

    Rcpp::List est_CV;
    est_CV = LM.estimate_CV( X_test, Y_test,
                             tau, omega_vect, N_iter, computation_norm, thres, non_zeros_pred,
                             lambda_estimation, verbose);
    return(est_CV);
};

