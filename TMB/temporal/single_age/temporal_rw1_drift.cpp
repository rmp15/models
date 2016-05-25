#include <TMB.hpp>

using namespace density;

template<class Type>
Type objective_function<Type>::operator() ()
{

// SIMULATED DATA FOR POISSON REGRESSION
// X ~ Po(deaths)
// log deaths = alpha_0 + beta_0 * t + pi_t + overdispersion

// DATA
DATA_MATRIX(log_counts);        // matrix of log of counts for single age group in multiple states, 
                                // with time across and states downwards
size_t T = log_counts.cols();   // number of time points
size_t N = log_counts.rows();   // number of states

// PARAMETERS
// intercepts
// PARAMETER(alpha_0);          //implicit in random walk
// slopes
PARAMETER(beta_0);              // global slope

// precisions
PARAMETER(log_tau_rw);          // log precision of rw1 
PARAMETER(log_tau_epsilon);     // log precision of overdispersion

// ESTIMATED OUTPUT
PARAMETER_MATRIX(log_counts_pred);      // estimated count
PARAMETER_MATRIX(pi); 

// INITIALISED NEGATIVE LOG-LIKELIHOOD
Type nll = Type(0.0);

// TRANSFORM PRECISIONS
Type log_sigma_rw       = (Type(-1) * log_tau_rw)       / Type(2) ;
Type log_sigma_epsilon  = (Type(-1) * log_tau_epsilon)  / Type(2) ;
   	
// RANDOM WALK
for (size_t n = 0; n < N; n++) {
        for (size_t t = 1; t < T; t++) {
                nll -= dnorm(pi(n,t), pi(n,t-1), exp(log_sigma_rw), TRUE);
        }
}

// PREDICTION
for (size_t n=0; n < N; n++) {
        for (size_t t=0; t < T; t++) {
                nll -= dnorm(log_counts_pred(n,t), beta_0 * (t + 1) + pi(n,t), exp(log_sigma_epsilon), TRUE);
        }
}

// data likelihood
for (size_t n=0; n < N; n++) {
        for (size_t t=0; t < T; t++) {
                nll -= dpois(log_counts(n,t), log_counts_pred(n,t), TRUE);
        }
}

return nll;
}
