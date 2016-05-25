#include <TMB.hpp>

using namespace density;

template<class Type>
Type objective_function<Type>::operator() ()
{

// SIMULATED DATA FOR POISSON REGRESSION
// X ~ Po(deaths)
// log deaths = alpha_0 + beta_0 * t

// DATA
DATA_MATRIX(log_counts);        // matrix of log of counts for single age group in multiple states, 
                                // with time across and states downwards
size_t T = log_counts.cols();   // number of time points
size_t N = log_counts.rows();   // number of states

// PARAMETERS
// intercepts
PARAMETER(alpha_0);             // global intercept
// slopes
PARAMETER(beta_0);              // global slope

// ESTIMATED OUTPUT
PARAMETER_MATRIX(log_counts_pred);      // estimated count

// INITIALISED NEGATIVE LOG-LIKELIHOOD
Type nll = Type(0.0);

// data likelihood
for (size_t n=1; n < N; n++) {
        for (size_t t=1; t < T; t++) {
                log_counts_pred(n,t) = alpha_0 + beta_0 * t;
                nll -= dpois(log_counts(n,t), log_counts_pred(n,t), TRUE);
        }
}

return nll;
}
