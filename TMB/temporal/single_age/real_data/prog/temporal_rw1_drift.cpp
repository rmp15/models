#include <TMB.hpp>

using namespace density;

template<class Type>
Type objective_function<Type>::operator() ()
{

// SIMULATED DATA FOR POISSON REGRESSION
// X ~ Po(deaths)
// log deaths = alpha_0 + beta_0 * t + pi_t + overdispersion

// DATA
DATA_MATRIX(deaths);        	// matrix of death counts for single age group in multiple states
DATA_MATRIX(E);		    	// matrix of population for single age group in multiple states
size_t T = deaths.cols();   	// number of time points
size_t N = deaths.rows();   	// number of states

// PARAMETERS
// intercepts
// PARAMETER(alpha_0);          // implicit in random walk
// slopes
PARAMETER(beta_0);              // global slope
// rates
PARAMETER_MATRIX(log_mu);	// matrix of log(deathrate) for single age groups in multiple states
// precisions
PARAMETER(log_prec_rw);         // log precision of rw1
PARAMETER(log_prec_epsilon);    // log precision of overdispersion
// random walk
PARAMETER_MATRIX(pi); 

// INITIALISE NEGATIVE LOG-LIKELIHOOD
Type nll = Type(0.0);
    
// ASSIGN HYPERPRIORS TO PRECISIONS
nll -= dlgamma(log_prec_rw, Type(1), Type(500000), TRUE);
nll -= dlgamma(log_prec_epsilon, Type(1), Type(500000), TRUE);

// TRANSFORM PRECISIONS
Type log_sigma_rw       = (Type(-1) * log_prec_rw)       / Type(2) ;
Type log_sigma_epsilon  = (Type(-1) * log_prec_epsilon)  / Type(2) ;
   	
// RANDOM WALK
for (size_t n = 0; n < N; n++) {
        for (size_t t = 1; t < T; t++) {
                nll -= dnorm(pi(n,t), pi(n,t-1), exp(log_sigma_rw), TRUE);
        }
}

// PREDICTION
for (size_t n=0; n < N; n++) {
        for (size_t t=0; t < T; t++) {
                nll -= dnorm(log_mu(n,t), beta_0 * (t + 1) + pi(n,t), exp(log_sigma_epsilon), TRUE);
        }
}

// data likelihood
for (size_t n=0; n < N; n++) {
        for (size_t t=0; t < T; t++) {
                nll -= dpois(deaths(n,t), exp(log_mu(n,t)+ log(E(n,t))), TRUE);
        }
}

return nll;
}
