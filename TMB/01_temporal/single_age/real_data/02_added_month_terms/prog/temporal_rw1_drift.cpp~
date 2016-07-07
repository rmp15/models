#include <TMB.hpp>
//#include <fenv.h>

using namespace density;

template<class Type>
Type objective_function<Type>::operator() ()
{

//feenableexcept(FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO | FE_UNDERFLOW);

// SIMULATED DATA FOR POISSON REGRESSION
// X ~ Po(deaths)
// log deaths = alpha_0 + alpha_m + beta_0 * t + pi_t + overdispersion

// DATA
DATA_MATRIX(deaths);        	// matrix of death counts for single age group in multiple states
DATA_MATRIX(E);                 // matrix of population for single age group in multiple states

size_t T = deaths.cols();   	// number of time points
size_t N = deaths.rows();   	// number of states

// PARAMETERS
// intercepts
//PARAMETER(alpha_0);		// global glope
PARAMETER_VECTOR(alpha_m);      // month specific intercept
// slopes
PARAMETER(beta_0);              // global slope
//PARAMETER_VECTOR(beta_m);       // month specific slope
// rates
PARAMETER_MATRIX(log_mu);       // matrix of log(deathrate) for single age groups in multiple states
// precisions
PARAMETER(log_prec_rw);         // log precision of rw1
PARAMETER(log_prec_epsilon);    // log precision of overdispersion
PARAMETER(log_prec_int_m);      // log precision of month intercepts
//PARAMETER(log_prec_slp_m);    // log precision of month slopes
// random walk
PARAMETER_MATRIX(pi); 

// variables to report
matrix<Type> mu(N, T);
matrix<Type> dth(N, T);
Type alpha_0;

// INITIALISE NEGATIVE LOG-LIKELIHOOD
Type nll = Type(0.0);

// ASSIGN PRIORS TO PARAMETERS
//nll -= dnorm(alpha_0, Type(0), Type(10), 1);
nll -= dnorm(beta_0, Type(0), Type(10), 1);

// ASSIGN PRIORS TO PRECISIONS
nll -= dlgamma(log_prec_rw, Type(1), Type(1000), TRUE);
nll -= dlgamma(log_prec_epsilon, Type(1), Type(1000), TRUE);
nll -= dlgamma(log_prec_int_m, Type(1), Type(1000), TRUE);
//nll -= dlgamma(log_prec_slp_m, Type(1), Type(1000), TRUE);

// TRANSFORM PRECISIONS
Type log_sigma_rw       = (Type(-1) * log_prec_rw)       / Type(2) ;
Type log_sigma_epsilon  = (Type(-1) * log_prec_epsilon)  / Type(2) ;
Type log_sigma_int_m 	= (Type(-1) * log_prec_int_m)  	 / Type(2) ;
//Type log_sigma_slp_m 	= (Type(-1) * log_prec_slp_m)  	 / Type(2) ;
   	
// RANDOM WALK OVER TIME
for (size_t n = 0; n < N; n++) {
        for (size_t t = 1; t < T; t++) {
                nll -= dnorm(pi(n,t), pi(n,t-1), exp(log_sigma_rw), TRUE);
        }
}

// RANDOM WALK FOR MONTH TERMS
for (size_t m = 1; m < 12; m++) {
        nll -= dnorm(alpha_m(m), m == 1 ? 0 : alpha_m(m-1), exp(log_sigma_int_m), TRUE);
        //nll -= dnorm(beta_m(m), beta_m(m-1), exp(log_sigma_slp_m), TRUE);
        //nll -= dnorm(beta_m(m), m == 1 ? 0 : beta_m(m-1), exp(log_sigma_slp_m), TRUE);
}

// PREDICTION
for (size_t n=0; n < N; n++) {
        for (size_t t=0; t < T; t++) {
		//nll -= dnorm(log_mu(n,t), (beta_m(t%12) + beta_0) * t + pi(n,t), exp(log_sigma_epsilon), TRUE);
        	nll -= dnorm(log_mu(n,t), alpha_m(t%12) + beta_0 * t + pi(n,t), exp(log_sigma_epsilon), TRUE);
		//nll -= dnorm(log_mu(n,t), beta_0 * t + pi(n,t), exp(log_sigma_epsilon), TRUE);
	}
}

// PREDICTION
for (size_t n=0; n < N; n++) {
        for (size_t t=0; t < T; t++) {
                mu(n, t) = exp(log_mu(n, t));
                dth(n, t) = mu(n, t) * E(n, t);
        }
}

// data likelihood
for (size_t n=0; n < N; n++) {
        for (size_t t=0; t < T; t++) {
                nll -= dpois(deaths(n,t), dth(n, t), TRUE);
        }
}

// report variables
alpha_0 = pi.sum() / T;

ADREPORT(alpha_0);
ADREPORT(mu);
ADREPORT(dth);
    
return nll;
}
