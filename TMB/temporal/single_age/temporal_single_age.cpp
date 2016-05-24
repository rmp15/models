#include <TMB.hpp>

using namespace density;

template<class Type>
Type objective_function<Type>::operator() ()
{

// tmb model to be compiled for a single age group over n states in t time points:
// 	model:	Poisson likelihood	
//  		log_mu = alpha_0 + alpha_m + alpha_s + alpha_ms
//			(beta_0 + beta_m + beta_s + beta_ms)*t
//			(pi_t)s + epsilon_st
  
// DATA
DATA_MATRIX(log_counts);        // matrix of log of counts for single age group in multiple states, 
				// with time across and states downwards
size_t T= log_counts.cols();	// number of time points
size_t N = log_counts.rows();	// number of states
  
// PARAMETERS
// intercepts
PARAMETER(alpha_0);             // global intercept
//PARAMETER_VECTOR(alpha_m);	// month-specific intercept
//PARAMETER_VECTOR(alpha_s);	// state-specific intercept
//PARAMETER_MATRIX(alpha_ms);	// month-state spatially correlated intercept

// slopes
PARAMETER(beta_0);		// global slope
//PARAMETER_VECTOR(beta_m);	// month-specific slope
//PARAMETER_VECTOR(beta_s);	// state-specific slope
//PARAMETER_MATRIX(beta_ms);	// month-state spatially correlated slope

// precisions
PARAMETER(log_tau_rw);      	// log precision of rw1 
PARAMETER(log_tau_epsilon);     // log precision of observations

// ESTIMATED OUTPUT
PARAMETER(log_counts_pred);	// estimated count

// INITIALISED NEGATIVE LOG-LIKELIHOOD
Type nll = Type(0.0);

// TRANSFORM PRECISIONS
Type log_sigma_rw 	= (Type(-1) * log_tau_rw) 	/ Type(2)
Type log_sigma_epsilon 	= (Type(-1) * log_epsilon_rw) 	/ Type(2)
   	
// RANDOM WALK
matrix<Type> pi(N,T);
for (size_t n = 1; n < N; n++) {
	for (size_t t = 0; t < T; t++) {
	//	if (t==0) {
	//		pi(n,t) =  Type(0.);
	//	}
	//	else {
	//		nll -= dnorm(pi(n,t), pi(n,t-1), exp(log_sigma_rw), TRUE);
	//	}
		nll -= dnorm(pi(n,t), pi(n,t-1), exp(log_sigma_rw), TRUE);
	}
}

// (RANDOM WALK EXAMPLE)
//for(int i = 1; i < n; i++){
//    Type m = a*u[i-1] ;    //linear model
//  nll -= dnorm(u[i], m, sigma_proc, true); //likelihood for random effects
//}

// PREDICTION
matrix<Type> log_counts_pred(N,T);
for (size_t t=0; t < T; t++) {
	for (size_t n=0; n < N; n++) {
		nll -= dnorm(log_counts_pred(n,t), alpha_0 + beta_0 * t + pi(n,t), exp(log_sigma_epsilon), TRUE);
	}
}

// data likelihood
for (size_t t=0; t < T; t++) {
	for (size_t n=0; n < N; n++) {
		nll -= dpois(log_counts(n,t), log_counts_pred(n,t), TRUE);
	}
}

return nll;
}
