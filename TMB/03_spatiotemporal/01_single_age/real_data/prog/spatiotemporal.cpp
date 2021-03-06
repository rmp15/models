// alpha + U + V; U ~ ICAR; V ~ Normal

#include <TMB.hpp>
//  using namespace density;

template<class Type>
Type objective_function<Type>::operator() () 
{
	
// DATA
DATA_MATRIX(deaths);		// matrix of death counts for single age group in multiple states
DATA_MATRIX(E);		    	// matrix of population for single age group in multiple states
DATA_SPARSE_MATRIX(P); 		// precision matrix for adjacency
size_t T = deaths.cols();   	// number of time points
size_t N = deaths.rows();   	// number of states

// PARAMETERS
PARAMETER(alpha);
PARAMETER(log_sigma2_V);
PARAMETER_MATRIX(V);

PARAMETER(log_sigma2_U);
PARAMETER_MATRIX(W);

// DEFINE DEATHS COUNTS DEATH RATE
matrix<Type> log_deaths_pred(N,T);
matrix<Type> deaths_pred = exp(log_deaths_pred);
matrix<Type> mu(N,T);

// INITIALISE NEGATIVE LOG-LIKELIHOOD
Type nll = 0;

// PRECISIONS
Type prec_V = 1 / exp(log_sigma2_V);
Type prec_U = 1 / exp(log_sigma2_U);

// ASSIGN HYPERPRIORS 
nll -= dnorm(alpha, Type(0), Type(10), 1);
nll -= dgamma(prec_V, Type(0.5), Type(2000), 1);
nll -= dgamma(prec_U, Type(0.5), Type(2000), 1);
nll -= dnorm(V, Type(0), exp(0.5 * log_sigma2_V), 1).sum();

// I DON'T UNDERSTAND THIS BIT ONWARDS...
matrix<Type> tmp = P * W;
nll -= -0.5 * (W * tmp).sum();
vector<Type> U = W * exp(0.5 * log_sigma2_U);

	nll -= dnorm(U.sum(), Type(0), Type(0.00001), 1);
	for (size_t i = 0; i < N; i++) log_deaths_pred(i) = log(E(i)) + alpha + V(i) + U(i);
	for (size_t i = 0; i < N; i++) nll -= dpois(deaths(i), exp(log_deaths_pred(i)), 1);
	for (size_t i = 0; i < N; i++) mu(i) = exp(alpha + V(i) + U(i));



	ADREPORT(U);
	ADREPORT(deaths_pred);
	ADREPORT(mu);

	return nll;
}
