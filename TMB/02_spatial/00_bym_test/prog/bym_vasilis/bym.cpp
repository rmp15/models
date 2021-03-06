// alpha + U + V; V ~ Normal; U ~ ICAR

#include <TMB.hpp>
//  using namespace density;

template<class Type>
Type objective_function<Type>::operator() () {
	DATA_VECTOR(E);
	DATA_VECTOR(O);
	DATA_SPARSE_MATRIX(P); // precision matrix

	PARAMETER(alpha);
	PARAMETER(log_sigma2_V);
	PARAMETER_VECTOR(V);

	PARAMETER(log_sigma2_U);
	PARAMETER_VECTOR(W);

	int N = E.size();

	vector<Type> logmu(N);
	vector<Type> RR(N);

	Type nll = 0;

	Type tau_V = 1 / exp(log_sigma2_V);
	Type tau_U = 1 / exp(log_sigma2_U);

	nll -= dnorm(alpha, Type(0), Type(10), 1);
	nll -= dgamma(tau_V, Type(0.5), Type(2000), 1);
	nll -= dgamma(tau_U, Type(0.5), Type(2000), 1);
	nll -= dnorm(V, Type(0), exp(0.5 * log_sigma2_V), 1).sum();

	vector<Type> tmp = P * W;
	nll -= -0.5 * (W * tmp).sum();
	vector<Type> U = W * exp(0.5 * log_sigma2_U);

	nll -= dnorm(U.sum(), Type(0), Type(0.00001), 1);
	for (size_t i = 0; i < N; i++) logmu(i) = log(E(i)) + alpha + V(i) + U(i);
	for (size_t i = 0; i < N; i++) nll -= dpois(O(i), exp(logmu(i)), 1);
	for (size_t i = 0; i < N; i++) RR(i) = exp(alpha + V(i) + U(i));

	ADREPORT(U);
	ADREPORT(logmu);
	ADREPORT(RR);

	return nll;
}
