// Laslett, Eveson and Polachek (2002) CJFAS | v0.4.2
// 2 randeff: Linf (Gaussian) and A (lognormal), indep
// model L1 and L2 jointly, cond. indep given Linf and A, observed T1 and T2
// option for uniform priors, user-supplied uniform bounds
#include <TMB.hpp>

template <class Type> 
Type square(Type x){
	return x*x;
}

template <class Type>
Type neglogdunif(Type x, Type a, Type b, int n){
	// neg log density of unif(a,b), a<b
	// n = sample size, so that largecst dominates other negloglik contrib
	// Type largecst = Type(10.0)*n; // must dominate log(b-a) and other nll contrib
	Type largecst = Type(100.0)*n+Type(100.0);
	// ^ 10*n was not enough for small sample sizes
	Type halfdensity = Type(0.5)*log(b-a);
	Type res1 = CppAD::CondExpGt(x, a, halfdensity, largecst);
	Type res2 = CppAD::CondExpLt(x, b, halfdensity, largecst);
	return res1+res2; // neg log density if a<x<b
}

template<class Type>
Type objective_function<Type>::operator() () {

	//--------------------------------------------------------------------------
	// Inputs
	//--------------------------------------------------------------------------

	// Data
	DATA_VECTOR(L1); // lengths at cap, dim n
	DATA_VECTOR(L2); // lengths at recap, dim n
	// DATA_VECTOR(T1); // time at cap, dim n
	// DATA_VECTOR(T2); // time at recap, dim n
	DATA_VECTOR(deltaT); // time diff betwee cap and recap, dim n
	// ^ as of v0.3: directly supply deltaT

	// Fixed parameters
	PARAMETER(logmuinf); // Gaussian mean von Bert Linf
	PARAMETER(logsigmainf); // Gaussian sd von Bert Linf
	PARAMETER(logK); // von Bert growth rate, common to all fish
	PARAMETER(logmuA); // lognormal mean (exp scale) age at cap A
	PARAMETER(logsigmaA); // lognormal sd (exp scale) age at cap A
	PARAMETER(logsigmaeps); // Gaussian sd obs error, iid for Tcap and Trecap

	// Random effects
	PARAMETER_VECTOR(logLinf); // iid fish-specific von Bert Linf, dim n
	PARAMETER_VECTOR(logA); // iid fish-specific age relative to T0, dim n

	// Hyperparameters
	DATA_VECTOR(lbubmuinf); // hyperparam unif prior muinf
	DATA_VECTOR(lbubK); // hyperparam unif prior K
	DATA_VECTOR(lbubsigmainf); // hyperparam unif prior sigmainf
	DATA_VECTOR(lbubmuA); // hyperparam unif prior muA
	DATA_VECTOR(lbubsigmaA); // hyperparam unif prior sigmaA
	DATA_VECTOR(lbubsigmaeps); // hyperparam unif prior sigmaeps

	// Misc
	DATA_INTEGER(enablepriors); // 0 = disabled, 1 = all enabled
	

	//--------------------------------------------------------------------------
	// Setup, procedures and init
	//--------------------------------------------------------------------------

	int n = L1.size(); // i = 1, ..., n

	Type muinf = exp(logmuinf);
	Type sigmainf = exp(logsigmainf);
	Type K = exp(logK);
	Type muA = exp(logmuA);
	Type sigmaA = exp(logsigmaA);
	Type sigmaeps = exp(logsigmaeps);

	vector<Type> Linf = exp(logLinf);
	vector<Type> A = exp(logA); // age relative to T0

	Type nll = 0.0; // init neg loglik


	//--------------------------------------------------------------------------
	// Priors
	//--------------------------------------------------------------------------

	if (enablepriors==1){ // only uniform
		nll += neglogdunif(muinf, lbubmuinf(0), lbubmuinf(1), n);
		nll += neglogdunif(sigmainf, lbubsigmainf(0), lbubsigmainf(1), n);
		nll += neglogdunif(K, lbubK(0), lbubK(1), n);
		nll += neglogdunif(muA, lbubmuA(0), lbubmuA(1), n);
		nll += neglogdunif(sigmaA, lbubsigmaA(0), lbubsigmaA(1), n);
		nll += neglogdunif(sigmaeps, lbubsigmaeps(0), lbubsigmaeps(1), n);
	}

	//--------------------------------------------------------------------------
	// Random effects
	//--------------------------------------------------------------------------

	// mean and sd on log scale, because both muA and sigmaA on exp scale
	Type meanlogA = log(square(muA)/sqrt(square(muA)+square(sigmaA)));
	Type sdlogA = sqrt(log(Type(1.0)+square(sigmaA/muA)));

	for (int i = 0; i < n; i++){
		nll -= dnorm(Linf(i), muinf, sigmainf, true);
		nll -= dnorm(logA(i), meanlogA, sdlogA, true);
	}


	//--------------------------------------------------------------------------
	// Observation equations
	//--------------------------------------------------------------------------

	// conditionally on Linf, L1 and L2 are independent, variance from epsilon
	for (int i = 0; i < n; i++){
		Type lT1 = Linf(i)*(Type(1.0)-exp(-K*A(i))); // vB at T1
		nll -= dnorm(L1(i), lT1, sigmaeps, true);
		// Type lT2 = Linf(i)*(Type(1.0)-exp(-K*(A(i)+T2(i)-T1(i)))); // vB at T2
		Type lT2 = Linf(i)*(Type(1.0)-exp(-K*(A(i)+deltaT(i)))); // vB at T2
		nll -= dnorm(L2(i), lT2, sigmaeps, true);
		// ^ vB reference point = T0
	}


	//--------------------------------------------------------------------------
	// Outputs
	//--------------------------------------------------------------------------

	REPORT(muinf); // to export MCMC post draw with tmbstan
	REPORT(K); // to export MCMC post draw with tmbstan

	REPORT(meanlogA); // expectation of age at cap, log scale
	REPORT(sdlogA); // sd of age at cap, log scale

	ADREPORT(muinf);
	ADREPORT(sigmainf);
	ADREPORT(K);
	ADREPORT(muA); // expectation of age at cap, original exp scale
	ADREPORT(sigmaA); // sd of age at cap, original exp scale
	ADREPORT(sigmaeps);

	ADREPORT(Linf);
	ADREPORT(A);

	return nll;
}
