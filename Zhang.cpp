// Model 1 of Zhang, Lessard and Campbell (2009) Fisheries Research | v0.4.2
// 3 randeff: Linf (Gaussian), K (Gausian) and A (gamma), all indep
// originally Bayesian
// option for uniform priors, user-supplied uniform bounds
// model L1 and L2 jointly, cond. indep give Linf and K, observed T1 and T2
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
	PARAMETER(logmuK); // Gaussian mean von Bert K, just positive
	PARAMETER(logsigmaK); //  Gaussian sd von Bert K
	PARAMETER(logmuA); // mean age A, gamma
	PARAMETER(logsigmaA); // sd age A, gamma
	PARAMETER(logsigmaeps); // Gaussian sd obs error

	// Random effects
	PARAMETER_VECTOR(logLinf); // iid fish-specific von Bert Linf, dim n
	PARAMETER_VECTOR(logK); // iid von Bert growth rate, dim n
	PARAMETER_VECTOR(logA); // iid fish-specific age, dim n

	// Hyperparameters
	DATA_VECTOR(lbubmuinf); // hyperparam unif prior muinf
	DATA_VECTOR(lbubK); // hyperparam unif prior K
	DATA_VECTOR(lbubgammashapeA); // hyperparam unif prior shape of A
	DATA_VECTOR(lbubgammarateA); // hyperparam unif prior rate of A
	DATA_VECTOR(lbubsigmainf); // hyperparam unif prior sigmainf
	DATA_VECTOR(lbubsigmaK); // hyperparam unif prior sigmaK
	DATA_VECTOR(lbubsigmaeps); // hyperparam unif prior sigmaeps

	// Misc
	DATA_INTEGER(enablepriors); // 0 = disabled, 1 = all enabled
	

	//--------------------------------------------------------------------------
	// Setup, procedures and init
	//--------------------------------------------------------------------------

	int n = L1.size(); // i = 1, ..., n

	Type muinf = exp(logmuinf);
	Type sigmainf = exp(logsigmainf);
	Type muK = exp(logmuK); // just positive
	Type sigmaK = exp(logsigmaK);
	Type muA = exp(logmuA);
	Type sigmaA = exp(logsigmaA);
	Type shapeA = square(muA/sigmaA); // gamma shape
	Type scaleA = square(sigmaA)/muA; // gamma scale
	Type sigmaeps = exp(logsigmaeps);

	vector<Type> Linf = exp(logLinf);
	vector<Type> K = exp(logK);
	vector<Type> A = exp(logA); // age relative to T0

	Type nll = 0.0; // init neg loglik


	//--------------------------------------------------------------------------
	// Priors
	//--------------------------------------------------------------------------
	
	if (enablepriors==1){
		// cf. Zhang et al. (2009, p. 291, top of right column)
		nll += neglogdunif(muinf, lbubmuinf(0), lbubmuinf(1), n);
		nll += neglogdunif(muK, lbubK(0), lbubK(1), n);
		Type rateA = Type(1.0)/scaleA; // unif prior on rate rather than scale
		nll += neglogdunif(shapeA, lbubgammashapeA(0), lbubgammashapeA(1), n);
		nll += neglogdunif(rateA, lbubgammarateA(0), lbubgammarateA(1), n);
		nll += neglogdunif(sigmainf, lbubsigmainf(0), lbubsigmainf(1), n);
		nll += neglogdunif(sigmaK, lbubsigmaK(0), lbubsigmaK(1), n);
		nll += neglogdunif(sigmaeps, lbubsigmaeps(0), lbubsigmaeps(1), n);
	}


	//--------------------------------------------------------------------------
	// Random effects
	//--------------------------------------------------------------------------

	// all randeff are iid (exchangeable)
	for (int i = 0; i < n; i++){
		nll -= dnorm(Linf(i), muinf, sigmainf, true);
		nll -= dnorm(K(i), muK, sigmaK, true);
		nll -= dgamma(A(i), shapeA, scaleA, true);
	}


	//--------------------------------------------------------------------------
	// Observation equations
	//--------------------------------------------------------------------------

	// conditionally on Linf and K, L1 and L2 are indep, variance from epsilon
	for (int i = 0; i < n; i++){
		Type lT1 = Linf(i)*(Type(1.0)-exp(-K(i)*A(i))); // vB at T1
		nll -= dnorm(L1(i), lT1, sigmaeps, true);
		// Type lT2 = Linf(i)*(Type(1.0)-exp(-K(i)*(A(i)+T2(i)-T1(i)))); // vB T2
		Type lT2 = Linf(i)*(Type(1.0)-exp(-K(i)*(A(i)+deltaT(i)))); // vB T2
		nll -= dnorm(L2(i), lT2, sigmaeps, true);
		// ^ vB reference point = T0
	}


	//--------------------------------------------------------------------------
	// Outputs
	//--------------------------------------------------------------------------

	REPORT(muinf); // to export MCMC post draw with tmbstan
	REPORT(muK); // to export MCMC post draw with tmbstan

	ADREPORT(muinf);
	ADREPORT(sigmainf);
	ADREPORT(muK);
	ADREPORT(sigmaK);
	ADREPORT(muA);
	ADREPORT(sigmaA);
	ADREPORT(sigmaeps);

	ADREPORT(Linf);
	ADREPORT(K);
	ADREPORT(A);

	return nll;
}
