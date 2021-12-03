// Bayesian formulation of Francis (1988) | v0.4.5
// Our take: not use reparametrization in terms of g_alpha and g_beta, only
//           consider extension of Fabens (1965) with various variance functions
//           added to basic error sd sigma.
// Gaussian likelihood assumed, as in Francis (1988) and Bayesian Fabens (1965)
// priors: choice on Linf, K, and sigma, hard-coded uniform for other par below
// model length at recap L2, given observed T1, T2 and L1
#include <TMB.hpp>

template <class Type>
Type neglogdunif(Type x, Type a, Type b, int n){
	// neg log density of unif(a,b), a<b
	// n = sample size, so that largecst dominates other negloglik contrib
	Type largecst = Type(100.0)*n+Type(100.0);
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
	DATA_VECTOR(deltaT); // time diff betwee cap and recap, dim n

	// Parameters
	PARAMETER(logLinf); // vB Linf
	PARAMETER(logK); // vB growth rate
	PARAMETER(logsigma); // additive Gaussian error sd
	PARAMETER(logsda); // varfunc extra par
	PARAMETER(logsdb); // varfunc extra par

	// Random effects: none

	// Hyperparameters
	DATA_VECTOR(hp_Linf); // hyperparam for prior on Linf, dim 2
	DATA_VECTOR(hp_K); // hyperparam for prior on K, dim 2
	DATA_VECTOR(hp_sigma); // hyperparam for prior on sigma, dim 2
	// ^ hp have different meaning depending on priordist:
	//   - priordist=0: unused
	//   - priordist=1: lower and upper bound of uniform density
	//   - priordist=2: mean and sd of Gaussian density
	//   - priordist=3: mean and sd of Gaussian density on log scale (lognormal)

	// Misc
	DATA_IVECTOR(priordist); // code for (Linf,K,sigma), dim 3
	// ^ 0 = no prior, 1 = uniform, 2 = Gaussian, 3 = lognormal
	DATA_INTEGER(priorsd); // 0 = no prior on sda/sdb, 1 = uniform (hard-coded)
	DATA_INTEGER(varfunc);
	// ^ code for variance function, with the following obs lkhd sd:
	//   - 0 = "constant" = sigma, i.e. equivalent to Fabens (1965)
	//   - 1 = "prop.dL" = sigma + a*E[deltaL] as in standard Francis (1988)
	//   - 2 = "prop.dT" = sigma + a*deltaT, cf. Francis (1988) remark on p. 43
	//   - 3 = "prop.L2" = sigma + a*E[L2], growth variability at later ages
	//   - 4 = "exp.dL" = sigma + a*(1-exp(b*E[deltaL])), Francis (1988) eq. (7)
	//   - 5 = "exp.L2" = sigma + a*(1-exp(b*E[L2])), some variant
	//   - 6 = "pow.dL" = sigma + a*E[deltaL]^b, Francis (1988) eq. (8)
	//   - 7 = "pow.L2" = sigma + a*E[L2]^b, some variant



	//--------------------------------------------------------------------------
	// Setup, procedures and init
	//--------------------------------------------------------------------------

	int n = L1.size(); // i = 1, ..., n

	Type Linf = exp(logLinf);
	Type K = exp(logK);
	Type sigma = exp(logsigma);
	Type sda = exp(logsda);
	Type sdb = exp(logsdb);
	
	Type nll = 0.0; // init neg loglik


	//--------------------------------------------------------------------------
	// Priors
	//--------------------------------------------------------------------------

	if (priordist(0)==1){ // uniform prior on Linf
		nll += neglogdunif(Linf, hp_Linf(0), hp_Linf(1), n);
	} else if (priordist(0)==2){ // Gaussian prior on Linf
		nll -= dnorm(Linf, hp_Linf(0), hp_Linf(1), true);
	} else if (priordist(0)==3){ // lognormal prior on Linf
		nll -= dnorm(logLinf, hp_Linf(0), hp_Linf(1), true) - logLinf;
		// ^ lognormal log-pdf evaluated at param on exp scale
	} // else no prior on Linf

	if (priordist(1)==1){ // uniform prior on K
		nll += neglogdunif(K, hp_K(0), hp_K(1), n);
	} else if (priordist(1)==2){ // Gaussian prior on K
		nll -= dnorm(K, hp_K(0), hp_K(1), true);
	} else if (priordist(1)==3){ // lognormal prior on K
		nll -= dnorm(logK, hp_K(0), hp_K(1), true) - logK;
		// ^ lognormal log-pdf evaluated at param on exp scale
	} // else no prior on K

	if (priordist(2)==1){ // uniform prior on sigma
		nll += neglogdunif(sigma, hp_sigma(0), hp_sigma(1), n);
	} else if (priordist(2)==2){ // Gaussian prior on sigma
		nll -= dnorm(sigma, hp_sigma(0), hp_sigma(1), true);
	} else if (priordist(2)==3){ // lognormal prior on sigma
		nll -= dnorm(logsigma, hp_sigma(0), hp_sigma(1), true) - logsigma;
		// ^ lognormal log-pdf evaluated at param on exp scale
	} // else no prior on sigma


	if (priorsd==1){ // then enable uniform prior on sda/sdb
		if (varfunc==1){ // prop.dL
			nll += neglogdunif(sda, Type(0.0), Type(1000.0), n);
		} else if (varfunc==2){ // prop.dT
			nll += neglogdunif(sda, Type(0.0), Type(1000.0), n);
		} else if (varfunc==3){ // prop.L2
			nll += neglogdunif(sda, Type(0.0), Type(1000.0), n);
		} else if (varfunc==4){ // exp.dL
			nll += neglogdunif(sda, Type(0.0), Type(1000.0), n);
			nll += neglogdunif(sdb, Type(0.0), Type(1000.0), n);
		} else if (varfunc==5){ // exp.L2
			nll += neglogdunif(sda, Type(0.0), Type(1000.0), n);
			nll += neglogdunif(sdb, Type(0.0), Type(1000.0), n);
		} else if (varfunc==6){ // pow.dL
			nll += neglogdunif(sda, Type(0.0), Type(1000.0), n);
			nll += neglogdunif(sdb, Type(0.0), Type(1000.0), n);
		} else if (varfunc==7){ // pow.L2
			nll += neglogdunif(sda, Type(0.0), Type(1000.0), n);
			nll += neglogdunif(sdb, Type(0.0), Type(1000.0), n);
		}
	}
	// ^ all hard-coded uniform are fairly uninformative



	//--------------------------------------------------------------------------
	// Random effects
	//--------------------------------------------------------------------------

	// none
	

	//--------------------------------------------------------------------------
	// Observation equations
	//--------------------------------------------------------------------------

	// different obs likelihood depending on variance function specification
	if (varfunc==0){ // constant
		for (int i = 0; i < n; i++){
			Type meanL2 = Linf-(Linf-L1(i))*exp(-K*deltaT(i)); // vB at recap
			nll -= dnorm(L2(i), meanL2, sigma, true); // Gaussian lkhd
			// ^ vB reference point = cap
		}
	} else if (varfunc==1){ // prop.dL
		for (int i = 0; i < n; i++){
			Type meandL = (Linf-L1(i))*(1.0-exp(-K*deltaT(i))); // vB E[deltaL]
			// ^ E[deltaL] = E[L2] - L1 since L1 fixed
			nll -= dnorm(L2(i)-L1(i), meandL, sigma+sda*meandL, true);
			// ^ Gaussian lkhd, vB reference point = cap
		}
	} else if (varfunc==2){ // prop.dT
		for (int i = 0; i < n; i++){
			Type meandL = (Linf-L1(i))*(1.0-exp(-K*deltaT(i))); // vB E[deltaL]
			// ^ E[deltaL] = E[L2] - L1 since L1 fixed
			nll -= dnorm(L2(i)-L1(i), meandL, sigma+sda*deltaT(i), true);
			// ^ Gaussian lkhd, vB reference point = cap
		}
	} else if (varfunc==3){ // prop.L2
		for (int i = 0; i < n; i++){
			Type meanL2 = Linf-(Linf-L1(i))*exp(-K*deltaT(i)); // vB at recap
			nll -= dnorm(L2(i), meanL2, sigma+sda*meanL2, true);
			// ^ Gaussian lkhd, vB reference point = cap
		}
	} else if (varfunc==4){ // exp.dL
		for (int i = 0; i < n; i++){
			Type meandL = (Linf-L1(i))*(1.0-exp(-K*deltaT(i))); // vB E[deltaL]
			// ^ E[deltaL] = E[L2] - L1 since L1 fixed
			Type sd_expdL = sigma + sdb*(1.0-exp(-sda*meandL));
			nll -= dnorm(L2(i)-L1(i), meandL, sd_expdL, true);
			// ^ Gaussian lkhd, vB reference point = cap
		}
	} else if (varfunc==5){ // exp.L2
		for (int i = 0; i < n; i++){
			Type meanL2 = Linf-(Linf-L1(i))*exp(-K*deltaT(i)); // vB at recap
			Type sd_expL2 = sigma + sdb*(1.0-exp(-sda*meanL2));
			nll -= dnorm(L2(i), meanL2, sd_expL2, true);
			// ^ Gaussian lkhd, vB reference point = cap
		}
	} else if (varfunc==6){ // pow.dL
		for (int i = 0; i < n; i++){
			Type meandL = (Linf-L1(i))*(1.0-exp(-K*deltaT(i))); // vB E[deltaL]
			// ^ E[deltaL] = E[L2] - L1 since L1 fixed
			nll -= dnorm(L2(i)-L1(i), meandL, sigma+sda*pow(meandL,sdb), true);
			// ^ Gaussian lkhd, vB reference point = cap
		}
	} else if (varfunc==7){ // pow.L2
		for (int i = 0; i < n; i++){
			Type meanL2 = Linf-(Linf-L1(i))*exp(-K*deltaT(i)); // vB at recap
			nll -= dnorm(L2(i), meanL2, sigma+sda*pow(meanL2,sdb), true);
			// ^ Gaussian lkhd, vB reference point = cap
		}
	}


	//--------------------------------------------------------------------------
	// Outputs
	//--------------------------------------------------------------------------

	REPORT(Linf);
	REPORT(K);
	REPORT(sigma);
	REPORT(sda);
	REPORT(sdb);

	ADREPORT(Linf);
	ADREPORT(K);
	// ADREPORT(sigma); // for later

	return nll;
}
