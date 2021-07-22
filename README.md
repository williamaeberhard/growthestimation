GrowthEstimation: R code for various estimation methods of von Bertalanffy growth parameters 
--------------------------------------------------------------------------------------------

These R scripts provide functions for estimating the von Bertalanffy growth parameters (Linf and K) from mark-recapture data. These allow to replicate the results in the simulation study of Dureuil et al. (2021+), including the CapRecapSim R function to simulate mark-recapture data under various scenarios.


### Contents

Files contained in this repository:
* GrowthEstimation_Main.r: the main R script to run, loading all required libraries, compiling the Template Model Builder (TMB) C++ templates, simulating data and estimating the von Bertalanffy growth parameters according to various methods as in Dureuil et al. (2020+);
* GrowthEstimation_Methods.r: an R script to source to create functions corresponding to all the estimation methods under study;
* GrowthEstimation_CapRecapSim.r: an R script to source to create the CapRecapSim function that simulated mark-recapture data according to various scenarios;
* FabensBayesian.cpp, Laslett.cpp, Zhang.cpp: C++ scripts to compile with the R package TMB;
* GrowthEstimation_Tests.r: an R script defining a likelihood ratio test (LRT) and a (Laplace-approximated) Bayes factor for comparing two independent populations in terms of their (Linf,K) parameters according to the (frequentist) Fabens (1965) method and its Bayesian formulation;
* FabensTwoPop.cpp: a C++ script to compile for the LRT comparing two populations;
* FabensTwoPopBayesian_M0.cpp and FabensTwoPopBayesian_M1.cpp: C++ scripts to compile for the Bayes factor comparing two populations;
* this README file.


### Version History

This is GrowthEstimation version 0.4.4. Changelog since last version:
* Adapted GrowthPriors for computation of Linf lognormal prior hyperparameters from supplied Lmax. New mandatory argument is lnorm.coef, a vector of length 2 with default c(1.2, 3), that defines narrow and wide lognormal priors. The lognormal mean parameter on the log scale is computed from the median (meanlog = log(median)), the latter set as Lmax/0.99 (the value for Linf reported so far). The lognormal sd parameter (also on log scale) is such that the 0.99 quantile macthes a given value, this being lnorm.coef[1]*median for the narrow case and lnorm.coef[2]*median for the wide case. The lognormal mean and sd on the original (exponential) scale are then computed from the parameters on the log scale. Additional outputs thus are:
  - lnorm.meanLinf.narrow: narrow lognormal distribution mean (original scale)
  - lnorm.sdLinf.narrow: narrow lognormal distribution sd (original scale)
  - lnorm.meanLinf.wide: wide lognormal distribution mean (original scale)
  - lnorm.sdLinf.wide: wide lognormal distribution sd (original scale)
* Removed all statistical tests (LRT and BF) stuff, since not pertaining to Dureuil et al. (2021+) paper. Thus deleted the following files:
  - GrowthEstimation_Tests.r
  - FabensTwoPop.cpp
  - FabensTwoPopBayesian_M0.cpp
  - FabensTwoPopBayesian_M1.cpp
* Removed methods that are not used in Dureuil et al. (2021+) paper in the end, i.e. la02, Bla02, zh09, and oldfr88. Thus deleted the following files:
  - Laslett.cpp
  - Zhang.cpp

Changelog v0.4.3:
* LRT_2pop_fa65:
  - now user supplies deltaT (=Trecap-Tcap) rather than the two separate (absolute) times as separate vector. FabensTwoPop.cpp modified too.
  - default values in case par not supplied changed from (0,0) to (1,0.5) and issuing a warning now.
* In GrowthEstimation_Tests.r, created BF_2pop_Bfa65: Bayesian analogue of LRT_2pop_fa65, it computes the (Laplace-approximated) Bayes factor between two competing models M1 and M0, where M0 sets the same Linf and K values between the two populations we compare (allowing different error sd parameters though) and M1 allows Linf and K to be different (with possibly different priors too).
* Adapted lines at the bottom of GrowthEstimation_Main.r to account for changes in LRT_2pop_fa65 and new function BF_2pop_Bfa65.

Changelog v0.4.2:
* Bfa65 and Bfr88: added DIC (both p_D and p_V versions) and WAIC as output.
* new wrapper function Bfr88.minIC: fits multiple Bfr88 models with same priors but different variance functions, and returns the best model in terms of min WAIC (if enablepriorsd=TRUE, all priors for Linf, K and sigma are specified, and onlyTMB=FALSE) or min AIC (if not all three conditions are met) along with character string specifying which sub-model it is (output $best.model).

Changelog v0.4.1:
* added new function Bfr88, our (Bayesian) take on Francis (1988).

Changelog v0.4:
* Francis (1988) estimation method fr88:
  - completely re-coded, now follows closely orginal specification in Francis (1988), notably the optional estimation of nu, m, and p
  - old fr88 function renamed oldfr88 for comparison, will likely be deprecated in future versions
  - new fr88 requires three new arguments par.nu, par.m, and par.p which are user-supplied initial values; if left NULL (default), then corresponding "component" disabled
  - new output AIC for model selection
  - new wrapper function fr88.minAIC: fits multiple fr88 models with increasing complexity (adding one param at a time among nu, m, and p) and returns the best model in terms of min AIC along with character string specifying which sub-model it is (output $best.model)

Changelog v0.3.2:
* Francis (1988) estimation method fr88:
  - estimation of nu param now on log scale
  - added my own initial values for g.alpha and g.beta, as good as multiple subsets yet faster than ini from fishmethods::grotag
  - many variance functions specification available through new sdfunc argument (with possible values "const", "prop.L2", "prop.dL", and "prop.dT").

Changelog v0.3.1:
* Bfa65: point estimate from MCMC is now posterior median (posterior dist of both Linf and K can be very skewed if small sample size, and very similar to posterior mean/mode if large sample) and point estimate of spread is now median absolute deviation about the median (MADAM, numerically most stable in small samples)
* all Bayesian estimation methods (zh09, Bfa65, and Bla02): default MCMC options are now 'nchains'=5, 'iter'=20000, and 'warmup'=10000.
* CapRecapSim: default values set now for low growth variability level (sd.Linf=0.017) and short times at liberty (scale.deltaT=3.5).

Changelog v0.3:
* for all Bayesian estimation methods (zh09, Bfa65, and Bla02): now rely on tmbstan
* for all estimation methods: now supply deltaT directly rather than T1 and T2 (they were never used individually anyway).
* CapRecapSim:
  - changed default value of sd.L from 0.3 to 0.2, better if simulated lengths at capture are slightly more concentrated
  - deleted output Tcap and Trecap, no use considering dates if true ages are expressed wrt birth = 0.
  - new mandatory argument Lbirth, so that we do not use Pauly (1979) equation for T0 anymore but rather use von Bertalanffy at birth to infer T0 from supplied Linf and K
  - use Lmax=0.99*Linf as in GrowthPriors
  - deleted warnings about ages at recapture exceeding 3/K, not useful since we use a max age based on vB at Lmax=0.99*Linf.
* GrowthPriors: drop Froese and Binohlan's (2000) relation between Lmax and Linf, very close to linear anyway, use now Linf=Lmax/0.99.
* Francis (1988) estimation method fr88: added try.many.ini boolean argument, if TRUE then try many starting values (based on data subsets as in fishmethods::grotag) in addition to the supplied one (optim highly sensitive to ini).

Changelog v0.2.2:
* Created GrowthEstimation_Tests.r where a likelihood ratio test (LRT) is implemented in the enw function LRT_2pop_fa65 to compare two populations in terms of (Linf,K) jointly, according to Fabens (1965) formulation. Added code at the end of GrowthEstimation_Main.r shows the usage on some simulated data.
* In Bfa65: now allow for different families for priors on Linf, K, and sigma.
* In Bfa65: for priordist="lognormal", user now provides mean and sd on original (exponential) scale, easier for interpretation.


### References

Dureuil, M., Aeberhard, W. H., Dowd, M., Froese, R., Pardo, S. A., Whoriskey, F. G., and Worm, B. (2021+) Reliable growth estimation from mark–recapture tagging data in elasmobranchs. In preparation


