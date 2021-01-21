GrowthEstimation: R code for various estimation methods of von Bertalanffy growth parameters 
--------------------------------------------------------------------------------------------

These R scripts provide functions for estimating the von Bertalanffy growth parameters (Linf and K) from mark-recapture data. These allow to replicate the results in the simulation study of Dureuil et al. (2021+), including the CapRecapSim R function to simulate mark-recapture data under various scenarios.


### Contents

Files contained in this repository:
* GrowthEstimation_Main.r: the main R script to run, loading all required libraries, compiling the Template Model Builder (TMB) C++ templates, simulating data and estimating the von Bertalanffy growth parameters according to various methods as in Dureuil et al. (2020+);
* GrowthEstimation_Methods.r: an R script to source to create functions corresponding to all the estimation methods under study;
* GrowthEstimation_CapRecapSim.r: an R script to source to create the CapRecapSim function that simulated mark-recapture data according to various scenarios;
* FabensBayesian.cpp, Laslett.cpp, Zhang.cpp: C++ scripts to compile with the R package TMB;
* GrowthEstimation_Tests.r: an R script defining a likelihood ratio test (LRT) for comparing two independent populations in terms of their (Linf,K) parameters according to the (frequentist) Fabens (1965) formulation;
* FabensTwoPop.cpp: a C++ script to compile for the LRT comparing two populations;
* this README file.


### Version History

This is GrowthEstimation version 0.3.1. Changelog since last version:
* Bfa65: point estimate from MCMC is now posterior median (posterior dist of both Linf and K can be very skewed if small sample size, and very similar to posterior mean/mode if large sample) and point estimate of spread is now median absolute deviation about the median (MADAM, numerically most stable in small samples)
* all Bayesian estimation methods (zh09, Bfa65, and Bla02): default MCMC options are now 'nchains'=5, 'iter'=20000, and 'warmup'=10000.
* CapRecapSim: default values set now for low growth variability level (sd.Linf=0.017) and short times at liberty (scale.deltaT=3.5).

Changelog v0.3:
* for all Bayesian estimation methods (zh09, Bfa65, and Bla02): now rely on tmbstan
* for all estimation methods: now supply deltaT directly rather than T1 and T2 (they were never used individually anyway).
* CapRecapSim:
  - changed default value of sd.L from 0.3 to 0.2, better if simulated lengths at capture are slightly more concentrated
  - deleted output Tcap and Trecap, no use considering dates if true ages are expressed wrt birth = 0.
  - new mandatory argument Lbirth, so that we do not use Pauly (1979) equation for T0 anymore but rather use von Bertlanffy at birth to infer T0 from supplied Linf and K
  - use Lmax=0.99*Linf as in GrowthPriors
  - deleted warnings about ages at recapture exceeding 3/K, not useful since we use a max age based on vB at Lmax=0.99*Linf.
* GrowthPriors: drop Froese and Binohlan's (2000) relation between Lmax and Linf, very close to linear anyway, use now Linf=Lmax/0.99.
* Francis (1988) estimation method fr88: added try.many.ini boolean argument, if TRUE then try many starting values (based on data subsets as in fishmethods::grotag) in addition to the supplied one (optim highly sensitive to ini).

Changelog v0.2.2:
* Created GrowthEstimation_Tests.r where a LRT is defined to compare two populations in terms of (Linf,K) jointly, according to Fabens (1965) formulation. Added code at the end of GrowthEstimation_Main.r shows the usage on some simulated data.
* In Bfa65: now allow for different families for priors on Linf, K, and sigma.
* In Bfa65: for priordist="lognormal", user now provides mean and sd on original (exponential) scale, easier for interpretation.


### References

Dureuil, M., Aeberhard, W. H., Dowd, M., Froese, R., Pardo, S. A., Whoriskey, F. G., and Worm, B. (2021+) Reliable growth estimation from mark–recapture tagging data in elasmobranchs. In preparation


