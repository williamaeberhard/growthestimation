GrowthEstimation: R code for various estimation methods of von Bertalanffy growth parameters 
--------------------------------------------------------------------------------------------

These R scripts provide functions for estimating the von Bertalanffy growth parameters (Linf and K) from mark-recapture data. These allow to replicate the results in the simulation study of Dureuil et al. (2020+), including the CapRecapSim R function to simulate mark-recapture data under various scenarios.

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

This is GrowthEstimation version 0.3. Changelog since last version:
* for all Bayesian estimation methods (zh09, Bfa65, and Bla02): now rely on tmbstan
* for all estimation methods: now supply deltaT directly rather than T1 and T2 (they were never used individually anyway).
* CapRecapSim:
  - changed default value of sd.L from 0.3 to 0.2, better if simulated lengths at capture are slightly more concentrated
  - deleted output Tcap and Trecap, no use considering dates if true ages are expressed wrt birth = 0.
  - new mandatory argument Lbirth, so that we do not use Pauly (1979) equation for T0 anymore but rather use von Bertlanffy at birth to infer T0 from supplied Linf and K

Changelog v0.2.2:
* Created GrowthEstimation_Tests.r where a LRT is defined to compare two populations in terms of (Linf,K) jointly, according to Fabens (1965) formulation. Added code at the end of GrowthEstimation_Main.r shows the usage on some simulated data.
* In Bfa65: now allow for different families for priors on Linf, K, and sigma.
* In Bfa65: for priordist="lognormal", user now provides mean and sd on original (exponential) scale, easier for interpretation.


### References

Dureuil, M., Aeberhard, W. H., Dowd, M., Froese, R., Pardo, S. A., Whoriskey, F. G., and Worm, B. (2020+) Reliable growth estimation from mark–recapture tagging data in elasmobranchs. In preparation


