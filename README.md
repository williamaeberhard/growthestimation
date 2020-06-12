GrowthEstimation: R code for various estimation methods of von Bertalanffy growth parameters 
--------------------------------------------------------------------------------------------

These R scripts provide functions for estimating the von Bertalanffy growth parameters (Linf and K) from mark-recapture data. These allow to replicate the results in the simulation study of Dureuil et al. (2020+), including the CapRecapSim R function to simulate mark-recapture data under various scenarios.

### Contents

Files contained in this repository:
* GrowthEstimation_Main.r: the main R script to run, loading all required libraries, compiling the Template Model Builder (TMB) C++ templates, simulating data and estimating the von Bertalanffy growth parameters according to various methods as in Dureuil et al. (2020+);
* GrowthEstimation_Methods.r: an R script to source to create functions corresponding to all the estimation methods under study;
* GrowthEstimation_CapRecapSim.r: an R script to source to create the CapRecapSim function that simulated mark-recapture data according to various scenarios;
* FabensBayesian.cpp, Laslett.cpp, Zhang.cpp: C++ scripts to compile with the R package TMB.
* this README file.

### Version History

This is GrowthEstimation version 0.2.1. Changelog since previous version:
* In Bfa65: now allow for different families for priors on Linf, K, and sigma.
* In Bfa65: for priordist="lognormal", user now provides mean and sd on original (exponential) scale, easier for interpretation.


### References

Dureuil, M., Aeberhard, W. H., Dowd, M., Froese, R., Pardo, S. A., Whoriskey, F. G., and Worm, B. (2020+) Re-evaluating methods to estimate growth from tagging data in elasmobranchs. In preparation


