#///////////////////////////////////////////////////////////////////
#### GrowthEstimation: compare methods on simulated data v0.4.4 ####
#///////////////////////////////////////////////////////////////////

# rm(list=ls())

### // setup ----
require(TMB)     # needed for Bayesian methods
require(tmbstan) # needed for Bayesian methods
# ^ for tmbstan function for MCMC on TMB obj, calls rstan

compile("FabensBayesian.cpp") # only need to run once
compile("FrancisBayesian.cpp") # only need to run once

dyn.load(dynlib("FabensBayesian")) # to run for every new R session
dyn.load(dynlib("FrancisBayesian")) # to run for every new R session

source('GrowthEstimation_CapRecapSim.r') # simulates lengths at cap and recap
source('GrowthEstimation_Methods.r') # loads TMB package and creates functions

# Available estimation methods are:
# * gh59: Gulland and Holt (1959)
# * fa65: Fabens (1965), least squares estimation
# * fr88: Francis (1988)
# * fr88.minAIC: wrapper for fr88 which tests sub-models corresponding to
#   combinations of the extra parameters nu, m, and p, and returns the best
#   sub-model in terms of smallest AIC
# * ja91: James (1991)
# * Bfa65: Bayesian version of Fabens (1965), with Gaussian likelihood
# * Bfr88: Bayesian version of Francis (1988), although without the extra
#   parameters m and p and with many variance functions to choose from
# * Bfr88.minIC: wrapper for Bfr88 which tests sub-models corresponding to
#   different variance functions, and returns the best sub-model in terms of
#   smallest WAIC

# All methods share the following mandatory arguments:
# * par: vector of length 2, starting values for Linf and K (in that order)
# * L1: vector of measured lengths at capture
# * L2: vector of measured lengths at recapture
# * deltaT: vector of time intervals between capture and recapture
# L1 and L2 must be in the same units, these implies the units of the Linf
# estimates. Similarly, the time units of deltaT imly the units of the K
# estimates, e.g. deltaT in years => K in years^{-1}

# In addition, the Bayesian methods Bfa65 and Bfr88 and their wrappers require:
# * priordist.Linf: a character string specifying the prior distributions on
#                   Linf, choice among "uniform", "normal" (equivalent to
#                   "Gaussian"), and "lognormal" (the default)
# * priordist.K: character string specifying the prior distributions on K, same
#                choice as for Linf with default "uniform"
# * priordist.sigma: character string specifying the prior distributions on
#                    sigma (the sd of the iid Gaussian error) same choice as for
#                    Linf with default "uniform"
# * hyperpar: a list of three numeric vectors, each of length 2, containing the
#             hyperparameters of the prior distribution on Linf, K, and sigma,
#             in that order. If priordist="uniform" then these are the lower and
#             upper bounds, if priordist="normal" then these are the mean and
#             sd, and if priordist="lognormal" then these are the mean and sd
#             on the original (exponential) scale.
# * mcmc.control: a list of three named scalars. nchains is the number of chains
#                 to run (potentially in parallel); iter is the total number of
#                 iterations, including burn-in; warmup is the number of
#                 iterations in the burn-in phase. Default is list(nchains=5,
#                 iter=20000, warmup=10000).

# Also, some specific arguments:
# * fr88 has three optional arguments par.nu, par.m, and par.p, which are
#   all NULL by default. Leaving any of these arguments NULL "disables" the
#   corresponding extra component and parameter. Supplying a numeric value to
#   any of these arguments "enables" the component and estimates the
#   corresponding  parameter, using the supplied numeric value as starting value
#   in the estimation. This means that by default (all three NULL), fr88 is a
#   maximum likelihood estimation version of fa65.
# * fr88.minAIC requires numeric values for par.nu, par.m, and par.p (NULL not
#   allowed), used as starting values in all the tested sub-models.
# * Bfr88 has the mandatory argument varfunc: this is a character string
#   specifying the variance function, choice among "constant" (equivalent to
#   Bfa65), "prop.dL" (standard Francis, 1988), "prop.dT" (proportional to
#   deltaT), "prop.L2" (proportional to lengths at recapture), "exp.dL" (eq. (7)
#   in Francis, 1988), "exp.L2" (alternative to exp.dL), "pow.dL" (eq. (8) in
#   Francis, 1988)), or "pow.L2" (alternative to pow.dL).
# * Bfr88.minIC also has the varfunc argument although here it is a character
#   vector of all variance functions to be tested. Any subset of the eight
#   listed above works. The default is varfunc="all" which includes all eight.

# The methods return different arguments, but they all have in common:
#  - $par: the vector of (point) estimates of Linf and K (in that order)
#  - $se: the vector of standard errors of Linf and K (in that order).


### // design, simulate data ----
n <- 100 # sample size, nb of capture-recapture pairs
trueLinf <- 123.5 # based on starry smoothhound Mustelus asterias
trueK <- 0.146 # based on starry smoothhound Mustelus asterias
theta <- c(trueLinf,trueK) # parameters of interest
length.theta <- length(theta)
names.theta <- c('Linf','K')
names(theta) <- names.theta

par.ini <- c(100, 0.3) # initial values for Linf and K, somewhat in ballpark

set.seed(1234) # for replicability
dat <- CapRecapSim(n=n,trueLinf=trueLinf,trueK=trueK,sizedist='norm')

str(dat)
# ^ what really matters is the "realistic" data: Lcap, Lrecap, and deltaT, which
#   include growth variability (GV) in terms of a Gaussian random intercept on
#   Linf and independent Gaussian measurement error (ME) on lengths
# All ages and time intervals are in years, so that K/Ki expressed in year^{-1}


### // for a given method, say fa65, see impact of GV and ME ----
cbind(c(trueLinf,trueK),
      fa65(par=par.ini,L1=dat$trueLcap,L2=dat$trueLrecap, # no GV, no ME
           deltaT=dat$truedeltaT,compute.se=F)$par,
      fa65(par=par.ini,L1=dat$gvLcap,L2=dat$gvLrecap, # GV, no ME
           deltaT=dat$deltaT,compute.se=F)$par,
      fa65(par=par.ini,L1=dat$eLcap,L2=dat$eLrecap, # no GV, ME
           deltaT=dat$truedeltaT,compute.se=F)$par,
      fa65(par=par.ini,L1=dat$Lcap,L2=dat$Lrecap, # both GV and ME
           deltaT=dat$deltaT,compute.se=F)$par
)
# 1st col = true values, 2nd col = est under no gv and no ME, etc.



### // compare estimates of (mean) Linf and K among methods ----
names.est <- c('gh59','fa65','fr88','ja91','Bfa65','Bfr88')
n.est <- length(names.est) # number of estimators to compare
est <- vector('list',n.est)

hplist <- list(c(120,40), # hyperpar for Linf: lognormal mean and sd (exp scale)
               c(0,2),    # hyperpar for K: uniform lower and upper bound
               c(0,50)    # hyperpar for sigma: uniform lower and upper bound
)
# ^ arbitrary reasonable values supplied here for the sake of the demo

Linf.meanlog <- 2*log(hplist[[1]][1])-log(hplist[[1]][2]^2+hplist[[1]][1]^2)/2
Linf.sdlog <- sqrt(-2*log(hplist[[1]][1]) +
                     + log(hplist[[1]][2]^2+hplist[[1]][1]^2))
# ^ lognormal mean and sd on log scale

plot(seq(0.5,300,0.5),
     dlnorm(seq(0.5,300,0.5),meanlog=Linf.meanlog,sdlog=Linf.sdlog),type='l')
abline(v=trueLinf,col='blue')
# ^ lognormal prior density on Linf and true value

mcmc.control <- list('nchains'=3,'iter'=4000,'warmup'=2000,
                     'adapt_delta'=0.8,'max_treedepth'=20)
# ^ nchains, iter, and warmup should be larger, just for the sake of the demo

# ## exact von Bertalanffy, no GV or ME
# L1 <- dat$trueLcap
# L2 <- dat$trueLrecap
# deltaT <- dat$truedeltaT
# ## GV, but no ME
# L1 <- dat$gvLcap
# L2 <- dat$gvLrecap
# deltaT <- dat$deltaT # deltaT modified with GV since neg correlated with L1
# ## no GV, but ME
# L1 <- dat$eLcap
# L2 <- dat$eLrecap
# deltaT <- dat$truedeltaT # with only ME on lengths, ages remain the same
##  with both GV and ME
L1 <- dat$Lcap
L2 <- dat$Lrecap
deltaT <- dat$deltaT # deltaT modified with GV since neg correlated with L1

# ^ different sets of data

system.time(est[[1]] <- gh59(par=par.ini,L1=L1,L2=L2,deltaT=deltaT))
system.time(est[[2]] <- fa65(par=par.ini,L1=L1,L2=L2,deltaT=deltaT))
system.time(est[[3]] <- fr88.minAIC(par=par.ini,L1=L1,L2=L2,deltaT=deltaT,
                                    par.nu=0.1,par.m=5,par.p=0.1))
system.time(est[[4]] <- ja91(par=par.ini,L1=L1,L2=L2,deltaT=deltaT))
system.time(est[[5]] <- Bfa65(par=par.ini,L1=L1,L2=L2,deltaT=deltaT,
                              priordist.Linf='lognormal',
                              priordist.K='uniform',
                              priordist.sigma='uniform',
                              hyperpar=hplist,
                              mcmc.control=mcmc.control))
# ^ takes some time because of MCMC
system.time(est[[6]] <- Bfr88.minIC(par=par.ini,L1=L1,L2=L2,deltaT=deltaT,
                                    priordist.Linf='lognormal',
                                    priordist.K='uniform',
                                    priordist.sigma='uniform',
                                    hyperpar=hplist,
                                    varfunc=c('constant'),
                                    mcmc.control=mcmc.control))
# ^ takes more time because of MCMC and fitting multiple sub-models

est.compare <- c(trueLinf,trueK)
for (j in 1:n.est){ # loop to look through all methods
  est.compare <- cbind(est.compare,est[[j]]$par)
}
dimnames(est.compare) <- list(c('Linf','K'),c('true',names.est[1:n.est]))
round(est.compare,4)


### // lognormal prior hyperparameters from likely value and upper bound ----
# rather than using GrowthPriors function, which relies on fishbase, compute
# lognormal hyperparameter from a given "most likely" value (considered as the 
# prior median) and a realistic upper bound (used as the prior 0.99 quantile)

hp.Linf <- hp.lognormal(median=120,upper.bound=200,plot=T)
# ^ most likely value and upper bound, assumed from species knowledge
hp.Linf
# ^ lognormal hyperparameters: mean and sd on original (exp) scale, to supply to
#   Bayesian estimation functions Bfa65, Bfr88, and Bfr88.minIC

hplist2 <- list(hp.Linf, # hyperpar for Linf: lognormal mean and sd (exp scale)
               c(0,2),   # hyperpar for K: uniform lower and upper bound
               c(0,50)   # hyperpar for sigma: uniform lower and upper bound
)

est5.alt <- Bfa65(par=par.ini,L1=L1,L2=L2,deltaT=deltaT,
                  priordist.Linf='lognormal',
                  priordist.K='uniform',
                  priordist.sigma='uniform',
                  hyperpar=hplist2,
                  mcmc.control=mcmc.control)

rbind(est5.alt$par, est[[5]]$par)
# ^ estimates (posterior medians) slightly different because different
#   hyperparameters for Linf lognormal prior

rbind(est5.alt$cred.int$Linf, est[[5]]$cred.int$Linf)
# ^ credible intervals for Linf also slightly different




# END GrowthEstimation_Main
