#///////////////////////////////////////////////////////////////////
#### GrowthEstimation: compare methods on simulated data v0.2.2 ####
#///////////////////////////////////////////////////////////////////

# rm(list=ls())

### setup
require(TMB) # needed for zh09, la02, Bfa65 and Bla02
require(adnuts) # needed for Bfa65, Bla02 and zh09

compile("FabensBayesian.cpp") # only need to run once
compile("Laslett.cpp") # only need to run once
compile("Zhang.cpp") # only need to run once

dyn.load(dynlib("FabensBayesian")) # to run for every new R session
dyn.load(dynlib("Laslett")) # to run for every new R session
dyn.load(dynlib("Zhang")) # to run for every new R session

source('GrowthEstimation_CapRecapSim.r') # simulates lengths at cap and recap
source('GrowthEstimation_Methods.r') # loads TMB package and creates functions

# Available functions are:
#  - gh59: Gulland and Holt (1959)
#  - fa65: Fabens (1965)
#  - fr88: Francis (1988)
#  - ja91: James (1991)
#  - la02: Laslett, Eveson and Polacheck (2002)
#  - zh09: Zhang, Lessard and Campbell (2009), Bayesian
#  - Bfa65: Bayesian version of Fabens (1965)
#  - Bla02: Bayesian version of Laslett, Eveson and Polacheck (2002)

# All methods take the exact same mandatory arguments:
#  - par: vector of length 2, starting values for Linf and K (in that order)
#  - L1: vector of lengths at cap
#  - L2: vector of lengths at recap
#  - T1: vector of dates at cap
#  - T2: vector of dates at recap

# In addition, the 2 Bayesian methods zh09 and Bla02 require:
#  - hyperpar: a list of two numeric vectors, both of length 2, containing the
#              lower and upper bounds (in that order) for uniform priors on Linf
#              and K.
# while the Bayesian method Bfa65 requires:
#  - priordist: a character string specifying the prior distributions on Linf, K,
#               and sigma (the sd of the iid Gaussian additive error on deltaL),
#               choice among "uniform", "normal" (equivalent to "Gaussian"), and
#               "lognormal"
#  - hyperpar: a list of three numeric vectors, each of length 2, containing the
#              hyperparameters of the prior distribution on Linf, K, and sigma,
#              respectively; if priordist="uniform" then these are the lower and
#              upper bounds, if priordist="normal" then these are the mean and
#              sd, and if priordist="lognormal" then these are the mean and sd
#              on the log scale.

# The methods return different arguments, but they all have:
#  - $par: the vector of estimates of Linf and K (in that order)
#  - $se: the vector of standard errors of Linf and K (in that order).


### design, simulate data
n <- 100 # sample size, nb of capture-recapture pairs
trueLinf <- 123.5 # based on starry smooth-hound Mustelus asterias
trueK <- 0.146 # based on starry smooth-hound Mustelus asterias
theta <- c(trueLinf,trueK) # parameters of interest
length.theta <- length(theta)
names.theta <- c('Linf','K')
names(theta) <- names.theta

par.ini <- c(100,0.3) # initial values for Linf and K, somewhat in ballpark

set.seed(1234) # for replicability
dat <- CapRecapSim(n=n,trueLinf=trueLinf,trueK=trueK,sizedist='norm')

str(dat) # what really matters is the "realistic" data: Lcap, Lrecap, Tcap and Trecap
# All dates are in years (numeric), with year 0 being "1970-01-01", R's default

### for a given method, say fa65, study impact of gv and measurement error
cbind(c(trueLinf,trueK),
      fa65(par=par.ini,L1=dat$trueLcap,L2=dat$trueLrecap, # no gv, no meas err
           T1=dat$Tcap,T2=dat$Trecap,compute.se=T)$par,
      fa65(par=par.ini,L1=dat$gvLcap,L2=dat$gvLrecap, # gv, no meas err
           T1=dat$Tcap,T2=dat$Trecap,compute.se=T)$par,
      fa65(par=par.ini,L1=dat$eLcap,L2=dat$eLrecap, # no gv, meas err
           T1=dat$Tcap,T2=dat$Trecap,compute.se=T)$par,
      fa65(par=par.ini,L1=dat$Lcap,L2=dat$Lrecap, # both gv and meas err
           T1=dat$Tcap,T2=dat$Trecap,compute.se=T)$par
)
# 1st col = true values, 2nd col = est under no gv and no meas err, etc.


### Compare estimates of (mean) Linf and K among methods
names.est <- c('gh59','fa65','fr88','ja91','la02','zh09','Bfa65','Bla02')
n.est <- length(names.est) # number of estimators to compare
est <- vector('list',n.est)

hp.unif <- list(c(0,500),c(0,2)) # hyperparam: lb and ub for Linf and K

mcmc.control <- list('nchains'=1,'iter'=1000,'warmup'=200)
# ^ all should be larger, just for the sake of the demo here

# L1 <- dat$trueLcap # pure vB, no gv or meas err
# L2 <- dat$trueLrecap # pure vB, no gv or meas err
# L1 <- dat$gvLcap # gv, but no meas err
# L2 <- dat$gvLrecap # gv, but no meas err
L1 <- dat$eLcap # no gv, but meas err
L2 <- dat$eLrecap # no gv, but meas err
# L1 <- dat$Lcap # with both gv and meas err
# L2 <- dat$Lrecap # with both gv and meas err
# ^ compare all methods on different sets of data

system.time(est[[1]] <- gh59(par=par.ini,L1=L1,L2=L2,T1=dat$Tcap,T2=dat$Trecap))
system.time(est[[2]] <- fa65(par=par.ini,L1=L1,L2=L2,T1=dat$Tcap,T2=dat$Trecap))
system.time(est[[3]] <- fr88(par=par.ini,L1=L1,L2=L2,T1=dat$Tcap,T2=dat$Trecap))
system.time(est[[4]] <- ja91(par=par.ini,L1=L1,L2=L2,T1=dat$Tcap,T2=dat$Trecap))
system.time(est[[5]] <- la02(par=par.ini,L1=L1,L2=L2,T1=dat$Tcap,T2=dat$Trecap,
                             lb.sd=1e-5)) # necessary for se with pure vB data
# ^ slower because of random effects
system.time(est[[6]] <- zh09(par=par.ini,L1=L1,L2=L2,T1=dat$Tcap,T2=dat$Trecap,
                             hyperpar=hp.unif,lb.sd=1e-5,
                             mcmc.control=mcmc.control))
# ^ even slower because of random effects and MCMC
system.time(est[[7]] <- Bfa65(par=par.ini,L1=L1,L2=L2,T1=dat$Tcap,T2=dat$Trecap,
                              priordist.Linf='uniform',
                              priordist.K='uniform',
                              priordist.sigma='uniform',
                              hyperpar=c(hp.unif,list(c(0,100))),
                              mcmc.control=mcmc.control))
# ^ somewhat slow because of MCMC, but no random effects
system.time(est[[8]] <- Bla02(par=par.ini,L1=L1,L2=L2,
                              T1=dat$Tcap,T2=dat$Trecap,
                              hyperpar=hp.unif,lb.sd=1e-5,
                              mcmc.control=mcmc.control))
# ^ also slow because of random effects and MCMC

est.compare <- c(trueLinf,trueK)
for (j in 1:n.est){ # loop to look through all methods
  est.compare <- cbind(est.compare,est[[j]]$par)
}
dimnames(est.compare) <- list(c('Linf','K'),c('true',names.est[1:n.est]))
round(est.compare,4)


### plot estimate and 95% CI as error bar
# Note: for Bayesian methods, more appropriate to compute credible intervals
#       from posterior draws, here CI just for the sake of the comparison

est.Linf <- unlist(lapply(est,FUN='[[','par'))[2*(1:n.est)-1]
est.K <- unlist(lapply(est,FUN='[[','par'))[2*(1:n.est)]

se.Linf <- unlist(lapply(est,FUN='[[','se'))[2*(1:n.est)-1]
se.K <- unlist(lapply(est,FUN='[[','se'))[2*(1:n.est)]

se.compare <- rbind(se.Linf,se.K)
dimnames(se.compare)[[2]] <- names.est
signif(se.compare,4)
# ^ sd.lb necessary for la02 under pure vB data, otherwise Hessian not pos def

par(mfrow=c(2,1))
# Linf
plot(1:n.est,est.Linf,pch=19,cex=0.5,ylim=c(50,200),
     xaxt='n',xlab='',ylab='',main='Linf')
axis(1,1:n.est,names.est[1:n.est])
arrows(x0=1:n.est,x1=1:n.est,angle=90,code=3,length=0.1,
       y0=est.Linf-2*se.Linf,y1=est.Linf+2*se.Linf)
abline(h=trueLinf,col='red')
# K
plot(1:n.est,est.K,pch=19,cex=0.5,ylim=c(0,1),
     xaxt='n',xlab='',ylab='',main='K')
axis(1,1:n.est,names.est[1:n.est])
arrows(x0=1:n.est,x1=1:n.est,angle=90,code=3,length=0.1,
       y0=est.K-2*se.K,y1=est.K+2*se.K)
abline(h=trueK,col='red')
par(mfrow=c(1,1))


#////////////////////////////////////////////////////////////////////////////
#### GrowthEstimation: test compare two population in terms of vB curves ####
#////////////////////////////////////////////////////////////////////////////

# rm(list=ls())

### setup
require(TMB) # needed for zh09, la02, Bfa65 and Bla02

compile("FabensTwoPop.cpp") # only need to run once
dyn.load(dynlib("FabensTwoPop")) # to run for every new R session

source('GrowthEstimation_Tests.r')


### sim data under exact Fabens formulation, but different pop1/pop2
n <- 100 # total sample size, nb of capture-recapture pairs
n1 <- 50 # sample size pop1
n2 <- n-n1 # sample size pop2

trueLinf1 <- 125 # close to starry smooth-hound Mustelus asterias
trueK1 <- 0.145 # close to starry smooth-hound Mustelus asterias
trueLinf2 <- 120 # pop2, slightly smaller than that of pop1
trueK2 <- 0.145 # pop2, same as pop1

sigmaeps1 <- 2 # errror sd, pop1
sigmaeps2 <- 2 # errror sd, pop2

par.ini <- c(100,0.2) # initial values for Linf and K, somewhat in ballpark

set.seed(1234) # for replicability

Lcap1 <- rnorm(n=n1,mean=60,sd=5) # pop1
Lcap2 <- rnorm(n=n2,mean=60,sd=5) # pop2

deltaT1 <- rgamma(n=n1,shape=1,scale=1) # in years
deltaT1 <- ifelse(deltaT1>3,3,deltaT1) # truncate at 3
deltaT2 <- rgamma(n=n2,shape=1,scale=1) # in years
deltaT2 <- ifelse(deltaT2>3,3,deltaT2) # truncate at 3

eps1 <- rnorm(n=n1,mean=0,sd=sigmaeps1) # pop1
eps2 <- rnorm(n=n2,mean=0,sd=sigmaeps2) # pop2

Lrecap1 <- trueLinf1 - (trueLinf1-Lcap1)*exp(-trueK1*deltaT1) + eps1 # pop 1
Lrecap2 <- trueLinf2 - (trueLinf2-Lcap2)*exp(-trueK2*deltaT2) + eps2 # pop 2
# ^ both follow the vB curve with some iid Gaussian error on lengths at recap


### compute likelihood ratio test (LRT) to compare the two pop
LRT_2pop_fa65(par=par.ini,alpha=0.05,
  L1.pop1=Lcap1,L2.pop1=Lrecap1,T1.pop1=rep(0,n1),T2.pop1=deltaT1,
  L1.pop2=Lcap2,L2.pop2=Lrecap2,T1.pop2=rep(0,n2),T2.pop2=deltaT2)
# ^ test stat > chi^2 critical value or equivalently p-value < alpha
#   => reject H0 at significance level alpha=0.05, the pair (Linf,K) differs
#      between the two pop



# END GrowthEstimation_Main
