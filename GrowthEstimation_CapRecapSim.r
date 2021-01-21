CapRecapSim <- function(n,trueLinf=123.5,trueK=0.146,Lbirth=30,
                        sizedist='norm',sd.L=0.2,
                        scale.deltaT=3.5,min.deltaT=0.25,
                        sd.Linf=0.017,sd.L1=0.033,sd.L2=0.067){
  
  #/////////////////////////////////////////////////////////////////////////////
  #### Documentation v0.3.1 #### 
  #/////////////////////////////////////////////////////////////////////////////
  
  # TODO: update doc
  #' @title Simulate fish lengths at capture and at recapture according to the
  #' von Bertalanffy (VB) model, with optional process and measurement error.
  #' @description
  #' TODO
  #' @param n positive integer, sample size.
  #' @param trueLinf positive scalar, VB asymptotic maximum length (average
  #' common to all fish). Default value (123.5 cm) is for starry smoothhound
  #' from Farrell et al. (2017).
  #' @param trueK positive scalar, VB rate (average common to all fish). Default
  #' value (0.146 year^(-1)) is for starry smoothhound from Farrell et al.
  #' (2017).
  #' @param Lbirth positive scalar, length at birth in same units as trueLinf.
  #' Default value (30 cm) for starry smoothhound from Farrell et al. (2010).
  #' @param sizedist distibution of the length at capture; choice among
  #' \code{"norm"} (Gaussian, the default), \code{"left"} (left-skewed),
  #' \code{"right"} (left-skewed) and \code{"bimod"} (bimodal, mixture of two
  #' normals).
  #' @param sd.L positive scalar, sd of generated lengths at capture.
  #' @param scale.deltaT positive scalar, scale parameter of gamma distribution
  #' for generating \code{deltaT}. The associated shape parameter is set so that
  #' \code{agemax-agecap} is the (1-0.01/n)th percentile of the gamma
  #' distribution, guaranteeing that generated ages at recapture
  #' (\code{agecap+deltaT}) rarely exceed \code{agemax}, where \code{agemax} is
  #' computed from the VB equation matching Lmax set as 95% of \code{trueLinf}.
  #' For long times at liberty set \code{scale.deltaT=1} (the default); for
  #' short times at liberty set \code{scale.deltaT=2}.
  #' @param min.deltaT positive scalar, minimum \code{deltaT} in years (added to
  #' gamma realizations). The default value (0.25) corresponds to 3 months.
  #' @param sd.Linf positive scalar, process error standard deviation (sd) of
  #' Linf (Gaussian random effect). The default value (0.067) ensures the
  #' individual Linf values will be within 20% of trueLinf.
  #' @param sd.L1 positive scalar, multiplicative coefficient for measurement
  #' error sd of L1 (Gaussian additive measurement noise). If set large and
  #' \code{deltaT} simulated small enough, then observed length at recapture
  #' can be smaller than observed length at capture (negative growth). Default
  #' value (0.033) yields roughly 10% relative error on average.
  #' @param sd.L2 positive scalar, multiplicative coefficient for measurement
  #' error sd of L2 (Gaussian additive measurement noise). If set large and
  #' \code{deltaT} simulated small enough, then observed length at recapture
  #' can be smaller than observed length at capture (negative growth). Default
  #' value (0.067) yields roughly 20% relative error on average.
  #' @details Randomness in: \code{Linf} (Gaussian, simulated); \code{K}
  #' (propagated from \code{Linf} through the slope of the growth performance
  #' function; length at capture \code{sizedist} (simulated); \code{deltaT}
  #' (lognormal, simulated, negatively correlated with age at capture); and
  #' measurement error in observed lengths at capture and recapture (both
  #' Gaussian). TODO: more details.
  #' @return A list with the following objects:
  #' \itemize{
  #'   \item{\code{Lcap}}{numeric vector of length}
  #'   \item{TODO}{TODO.}
  #' }
  
  # help(document) # document::document
  
  
  #/////////////////////////////////////////////////////////////////////////////
  #### Define true parameters and setup #### 
  #/////////////////////////////////////////////////////////////////////////////
  
  # t0 <- -10^(-0.3922-0.2752*log10(trueLinf)-1.038*log10(trueK)) # Pauly (1979)
  # # ^ birth assumed at t=0 so that t0 is negative
  # ^ as of v0.3: don't use Pauly (1979) eq anymore, rather start from Lbirth
  
  ### lengths at birth, max and mean
  # Lbirth <- trueLinf*(1-exp(trueK*t0)) # length at birth
  T0 <- log(1-Lbirth/trueLinf)/trueK # as of v0.3: T0 from supplied Lbirth, vB
  # ^ time unit from supplied trueK
  
  trueLmax <- 0.99*trueLinf # as of v0.3: assumption about observed max length
  trueLmean <- (trueLmax+Lbirth)/2
  
  ### maximum age
  # trueTmax <- 3/trueK
  trueTmax <- -log(1-trueLmax/trueLinf)/trueK + T0
  # ^ inverted vB for trueLmax with T0 as ref point

  ### number of days in a year, for years-days conversion
  ndays <- 365.2425
  
  # ### set date of capture T1 (not random), same for all by default
  # if (is.null(T1)){T1 <- rep(as.Date("1990-01-01"),n)}
  # ^ as of v0.3: deprecated, T1=agecap inferred from vB at cap in STEP 2
  
  
  #/////////////////////////////////////////////////////////////////////////////
  ### STEP 1: sim lengths at cap ####
  #/////////////////////////////////////////////////////////////////////////////
  
  ### distribution for length at cap
  if (sizedist=="left"){
    trueLlarge <- (trueLmax+trueLmean)/2 # peak for left-skewed dist
    L1 <- rnorm(n,mean=trueLlarge^2,sd=trueLlarge^2*sd.L)
    L1 <- ifelse(L1<=Lbirth^2,(1.01*Lbirth)^2,L1)
    L1 <- ifelse(L1>=trueLmax^2,(0.99*trueLmax)^2,L1)
    L1 <- sqrt(L1) # sqrt(Gaussian) to create left skewness
  } else if (sizedist=="right"){
    trueLsmall <- (Lbirth+trueLmean)/2 # peak for right-skewed dist
    shape.gamma <- trueLsmall/4+1 # scale=4
    L1 <- rgamma(n,shape=shape.gamma,scale=4)
    L1 <- ifelse(L1<=Lbirth,1.01*Lbirth,L1)
    L1 <- ifelse(L1>=trueLmax,0.99*trueLmax,L1)
  } else if (sizedist=="norm"){
    L1 <- rnorm(n,mean=trueLmean,sd=trueLmean*sd.L)
    L1 <- ifelse(L1<=Lbirth,1.01*Lbirth,L1)
    L1 <- ifelse(L1>=trueLmax,0.99*trueLmax,L1)
  } else if (sizedist=="bimod"){
    trueLsmall <- (Lbirth+trueLmean)/2 # peak for first peak
    trueLlarge <- (trueLmax+trueLmean)/2 # peak for second peak
    L1 <- c(rnorm(floor(n/2),mean=trueLsmall,sd=trueLmean*sd.L/2),
                   rnorm(ceiling(n/2),mean=trueLlarge,sd=trueLmean*sd.L/2))
    L1 <- ifelse(L1<=Lbirth,1.01*Lbirth,L1)
    L1 <- ifelse(L1>=trueLmax,0.99*trueLmax,L1)
  } else {stop('sizedist can only be "left", "right", "norm" or "bimod".')}
  
  # summary(L1)
  # hist(L1,freq=F,xlim=c(0,200))
  # abline(v=c(Lbirth,trueLmax))
  # abline(v=trueLinf,col='blue')
  # lines(density(L1),col='red')
  
  # if (any(L1<Lbirth) | any(L1>trueLmax)){
  #   stop('Lengths at capture generated below Lbirth or above trueLmax.')
  # }
  
  ### standardize length at cap within [0,1] to define L1
  # min.Lcap <- min(L1)
  # max.Lcap <- max(L1)
  # props <- (L1 - min.Lcap)/(max.Lcap - min.Lcap)
  # L1 <- (trueLmax*(props))/(max(props)) # bad, entirely spans [0,trueLmax]
  
  
  #/////////////////////////////////////////////////////////////////////////////
  #### STEP 2: age at cap, sim deltaT ####
  #/////////////////////////////////////////////////////////////////////////////
  
  ### calculate age at cap
  agecap <- -log(1-L1/trueLinf)/trueK + T0 # ages relative to birth
  # ^ inverted vB at cap with T0 as ref point, agecap=T1 if Tbirth=0
  
  ### simulate gamma-distributed deltaT (in years)
  # negatively correlated with agecap
  # truncate gamma dist at (1-0.01/n)th max quantile, so never exceeds trueTmax
  # only happens in 1 sample every 100 indep samples of size n, on average
  dTmaxq <- trueTmax-agecap # max quantile allowed = upper bound
  dTshape <- sapply(X=dTmaxq,FUN=function(q99vec,scale.deltaT){
    uniroot(f=function(mu,q99,scale.deltaT){
      qgamma(p=(1-0.01/n),scale=scale.deltaT,shape=mu)-q99
    }, # find shape param so that gamma quantile matches dTmaxq
    q99=q99vec,scale.deltaT=scale.deltaT,interval=c(1e-10,100))$root
    # ^ TODO: better bounds? something at user level?
  },scale.deltaT=scale.deltaT)
  deltaT <- rgamma(n,shape=dTshape,scale=scale.deltaT)+min.deltaT
  deltaT <- ifelse(deltaT>dTmaxq,dTmaxq,deltaT) # truncate at (trueTmax-agecap)

  ### compute date of recap T2 and age at recap (years)
  # T2 <- T1+deltaT*ndays # arithmetic in days # as of v0.3: deleted, no use
  agerecap <- agecap+deltaT # in years
  # if (any((3/trueK)<agerecap)){
  #   warning('Some ages at recapture exceed 3/trueK.')
  # }
  # # ^ as of v0.3: warning not useful, we now use trueTmax from vB at Lmax

  # cbind(agecap,deltaT,agerecap)[agecap%in%sort(agecap)[c(1,floor(n/2),n)],]
  # # ^ check neg corr between agecap and deltaT so that agerecap not too high
  
  
  #/////////////////////////////////////////////////////////////////////////////
  #### STEP 3: sim lengths at recap (pure vB) ####
  #/////////////////////////////////////////////////////////////////////////////
  
  ### compute L2 without measurement error or individual growth variability
  L2 <- trueLinf*(1-exp(-trueK*(agerecap-T0))) # vB, ref point is T0
  
  # plot(x=c(rep(T0,n),rep(0,n),agecap,agerecap),
  #      y=c(rep(0,n),rep(Lbirth,n),L1,L2))
  # # ^ 4 time points for each fish: T0, birth, cap, recap. All aligned on same vB
  # #   curve because all have same Linf and K param, same Lbirth, and shifted so
  # #   that birth=0.

  
  #/////////////////////////////////////////////////////////////////////////////
  #### STEP 4: sim growth variability in Linf, propagate ####
  #/////////////////////////////////////////////////////////////////////////////
  
  ### introduce individual growth variability via the growth performance index
  Linfi <- trueLinf + rnorm(n,mean=0,sd=sd.Linf*trueLinf)
  truephi <- log10(trueK)+2*log10(trueLinf)
  # ^ growth performance index linking Linf to K, Pauly
  
  ### individual Ki values based on random effect Linfi and single truephi
  # relation: 10^phi = K*Linf^2 <=> K = 10^(phi-2*log10(Linf))
  Ki <- 10^(truephi-2*log10(Linfi)) # single source of randomness through Linfi
  
  ### compute negative time diff (years) between birth and date of length 0
  # T0i <- -10^(-0.3922-0.2752*log10(Linfi)-1.038*log10(Ki)) # Pauly (1979)
  T0i <- log(1-Lbirth/Linfi)/Ki # as of v0.3: T0i from Lbirth, vB
  
  ### max lengths and age given Linfi and Ki
  Lmaxi <- 0.99*Linfi # as of v0.3
  trueTmaxi <- -log(1-Lmaxi/Linfi)/Ki + T0i # vB for Lmaxi, T0i as ref point
  # trueTmaxi <- 3/Ki
  
  ### recompute deltaT to ensure no fish exceed max age given Linfi and Ki
  agecapgv <- agecap
  ind.toolarge <- agecapgv>=(trueTmaxi-min.deltaT)
  agecapgv[ind.toolarge] <- trueTmaxi[ind.toolarge]-min.deltaT
  dTmaxqi <- trueTmaxi-agecapgv # max quantile allowed = upper bound
  dTshapei <- sapply(X=dTmaxqi,FUN=function(q99vec,scale.deltaT){
    uniroot(f=function(mu,q99,scale.deltaT){
      qgamma(p=(1-0.01/n),scale=scale.deltaT,shape=mu)-q99
    },
    q99=q99vec,scale.deltaT=scale.deltaT,interval=c(1e-10,100))$root
    # ^ TODO: better bounds? something at user level?
  },scale.deltaT=scale.deltaT)
  
  # head(cbind(trueTmaxi,dTshapei,dTshape))

  deltaTi <- rgamma(n,shape=dTshapei,scale=scale.deltaT)
  deltaTi <- ifelse(deltaTi>dTmaxqi,dTmaxqi,deltaTi)
  # ^ truncate ages at recap at trueTmaxi
  deltaTi <- ifelse(deltaTi<min.deltaT,min.deltaT,deltaTi)
  
  # plot(deltaT,deltaTi,xlim=c(0,20),ylim=c(0,20))
  # abline(0,1,col='red')
  # 
  # par(mfrow=c(2,1))
  # hist(deltaT,freq=F,ylim=c(0,0.5),xlim=c(0,20))
  # lines(density(deltaT),col='red')
  # hist(deltaTi,freq=F,ylim=c(0,0.5),xlim=c(0,20))
  # lines(density(deltaTi),col='red')
  # par(mfrow=c(1,1))
  
  ### compute L1 and L2 with growth variability, based on ages with gv
  agerecapgv <- agecapgv+deltaTi # in years
  # if (any((3/Ki)<agerecapgv)){
  #   warning('Some ages at recapture exceed 3/Ki.')
  # }
  # ^ as of v0.3: warning not useful, we now use trueTmaxi from vB at Lmax
  
  L1gv <- Linfi*(1-exp(-Ki*(agecapgv-T0i)))
  L2gv <- Linfi*(1-exp(-Ki*(agerecapgv-T0i)))
  # ^ vB at cap and recap, randomness in Linf propagated everywhere

  # plot(x=c(T0i,rep(0,n),agecapgv,agerecapgv),
  #      y=c(rep(0,n),rep(Lbirth,n),L1gv,L2gv))
  # # ^ 4 time points for each fish: T0, birth, cap, recap. Different vB curves
  # #   now because different Linf and K param, but same Lbirth and birth=0 so
  # #   all T0i are different.
  
  ### set L2gv as NA if negative growth, keep same lengths as L1 and L2
  L2gv <- ifelse(L2gv<L1gv,NA,L2gv) # set negative growth as NA
  # ^ rare, but can happen if large sd.Linf and small deltaTi
  
  
  #/////////////////////////////////////////////////////////////////////////////
  #### STEP 5: sim measurement error ####
  #/////////////////////////////////////////////////////////////////////////////
  
  L1.error <- rnorm(n,mean=0,sd=sd.L1*L1)
  L1e <- L1 + L1.error
  L2.error <- rnorm(n,mean=0,sd=sd.L2*L2)
  L2e <- L2 + L2.error
  
  ### set L2e as NA if negative growth, keep same lengths as L1 and L2
  # L2e <- L1e+abs(L2e-L1e) # overwrite, force Lrecap > Lcap
  L1e.clean <- L1e # simpler to leave it full
  ind.ng <- which(L2e<L1e) # index of obs with negative growth
  if (length(ind.ng)==0L){ # no negative growth, index is a length zero vector
    L2e.clean <- L2e
  } else {
    L2e.clean <- ifelse(L2e<L1e,NA,L2e) # L2e[-ind.ng]
  }
  
  
  #/////////////////////////////////////////////////////////////////////////////
  #### STEP 6: add measurement error to lengths with growth variability ####
  #/////////////////////////////////////////////////////////////////////////////
  
  ### compute L1 and L2 with both growth variability and measurement error
  L1sim <- L1gv + L1.error # same meas err as in L1e
  L2sim <- L2gv + L2.error # same meas err as in L2e, NAs propagated
  # L2sim <- L1sim+abs(L2sim-L1sim) # overwrite, force Lrecap > Lcap
  
  ### set L2sim as NA if negative growth, keep same lengths as L1 and L2
  L1sim.clean <- L1sim
  ind.ngsim <- which(L2sim<L1sim) # index of obs with negative growth
  if (length(ind.ngsim)==0L){
    L2sim.clean <- L2sim
  } else {
    L2sim.clean <- ifelse(L2sim<L1sim,NA,L2sim) # L2sim[-ind.ngsim]
  }
  
  
  #/////////////////////////////////////////////////////////////////////////////
  #### Output ####
  #/////////////////////////////////////////////////////////////////////////////
  
  return(list('Lcap'=L1sim,'Lrecap'=L2sim, # gv+me
              'Lcap.clean'=L1sim.clean,'Lrecap.clean'=L2sim.clean, # gv+me, NAs
              'eLcap'=L1e,'eLrecap'=L2e, # no gv, me
              'eLcap.clean'=L1e.clean,'eLrecap.clean'=L2e.clean, # no gv, me, NAs
              'gvLcap'=L1gv,'gvLrecap'=L2gv, # gv, no me
              'trueLcap'=L1,'trueLrecap'=L2, # neither gv nor me
              'trueLinf'=trueLinf,'trueK'=trueK,'phi'=truephi, # scalars
              'Linfi'=Linfi,'Ki'=Ki, # random effects values, both vectors
              # 'Tcap'=as.numeric(T1)/ndays,'Trecap'=as.numeric(T2)/ndays,
              'deltaT'=deltaTi,'agecap'=agecapgv,'agerecap'=agerecapgv,
              'truedeltaT'=deltaT,'trueagecap'=agecap,'trueagerecap'=agerecap,
              'T0'=T0,'T0i'=T0i,'Lbirth'=Lbirth
              # ^ TODO: document output
  ))
}

# END CapRecapSim
