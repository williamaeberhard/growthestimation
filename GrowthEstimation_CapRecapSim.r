CapRecapSim <- function(n,trueLinf=123.5,trueK=0.146,
                        sizedist='norm',sd.L=0.3,
                        T1=NULL,scale.deltaT=1,min.deltaT=0.25,
                        sd.Linf=0.067,sd.L1=0.033,sd.L2=0.067){
  
  #/////////////////////////////////////////////////////////////////////////////
  #### Documentation v0.2.1 #### 
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
  #' @param sizedist distibution of the length at capture; choice among
  #' \code{"norm"} (Gaussian, the default), \code{"left"} (left-skewed),
  #' \code{"right"} (left-skewed) and \code{"bimod"} (bimodal, mixture of two
  #' normals).
  #' @param sd.L positive scalar, sd of generated lengths at capture.
  #' @param T1 vector of dates of capture (length \code{n}); if left
  #' \code{NULL}, then arbitrarily set to January 1st, 1990 for all fish.
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
  #' value (0.033) yields 10% relative error on average.
  #' @param sd.L2 positive scalar, multiplicative coefficient for measurement
  #' error sd of L2 (Gaussian additive measurement noise). If set large and
  #' \code{deltaT} simulated small enough, then observed length at recapture
  #' can be smaller than observed length at capture (negative growth). Default
  #' value (0.067) yields 20% relative error on average.
  #' @details Randomness in: \code{Linf} (Gaussian); \code{K} (Gaussian)
  #' thorugh randomness in the slope of the growth performance function; 
  #' length at capture (\code{sizedist}); \code{deltaT} (log-normal, negatively
  #' correlated with age at capture); and measurement error in observed lengths
  #' at capture and recapture (Gaussian). TODO: more details.
  #' @return A list with the following objects:
  #' \itemize{
  #'   \item{\code{L1}}{numeric vector of length}
  #'   \item{TODO}{TODO.}
  #' }
  
  # help(document) # document::document
  
  #/////////////////////////////////////////////////////////////////////////////
  #### Define true parameters and setup #### 
  #/////////////////////////////////////////////////////////////////////////////
  
  ### compute negative time diff (in years) between birth and date of length 0
  # if (is.null(t0)){
  #   t0 <- -10^(-0.3922-0.2752*log10(trueLinf)-1.038*log10(trueK)) # Pauly (1979)
  # } else {
  #   t0 <- -abs(t0) # ensure it is negative even if supplied t0 is positive
  # } # ^ wasn't working with t0i, back to Pauly and no user-level argument t0
  t0 <- -10^(-0.3922-0.2752*log10(trueLinf)-1.038*log10(trueK)) # Pauly (1979)
  # ^ birth assumed at t=0 so that t0 is negative
  
  ### lengths at birth, max and mean
  trueLB <- trueLinf*(1-exp(trueK*t0)) # length at birth
  trueLmax <- 0.95*trueLinf # trueLinf*(1-exp(-trueK*(truetmax-t0)))
  trueLmean <- (trueLmax+trueLB)/2
  
  ### maximum age
  truetmax <- t0-log(1-trueLmax/trueLinf)/trueK # VB for trueLmax
  # truetmax <- 3/trueK
  
  ### number of days in a year, for years-days conversion
  ndays <- 365.2425
  
  ### set date of capture T1 (not random), same for all by default
  if (is.null(T1)){T1 <- rep(as.Date("1990-01-01"),n)}
  # ^ TODO: if T1 supplied, make some checks on date format
  
  #/////////////////////////////////////////////////////////////////////////////
  ### STEP 1 ####
  #/////////////////////////////////////////////////////////////////////////////
  
  ### distribution for length at cap
  if (sizedist=="left"){
    trueLlarge <- (trueLmax+trueLmean)/2 # peak for second peak
    Lcap.dist <- rnorm(n,mean=trueLlarge^2,sd=trueLlarge^2*sd.L)
    Lcap.dist <- ifelse(Lcap.dist<=trueLB^2,1.01*trueLB^2,Lcap.dist)
    Lcap.dist <- ifelse(Lcap.dist>=trueLmax^2,0.99*trueLmax^2,Lcap.dist)
    Lcap.dist <- sqrt(Lcap.dist)
  } else if (sizedist=="right"){
    trueLsmall <- (trueLB+trueLmean)/2 # peak for right-skewed dist
    shape.gamma <- trueLsmall/4+1 # scale=4
    Lcap.dist <- rgamma(n,shape=shape.gamma,scale=4)
    Lcap.dist <- ifelse(Lcap.dist<=trueLB,1.01*trueLB,Lcap.dist)
    Lcap.dist <- ifelse(Lcap.dist>=trueLmax,0.99*trueLmax,Lcap.dist)
  } else if (sizedist=="norm"){
    Lcap.dist <- rnorm(n,mean=trueLmean,sd=trueLmean*sd.L)
    Lcap.dist <- ifelse(Lcap.dist<=trueLB,1.01*trueLB,Lcap.dist)
    Lcap.dist <- ifelse(Lcap.dist>=trueLmax,0.99*trueLmax,Lcap.dist)
  } else if (sizedist=="bimod"){
    trueLsmall <- (trueLB+trueLmean)/2 # peak for first peak
    trueLlarge <- (trueLmax+trueLmean)/2 # peak for second peak
    Lcap.dist <- c(rnorm(floor(n/2),mean=trueLsmall,sd=trueLmean*sd.L/2),
                   rnorm(ceiling(n/2),mean=trueLlarge,sd=trueLmean*sd.L/2))
    Lcap.dist <- ifelse(Lcap.dist<=trueLB,1.01*trueLB,Lcap.dist)
    Lcap.dist <- ifelse(Lcap.dist>=trueLmax,0.99*trueLmax,Lcap.dist)
    # summary(Lcap.dist)
    # hist(Lcap.dist,freq=F,xlim=c(0,150))
    # abline(v=c(trueLB,trueLmax))
    # lines(density(Lcap.dist),col='red')
  } else {stop('sizedist can only be "left", "right", "norm" or "bimod".')}
  
  # if (any(Lcap.dist<trueLB) | any(Lcap.dist>trueLmax)){
  #   stop('Lengths at capture generated below trueLB or above trueLmax.')
  # }
  
  ### standardize length at cap within [0,1] to define L1
  # min.Lcap <- min(Lcap.dist)
  # max.Lcap <- max(Lcap.dist)
  # props <- (Lcap.dist - min.Lcap)/(max.Lcap - min.Lcap)
  # L1 <- (trueLmax*(props))/(max(props)) # bad, entirely spans [0,trueLmax]
  L1 <- Lcap.dist # TODO: rename from the beginning, with notation cleanup
  
  
  #/////////////////////////////////////////////////////////////////////////////
  #### STEP 2 ####
  #/////////////////////////////////////////////////////////////////////////////
  
  ### calculate age at cap
  agecap <- -(log(1-L1/trueLinf))/trueK+t0 # positive ages
  
  ### simulate gamma-distributed deltaT (in years)
  # negatively correlated with agecap
  # truncate gamma dist at (1-0.01/n)th max quantile, so never exceeds truetmax
  # only happens in 1 sample every 100 indep samples of size n, on average
  dTmaxq <- truetmax-agecap # max quantile allowed = upper bound
  dTshape <- sapply(X=dTmaxq,FUN=function(q99vec,scale.deltaT){
    uniroot(f=function(mu,q99,scale.deltaT){
      qgamma(p=(1-0.01/n),scale=scale.deltaT,shape=mu)-q99
    }, # find shape param so that gamma quantile matches dTmaxq
    q99=q99vec,scale.deltaT=scale.deltaT,interval=c(1e-10,100))$root
    # ^ TODO: better bounds? something at user level?
  },scale.deltaT=scale.deltaT)
  deltaT <- rgamma(n,shape=dTshape,scale=scale.deltaT)+min.deltaT
  deltaT <- ifelse(deltaT>dTmaxq,dTmaxq,deltaT) # truncate at (truetmax-agecap)

  ### compute date of recap T2 and age at recap (years)
  T2 <- T1+deltaT*ndays # arithmetic in days
  agerecap <- agecap+deltaT # in years
  if (any((3/trueK)<agerecap)){
    warning('Some ages at recapture exceed 3/trueK.')
  }
  
  #/////////////////////////////////////////////////////////////////////////////
  #### STEP 3 ####
  #/////////////////////////////////////////////////////////////////////////////
  
  ### compute L2 without error or indivifual growth variability
  L2 <- trueLinf*(1-exp(-trueK*(agerecap-t0)))
  
  #/////////////////////////////////////////////////////////////////////////////
  #### STEP 4 ####
  #/////////////////////////////////////////////////////////////////////////////
  
  ### introduce individual growth variability via the growth performance index
  Linfi <- trueLinf+rnorm(n,mean=0,sd=sd.Linf*trueLinf)
  truephi <- log10(trueK)+2*log10(trueLinf) # true growth performance index, Pauly
  
  ### individual Ki values based on random effect Linfi and single truephi
  # relation: 10^phi = K*Linf^2 <=> K = 10^(phi-2*log10(Linf))
  Ki <- 10^(truephi-2*log10(Linfi)) # single source of randomness through Linfi
  
  ### compute negative time diff (years) between birth and date of length 0
  t0i <- -10^(-0.3922-0.2752*log10(Linfi)-1.038*log10(Ki)) # Pauly (1979)
  
  ### max lengths and age given Linfi and Ki
  Lmaxi <- 0.95*Linfi
  truetmaxi <- t0i-log(1-Lmaxi/Linfi)/Ki # VB for Lmaxi
  # truetmaxi <- 3/Ki
  
  ### recompute deltaT to ensure no fish exceed max age given Linfi and Ki
  agecapgv <- agecap
  ind.toolarge <- agecapgv>=(truetmaxi-min.deltaT)
  agecapgv[ind.toolarge] <- truetmaxi[ind.toolarge]-min.deltaT
  dTmaxqi <- truetmaxi-agecapgv # max quantile allowed = upper bound
  dTshapei <- sapply(X=dTmaxqi,FUN=function(q99vec,scale.deltaT){
    uniroot(f=function(mu,q99,scale.deltaT){
      qgamma(p=(1-0.01/n),scale=scale.deltaT,shape=mu)-q99
    },
    q99=q99vec,scale.deltaT=scale.deltaT,interval=c(1e-10,100))$root
    # ^ TODO: better bounds? something at user level?
  },scale.deltaT=scale.deltaT)

  deltaTi <- rgamma(n,shape=dTshapei,scale=scale.deltaT)
  deltaTi <- ifelse(deltaTi>dTmaxqi,dTmaxqi,deltaTi)
  # ^ truncate ages at recap at truetmaxi
  deltaTi <- ifelse(deltaTi<min.deltaT,min.deltaT,deltaTi)
  
  # deltaT <- rgamma(n,shape=dTshape,scale=scale.deltaT)+min.deltaT
  # deltaTi <- rgamma(n,shape=dTshapei,scale=scale.deltaT)+min.deltaT
  # plot(deltaT,deltaTi)
  # par(mfrow=c(2,1))
  # hist(deltaT,freq=F,ylim=c(0,0.5),xlim=c(0,20))
  # lines(density(deltaT))
  # hist(deltaTi,freq=F,ylim=c(0,0.5),xlim=c(0,20))
  # lines(density(deltaTi))
  # par(mfrow=c(1,1))
  
  
  ### compute L1 and L2 with growth variability, based on agerecapgv
  agerecapgv <- agecapgv+deltaTi # in years
  if (any((3/Ki)<agerecapgv)){
    stop('Some ages at recapture exceed 3/Ki.')
  }
  L1gv <- Linfi*(1-exp(-Ki*(agecapgv-t0i)))
  L2gv <- Linfi*(1-exp(-Ki*(agerecapgv-t0i)))
  
  ### set L2gv as NA if negative growth, keep same lengths as L1 and L2
  L2gv <- ifelse(L2gv<L1gv,NA,L2gv) # set negative growth as NA
  # ^ rare, but can happen if large sd.Linf and small deltaTi
  
  
  #/////////////////////////////////////////////////////////////////////////////
  #### STEP 5 ####
  #/////////////////////////////////////////////////////////////////////////////
  
  ### compute L1 and L2 with measurement error
  L1ee <- rnorm(n,mean=0,sd=sd.L1*L1)
  L1e <- L1+L1ee
  L2ee <- rnorm(n,mean=0,sd=sd.L2*L2)
  L2e <- L2+L2ee
  
  ### set L2e as NA if negative growth, keep same lengths as L1 and L2
  # L2e <- L1e+abs(L2e-L1e) # overwrite, force Lrecap > Lcap
  L1e.clean <- L1e # simpler to leave it full
  ng.ind <- which(L2e<L1e) # index of obs with negative growth
  if (length(ng.ind)==0L){ # no negative growth, index is a length zero vector
    L2e.clean <- L2e
  } else {
    L2e.clean <- ifelse(L2e<L1e,NA,L2e) # L2e[-ng.ind]
  }
  
  
  #/////////////////////////////////////////////////////////////////////////////
  #### STEP 6 ####
  #/////////////////////////////////////////////////////////////////////////////
  
  ### compute L1 and L2 with both growth variability and measurement error
  # L1sim <- L1gv+sample(c(1,-1),n,replace=TRUE)*L1ee # errors can be pos or neg
  L1sim <- L1gv+L1ee
  L2sim <- L2gv+L2ee # NAs should be propagated
  # L2sim <- L1sim+abs(L2sim-L1sim) # overwrite, force Lrecap > Lcap
  
  ### set L2sim as NA if negative growth, keep same lengths as L1 and L2
  L1sim.clean <- L1sim
  ngsim.ind <- which(L2sim<L1sim) # index of obs with negative growth
  if (length(ngsim.ind)==0L){
    L2sim.clean <- L2sim
  } else {
    L2sim.clean <- ifelse(L2sim<L1sim,NA,L2sim) # L2sim[-ngsim.ind]
  }
  
  ### extra: length at birth with gv, date of birth TB, date of 0 length T0
  # LBi <- Linfi*(1-exp(Ki*t0i))
  # TB <- T1-agecap*ndays # arithmetic in days
  # T0 <- TB-(t0i*ndays) # arithmetic in days 
  
  #/////////////////////////////////////////////////////////////////////////////
  #### Output ####
  #/////////////////////////////////////////////////////////////////////////////
  
  return(list('Lcap'=L1sim,'Lrecap'=L2sim, # gv and meas err
              'Lcap.clean'=L1sim.clean,'Lrecap.clean'=L2sim.clean, # gv and meas err
              'eLcap'=L1e,'eLrecap'=L2e, # no gv but meas err
              'eLcap.clean'=L1e.clean,'eLrecap.clean'=L2e.clean, # no gv, meas err
              'gvLcap'=L1gv,'gvLrecap'=L2gv, # gv, no meas err
              'trueLcap'=L1,'trueLrecap'=L2, # neither gv nor meas err
              'trueLinf'=trueLinf,'trueK'=trueK,'theta'=truephi, # scalars
              'Linfi'=Linfi,'Ki'=Ki, # random effects values, both vectors
              'Tcap'=as.numeric(T1)/ndays,'Trecap'=as.numeric(T2)/ndays,
              'deltaT'=deltaTi,'agecap'=agecapgv,'agerecap'=agerecapgv,
              'truedeltaT'=deltaT,'trueagecap'=agecap,'trueagerecap'=agerecap,
              'T0'=t0,'LB'=trueLB
              # ^ TODO: cleanup output and document
              # 'Tbirth'=as.numeric(TB)/ndays,'T0'=as.numeric(T0)/ndays)
  ))
}
# END CapRecapSim
